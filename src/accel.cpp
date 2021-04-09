/*
Rodrigo Ratto
*/

#include <nori/accel.h>
#include <Eigen/Geometry>
#include <algorithm>
#include <stack>
#include <array>

const int MAX_DEPTH = 2;
int nodeNum = 0;

NORI_NAMESPACE_BEGIN


void Accel::addMesh(Mesh* mesh) {
    if (m_mesh)
        throw NoriException("Accel: only a single mesh is supported!");
    m_mesh = mesh;
    m_bbox = m_mesh->getBoundingBox();
}

Knode* kd_root;
float Kt = 1.0; //traversal step
float Ki = 2.0; // intersection step
int numberT = 0;

float lambda(int TL, int TR, float probL, float probR) {
    
    if ((TL == 0 || TR == 0) && !(probL == 1 || probR == 1)) // if probability == 1 go with 1.0 value, not a good split
        return 0.8f;
    return 1.0f;
}

float cost(float PL, float PR, int NL, int NR) {
    return(lambda(NL, NR, PL, PR) * (Ki + Kt * (PL * NL + PR * NR)));
}

void Accel::splitBox(const BoundingBox3f& bbox, const SplitPlane& p, BoundingBox3f& boxL, BoundingBox3f& boxR) const {
    boxL = boxR = bbox;
    boxL.max[p.axis] = p.pos;
    boxR.min[p.axis] = p.pos;
}

void Accel::SAH(const SplitPlane& p, const BoundingBox3f& bbox, int NL, int NR, int NP, float& Cp, int& pside) const {

    BoundingBox3f boxL, boxR;
    splitBox(bbox, p, boxL, boxR);

    float PL = boxL.getSurfaceArea() / bbox.getSurfaceArea();
    float PR = boxR.getSurfaceArea() / bbox.getSurfaceArea();
    float costL, costR;

    costL = cost(PL, PR, NL + NP, NR);
    costR = cost(PL, PR, NL, NP + NR);
    if (costL < costR) {
        Cp = costL;
        pside = 1; //left == 1;
    }
    else {
        Cp = costR;
        pside = 2; //right == 2;
    }
}

bool sortEvents(Event& lhs, Event& rhs) {

    return ((lhs.splitPlane.pos < rhs.splitPlane.pos) || (lhs.splitPlane.pos == rhs.splitPlane.pos && lhs.event < rhs.event));
}

void Accel::findPlane(const std::vector<int>& triangles, const BoundingBox3f& bbox, SplitPlane& pFinal, float& Cfinal, int& psideFinal) const {
 
    Cfinal = INFINITY;

    for (int k = 0; k < 3; ++k) {

        std::vector<Event> events;

        for (int i = 0; i < triangles.size(); i++) {

            BoundingBox3f tribox = m_mesh->getBoundingBox(triangles[i]);
            tribox.clip(bbox);

            Point3f min = tribox.getMin();
            Point3f max = tribox.getMax();
            float x = max[0] - min[0];
            float y = max[1] - min[1];
            float z = max[2] - min[2];
            if (x <= 0.1 || y <= 0.1 || z <= 0.1) { //test for planar
                events.push_back(Event(triangles[i], tribox.min[k], k , 1));
            }
            else {
                events.push_back(Event(triangles[i], tribox.min[k], k, 2));
                events.push_back(Event(triangles[i], tribox.max[k], k, 1));
            }
        }
        sort(events.begin(), events.end(), sortEvents);
        int NL = 0, NP = 0, NR = triangles.size();

        for (std::vector<Event>::size_type event = 0; event < events.size(); ++event) {

            SplitPlane& p = events[event].splitPlane;
            int pPlus = 0, pMinus = 0, pPlanar = 0;
            while (event < events.size() && events[event].splitPlane.pos == p.pos && events[event].event == 0) {
                ++pMinus;
                event++;
            }
            while (event < events.size() && events[event].splitPlane.pos == p.pos && events[event].event == 1) {
                ++pPlanar;
                event++;
            }
            while (event < events.size() && events[event].splitPlane.pos == p.pos && events[event].event== 2) {
                ++pPlus;
                event++;
            }
            NP = pPlanar;
            NR -= pPlanar;
            NR -= pMinus;

            float C;
            int pside = 0;

            SAH(p, bbox, NL, NR, NP, C, pside);
            if (C < Cfinal) {
                Cfinal = C;
                pFinal = p;
                psideFinal = pside;
            }
            NL += pPlus;
            NL += pPlanar;
            NP = 0;
        }
    }
}

Knode* Accel::KDBuild(std::vector<int> triangles, const BoundingBox3f& bbox, int depth, const SplitPlane& lastP) {

    float Cp;
    SplitPlane p;
    int pside;
    findPlane(triangles, bbox, p, Cp, pside); //find plane before terminating if to compare next plane with previous

    if (p == lastP) //delimiter to help with terminating 
        depth++;
     
    if (Cp > Ki * triangles.size() || depth == 1) {
  
        Knode* node = new Knode();
        node->tri = triangles;
        node->isLeaf = true;
        numberT += triangles.size();

        //std::cout << "Number of triangles: " << triangles.size() << "  total triangles: " << numberT << std::endl
        //if (depth > 1)
            //std::cout << " Limiter hit " << std::endl;
        return node;
    }
    BoundingBox3f VL, VR;
    std::vector<int> TL, TR;
    Knode* node = new Knode();
    splitBox(bbox, p, VL, VR); 

    for (int i = 0; i < triangles.size(); i++) {
        BoundingBox3f tribox = m_mesh->getBoundingBox(triangles[i]);
        if (tribox.overlaps(VL))
            TL.push_back(triangles[i]);
        if (tribox.overlaps(VR))
            TR.push_back(triangles[i]);
    }
    
    node->splitPlane = p;
    node->left = KDBuild(TL, VL, depth , p);
    node->right = KDBuild(TR, VR, depth , p);
    //nodeNum++;
    //nodeNum++;
    return node;
}

void Accel::build() {

    std::vector<int> triangles;
    for (uint32_t idx = 0; idx < m_mesh->getTriangleCount(); ++idx) {
        triangles.push_back(idx);
    }
    SplitPlane p;
    kd_root = this->KDBuild(triangles, m_bbox, 0, p);
    //buffer << "Number of nodes: " << nodeNum << std::endl;
}

bool Accel::rayIntersect(const Ray3f& ray_, Intersection& its, bool shadowRay) {
    bool foundIntersection = false;  // Was an intersection found so far?
    uint32_t f = (uint32_t)-1;      // Triangle index of the closest intersection
    Ray3f ray(ray_); /// Make a copy of the ray (we will need to update its '.maxt' value)

    float t;
    std::stack<Node*> s;;
    float tMin;
    float tMax;
    m_bbox.rayIntersect(ray, tMin, tMax);
    Knode* currentNode = kd_root;
    bool notDone = true;

    while (notDone) {

        while (!currentNode->isLeaf) {
            int axis = currentNode->splitPlane.axis;
            float diff = currentNode->splitPlane.pos - ray.o[axis];
            float tDist = diff / ray.d[axis];

            Knode* near, * far;
            if (ray.o[currentNode->splitPlane.axis] < currentNode->splitPlane.pos) {
                near = currentNode->left;
                far = currentNode->right;
            }
            else {
                near = currentNode->right;
                far = currentNode->left;
            }

            if (tDist < 0 || tDist >= tMax)
                currentNode = near;
            else {
                
                if (tDist <= tMin)
                    currentNode = far;

                else {
                    Node* node = new Node();
                    node->node = far;
                    node->t = tMax;
                    s.push(node);
                    currentNode = near;
                    tMax = tDist;
                }
            }
        }

        if (currentNode != nullptr && currentNode->tri.size() > 0)
            for (uint32_t idx = 0; idx < currentNode->tri.size(); ++idx) {
                float u, v, t;
                if (m_mesh->rayIntersect(currentNode->tri[idx], ray, u, v, t)) {
                    /* An intersection was found! Can terminate
                    immediately if this is a shadow ray query */
                    if (shadowRay)
                        return true;
                    ray.maxt = its.t = t;
                    its.uv = Point2f(u, v);
                    its.mesh = m_mesh;
                    f = currentNode->tri[idx];
                    foundIntersection = true;
                    notDone = false;
                }
            }

        if (s.empty())
            break;

        tMin = tMax;
        currentNode = s.top()->node;
        tMax = s.top()->t;
        s.pop();
    }

    if (foundIntersection) {
        /* At this point, we now know that there is an intersection,
           and we know the triangle index of the closest such intersection.

           The following computes a number of additional properties which
           characterize the intersection (normals, texture coordinates, etc..)
        */

        /* Find the barycentric coordinates */
        Vector3f bary;
        bary << 1 - its.uv.sum(), its.uv;

        /* References to all relevant mesh buffers */
        const Mesh* mesh = its.mesh;
        const MatrixXf& V = mesh->getVertexPositions();
        const MatrixXf& N = mesh->getVertexNormals();
        const MatrixXf& UV = mesh->getVertexTexCoords();
        const MatrixXu& F = mesh->getIndices();

        /* Vertex indices of the triangle */
        uint32_t idx0 = F(0, f), idx1 = F(1, f), idx2 = F(2, f);

        Point3f p0 = V.col(idx0), p1 = V.col(idx1), p2 = V.col(idx2);

        /* Compute the intersection positon accurately
           using barycentric coordinates */
        its.p = bary.x() * p0 + bary.y() * p1 + bary.z() * p2;

        /* Compute proper texture coordinates if provided by the mesh */
        if (UV.size() > 0)
            its.uv = bary.x() * UV.col(idx0) +
            bary.y() * UV.col(idx1) +
            bary.z() * UV.col(idx2);

        /* Compute the geometry frame */
        its.geoFrame = Frame((p1 - p0).cross(p2 - p0).normalized());

        if (N.size() > 0) {
            /* Compute the shading frame. Note that for simplicity,
               the current implementation doesn't attempt to provide
               tangents that are continuous across the surface. That
               means that this code will need to be modified to be able
               use anisotropic BRDFs, which need tangent continuity */

            its.shFrame = Frame(
                (bary.x() * N.col(idx0) +
                    bary.y() * N.col(idx1) +
                    bary.z() * N.col(idx2)).normalized());
        }
        else {
            its.shFrame = its.geoFrame;
        }
    }

    return foundIntersection;
}

NORI_NAMESPACE_END

