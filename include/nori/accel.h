/*
Name: Rodrigo Ratto
Student ID: 89526973
*/

#pragma once

#include <nori/mesh.h>
#include <stack>

NORI_NAMESPACE_BEGIN


/**
 * \brief Acceleration data structure for ray intersection queries
 *
 * The current implementation falls back to a brute force loop
 * through the geometry.
 */

struct SplitPlane {
    SplitPlane() {};
    SplitPlane(const int axis, const float pos) : axis(axis), pos(pos) {};

    int axis;
    float pos;

    bool operator==(const SplitPlane& sp) {
        return(axis == sp.axis && pos == sp.pos);
    }
};

struct Knode {
    Knode() {};
    Knode(const bool isLeaf, Knode* left,  Knode* right, std::vector<int> tri, SplitPlane splitPlane) :
        isLeaf(isLeaf), left(left), right(right), tri(tri), splitPlane(splitPlane) {};

    bool isLeaf = false;
    Knode* left;
    Knode* right;
    std::vector<int> tri;
    SplitPlane splitPlane;
};

struct Node { //Node for stack
    Knode* node;
    float t;
};

struct Event {
    
    SplitPlane splitPlane;
    int event; // end event =0, planar = 1, start event = 2
    int triangle;

    Event(int tri, float pos, int axis, int event) : triangle(tri), event(event), splitPlane(SplitPlane(axis, pos)) {}
};

class Accel {
public:
    /**
     * \brief Register a triangle mesh for inclusion in the acceleration
     * data structure
     *
     * This function can only be used before \ref build() is called
     */
    void addMesh(Mesh* mesh);

    /// Build the acceleration data structure (currently a no-op)
    //std::stringstream buffer;
    //std::stringstream errorb;

    void build();

    Knode* KDBuild(std::vector<int> triangles, const BoundingBox3f& bbox, int depth, const SplitPlane& lastP);
    void splitBox(const BoundingBox3f& V, const SplitPlane& p, BoundingBox3f& VL, BoundingBox3f& VR) const;
    void SAH(const SplitPlane& p, const BoundingBox3f& V, int NL, int NR, int NP, float& CP, int& pside) const;
    void findPlane(const std::vector<int>& triangles, const BoundingBox3f& bbox, 
        SplitPlane& pFinal, float& Cfinal, int& psideFinal) const;

    /// Return an axis-aligned box that bounds the scene
    const BoundingBox3f& getBoundingBox() const { return m_bbox; }

    /**
     * \brief Intersect a ray against all triangles stored in the scene and
     * return detailed intersection information
     *
     * \param ray
     *    A 3-dimensional ray data structure with minimum/maximum extent
     *    information
     *
     * \param its
     *    A detailed intersection record, which will be filled by the
     *    intersection query
     *
     * \param shadowRay
     *    \c true if this is a shadow ray query, i.e. a query that only aims to
     *    find out whether the ray is blocked or not without returning detailed
     *    intersection information.
     *
     * \return \c true if an intersection was found
     */
    bool rayIntersect(const Ray3f& ray, Intersection& its, bool shadowRay);

private:
    Mesh* m_mesh = nullptr; ///< Mesh (only a single one for now)
    BoundingBox3f m_bbox;           ///< Bounding box of the entire scene
};

NORI_NAMESPACE_END
