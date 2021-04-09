This repository contains my KD-tree implementation utilizing the Nori codebase described bellow. The KD-tree construction utilizes a SAH heuristic and largely follows the
algorithms provided by the paper “On building fast kd-Trees for Ray Tracing, and on doing that
in O(N log N)” which will be linked at the bottom of this report. Similarly my traversal utilizes a
standard kd traversal given by the paper entitled “Understanding the Efficiency of kD-tree
Ray-Traversal Techniques over a GPGPU Architecture”.

http://www.irisa.fr/prive/kadi/Sujets_CTR/kadi/Kadi_sujet2_article_Kdtree.pdf
https://www.cin.ufpe.br/~als3/saap/ArturLiraDosSantos-ArtigoIJPP.pdf


[![CS440 Banner](https://rgl.s3.eu-central-1.amazonaws.com/media/uploads/wjakob/2017/02/16/cs440-logo_web.jpg)](https://rgl.s3.eu-central-1.amazonaws.com/media/uploads/wjakob/2017/02/20/cs440-rgl.jpg)

## Nori Version 2
![Build status](https://github.com/wjakob/nori/workflows/Build/badge.svg)

Nori is a simple ray tracer written in C++. It runs on Windows, Linux, and
Mac OS and provides basic functionality that is required to complete the
assignments in the course Advanced Computer Graphics taught at EPFL.

### Course information and framework documentation

For access to course information including slides and reading material, visit the main [Advanced Computer Graphics website](https://rgl.epfl.ch/courses/ACG17). The Nori 2 framework and coding assignments will be described on the [Nori website](https://wjakob.github.io/nori).

### Note to researchers and students from other institutions

Last year's version of Nori including a full set of assignment descriptions can
be found in the following [archive](https://github.com/wjakob/nori-old).


### Known Issues
There is a known issue with the NanoGUI version that Nori uses: on Linux systems with an integrated Intel GPU, a bug in the Mesa graphics drivers causes the GUI to freeze on startup. A workaround is to temporarily switch to an older Mesa driver to run Nori. This can be done by running
```
export MESA_LOADER_DRIVER_OVERRIDE=i965
```
