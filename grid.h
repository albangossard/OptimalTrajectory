#ifndef H_GRID
#define H_GRID

#include <iostream>
#include <list>
#include <math.h>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>

using namespace std;


struct Point{
    double x,y;
};
struct PointPos{
    double x,y,t;
};

class Grid{
private:
    int Nx, Ny;
    double Dx, Dy;
    struct ProjectionPoint{
        double val;
        double t;
    };
    double* u;
    double* r0;
    int VERBOSE;
public:
    Grid(int Nx, int Ny, double* density, int VERBOSE);
    ~Grid();
    void intersectSegment(Point* pt1, Point* pt2, list<PointPos>& intersectPt);
    double computeCost(std::list<Point>& pts);
    void print_r0();
};

#endif