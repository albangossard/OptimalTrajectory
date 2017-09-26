#ifndef H_SET
#define H_SET

#include <iostream>
#include <list>
#include <fstream>

using namespace std;


const double PI = 3.14159265358979323846264338327950288419716939937510582097494459231;

struct Integration_result{
    double cost;
};

class Intersection{
    private:
        int VERBOSE;
    public:
        list<Point> points;
        list<Point> mesh;
        Intersection(list<Point>& pts, int VERBOSE);
        void generate_mesh(double epsPath);
        void print_mesh();
};

// class PathFinding{
//     public:
//         PathFinding(list<Point>& pts);
//         ~PathFinding();
//         Intersection * intersec;
// };

#endif