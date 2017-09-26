#include "grid.h"
#include "set.h"


Intersection::Intersection(list<Point>& pts, int VERBOSE){
    this->VERBOSE=VERBOSE;
    this->points=pts;
};
void Intersection::generate_mesh(double epsPath){
    double epaisseur_trait_parcours=epsPath;

    list<Point> listForward;
    list<Point> listBackward;
    Point tmp;
    list<Point>::iterator it;
    list<Point>::iterator itPrev;
    list<Point>::iterator itNext;
    itPrev=this->points.begin();
    itNext=this->points.begin();
    for (it=this->points.begin(); it!=--this->points.end(); ++it){
        ++itNext;
        if (it!=this->points.begin()){
            if (it!=++this->points.begin()){
                ++itPrev;
            }
            double aX=it->x-itPrev->x;
            double aY=it->y-itPrev->y;
            double bX=-(itNext->x-it->x);
            double bY=-(itNext->y-it->y);
            aX=aX/sqrt(aX*aX+aY*aY);
            aY=aY/sqrt(aX*aX+aY*aY);
            bX=bX/sqrt(bX*bX+bY*bY);
            bY=bY/sqrt(bX*bX+bY*bY);
            double nX=aX+bX;
            double nY=aY+bY;
            nX=nX/sqrt(nX*nX+nY*nY);
            nY=nY/sqrt(nX*nX+nY*nY);
            double alpha = atan(aY/aX)-atan(bY/bX);
            // cout<<"alpha="<<alpha<<endl;
            if (alpha>PI/2.0){
                // cout<<"reverse"<<endl;
                tmp.x=nX*epaisseur_trait_parcours+it->x;
                tmp.y=nY*epaisseur_trait_parcours+it->y;
                listBackward.push_front(tmp);
                tmp.x=-nX*epaisseur_trait_parcours+it->x;
                tmp.y=-nY*epaisseur_trait_parcours+it->y;
                listForward.push_back(tmp);
            }else{//situation d'origine
                tmp.x=nX*epaisseur_trait_parcours+it->x;
                tmp.y=nY*epaisseur_trait_parcours+it->y;
                listForward.push_back(tmp);
                tmp.x=-nX*epaisseur_trait_parcours+it->x;
                tmp.y=-nY*epaisseur_trait_parcours+it->y;
                listBackward.push_front(tmp);
            }
        }
    }
    tmp.x=this->points.front().x;
    tmp.y=this->points.front().y;
    this->mesh.push_back(tmp);
    for (it=listForward.begin(); it!=listForward.end(); ++it){
        this->mesh.push_back(*it);
    }
    itPrev=--this->points.end();
    tmp.x=itPrev->x;
    tmp.y=itPrev->y;
    this->mesh.push_back(tmp);
    for (it=listBackward.begin(); it!=listBackward.end(); ++it){
        this->mesh.push_back(*it);
    }
}
void Intersection::print_mesh(){
    ofstream meshfile;
    meshfile.open("mesh.txt");
    list<Point>::iterator it;
    for (it=this->mesh.begin(); it!=this->mesh.end(); ++it){
        meshfile<<it->x<<" "<<it->y<<endl;
    }
    meshfile.close();
}


// PathFinding::PathFinding(list<Point>& pts){
//     this->intersec = new Intersection(pts);
// };
// PathFinding::~PathFinding(){
//     delete this->intersec;
// };