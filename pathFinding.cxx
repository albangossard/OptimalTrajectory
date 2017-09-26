#include "grid.h"
#include "set.h"
#include <iostream>
#include <list>
using namespace std;

bool checkPositionPt(list<Point> &pts){
    list<Point>::iterator it;
    for (it=pts.begin(); it!=pts.end(); ++it){
        if (it->x<0.0 || it->x>1.0 || it->y<0.0 || it->y>1.0){
            return false;
        }
    }
    return true;
}

int main(){
    int Nx=5;
    int Ny=5;
    double* density=(double*)calloc((Nx-1)*(Ny-1), sizeof(double));
    for (int i=0; i<(Nx-1)/2; i++){
        for (int j=0; j<Ny-1; j++){
            density[i*(Nx-1)+j]=1.0;
        }
    }
    // density[1*(Nx-1)+1]=1.0;
    cout<<">>> Map"<<endl;
    for (int j=Ny-2; j>=0; j--){
        for (int i=0; i<Nx-1; i++){
            cout<<" "<<density[i*(Nx-1)+j];
        }
        cout<<endl;
    }
    cout<<"<<<<<<<<<<<<<<<"<<endl<<endl;
    
    list<Point> pts;
    Point tmp;
    double epsPath;
    int VERBOSE=0;


    /*** INTEGRATION TESTS ***/
    // Intersection * intersecTest;
    // Grid * MapTest;
    // double eps=3e-1*0+3e-2;
    // epsPath=0.02*0+0.01;
    // /*
    // --
    //  |
    // */
    // tmp.x=0.0+eps; tmp.y=1.0-eps; pts.push_back(tmp);
    // tmp.x=1.0-eps; tmp.y=1.0-eps; pts.push_back(tmp);
    // tmp.x=1.0-eps; tmp.y=0.0+eps; pts.push_back(tmp);
    // /*
    // --
    // |
    // */
    // tmp.x=0.0+eps; tmp.y=0.0+eps; pts.push_back(tmp);
    // tmp.x=0.0+eps; tmp.y=1.0-eps; pts.push_back(tmp);
    // tmp.x=1.0-eps; tmp.y=1.0-eps; pts.push_back(tmp);
    // // tmp.x=1.0-eps; tmp.y=0.0+eps; pts.push_back(tmp);

    // intersecTest = new Intersection(pts, VERBOSE);
    // intersecTest->generate_mesh(epsPath);
    // MapTest = new Grid(Nx,Ny,density,VERBOSE);
    // Map.print_r0();
    // double testCost=MapTest->computeCost(intersecTest->mesh);
    // cout<<"testCost="<<testCost<<endl;
    // delete intersecTest;
    // delete MapTest;
    // return 0;
    /*************************/
    
    // tmp.x=0.2; tmp.y=0.2; pts.push_back(tmp);
    // tmp.x=0.4; tmp.y=0.85; pts.push_back(tmp);
    // tmp.x=0.8; tmp.y=0.55; pts.push_back(tmp);
    // tmp.x=0.6; tmp.y=0.0; pts.push_back(tmp);

    // Diagonale
    // tmp.x=0.2; tmp.y=0.2; pts.push_back(tmp);
    // tmp.x=0.25; tmp.y=0.85; pts.push_back(tmp);
    // tmp.x=0.8; tmp.y=0.8; pts.push_back(tmp);

    // Horizontale
    tmp.x=0.2; tmp.y=0.5; pts.push_back(tmp);
    tmp.x=0.5*0+0.7; tmp.y=0.90*0+0.5; pts.push_back(tmp);
    tmp.x=0.8; tmp.y=0.5+0.6*0; pts.push_back(tmp);
    
    // double eps=1e-7;
    // tmp.x=0; tmp.y=0; pts.push_back(tmp);
    // tmp.x=0; tmp.y=1-eps; pts.push_back(tmp);
    // tmp.x=1-eps; tmp.y=1-eps; pts.push_back(tmp);
    // tmp.x=1-eps; tmp.y=0; pts.push_back(tmp);

    epsPath=0.02;

    double epsGrad = 1e-7*0+1e-3;
    int nIterMax = 10*0+2500;
    double tolDiffCost = 1e-10;
    double stepGrad=0.5*0+0.2*0+0.05+0.001*0;
    double dir[2];
    list<Point> ptsTmp;
    int numPt;
    double deltaX, deltaY;
    double theta1, theta2, theta;
    double differentialCost = tolDiffCost*2.0;
    list<Point>::iterator iterPt2;
    Intersection * intersec;
    Grid * Map;
    list<Point> ptsNew;
    list<Point>::iterator iterPtNew;

    cout<<" Parameters :"<<endl<<"epsGrad="<<epsGrad<<endl;

    ptsNew = pts;

    int nIter=0;
    while (nIter<nIterMax && abs(differentialCost)>tolDiffCost){
        nIter++;
        cout<<endl<<"### Iter "<<nIter<<" ###"<<endl;
        list<Point>::iterator iterPt;
        list<Point>::iterator iterPrev;
        list<Point>::iterator iterNext;
        numPt=1;
        for (iterPt=++pts.begin(); iterPt!=--pts.end(); ++iterPt){
            iterPrev=iterPt; iterPrev--;
            iterNext=iterPt; iterNext++;
            // cout<<"iterPrev "<<iterPrev->x<<" "<<iterPrev->y<<endl;
            // cout<<"iterPt "<<iterPt->x<<" "<<iterPt->y<<endl;
            // cout<<"iterNext "<<iterNext->x<<" "<<iterNext->y<<endl;
            differentialCost=0.0;
            // dir[0]=(iterPt->x-iterPrev->x)-(iterNext->x-iterPt->x);
            // dir[1]=(iterPt->y-iterPrev->y)-(iterNext->y-iterPt->y);
            deltaX=(iterPrev->x-iterPt->x);
            deltaY=(iterPrev->y-iterPt->y);
            if (deltaX>=0.0){
                theta1=atan(deltaY/deltaX);
            }else{
                theta1=PI+atan(deltaY/deltaX);
            }
            deltaX=(iterNext->x-iterPt->x);
            deltaY=(iterNext->y-iterPt->y);
            if (deltaX>=0.0){
                theta2=atan(deltaY/deltaX);
            }else{
                theta2=PI+atan(deltaY/deltaX);
            }
            theta=(theta1+theta2)/2.0;
            cout<<"theta1="<<theta1<<" theta2="<<theta2<<" theta="<<theta<<endl;
            dir[0]=cos(theta);
            dir[1]=sin(theta);
            cout<<"Direction ("<<dir[0]<<", "<<dir[1]<<")"<<endl;

            ptsTmp = pts;
            iterPt2=ptsTmp.begin();
            for (int i=0; i<numPt; i++){
                ++iterPt2;
            }
            iterPt2->x+=epsGrad*dir[0];
            iterPt2->y+=epsGrad*dir[1];
            intersec = new Intersection(ptsTmp, VERBOSE);
            intersec->generate_mesh(epsPath);
            if (!checkPositionPt(ptsTmp)){
                cout<<"position of the points out of range !"<<endl;
                exit(-1);
            }
            Map = new Grid(Nx,Ny,density,VERBOSE);
            differentialCost+=Map->computeCost(intersec->mesh);
            delete intersec;
            delete Map;
            
            ptsTmp = pts;
            iterPt2=ptsTmp.begin();
            for (int i=0; i<numPt; i++){
                ++iterPt2;
            }
            iterPt2->x-=epsGrad*dir[0];
            iterPt2->y-=epsGrad*dir[1];
            intersec = new Intersection(ptsTmp, VERBOSE);
            intersec->generate_mesh(epsPath);
            if (!checkPositionPt(ptsTmp)){
                cout<<"position of the points out of range !"<<endl;
                exit(-1);
            }
            Map = new Grid(Nx,Ny,density,VERBOSE);
            differentialCost-=Map->computeCost(intersec->mesh);
            delete intersec;
            delete Map;

            differentialCost=differentialCost/(2.0*epsGrad);

            cout<<"differentialCost="<<differentialCost<<endl;

            iterPtNew=ptsNew.begin();
            for (int i=0; i<numPt; i++){
                ++iterPtNew;
            }
            double computedStepX = differentialCost*stepGrad*dir[0];
            double computedStepY = differentialCost*stepGrad*dir[1];
            cout<<"computedStep="<<computedStepX<<" "<<computedStepY<<endl;
            iterPtNew->x-=computedStepX;
            iterPtNew->y-=computedStepY;

            cout<<"New point ("<<iterPtNew->x<<", "<<iterPtNew->y<<")"<<endl;

            numPt++;
        }
        pts=ptsNew;
    }

    return 0;
}

