#include "grid.h"


Grid::Grid(int Nx, int Ny, double* density, int VERBOSE){
    this->VERBOSE=VERBOSE;
    this->Nx=Nx;
    this->Ny=Ny;
    this->Dx=1.0/((double)(Nx-1));
    this->Dy=1.0/((double)(Ny-1));
    // cout<<">>> Grid parameters"<<endl<<"Dx="<<this->Dx<<"   Dy="<<this->Dy<<endl<<">>>>>>>>>>>>>"<<endl<<endl;
    this->u=density;
    this->r0=(double*)malloc(sizeof(double)*(Nx-1)*(Ny-1));
    for (int j=0; j<Ny-1; j++){
        this->r0[j]=0.0;
        for (int i=0; i<Nx-2; i++){
            this->r0[(i+1)*(Nx-1)+j]=this->r0[i*(Nx-1)+j]+this->u[i*(Nx-1)+j];
        }
    }
}
Grid::~Grid(){
    free(this->r0);
}
void Grid::print_r0(){
    cout<<endl<<"i j val"<<endl;
    for (int i=0; i<this->Nx-1; i++){
        for (int j=0; j<this->Ny-1; j++){
            cout<<i<<" "<<j<<" "<<this->r0[i*(Nx-1)+j]<<endl;
        }
    }
    cout<<endl;
}
void Grid::intersectSegment(Point* pt1, Point* pt2, list<PointPos>& intersectPt){
    list<ProjectionPoint> listX;
    list<ProjectionPoint> listY;
    int floorX1 = floor(pt1->x*(this->Nx-1));
    int floorX2 = floor(pt2->x*(this->Nx-1));
    int floorY1 = floor(pt1->y*(this->Ny-1));
    int floorY2 = floor(pt2->y*(this->Ny-1));
    double lenX = abs(pt1->x-pt2->x);
    double lenY = abs(pt1->y-pt2->y);
    ProjectionPoint tmp;
    if (floorX1<floorX2){
        for (int indX=floorX1+1; indX<=floorX2; indX++){
            tmp.val=indX*this->Dx;
            tmp.t=abs(tmp.val-pt1->x)/lenX;
            listX.push_back(tmp);
        }
    }else{
        for (int indX=floorX1; indX>=floorX2+1; indX--){
            tmp.val=indX*this->Dx;
            tmp.t=abs(tmp.val-pt1->x)/lenX;
            listX.push_back(tmp);
        }
    }
    if (floorY1<floorY2){
        for (int indY=floorY1+1; indY<=floorY2; indY++){
            tmp.val=indY*this->Dy;
            tmp.t=abs(tmp.val-pt1->y)/lenY;
            listY.push_back(tmp);
        }
    }else{
        for (int indY=floorY1; indY>=floorY2+1; indY--){
            tmp.val=indY*this->Dy;
            tmp.t=abs(tmp.val-pt1->y)/lenY;
            listY.push_back(tmp);
        }
    }

    double t;
    PointPos tmpPt;
    tmpPt.x=pt1->x;
    tmpPt.y=pt1->y;
    tmpPt.t=0.0;
    intersectPt.push_back(tmpPt);
    while (!listX.empty() || !listY.empty()){
        if (!listX.empty() && !listY.empty()){
            if (listX.front().t==listY.front().t){
                tmpPt.x=listX.front().val;
                tmpPt.y=listY.front().val;
                t=listX.front().t;
                tmpPt.t=t;
                // cout<<"corner : "<<tmpPt.x<<" "<<tmpPt.y<<endl;
                listX.pop_front();
                listY.pop_front();
            }else if (listX.front().t<listY.front().t){
                tmpPt.x=listX.front().val;
                t=listX.front().t;
                tmpPt.y=(1.0-t)*pt1->y+t*pt2->y;
                tmpPt.t=t;
                listX.pop_front();
            }else{
                tmpPt.y=listY.front().val;
                t=listY.front().t;
                tmpPt.x=(1.0-t)*pt1->x+t*pt2->x;
                tmpPt.t=t;
                listY.pop_front();
            }
        }else if(!listX.empty()){
            tmpPt.x=listX.front().val;
            t=listX.front().t;
            tmpPt.y=(1.0-t)*pt1->y+t*pt2->y;
            tmpPt.t=t;
            listX.pop_front();
        }else{
            tmpPt.y=listY.front().val;
            t=listY.front().t;
            tmpPt.x=(1.0-t)*pt1->x+t*pt2->x;
            tmpPt.t=t;
            listY.pop_front();
        }
        intersectPt.push_back(tmpPt);
    }
    tmpPt.x=pt2->x;
    tmpPt.y=pt2->y;
    tmpPt.t=1.0;
    intersectPt.push_back(tmpPt);
}

double Grid::computeCost(std::list<Point>& pts){
    double cost = 0.0;
    double localCost;
    double dir[2];
    double alpha;
    double normDir;
    int i,j;
    Point current, next;
    list<Point>::iterator itNode;
    list<Point>::iterator itNodeNext;
    itNodeNext=pts.begin();
    for (itNode=pts.begin(); itNode!=pts.end(); ++itNode){
        if (itNodeNext!=pts.end()){
            ++itNodeNext;
        }
        current = *itNode;
        if (itNode==--pts.end()){
            next = *pts.begin();
        }else{
            next = *itNodeNext;
        }
        // cout<<"new segment starting from ("<<current.x<<" "<<current.y<<") to ("<<next.x<<" "<<next.y<<")"<<endl;
        dir[0]=(next.x-current.x);
        dir[1]=(next.y-current.y);
        alpha = -dir[1]/sqrt(pow(dir[0],2.0)+pow(dir[1],2.0));
        normDir = sqrt(pow(next.x-current.x,2.0)+pow(next.y-current.y,2.0));
        // A VERIFIER !
        // normDir = 1.0/sqrt(pow(next.x-current.x,2.0)+pow(next.y-current.y,2.0));
        if (this->VERBOSE>=2)
            cout<<endl<<"Integrating over part ("<<current.x<<" "<<current.y<<") to ("<<next.x<<" "<<next.y<<")"<<endl<<"    alpha="<<alpha<<"\t normDir="<<normDir<<endl;

        list<PointPos> intersect;
        this->intersectSegment(&current,&next,intersect);
        list<PointPos>::iterator itPt;
        list<PointPos>::iterator itPtNext;
        itPtNext=intersect.begin();
        // ++itPtNext;
        for (itPt=intersect.begin(); itPt!=--intersect.end(); ++itPt){
            // if (itPtNext!=intersect.end()){
                ++itPtNext;
                // }
            i = floor(itPt->x*(this->Nx-1));
            j = floor(itPt->y*(this->Ny-1));
            if (this->VERBOSE>=2){
                cout<<"   Integrating over segment ("<<itPt->x<<" "<<itPt->y<<") to ("<<itPtNext->x<<" "<<itPtNext->y<<")"<<endl;
                cout<<"     In cell ("<<i<<" "<<j<<")";
                cout<<"\t u[(i,j)]="<<this->u[i*(this->Nx-1)+j];
                cout<<"\t r0[(i,j)]="<<this->r0[i*(this->Nx-1)+j]<<endl;
            }
            // cout<< "r0("<<i<<" "<<j<<")="<< this->r0[i*(this->Nx-1)+j] <<endl;
            localCost = this->Dx*this->u[i*(this->Nx-1)+j]*alpha*(itPtNext->t-itPt->t)*normDir*(itPt->x+itPtNext->x)/2.0
                        + this->Dx*this->r0[i*(this->Nx-1)+j]*alpha*(itPtNext->t-itPt->t)*normDir;
            // cout<<"localCost="<<localCost<<endl;
            cost += localCost;
        }
    }
    if (this->VERBOSE>=2){
        cout<<"cost="<<cost<<endl;
    }
    return cost;
}