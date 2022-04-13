#ifndef CKERNEL_HPP
#define CKERNEL_HPP

struct pts2D {
        double x,y;
};

const double PI = 3.1415926535897932384;
double kappa = 1.0;
int Qchoice = 0;

#include "HODLR_Matrix.hpp"
#include "HODLR.hpp"
#include <vector>
#include "KERNEL_FUNCTIONS.hpp"

class FMM2DTree
{
public:
std::vector<pts2D> gridPoints;
FMM2DTree(std::vector<pts2D> data)
{
        this->gridPoints = data;
}
dtype_base getMatrixEntry(int i,int j)
{
        dtype_base output = dtype_base(0.0);
        output = Ker::KernelFunc(gridPoints[i].x,gridPoints[i].y,gridPoints[j].x,gridPoints[j].y);
        return output;
}
~FMM2DTree(){
}
};

class points_dt
{
int N;

VectorXd Xdir,Ydir;
public:
FMM2DTree* F;
std::vector<pts2D> gridPoints;     // location of particles in the domain
points_dt(int nChebNodes, dtype_base xmax, dtype_base xmin, dtype_base ymax,dtype_base ymin)
{
        std::cout << "Cheb Points"<<std::endl;
        this->Xdir = cheb_nodes(xmin,xmax,nChebNodes);
        this->Ydir = cheb_nodes(ymin,ymax,nChebNodes);
        this->N = nChebNodes * nChebNodes;
        for(size_t i=0; i<nChebNodes; i++)
                for(size_t j=0; j<nChebNodes; j++)
                {
                        pts2D temp;
                        temp.y = Ydir[i];
                        temp.x = Xdir[j];
                        gridPoints.push_back(temp);
                }
        this->F = new FMM2DTree(gridPoints);
        Ker::Kii = N;
}
Eigen::VectorXd cheb_nodes(double a,double b,int n)
{
        Eigen::VectorXd X(n);
        double l,l1,param;
        l = 0.5*(a+b);
        l1 = 0.5*(b-a);
        for(int k=0; k<n; k++)
        {
                param = (double)(k+0.5)/n;
                X(k) = l - l1*cos(param*M_PI);
        }
        return X;
}
int len()
{
        return N;
}
~points_dt(){
        delete F;
}
};

class Kernel : public HODLR_Matrix
{
points_dt* FKerneldat;
int Ndt;
std::vector<int>* src_idx;
std::vector<int>* dst_idx;
double start, end;

public:
double Matrix_TIME,MAT_VEC_TIME;
Vec x,rhs;
Kernel(points_dt*& data,std::vector<int>*& src,std::vector<int>*& dst,int N) : HODLR_Matrix(N)
{
        FKerneldat = data;
        src_idx = src;
        dst_idx = dst;
        Ndt = FKerneldat->len();
}

Kernel(points_dt*& data,std::vector<int>*& src,int N) : HODLR_Matrix(N)
{
        FKerneldat = data;
        src_idx = src;
        dst_idx = src;
        Ndt = FKerneldat->len();
}

dtype getMatrixEntry(int i, int j)
{
        dtype output;
        output = FKerneldat->F->getMatrixEntry(src_idx->at(i), dst_idx->at(j));
        return output;
}

void set_test(bool doComputeAll=true)
{
        x = Vec::Zero(N);
        rhs = Vec::Zero(N);
        rhs  = getCol(0, 0, N);
        x(0) = 1;
}

~Kernel()
{
}
};
#endif

/*
   Pa Pa Ga Sa Dha Ni Pa
   Ma Pa Ga Ma Ri Sa Sa
   Pa Ma Ri Sa Dha Ni Pa
   Ma Pa Ga Ma Pa

   Pa Sa Ni Ma PA Ni Sa
   Ga GA Ma Ri Sa
   Sa Ni Sa Pa  Ni Ma Pa
   Ma Ri Sa Ni Sa Ri Ma
 */
