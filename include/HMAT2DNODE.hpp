#ifndef HMAT2DNODE_HPP
#define HMAT2DNODE_HPP

#include "myHeaders.hpp"
#include "HODLR_Matrix.hpp"
#include "LowRank.hpp"
#include "CKERNEL.hpp"


const static double tol  = pow(10,-12);
size_t MAX_RANK_AFMM = 0;

class Node
{
int north,south,east,west;

public:
points_dt *charges;
Node *parent;
long int wload;
double fload;
std::vector<Node*> Child;
std::vector<size_t> crank;
std::vector<Node*> my_neighbour_addr;
std::vector<Node*> my_intr_list_addr;
std::vector<int> *idx = new std::vector<int>;
std::vector<int> *holdMyData = new std::vector<int>;

unsigned level_num,self_num,self_block;
double x_start,y_start,x_end,y_end;
double cx,cy;

int charge_start,charge_end;
int n_particles,n_neighbours,n_intraction;
bool isleaf = false;
bool isroot;
double eta_adm = 1.0;

std::vector<int> DLR,AFLR; // Diagonal Low Rank, Adjacent Full Rank Block


Mat P2P_self;
Mat *P2P,*P2M, *L2P;

Node(double xs,double xe,double ys,double ye,points_dt*& charges,double eta_) // Parent Node
{
        level_num = 0;
        self_num = 0;
        north = -1;
        south = -1;
        east = -1;
        west = -1;
        x_start = xs;
        x_end = xe;
        y_start = ys;
        y_end = ye;
        isroot = true;
        cx = 0.5*(xe + xs);
        cy = 0.5*(ye + ys);
        //parent = NULL;
        this->charges = charges;
        eta_adm =eta_;
        n_neighbours = 0;
        n_intraction = 0;
        wload = 0;
        std::cout << "Tol - " << tol << std::endl;
}

Node(unsigned level_num,unsigned self_num,Node*& obj,points_dt*& charges,double eta_)
{
        parent = obj;
        eta_adm =eta_;
        this->charges = charges;
        isroot = false;
        this->level_num = level_num;
        this->self_num = self_num;
        self_block = self_num%4;
        get_domain();
}

// Gets the domain limits for a Node
void get_domain();
// FMM Intraction list
void set_intraction_list(dtype_base rval);


// Initalizes the IFMM operators
void Initialize_node(bool isPrecond = false);
void mat_vec(const Vec& x,Vec& b);
// Utilities
void print_self();
void print_self(string fname,bool print_flag);

void myload();
double myload_full();
void particle_to_child();

void mark_leaf();

~Node(){
        delete idx;
}

};

void Node::mark_leaf()
{
        this->isleaf = true;
        for(int i=0; i<holdMyData->size(); i++)
                this->idx->push_back(holdMyData->at(i));
        n_particles = this->idx->size();
        delete holdMyData;
}
/*
   Sub Domain  Numbered as follows
   ____________________
 |         |        |
 |  3      |    2   |
 |_________|________|
 |         |        |
 |  0      |    1   |
 |_________|________|

 */
void Node::get_domain()
{
        if(self_block == 0)
        {
                x_start = parent->x_start;
                x_end = parent->cx;
                y_start = parent->y_start;
                y_end = parent->cy;
        }
        if(self_block == 1)
        {
                x_start =  parent->cx;
                x_end = parent->x_end;
                y_start = parent->y_start;
                y_end = parent->cy;
        }
        if(self_block == 2)
        {
                x_start = parent->cx;
                x_end = parent->x_end;
                y_start = parent->cy;
                y_end = parent->y_end;
        }
        if(self_block == 3)
        {
                x_start = parent->x_start;
                x_end = parent->cx;
                y_start = parent->cy;
                y_end = parent->y_end;
        }
        cx = 0.5*(x_end + x_start);
        cy = 0.5*(y_end + y_start);
}

void Node::set_intraction_list(dtype_base Lmax)
{
        dtype_base rval = Lmax/pow(2,level_num-1);
        if(!isroot)
        {
                // Push brothers in neighbors list
                for(int k=0; k<4; k++)
                {
                        Node *KD;
                        KD = parent->Child[k];
                        double dist = Calc_dist(cx,cy,KD->cx,KD->cy);
                        if(self_num != KD->self_num)
                        {
                                if(dist <= eta_adm * rval)
                                {
                                        AFLR.push_back(KD->self_num);
                                        my_neighbour_addr.push_back(KD);
                                }
                                else
                                {
                                        DLR.push_back(KD->self_num);
                                        my_intr_list_addr.push_back(KD);
                                }
                        }
                }
                // Search in Parents Neighbors for Intraction list
                for(size_t i=0; i< parent->AFLR.size(); i++)
                {
                        Node *KD;
                        for(int k=0; k<4; k++)
                        {
                                KD = parent->my_neighbour_addr[i]->Child[k];
                                double dist = Calc_dist(cx,cy,KD->cx,KD->cy);
                                if(dist <= eta_adm * rval)
                                {
                                        AFLR.push_back(KD->self_num);
                                        my_neighbour_addr.push_back(KD);
                                }
                                else
                                {
                                        DLR.push_back(KD->self_num);
                                        my_intr_list_addr.push_back(KD);
                                }
                        }
                }
        }
        n_neighbours = AFLR.size();  // AFLR - Neighbours list
        n_intraction = DLR.size();  // DLR - Interaction list
}


void Node::Initialize_node(bool isPrecond)
{
        if(n_particles != 0)
        {
                P2P = new Mat[n_neighbours];
                L2P = new Mat[n_intraction];
                P2M = new Mat[n_intraction];
                if(isleaf)
                {
                        Kernel* K_self  = new Kernel(charges,idx,n_particles);
                        P2P_self = Mat::Zero(n_particles,n_particles);
                        P2P_self = K_self->getMatrix(0,
                                                     0,
                                                     n_particles,
                                                     n_particles);
                        delete K_self;
                        for(int i=0; i<n_neighbours; i++)
                        {
                                if(my_neighbour_addr[i]->n_particles != 0)
                                {
                                        Kernel* K_neighbors  = new Kernel(charges,idx,my_neighbour_addr[i]->idx,n_particles);
                                        P2P[i] = Mat::Zero(n_particles,my_neighbour_addr[i]->n_particles);
                                        P2P[i] = K_neighbors->getMatrix(0,0,
                                                                        n_particles,
                                                                        my_neighbour_addr[i]->n_particles);
                                        delete K_neighbors;
                                }
                        }

                }
                for(int i=0; i<n_intraction; i++)
                {
                        if(my_intr_list_addr[i]->n_particles != 0)
                        {
                                Kernel* K_intraction   = new Kernel(charges,idx,my_intr_list_addr[i]->idx,n_particles);
                                LowRank* lr = new LowRank(K_intraction, "rookPivoting");
                                lr->getFactorization(L2P[i], P2M[i], tol,0,
                                                     0,
                                                     n_particles,
                                                     my_intr_list_addr[i]->n_particles);
                                crank.push_back(L2P[i].cols());
                                if(MAX_RANK_AFMM < L2P[i].cols())
                                        MAX_RANK_AFMM = L2P[i].cols();

                                //std::cout << std::endl<< "L2P_"<<self_num<<"_"<<my_intr_list_addr[i]->self_num  << " : "<< L2P[i].cols()<< std::endl;
                                delete K_intraction;
                        }

                }
        }
}

// Performs individual Mat-Vec
void Node::mat_vec(const Vec& x,Vec& b)
{
        if(n_particles!=0)
        {
                Vec node_charge = Vec::Zero(n_particles);
                if(isleaf)
                {
                        node_charge += P2P_self * x.segment(charge_start, n_particles);

                        for(int i=0; i<n_neighbours; i++)
                        {
                                if(my_neighbour_addr[i]->n_particles!=0)
                                {
                                        node_charge +=
                                                P2P[i] * x.segment(my_neighbour_addr[i]->charge_start, my_neighbour_addr[i]->n_particles);

                                }

                        }
                }
                for(int i=0; i<n_intraction; i++)
                {
                        if(my_intr_list_addr[i]->n_particles!=0)
                        {
                                node_charge +=
                                        (L2P[i] * (P2M[i].transpose() * x.segment(my_intr_list_addr[i]->charge_start, my_intr_list_addr[i]->n_particles)));

                        }
                }
          #pragma omp critical
                {
                        b.segment(charge_start, n_particles) += node_charge;
                }
        }
}

void Node::print_self()
{
        std::cout<<"Level :"<<level_num<<" Self ID : " <<self_num<<std::endl;
        std::cout << "Neighbour incl self : ";
        for(int i=0; i<n_neighbours; i++)
                std::cout << AFLR[i] << " ";
        std::cout<<std::endl;
        std::cout << "DALR : ";
        for(int i=0; i<n_intraction; i++)
                std::cout << DLR[i] << " ";
        std::cout<<std::endl;
        std::cout << std::endl;
        for(int i=0; i<crank.size(); i++)
                std::cout << crank[i] << " ";
        std::cout << std::endl;
        std::cout<<"Charge Index : ["<<charge_start<<","<<charge_end<<"]"<<std::endl;
        std::cout<<"X : ["<<x_start<<","<<x_end<<"]"<<std::endl;
        std::cout<<"Y : ["<<y_start<<","<<y_end<<"]"<<std::endl;
        std::cout << "N Particles : " << n_particles <<std::endl;
}

void Node::print_self(string fname,bool print_flag)
{
        std::ofstream fout;
        fout.open(fname,std::ios::app);
        fout << "Level : " << level_num << " Self ID : " << self_num << std::endl;
        fout << "Neighbours : ";
        for(int i=0; i<n_neighbours; i++)
                fout << AFLR[i] << " ";
        fout << std::endl;
        fout << "DALR : ";
        for(int i=0; i<n_intraction; i++)
                fout << DLR[i] << " ";
        fout << std::endl;
        for(int i=0; i<crank.size(); i++)
                fout << crank[i] << " ";
        fout<<std::endl;
        fout<<"X : ["<<x_start<<","<<x_end<<"]"<<std::endl;
        fout<<"Y : ["<<y_start<<","<<y_end<<"]"<<std::endl;
        fout << "N Particles : " << n_particles <<std::endl;
        fout.close();
}


double Node::myload_full()
{
        fload = 0;
        if(n_particles != 0)
        {
                if(isleaf)
                {
                        fload += n_particles * n_particles;
                        for(int i=0; i<n_neighbours; i++)
                                fload += my_neighbour_addr[i]->n_particles * this->n_particles;
                }
                for(int i=0; i<n_intraction; i++)
                {
                        fload += my_intr_list_addr[i]->n_particles * P2M[i].cols();
                        fload += L2P[i].cols() * this->n_particles;
                }
        }
        return fload;
}

void Node::particle_to_child()
{
        for(size_t iter=0; iter<holdMyData->size(); iter++)
        {
                double hx = 0.0,hy = 0.0;
                int id = 0;
                hx = charges->gridPoints[holdMyData->at(iter)].x - this->cx;
                hy = charges->gridPoints[holdMyData->at(iter)].y - this->cy;
                if(hx < 0.0 && hy <= 0.0)
                        id = 0;
                if(hx >= 0.0 && hy < 0.0)
                        id = 1;
                if(hx > 0.0 && hy >= 0.0)
                        id = 2;
                if(hx <= 0.0 && hy > 0.0)
                        id = 3;
                Child[id]->holdMyData->push_back(holdMyData->at(iter));
        }
        delete holdMyData;
}

/*cout<<"North ..."<<north<<endl;
   cout<<"South ..."<<south<<endl;
   cout<<"East  ..."<<east<<endl;
   cout<<"West  ..."<<west<<endl;*/


#endif
