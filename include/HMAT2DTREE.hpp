#ifndef HMAT2DTREE_HPP
#define HMAT2DTREE_HPP

#include "GMRES.hpp"
#include "myHeaders.hpp"
#include "HMAT2DNODE.hpp"

const static size_t max_particles_leaf  = 1000;

class Tree : public iterSolver
{
int level,n,N,fam_count,llnum;        // n - number in xdir, N - Number of charges
double etaAdm;
int Hchoice = 1;
Vec testx,testrhs;
Vec solx,solrhs; // For Backward Error
HODLR *LSHodlr;
Kernel *KHodlr;

double wload_eval = 0;
public:
std::vector<int> *idx = new std::vector<int>;
std::vector<std::vector<Node *> > obj_arr;
std::vector<Node *>  obj_arr_lin;
Node *root;
points_dt *charges;

// ------ Time Profile -----
double start, end;
double startp, endp;
double start_iter, end_iter;
double TREE_FORMATION_TIME,TREE_INIT_TIME,MAT_VEC_TIME_PER_ITER,SOL_TIME,MAT_VEC_TIME;
double RESIDUAL;
dtype_base hlen;
double BE_x, BE_rhs;


Tree(points_dt*&charges,int hchoice,double xs,double xe,double ys,double ye,int a,int b,double c) : iterSolver(a,b,c)
{
        hlen = xe;
        this->charges = charges;

        N = charges->len();
        // HODLR
        if (hchoice == 0)
        {
                etaAdm = 0.5;
                this->Hchoice = 0;
                std::cout << "HODLR" << std::endl;
        }
        // HODLR2D
        else if(hchoice == 1)
        {
                etaAdm = 1.0;
                this->Hchoice = 1;
                std::cout << "HODLR2D" << std::endl;
        }
        // Hmatrix with standard admissibility
        else
        {
                etaAdm = 1.5;
                this->Hchoice = 2;
                std::cout << "HMATRIX" << std::endl;
        }
        start = omp_get_wtime();
        create_tree(xs, xe, ys, ye);
        end = omp_get_wtime();
        TREE_FORMATION_TIME = end -start;
        KHodlr  = new Kernel(charges,idx,idx,N);
        KHodlr->set_test();
}

Vec matvec(Vec& x);
void set_rhs();

Vec getX()
{
        return testx;
}
Vec getRhs()
{
        return testrhs;
}
void create_tree(double xs,double xe,double ys,double ye);
int recursive_find_domain(Node*& node, dtype_base x,dtype_base y);

// Gets the indices of charges inside each Node at Leaf level
void get_points_inside_box();

// Updates the indices for nodes in the tree from leaf to root
void update_tree();
void update_child_info();

void Initialize_tree();
Vec solve(Vec& x);
void print_tree();
void wload_eval_fun();
// Added at V8
void transfer_particles();
bool do_partition();
void get_child_info();


~Tree()
{
        for(int i = 0; i < level; i++)
        {
      #pragma omp parallel for
                for(size_t j = 0; j < (size_t)pow(4,i); j++)
                {
                        delete obj_arr[i][j];
                }
        }

        //delete charges;
        delete idx;
        std::cout << "Deleted Data" << std::endl;
}
};

void Tree::create_tree(double xs,double xe,double ys,double ye)
{
        llnum = 0;
        level = llnum + 1;

        fam_count = int((pow(4,level)-4)/3) + 1;
        std::cout << "Started" << std::endl;
        obj_arr.resize(level);

        // Set Root of the Tree
        Node* A = new Node(xs,xe,ys,ye,charges,etaAdm);
        root = A;

        for(int iter=0; iter<N; iter++)
                A->holdMyData->push_back(iter);

        obj_arr[0].push_back(A);

        while(do_partition())
        {
                obj_arr.resize(level+1);
                // Create Leaf Nodes
                for(int j=0; j<pow(4,level); j++)
                {
                        Node* temp= new Node(level,j,obj_arr[level-1][j/4],charges,etaAdm);
                        obj_arr[level].push_back(temp);
                }
                // Child Nodes linked to Parent
                get_child_info();
                // Update Child With their parents particles
                transfer_particles();
                for(size_t j=0; j<obj_arr[level].size(); j++)
                        obj_arr[level][j]->set_intraction_list(hlen);
                llnum++;
                level++;
        }
        std::cout << "Tree Created and points assigned to Leaf" << std::endl;
        get_points_inside_box();
        std::cout << "Points distributed till Root" << std::endl;
        update_tree();
        std::cout << std::endl <<"Charges Created : " << charges->len() << std::endl;
        for(int k = 0; k<level; k++)
                for(size_t i=0; i<obj_arr[k].size(); i++)
                        obj_arr_lin.push_back(obj_arr[k][i]);
        std::cout << "RHS set" << std::endl;
}

bool Tree::do_partition()
{
        // Allows to perform partition and if condition fails it marks the nodes a leaf
        size_t max_particles = obj_arr[llnum][0]->holdMyData->size();
        for(size_t j=1; j<obj_arr[llnum].size(); j++)
        {
                if(obj_arr[llnum][j]->holdMyData->size() > max_particles)
                        max_particles = obj_arr[llnum][j]->holdMyData->size();
        }
        if(max_particles_leaf < max_particles)
                return true;
        else
        {
                for(size_t j=0; j<obj_arr[llnum].size(); j++)
                        obj_arr[llnum][j]->mark_leaf();
                return false;
        }
}

void Tree::get_child_info()
{
// Assumes that the domain is subdivided and Nodes are formed and done before changing the Leaf
        for(size_t j = 0; j < obj_arr[llnum].size(); j++)
        {
                for (int c=0; c<4; c++)
                {
                        int tmp = 4 * obj_arr[llnum][j]->self_num + c;
                        obj_arr[llnum][j]->Child.push_back(obj_arr[level][tmp]);
                }
        }
        std::cout <<"Child Updated for Level " << llnum <<std::endl;
}

// Transfer the information from root to leaf
void Tree::transfer_particles()
{
        for(size_t j = 0; j < obj_arr[llnum].size(); j++)
                obj_arr[llnum][j]->particle_to_child();
}


// Provides the charges for each node
void Tree::get_points_inside_box()
{
        for(size_t i = 0; i<obj_arr[llnum].size(); i++)
        {
                obj_arr[llnum][i]->charge_start = this->idx->size();
                for(int j=0; j<obj_arr[llnum][i]->n_particles; j++)
                        this->idx->push_back(obj_arr[llnum][i]->idx->at(j));
                obj_arr[llnum][i]->charge_end = this->idx->size()-1;
        }
}

// Updates the domain of each node in the tree from leaf to root
void Tree::update_tree()
{
        for(int i = llnum; i > 0; --i)
        {
                for(size_t j = 0; j < obj_arr[i].size(); j++)
                {
                        if(obj_arr[i][j]->self_block == 0)
                                obj_arr[i][j]->parent->charge_start = obj_arr[i][j]->charge_start;
                        if(obj_arr[i][j]->self_block == 3)
                                obj_arr[i][j]->parent->charge_end = obj_arr[i][j]->charge_end;
                        for (std::vector<int>::iterator ptr = obj_arr[i][j]->AFLR.begin(); ptr != obj_arr[i][j]->AFLR.end(); ptr++)
                                obj_arr[i][j]->my_neighbour_addr.push_back(obj_arr[i][*ptr]);
                        for (std::vector<int>::iterator ptr = obj_arr[i][j]->DLR.begin(); ptr != obj_arr[i][j]->DLR.end(); ptr++)
                                obj_arr[i][j]->my_intr_list_addr.push_back(obj_arr[i][*ptr]);
                        if(!obj_arr[i][j]->isleaf)
                        {
                                for(int c=0; c<4; c++)
                                {
                                        for(size_t k=0; k<obj_arr[i][j]->Child[c]->idx->size(); k++)
                                                obj_arr[i][j]->idx->push_back(obj_arr[i][j]->Child[c]->idx->at(k));
                                }
                                obj_arr[i][j]->n_particles = obj_arr[i][j]->idx->size();
                        }
                }
        }
        obj_arr[0][0]->n_particles = this->idx->size();
        if(obj_arr[0][0]->Child.size()!=0)
                for(int c=0; c<4; c++)
                {
                        for(size_t k=0; k<obj_arr[0][0]->Child[c]->idx->size(); k++)
                                obj_arr[0][0]->idx->push_back(obj_arr[0][0]->Child[c]->idx->at(k));
                }
}

// Initialization of Nodes and formation of tree
void Tree::Initialize_tree()
{

        std::cout << "Initialising : "  <<std::endl;
        start = omp_get_wtime();
        /*  for(int k = 0; k<level; k++)
                  for(size_t i=0; i<obj_arr[k].size(); i++)
                  {
                          std::cout << "Initialising : " << obj_arr[k][i]->self_num << std::endl;
                          obj_arr[k][i]->Initialize_node();
                  }*/
    #pragma omp parallel for
        for(size_t i=0; i<obj_arr_lin.size(); i++)
                obj_arr_lin[i]->Initialize_node();
        end = omp_get_wtime();
        TREE_INIT_TIME = end -start;
        std::cout << std::endl;
        std::cout << "Initialized .... "<< std::endl;
        std::cout << std::endl;
        set_rhs();
}

//Mat vec, Also an explicit for iterSolver
Vec Tree::matvec(Vec& x)
{
        Vec b = Vec::Zero(x.size());
        /*for(int k = 0; k<level; k++)
                for(size_t i=0; i<obj_arr[k].size(); i++)
                        obj_arr[k][i]->mat_vec(x,b);*/
    #pragma omp parallel for
        for(size_t i=0; i<obj_arr_lin.size(); i++)
                obj_arr_lin[i]->mat_vec(x,b);
        return b;
}

//Public routine called form driver function that provides the solution
Vec Tree::solve(Vec& rhs)
{
        Vec x = Vec::Zero(rhs.size());
        iterSolver* Solver;
        Solver = this;
        x = rhs;
        string fname;
        fname = "gmres_out_" + std::to_string((int)Ker::Kchoice) + "_" + std::to_string(N) + "_" + std::to_string(level) + ".txt";
        Solver->set_output_file(fname);
        start = omp_get_wtime();
        int flag = this->GMRES(x,rhs);
        end = omp_get_wtime();
        SOL_TIME = end -start;
        RESIDUAL = getResidual();
        int tmp = getMaxIterations();
        if(tmp !=0)
                MAT_VEC_TIME_PER_ITER = SOL_TIME/tmp;

        this->solx = Vec::Zero(rhs.size());
        this->solx = x;
        BE_x = Vec_ops::relative_error(testx,solx);
        return x;
}

void Tree::set_rhs()
{
        this->testrhs = Vec::Zero(N);
        this->testx = Vec::Zero(N);
        this->testrhs = KHodlr->rhs;
        this->testx = KHodlr->x;
}

// This Evaluates number of doubles in physical memory
void Tree::wload_eval_fun()
{
        wload_eval = 0;
        for(int i=1; i<level; i++)
                for(size_t j=0; j<obj_arr[i].size(); j++)
                {
                        wload_eval += obj_arr[i][j]->myload_full();
                }
}

// Utility to check the correctness of various node in Tree
void Tree::print_tree()
{
        string fname;
        fname = "hmat_" + std::to_string((int)Ker::Kchoice) + "_" + std::to_string(N) + "_" + std::to_string(level) + ".txt";
        std::ofstream fout;
        fout.open(fname);        //,std::ios::app);
        fout << (std::string)getTimeStamp() << std::endl;
        fout << "Residual : " << RESIDUAL << std::endl;
        fout << "Backward Error X : " << BE_x << std::endl;
        fout << "Tree Formation Time : " <<TREE_FORMATION_TIME << "s" << std::endl;        // .count()
        fout << "Tree Initialization Time : " <<TREE_INIT_TIME << "s" << std::endl;        // .count()
        fout << "Mat-Vec Per Iteration Time : " <<MAT_VEC_TIME_PER_ITER << "s" << std::endl;        // .count()
        fout << "Time to Solution : " << SOL_TIME << "s" << std::endl;        // .count()
        wload_eval_fun();
        fout << "Number of FLOP : " << wload_eval << std::endl;
        double N2 = double(N) * double(N);
        fout << "Compression : " << wload_eval/N2 << std::endl;
        fout << "MAX Rank across Tree : " << MAX_RANK_AFMM <<std::endl;
        fout << std::endl << "==============================================================================" << std::endl;
        fout.close();
        for(int i=0; i<level; i++)
                for(size_t j=0; j<obj_arr[i].size(); j++)
                {
                        obj_arr[i][j]->print_self();
                        obj_arr[i][j]->print_self(fname, true);
                        fout.open(fname,std::ios::app);
                        if(obj_arr[i][j]->isleaf)
                        {
                                std::cout << std::endl << "I am Leaf" << std::endl;
                                fout << "I am Leaf" << std::endl;
                        }
                        else
                        {
                                std::cout << std::endl <<"I am Not a Leaf" << std::endl;
                                fout <<"I am Not a Leaf" << std::endl;
                        }
                        fout << std::endl << "==============================================================================" << std::endl;
                        fout.close();
                }
}
#endif
