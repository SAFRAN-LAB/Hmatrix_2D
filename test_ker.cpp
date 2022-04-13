#include "HMAT2DTREE.hpp"

int main(int argc, char* argv[])
{
        Vec x,rhs;
        int nChebNodes    = atoi(argv[1]);  //  Number of Cheb nodes in single dimension
        double L          = atof(argv[2]);  // Half Length of domain centered at (0,0,..)
        Ker::Kchoice      = atoi(argv[3]);  // Type of Kernel refer KERNEL_FUNCTIONS
        int Hchoice       = atoi(argv[4]);  // Type of Hmatrix

        double xmin = -L,xmax = L;
        double ymin = -L,ymax = L;
        // creating N = nChebNodes*nChebNodes points [xmin,xmax]x[ymin,ymax]
        points_dt *charges;
        charges = new points_dt(nChebNodes,xmax,xmin,ymax,ymin);
        int N = charges->len();
        // Create Hierarchical quad tree
        Tree *At;
        At = new Tree(charges,  //  Cheb points
                      Hchoice,  // Hmatrix choice
                      xmin,xmax, // Domain
                      ymin,ymax,
                      100, 100, (double)pow(10,-10)); // GMRES(restart,Max_iter,tolerance)
        std::cout << "Tree Formed" << std::endl;
        // Initializes Tree
        At->Initialize_tree();
        std::cout << "H-matrix Initialized..." << std::endl;
        //At.set_rhs();
        Vec b;
        b = At->getRhs();
        std::cout << "RHS" << std::endl;
        Vec error;
        // GMRES based Iterative Solver
        error = At->solve(b);
        std::cout << "Err" << std::endl;
        At->print_tree();
        delete At;
        delete charges;
        return 0;
}
