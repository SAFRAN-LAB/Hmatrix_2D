#ifndef GMRES_HPP
#define GMRES_HPP

//*****************************************************************
// Iterative template routine -- GMRES
//
// Remodified from : https://github.com/amiraa127/Sparse_MultiFrontal/blob/master/IML/include/gmres.h
// GMRES solves the unsymmetric linear system Ax = b using the
// Generalized Minimum Residual method
//
// GMRES follows the algorithm described on p. 20 of the
// SIAM Templates book.
//
// The return value indicates convergence within max_iter (input)
// iterations (0), or no convergence within max_iter iterations (1).
//
// Upon successful return, output arguments have the following values:
//
//        x  --  approximate solution to Ax = b
// max_iter  --  the number of iterations performed before the
//               tolerance was reached
//      tol  --  the residual after the final iteration
// Remodification: (Functionality remains the same)
//      Preconditioned Mat-Vec can be user defined
//      Modified for complex
//      myHeader has macro for different data types
//*****************************************************************

#include "myHeaders.hpp"

class iterSolver
{
int restart,max_iter;
string fname;
double tol;
dtype_base resid;
bool ptofile = false;
std::ofstream fout;
void GeneratePlaneRotation(dtype& dx, dtype& dy,dtype& cs, dtype& sn);
void ApplyPlaneRotation(dtype& dx, dtype& dy, dtype& cs, dtype& sn);
void Update(Vec& x, int k, Mat &h, Vec& s, Vec*& v);
public:
virtual Vec matvec(Vec& x)
{
        x.setZero(x.size());
        std::cout <<"Sanity Check .. Define Mat-Vec"<<std::endl;
        return x;
}
iterSolver(int restart,int max_iter,double tol);
iterSolver(int max_iter,double tol);
int GMRES(Vec &x, const Vec &b);
void set_output_file(string fname);
double getResidual();
int getMaxIterations();
~iterSolver() {
}
};

// Constructor with restart, MAX Iterations and tolerance
iterSolver::iterSolver(int restart,int max_iter,double tol)
{
        this->restart = restart;
        this->max_iter = max_iter;
        this->tol =  tol;
        resid = 0.0;
}

// Constructor with MAX Iterations and tolerance Restart set to MAX iterations
iterSolver::iterSolver(int max_iter,double tol)
{
        this->restart = max_iter;
        this->max_iter = max_iter;
        this->tol =  tol;
        resid = 0.0;
}


void iterSolver::set_output_file(string fname)
{
        this->fname = fname;
        ptofile = true;
}

void iterSolver::GeneratePlaneRotation(dtype& dx, dtype& dy,dtype& cs, dtype& sn)
{
        if (dy == 0.0)
        {
                cs = 1.0;
                sn = 0.0;
        }
        else if (abs(dy) > abs(dx))
        {
                dtype temp = dx / dy;
                sn = 1.0 / sqrt( 1.0 + temp*temp );
                cs = temp * sn;
        }
        else
        {
                dtype temp = dy / dx;
                cs = 1.0 / sqrt( 1.0 + temp*temp );
                sn = temp * cs;
        }
}


void iterSolver::ApplyPlaneRotation(dtype& dx, dtype& dy, dtype& cs, dtype& sn)
{
        dtype temp  =  cs * dx + sn * dy;
        dy = (-conj_(sn)) * dx + cs * dy;
        dx = temp;
}

void iterSolver::Update(Vec& x, int k, Mat &h, Vec& s, Vec*& v)
{
        Vec y = s;
        // Backsolve:
        for (int i = k; i >= 0; i--)
        {
                y(i) /= h(i,i);
                for (int j = i - 1; j >= 0; j--)
                        y(j) -= h(j,i) * y(i);
        }

        for (int j = 0; j <= k; j++)
                x += v[j] * y(j);
}

dtype_base iterSolver::getResidual()
{
        return resid;
}

int iterSolver::getMaxIterations()
{
        return max_iter;
}

int iterSolver::GMRES(Vec &x, const Vec &b)
{
        int counter = 0;
        double start,end;
        double MatVecTIME = 0.0;
        std::cout << "Restarted GMRES Parameters " << std::endl;
        std::cout << "Maximum Iterations : " << max_iter << std::endl;
        std::cout << "Tolerance : " << tol << std::endl;
        std::cout << "Restart Iteration at " << restart << std::endl;
        if(ptofile)
        {
                fout.open(fname);//,std::ios::app);
                fout << (std::string)getTimeStamp() << std::endl;
                fout << "Restarted GMRES Parameters " << std::endl;
                fout << "Maximum Iterations : " << max_iter << std::endl;
                fout << "Tolerance :" << tol << std::endl;
                fout << "Restart Iteration at " << restart << std::endl;
        }

        int m = restart;
        // It is expected that you apply the precondtioner to the rhs and supply
        // The Upper Hessenberg Matrix
        Mat H = Mat::Zero(restart+1,restart);

        int i, j = 1, k;
        // Vectors for Rotation

        Vec s = Vec::Zero(restart+1),
            cs = Vec::Zero(restart+1),
            sn = Vec::Zero(restart+1), w;

        dtype_base normb = b.norm();

        Vec r = b - matvec(x); //Explict method to perform the MatVec

        dtype_base beta = r.norm();

        if (normb == 0.0)
                normb = 1;
        resid = beta/normb;
        if(resid <= tol)
        {
                tol = resid;
                max_iter = 0;
                std::cout << "Exited before any Progress" <<std::endl;
                if(ptofile)
                {
                        fout << "Exited before any Progress" <<std::endl;
                        fout << "Residual : " <<resid  << std::endl;
                        fout.close();
                }
                return 0;
        }

        Vec *v = new Vec[restart+1];

        while (j <= max_iter)
        {

                v[0] = r * (1.0 / beta);
                s = Vec::Zero(restart+1);
                s(0) = beta;

                for (i = 0; i < m && j <= max_iter; i++, j++)
                {
                        std::cout << "Iter["<<j<<"] : "<<std::setprecision(4)<<resid<<std::endl;
                        if(ptofile)
                                fout << "Iter["<<j<<"] : "<<std::setprecision(8)<<resid<<std::endl;

                        w = matvec(v[i]);
                        MatVecTIME += end - start;

                        for (k = 0; k <= i; k++)
                        {
                                H(k, i) = w.dot(v[k]);
                                w -= H(k, i) * v[k];
                        }
                        H(i+1, i) = w.norm();
                        v[i+1] = (1.0 / H(i+1, i)) * w;   // ??? w / H(i+1, i)

                        for (k = 0; k < i; k++)
                                ApplyPlaneRotation(H(k,i), H(k+1,i), cs(k), sn(k));

                        GeneratePlaneRotation(H(i,i), H(i+1,i), cs(i), sn(i));
                        ApplyPlaneRotation(H(i,i), H(i+1,i), cs(i), sn(i));
                        ApplyPlaneRotation(s(i), s(i+1), cs(i), sn(i));
                        resid = abs_(s(i+1)) / normb;
                        if (resid < tol)
                        {
                                Update(x, i, H, s, v);
                                tol = resid;
                                max_iter = j;
                                std::cout << "Reached Solution before Max_Iterations" <<std::endl;
                                if(ptofile)
                                {
                                        fout << "Reached Solution before Max_Iterations" <<std::endl;
                                        fout << "Residual : " <<resid  << std::endl;
                                        fout.close();
                                }
                                delete [] v;
                                return 0;
                        }
                }
                std::cout<<"Restarted"<<std::endl;
                Update(x, m - 1, H, s, v);

                r = b - matvec(x);

                beta = r.norm();
                resid = beta / normb;
                if(resid < tol)
                {
                        max_iter = j;
                        std::cout << "Reached Solution before Max_Iterations" <<std::endl;
                        if(ptofile)
                        {
                                fout << "Reached Solution before Max_Iterations" <<std::endl;
                                fout << "Residual : " <<resid  << std::endl;
                                fout.close();
                        }
                        delete [] v;
                        return 0;
                }
                //std::cout << "Iter["<<j<<"] : "<<std::setprecision(4)<<".."<<resid<<std::endl;
        }
        if(ptofile)
                fout.close();
        delete [] v;
        return 1;
}


#endif

// The Verion of Plane rotation from the source
/*
   void iterSolver::GeneratePlaneRotation(dtype& dx, dtype& dy,dtype& cs, dtype& sn)
   {
        if (dy == 0.0)
        {
                cs = 1.0;
                sn = 0.0;
        }
        else if (abs(dy) > abs(dx))
        {
                dtype temp = dx / dy;
                sn = 1.0 / sqrt( 1.0 + temp*temp );
                cs = temp * sn;
        }
        else
        {
                dtype temp = dy / dx;
                cs = 1.0 / sqrt( 1.0 + temp*temp );
                sn = temp * cs;
        }
   }

   // FOR COMPLEX MATRICES
   void iterSolver::GeneratePlaneRotation(dtype& dx, dtype& dy,dtype& cs, dtype& sn)
   {
   if (dx == dtype(0))
   {
           cs = dtype(0);
           sn = dtype(1);
   }
   else
   {
           dtype_base scale = abs_(dx) + abs_(dy);
           dtype_base norm = scale * std::sqrt(abs_(dx / scale) * abs_(dx / scale) +
                                               abs_(dy / scale) * abs_(dy / scale));
           dtype alpha = dx / abs_(dx);
           cs = abs_(dx) / norm;
           sn = alpha * conj_(dy) / norm;
   }
      }
 */
