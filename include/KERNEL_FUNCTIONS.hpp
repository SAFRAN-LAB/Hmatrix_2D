#ifndef KERNEL_FUNCTIONS_HPP
#define KERNEL_FUNCTIONS_HPP

#include "myHeaders.hpp"

namespace Ker
{
dtype_base a = 0.001;
dtype_base Kchoice = 0;
dtype_base Kii = 100;


// Singular Kernels
dtype_base Ker1(const dtype_base r)
{
        return (r >= a) ? log(r)/log(a) : (r*(log10(r)-1))/(a*(log10(a)-1));
}

dtype_base Ker2(const dtype_base r)
{
        return (r >= a) ? (a/r) : (r/a);
}

dtype_base Ker3(const dtype_base r)
{
        return (r >= a) ? (r*r*log(r))/(a*a*log(a)) : (pow(r,3)*(3*log10(r)-1))/(pow(a,3)*(3*log10(a)-1));
}
// 1/r Kernel
dtype_base Ker4(const dtype_base r)
{
        return 1/r;
}

dtype_base inline KernelFunc(const dtype_base rowX, const dtype_base rowY,const dtype_base colX,const dtype_base colY)
{
        dtype_base r = dtype_base(0);
        r = (rowX - colX) * (rowX - colX) + (rowY - colY) * (rowY - colY);
        // Select the Kernel function
        if(r == dtype_base(0))
                r = Kii;
        else
        {
                r = sqrt(r);
                if (Kchoice == 1)
                        r = Ker1(r);
                else if (Kchoice == 2)
                        r = Ker2(r);
                else if (Kchoice == 3)
                        r = Ker3(r);
                else
                        r = Ker4(r);
        }
        return r;
}
}

#endif
