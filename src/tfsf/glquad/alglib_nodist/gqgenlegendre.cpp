/*************************************************************************
Copyright (c) 2005-2007, Sergey Bochkanov (ALGLIB project).

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

- Redistributions of source code must retain the above copyright
  notice, this list of conditions and the following disclaimer.

- Redistributions in binary form must reproduce the above copyright
  notice, this list of conditions and the following disclaimer listed
  in this license in the documentation and/or other materials
  provided with the distribution.

- Neither the name of the copyright holders nor the names of its
  contributors may be used to endorse or promote products derived from
  this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*************************************************************************/

// #include <stdafx.h>
#include "gqgenlegendre.h"

/*************************************************************************
Computation of nodes and weights for a Gauss-Legendre quadrature formula

The  algorithm  calculates  the  nodes  and  weights of the Gauss-Legendre
quadrature formula on domain [-1, 1].

Input parameters:
    n   �   a required number of nodes.
            n>=1.

Output parameters:
    x   -   array of nodes.
            Array whose index ranges from 0 to N-1.
    w   -   array of weighting coefficients.
            Array whose index ranges from 0 to N-1.

The algorithm was designed by using information from the QUADRULE library.
*************************************************************************/
void buildgausslegendrequadrature(int n,
     ap::real_1d_array& x,
     ap::real_1d_array& w)
{
    int i;
    int j;
    double r;
    double r1;
    double p1;
    double p2;
    double p3;
    double dp3;

    x.setbounds(0, n-1);
    w.setbounds(0, n-1);
    for(i = 0; i <= (n+1)/2-1; i++)
    {
        r = cos(ap::pi()*(4*i+3)/(4*n+2));
        do
        {
            p2 = 0;
            p3 = 1;
            for(j = 0; j <= n-1; j++)
            {
                p1 = p2;
                p2 = p3;
                p3 = ((2*j+1)*r*p2-j*p1)/(j+1);
            }
            dp3 = n*(r*p3-p2)/(r*r-1);
            r1 = r;
            r = r-p3/dp3;
        }
        while(fabs(r-r1)>=ap::machineepsilon*(1+fabs(r))*100);
        x(i) = r;
        x(n-1-i) = -r;
        w(i) = 2/((1-r*r)*dp3*dp3);
        w(n-1-i) = 2/((1-r*r)*dp3*dp3);
    }
}



