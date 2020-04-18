//blitz++ wrapper around the Gaussian quadrature rule generator "gqgenlegendre.cpp"

#include "gaussquad.h"

#include "gqgenlegendre.h"

void gaussquadrule(const int& n, Array<double,1>& x, Array<double,1>& w)
{
	//define the custom gqgenlegendre arrays and fill them
	ap::real_1d_array X,W;
	buildgausslegendrequadrature(n,X,W);
	//then copy the values to blitz++ arrays
	x.resize(n);
	w.resize(n);
	for (int i=0;i<n;i++)
	{
		x(i)=X(i);
		w(i)=W(i);
	}
}
