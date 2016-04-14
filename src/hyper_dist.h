/*###########################################
# Descriptions   : hyper_dist
# Usage          : 
# Parameters	 : none
# Sample Input   : 
# Sample Output  : 
# Depedency      : none
# Temp File      : none
# Comments       : none
# See Also       : none
# Data           : 
# Template       : Last modified data 08/14/18
# Author         : setupX
############################################ */

#ifndef _hyper_dist_setupx
#define _hyper_dist_setupx
#include "statistics.h"

NUM dhyper_lanczos_imp(unsigned x, unsigned r, unsigned n, unsigned N)
{
   return exp(
      lgamma( (NUM)(n + 1) ) 
			+ lgamma( (NUM) (r + 1) )
      + lgamma( (NUM) (N - n + 1) ) 
			+ lgamma( (NUM) (N - r + 1) )
      - lgamma( (NUM) (N + 1) )
      - lgamma( (NUM) (x + 1) )
      - lgamma( (NUM) (n - x + 1) )
      - lgamma( (NUM) (r - x + 1) )
      - lgamma( (NUM) (N - n - r + x + 1) ) );
}

double dhyper_factorial_imp(unsigned x, unsigned r, unsigned n, unsigned N)
{
   double result = factorial(n);
   double num[3] = {
      factorial(r),
      factorial(N - n),
      factorial(N - r)
   };
   double denom[5] = {
      factorial(N),
      factorial(x),
      factorial(n - x),
      factorial(r - x),
      factorial(N - n - r + x)
   };
   int i = 0;
   int j = 0;
   while((i < 3) || (j < 5))
   {
      while((j < 5) && ((result >= 1) || (i >= 3)))
      {
         result /= denom[j];
         ++j;
      }
      while((i < 3) && ((result <= 1) || (j >= 5)))
      {
         result *= num[i];
         ++i;
      }
   }
   return result;
}

// choose(n,x)choose(N-n,r-x)/choose(N,r)
// var r * (n / N) * (1 - n / N) * (N - r) / (N - 1);
// mode floor((r + 1) * (n + 1) / (N + 2))
// skewness (N - 2 * n) * sqrt(N - 1) * (N - 2 * r) / (sqrt(n * r * (N - n) *
// (N - r)) * (N - 2))
// kurtosis_excess
// NUM t1 = N * N * (N - 1) / (r * (N - 2) * (N - 3) * (N - r));
// NUM t2 = (N * (N + 1) - 6 * N * (N - r)) / (n * (N - n))
//                + 3 * r * (N - r) * (N + 6) / (N * N) - 6;
// return t1 * t2;
//	kurtosis 
//	kurtosis_excess(dist) + 3;
NUM dhyper(unsigned x, unsigned r, unsigned n, unsigned N){
	unsigned lower_limit = MAX(0, (int)(n+r)-(int)(N));
	unsigned u = MIN(r, n);
	NUM result;
	if(r > N || n > N 
			|| x > r || x > n
			|| x < 0 || x < lower_limit)
	{
		_error(1);
		return nanf();
	}
	if(N <= MAX_FACTORIAL)
	{
		// If N is small enough then we can evaluate the PDF via the factorials
		// directly: table lookup of the factorials gives the best performance
		// of the methods available:
		result = dhyper_factorial_imp(x, r, n, N);
	}
	else
	{
		// Catch all case - use the lanczos approximation - where available - 
		// to evaluate the ratio of factorials.  This is reasonably fast
		// (almost as quick as using logarithmic evaluation in terms of lgamma)
		// but only a few digits better in accuracy than using lgamma:
		result = dhyper_lanczos_imp(x, r, n, N);
	}
	if(result > 1)
	{
		result = 1;
	}
	if(result < 0)
	{
		result = 0;
	}
	return result;
}


NUM phyper(unsigned x, unsigned r, unsigned n, unsigned N){
	int mode = floor( (int)(r + 1)*(int)(n + 1) / (N + 2) );
	NUM result = 0.0;
	NUM diff;
	unsigned lower_limit;
	if(x < mode)
	{
		result = dhyper(x, r, n, N);
		diff = result;
		lower_limit = MAX(0, (int)(n+r)-(int)(N));
		while(diff > (result * EPS))
		{
			diff = (int)x * (int)((N + x) - n - r) * diff 
				/ ((int)(1 + n - x) * (int)(1 + r - x));
			result += diff;
			if(x == lower_limit)
				break;            
			--x;        
		}
	}
	else
	{
		unsigned upper_limit = MIN(r, n);
		if(x != upper_limit)
		{
			++x;
			result = dhyper(x, r, n, N);
			diff = result;
			while( (x <= upper_limit) && (diff > EPS) )
			{
				diff = (int)(n - x) * (int)(r - x) * diff 
					/ ((int)(x + 1) * (int)((N + x + 1) - n - r));
				result += diff;
				++x;
			}
		}
		result=1-result;
	}
	if(result>1)
		result = 1;
	if(result < 0)
		result = 0;
	return result;
}

#endif //_hyper_dist_setupx

