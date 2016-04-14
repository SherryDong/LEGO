/* Created by Sherry 
 * Network start with Node:1
 * self-loop will be calculated
 * */

#include "hyper_dist.h"
#include "statistics.h"

int main (int argc, char *argv[]) {	
// argc	
	int a,b,c,d;
	int q,m,n,k; // phyper in R
	double pval;
	if (argc <= 2) {
		fprintf(stderr, "Usage: %s <D> <B> <C> <A>\nsame as fisher.test(matrix(c(D,B,C,A),2,2),\"greater\") in R\n",argv[0]); 
		exit(0); 
	} 
// argc --> par	
	a = atoi(argv[1]);
	b = atoi(argv[2]);
	c = atoi(argv[3]);
	d = atoi(argv[4]);
	// R: phyper
	q = a-1;
	m = a+c;
	n = b+d;
	k = a+b;
	pval=1-phyper(q,k,m,m+n);
	printf("%e\n",pval);
	return(0);
}

