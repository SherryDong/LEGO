/*###########################################
# Descriptions   : statistics
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

#ifndef _statistics_setupx
#define _statistics_setupx
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "io_filter.h"


#define PI 3.14159265359
//	gamma function
#define ITMAX 100
#define EPS 3.0e-9
#define FPMIN 1.0e-30 //Number near the smallest representable
					  //floating-point number.
void _error(int error_id){
	printf("Error ");
	switch(error_id){
		case 1: printf("%d: input error\n",error_id);
		case 2: printf("%d: ITMAX too small\n",error_id);
	}

}
NUM max(NUM *in, int len){
	int i;
	NUM max=in[0];
	NUM ini;
	for (i = 1; i < len; i++){
		ini=in[i];
		if(ini > max){
			max = ini;
		}

	}
	return max;
}
NUM min(NUM *in, int len){
	int i;
	NUM min=in[0];
	NUM ini;
	for (i = 1; i < len; i++){
		ini=in[i];
		if(ini < min){
			min = ini;
		}

	}
	return min;
}
int rank(NUM *in, int len, NUM **out){
	return 0;
}
int order(NUM *in, int len, NUM **out){
	return 0;
}
NUM psum(int power, NUM *in, int len){
	int i;
	NUM psum=0;
	for (i = 0; i < len; i++){

		psum += pow(in[i],power);
	}

	return psum;
}
NUM sumsq(NUM *in, int len){
	int i;
	NUM sumsq=0;
	for (i = 0; i < len; i++){
		sumsq += in[i]*in[i];
	}

	return sumsq;
}
NUM sum(NUM *in, int len){
	int i;
	NUM sum=0;
	for (i = 0; i < len; i++){

		sum += in[i];
	}

	return sum;
}
NUM avg(NUM *in, int len){
	int i;
	NUM avg=0.0;
	for (i = 0; i < len; i++){
		avg += in[i];
	}
	avg /= len;
	return avg;
}
NUM mean(NUM *in, int len){
	return avg(in, len);
}
//	standard power sum
NUM spsum(int power, NUM *in, int len){
	int i;
	NUM mean;
	mean = avg(in, len);
	NUM spsum=0;
	for (i = 0; i < len; i++){
		spsum+=pow(in[i]-mean,power);
	}
	return spsum;
}
NUM moment(int power, NUM *in, int len){
	return spsum(power, in, len)/len;
}
NUM var(NUM *in, int len){
	NUM s1=psum(2, in, len);
	NUM s2=sum(in, len);
	return s1/(len-1)-(s2*s2)/(len*(len-1));
}
NUM sd(NUM *in, int len){
	return sqrt(var(in, len));
}
NUM sd_sum(NUM *in, int len, NUM sum){
	return sqrt(sumsq(in, len)/(len-1)-(sum*sum)/(len*(len-1)));
}
NUM inner_prod(NUM *in1, NUM *in2, int len){
	int i;
	NUM prod=0;
	for (i = 0; i < len; i++){
		prod+=in1[i]*in2[i];
	}
	return prod;
}
NUM cov(NUM *in1, NUM *in2, int len){
	int i;
	NUM sum=0;
	NUM mean1=avg(in1, len);
	NUM mean2=avg(in2, len);
	for (i = 0; i < len; i++){
		sum+=(in1[i]-mean1)*(in2[i]-mean2);
	}
	return sum/(len -1);
}
NUM cor(NUM *in1, NUM *in2, int len){
	NUM pcd;
	pcd=cov(in1, in2, len)/(sd(in1, len) * sd(in2, len));
	return pcd;
}

void Z_transform(NUM **in, int len){
	int i;
	NUM vsum=sum(*in, len);
	NUM vavg=vsum/(NUM)len;
	NUM vsd=sd_sum(*in, len, vsum);
	for(i=0; i<len; i++){
		(*in)[i] = ( (*in)[i] - vavg ) / vsd;
	}
}

void unit_magnitude(NUM **in, int len){
	int i;
	NUM mag=sqrt(sumsq(*in, len));
	for(i=0; i<len; i++){
		(*in)[i] =  (*in)[i] / mag;
	}
}

//	Gamma Function
NUM lgamma(NUM in)
//	Returns the value ln[.(in)] for in > 0.
{
//	Internal arithmetic will be done in double precision, a nicety that you can omit if five-figure
//	accuracy is good enough.
	NUM x,y,tmp,ser;
	static NUM cof[6]={76.18009172947146,
		-86.50532032941677,
		24.01409824083091,
		-1.231739572450155,
		0.1208650973866179e-2,
		-0.5395239384953e-5
	};
	int j;
//	Stirling's approximation
	y=x=in;
	tmp=x+5.5;
	tmp -= (x+0.5)*log(tmp);

	ser=1.000000000190015;
	for (j=0;j<=5;j++) ser += cof[j]/++y;

	return -tmp+log(2.5066282746310005*ser/x);
}
NUM lfactorial(int in){
	return lgamma(in+1);
}
#ifndef MAX_FACTORIAL
#define MAX_FACTORIAL 34
#endif
NUM factorial(int in)
//	Returns the value in!.
{
	static int ntop=4;
//Fill in table only as required.
	static double a[35]={
      1.0F,
      1.0F,
      2.0F,
      6.0F,
      24.0F,
      120.0F,
      720.0F,
      5040.0F,
      40320.0F,
      362880.0F,
      3628800.0F,
      39916800.0F,
      479001600.0F,
      6227020800.0F,
      87178291200.0F,
      1307674368000.0F,
      20922789888000.0F,
      355687428096000.0F,
      6402373705728000.0F,
      121645100408832000.0F,
      0.243290200817664e19F,
      0.5109094217170944e20F,
      0.112400072777760768e22F,
      0.2585201673888497664e23F,
      0.62044840173323943936e24F,
      0.15511210043330985984e26F,
      0.403291461126605635584e27F,
      0.10888869450418352160768e29F,
      0.304888344611713860501504e30F,
      0.8841761993739701954543616e31F,
      0.26525285981219105863630848e33F,
      0.822283865417792281772556288e34F,
      0.26313083693369353016721801216e36F,
      0.868331761881188649551819440128e37F,
      0.29523279903960414084761860964352e39F
	};
	int j;
	if (in < 0) _error(1);
	if (in > 34) return exp(lfactorial(in));
//	Larger value than size of table is required. Actually, this big a value is going to over.ow
//	on many computers, but no harm in trying.
	return a[in];
}
NUM lchoose(int n,int k){
	return lfactorial(n)-lfactorial(k)-lfactorial(n-k);
}
NUM choose(int n, int k){
	return floor(0.5 + exp(lchoose(n,k)));
}
NUM lbeta(NUM z, NUM w){
	return lgamma(z)+lgamma(w)-lgamma(z+w);
}

NUM _gamma_ser(NUM a, NUM x){
//	Returns the incomplete gamma function P(a; x) evaluated by its series representation as gamser.
	int n;
	NUM sum,del,ap;
	
	NUM gln=lgamma(a);
	if (x <= 0.0) {
		if (x < 0.0) error(1);
		return 0;
	} else {
		ap=a;
		del=sum=1.0/a;
		for (n=1;n<=ITMAX;n++) {
			++ap;
			del *= x/ap;
			sum += del;
			if (fabs(del) < fabs(sum)*EPS) {
				return sum*exp(-x+a*log(x)-(gln));
			}
		}
		_error(2);
		return 0.0;
	}
}
NUM _gamma_cf(NUM a, NUM x){
// Returns the incomplete gamma function Q(a; x) evaluated by its continued fraction representation
	int i;
	NUM an,b,c,d,del,h;
	NUM gln=lgamma(a);
	b=x+1.0-a;
	c=1.0/FPMIN;
	d=1.0/b;
	h=d;
	for (i=1;i<=ITMAX;i++) {
		an = -i*(i-a);
		b += 2.0;
		d=an*d+b;
		if (fabs(d) < FPMIN) d=FPMIN;
		c=b+an/c;
		if (fabs(c) < FPMIN) c=FPMIN;
		d=1.0/d;
		del=d*c;
		h *= del;
		if (fabs(del-1.0) < EPS) break;
	}
	if (i > ITMAX) _error(2);
	return exp(-x+a*log(x)-(gln))*h;
}
NUM pgamma(NUM a, NUM x){
	if (x < 0.0 || a <= 0.0) _error(1);

	if (x < (a+1.0)) { 	//Use the series representation.
		return _gamma_ser(a,x);
	} else { 			//Use the continued fraction representation
		return 1.0-_gamma_cf(a,x); //and take its complement.
	}
}
NUM pchisq(NUM chisq, NUM df){
	return pgamma(df/2, chisq/2);
}
NUM _beta_cf(NUM a, NUM b, NUM x){
//	Evaluates continued fraction for incomplete beta function by modified Lentz's method 
	int m,m2;
	NUM aa,c,d,del,h,qab,qam,qap;
	qab=a+b; //These q's will be used in factors that occur in the coecients. 
	qap=a+1.0;
	qam=a-1.0;
	c=1.0; //First step of Lentz's method.
	d=1.0-qab*x/qap;
	if (fabs(d) < FPMIN) d=FPMIN;
	d=1.0/d;
	h=d;
	for (m=1;m<=ITMAX;m++) {
		m2=2*m;
		aa=m*(b-m)*x/((qam+m2)*(a+m2));
		d=1.0+aa*d; //One step (the even one) of the recurrence.
		if (fabs(d) < FPMIN) d=FPMIN;
		c=1.0+aa/c;
		if (fabs(c) < FPMIN) c=FPMIN;
		d=1.0/d;
		h *= d*c;
		aa = -(a+m)*(qab+m)*x/((a+m2)*(qap+m2));
		d=1.0+aa*d; //Next step of the recurrence (the odd one).
		if (fabs(d) < FPMIN) d=FPMIN;
		c=1.0+aa/c;
		if (fabs(c) < FPMIN) c=FPMIN;
		d=1.0/d;
		del=d*c;
		h *= del;
		if (fabs(del-1.0) < EPS) break; //Are we done?
	}
	if (m > ITMAX) _error(2);
	return h;
}
NUM pbeta(NUM a, NUM b, NUM x){
	NUM bt;
	if (x < 0.0 || x > 1.0) _error(1);
	if (x == 0.0 || x == 1.0) 
		bt=0.0;
	else 	//	Factors in front of the continued fraction.
		bt=exp(lgamma(a+b)-lgamma(a)-lgamma(b)+a*log(x)+b*log(1.0-x));

	if (x < (a+1.0)/(a+b+2.0)) //Use continued fraction directly.
		return bt*_beta_cf(a,b,x)/a;
	else 	//	Use continued fraction after making the symmetry transformation. 
		return 1.0-bt*_beta_cf(b,a,1.0-x)/b;
}
NUM pt(NUM t, NUM df){
	NUM sq=sqrt(pow(t,2)+df);
	return pbeta( df/2, df/2, (t+sq)/(2*sq) );
}
NUM pf(NUM f, NUM df1, NUM df2){
	return 1-pbeta( df2/2, df1/2, df2/(df2+df1*f) );
}
NUM pbinom(NUM k, int n, NUM p){
	return pbeta( k, n-k+1, p);
}
int freqs(NUM *in, NUM **pout, int len){
	NUM s=sum(in, len);
	int i;
	for(i=0; i<len; i++){
		(*pout)[i]=in[i]/s;
	}
	return 0;
}
NUM entropy(NUM *in, int len){
	NUM H=0.0;
	NUM f[len]; //store freqs
	NUM *pf=f;
	int i;
	freqs(in, &pf, len);
	for(i=0; i<len; i++){
		H-=f[i]*log(f[i]);
	}
	return H;
}
NUM mb(NUM **x, int xlen, int ylen){
	return 0;
}

NUM mi(int xlen, int ylen, NUM x[][ylen]){
	int xylen=xlen*ylen;
	NUM xy[xylen];
	NUM xf[ylen],yf[xlen];
	int i,j;
	for (i=0; i<ylen; i++){
		yf[i]=0;
	}
	for (j=0; j<xlen; j++){
		xf[j]=0;
	}
	for(i=0; i<xlen; i++){
		for(j=0; j<ylen; j++){
			xy[i+j*xlen]=x[i][j];
			xf[j]+=x[i][j];
			yf[i]+=x[i][j];
		}
	}
	return entropy(xf,ylen)+entropy(yf,xlen)-entropy(xy,xylen);
}
// two dimension MI calculate
NUM mi2d(NUM *x, NUM *y, int len){
	NUM mi2d[2][len];
	int i;
	for(i=0; i<len; i++){
		mi2d[0][i]=x[i];
		mi2d[1][i]=y[i];
	}
	return mi(2,len,mi2d);
}
#endif //_statistics_setupx

