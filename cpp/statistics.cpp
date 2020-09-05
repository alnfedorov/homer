/*
 * This is a file devoted to statistics
 */
#include "statistics.h"

double logPoisson(int x, double lambda) {
	double rv = ((double)x)*log(lambda) - lambda - factln(x);
	return rv;
}

double ilogCumulativePoisson(int x, double lambda) {
	double P = 0.0;
	double diff = 1;
	int limit = (int)lambda*10+x;
	for (int i=x;i<=limit;i++) {
		double p = logPoisson(i, lambda);
		if (i==x) {
			P = p;
		} else {
			diff = log(1+exp(p-P));
			P += diff;
		}
		if (diff <= 0) break;
	}
	return P;
}


double logbetaiD(double a, double b, double x) {
    if (x <= 0.0 || x >= 1.0) {
        fprintf(stderr, "X is not right in betai (%lf)\n",x);
        //exit(0);
    }
    //if (x==0.0 || x == 1.0) bt = -1*FLT_MAX;
    double bt = gammlnD(a+b)-gammlnD(a)-gammlnD(b)+a*log(x)+b*log(1.0-x);
    if (x < ((a+1.0)/(a+b+2.0))) {
        return bt + log(betacfD(a,b,x))-log(a);
	} else {
        return log(1.0-exp(bt)*betacfD(b,a,1.0-x)/b);
	}
}


#define MAXIT 1000
#define EPS 3.0e-7
#define FPMIN 1.0e-30

//#define EPS 3.0e-7
#define FPMIND 1.0e-30
double betacfD(double a, double b, double x) {
    int m,m2;
    double aa, c,d,del,h,qab,qam,qap;
    qab=a+b;
    qap=a+1.0;
    qam=a-1.0;
    c=1.0;
    d=1.0-qab*x/qap;
    if (fabs(d) < FPMIND) d=FPMIND;
    d=1.0/d;
    h=d;
    for (m=1;m<=MAXIT;m++) {
        m2=2*m;
        aa=m*(b-m)*x/((qam+m2)*(a+m2));
        d=1.0+aa*d;
        if (fabs(d) < FPMIND) d=FPMIND;
        c=1.0+aa/c;
        if (fabs(c) < FPMIND) c=FPMIND;
        d=1.0/d;
        h *= d*c;
        aa = -(a+m)*(qab+m)*x/((a+m2)*(qap+m2));
        d=1.0+aa*d;
        if (fabs(d) < FPMIND) d=FPMIND;
        c=1.0+aa/c;
        if (fabs(c) < FPMIND) c=FPMIND;
        d=1.0/d;
        del=d*c;
        h*=del;
        if (fabs(del-1.0) < EPS) break;
    }
    if (m > MAXIT) {
        fprintf(stderr, "Error calculating betacf\n");
        //exit(0);
    }
    return h;
}


float gammln(float xx) {
    double x,y,tmp,ser;
    static double cof[6]={76.18009172947146,-86.50532032941677,
            24.01409824083091,-1.231739572450155,
            0.1208650973866179e-2,-0.5396239384953e-5};
    int j;
    y=x=xx;
    tmp=x+5.5;
    tmp -= (x+0.5)*log(tmp);
    ser=1.000000000190015;
    for (j=0;j<=5;j++) ser += cof[j]/++y;
    return -tmp+log(2.5066282746310005*ser/x);
}

//float bicoln(int n, int k);
//float factln(int n);
//float gammln(float xx);
//float hypergeo(int N, int n1, int n2, int n);


// double routines


float factln(int nn) {
    static float a[HYPERGEO_CACHE_SIZE+1];
    if (nn<0) { 
		fprintf(stderr,"negative factorial %d\n", nn);
		//fprintf(stdout,"[%d %d %d %d]",N,n1,n2,n);
		//fprintf(stdout,"1.0");
		//exit(0);
		return 0.0;
    } else if (nn <= 1) {
		return 0.0;
    } else if (nn <= HYPERGEO_CACHE_SIZE) {
		return a[nn] ? a[nn] : (a[nn]=gammln(nn+1.0));
    } else {
		return gammln(nn+1.0);
    }
}

double gammlnD(double xx) {
    double x,y,tmp,ser;
    static double cof[6]={76.18009172947146,-86.50532032941677,
            24.01409824083091,-1.231739572450155,
            0.1208650973866179e-2,-0.5396239384953e-5};
    unsigned int j;
    y=x=xx;
    tmp=x+5.5;
    tmp -= (x+0.5)*log(tmp);
    ser=1.000000000190015;
    for (j=0;j<=5;j++) ser += cof[j]/++y;
    return -tmp+log(2.5066282746310005*ser/x);
}

float gammp(float a, float x) {
	float gamser=0.0,gammcf=0.0,gln=0.0;
	if (x < 0.0 || a <= 0.0) {
		return 0;
	}
	if (x < (a+1.0)) {
		gser( a, x, gamser, gln);
		return gamser;
	} else {
		gcf(a, x, gammcf, gln);
		return 1.0-gammcf;
	}
}
float gammq(float a, float x) {
	float gamser = 0.0;
	float gln = 0.0;
	if (x< 0 || a <= 0) {
		fprintf(stderr, "Bad inputs (%f,%f) into gammq\n", a, x);
		return 0.0;
	}
	if (x < a+1.0) {
		gser(a,x,gamser,gln);
		return 1.0-gamser;
	} else {
		gcf(a,x,gamser,gln);
		return gamser;
	}
}

void gser(float a, float x, float &gamser, float &gln) {
	gln = gammln(a);
	gamser =0;
	int itmax = 100;
	float eps = 3.0e-10;
	if (x <= 0.0) {
		gamser = 0;
		return;
	} else {
		float ap = a;
		float del = 1.0/a;
		float sum = del;
		for (int i=1;i<itmax;i++) {
			ap++;
			del *= x/ap;
			sum += del;
			if (fabs(del) < fabs(sum)*eps) {
				gamser = sum*exp(-1*x+a*log(x)-gln);
				return;
			}
		}
		fprintf(stderr, "reached max iterations in gser\n");
		return;
	}
}

void gcf(float a, float x, float &gammcf, float& gln) {
	int itmax = 100;
	float eps = 2.0e-8;
	gammcf = 0;
	gln = gammln(a);
	float b = x+1.0-a;
	float fpmin = 1e-30;
	float c = 1.0/fpmin;
	float d = 1.0/b;
	float h = d;
	int i=0;
	float del;
	for (i=1;i<=itmax;i++) {
		float an = i*(i-a);
		b += 2.0;
		d = an*d+b;
		if (fabs(d) < fpmin) {
			d = fpmin;
		}
		c = b+an/c;
		if (fabs(c) < fpmin) {
			c = fpmin;
		}
		d = 1.0/d;
		del = d*c;
		h*=del;
		if (fabs(del-1.0) < eps) {
			break;
		}
	}
	if (i>itmax) {
		//fprintf(stderr, "reaced maxit in gcf\n");
	}
	gammcf = exp(-1*x+a*log(x)-gln)*h;
}

float chi2pvaluelog(float chi2, int df) {
	//return gammqlog(((float)df)/2.0,chi2/2.0);
	return 0;
}



