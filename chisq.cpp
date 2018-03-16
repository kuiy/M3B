//Chi-square distribution.
//From Numerical Recipes at http://www.library.cornell.edu/nr/bookcpdf.html.

#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>

#define ITMAX 100 //Maximum allowed number of iterations.
#define EPS 3.0e-7 //Relative accuracy.
#define FPMIN 1.0e-30 //Number near the smallest representable floating-point number.

void gser(double *gamser, double a, double x, double *gln);
void gcf(double *gammcf, double a, double x, double *gln);
double gammln(double xx);

double gammq(double a, double x)

//Returns the incomplete gamma function Q(a, x) = 1 - P(a, x).

{double gamser,gammcf,gln;

if (x < 0.0 || a <= 0.0) 
	{printf("Invalid arguments in routine gammp df=%.2f dep=%.2f", a, x);
	 exit(0);
	}

if (x < (a+1.0)) //Use the series representation and take its complement.
	{gser(&gamser,a,x,&gln);
	 return 1.0-gamser; 
	} 
else //Use the continued fraction representation.
	{gcf(&gammcf,a,x,&gln);
	 return gammcf;
	}
}

void gser(double *gamser, double a, double x, double *gln)

//Returns the incomplete gamma function P(a, x) evaluated by its series representation as gamser.
//Also returns ln gamma(a) as gln.

{int n;
 double sum,del,ap;

*gln=gammln(a);
if (x <= 0.0) 
	{if (x < 0.0) 
		{printf("x less than 0 in routine gser");
		 exit(0);
		}

	 *gamser=0.0;
	 return;
	} 
else 
	{ap=a;
	 del=sum=1.0/a;
	 for (n=1;n<=ITMAX;n++) 
		{++ap;
		 del *= x/ap;
		 sum += del;
		 if (fabs(del) < fabs(sum)*EPS) 
			{*gamser=sum*exp(-x+a*log(x)-(*gln));
			 return;
			}
		}
	 printf("a too large, ITMAX too small in routine gser");
	 exit(0);

	 return;
	}
}

void gcf(double *gammcf, double a, double x, double *gln)

//Returns the incomplete gamma function Q(a, x) evaluated by its continued fraction representation as gammcf. 
//Also returns ln gamma(a) as gln.

{int i;
 double an,b,c,d,del,h;

*gln=gammln(a);
b=x+1.0-a; //Set up for evaluating continued fraction by modified Lentz’s method (5.2) with b0 = 0.
c=1.0/FPMIN;
d=1.0/b;
h=d;
for (i=1;i<=ITMAX;i++) //Iterate to convergence.
	{an = -i*(i-a);
	 b += 2.0;
	 d=an*d+b;
	 if (fabs(d) < FPMIN) 
		 d=FPMIN;
	 c=b+an/c;
	 if (fabs(c) < FPMIN) 
		 c=FPMIN;
	 d=1.0/d;
	 del=d*c;
	 h *= del;
	 if (fabs(del-1.0) < EPS) 
		 break;
	}
if (i > ITMAX) 
	{printf("a too large, ITMAX too small in gcf");
	 exit(0);
	}

*gammcf=exp(-x+a*log(x)-(*gln))*h; //Put factors in front.
}

double gammln(double xx)

//Returns the value ln[gamma(xx)] for xx > 0.

{double x,y,tmp,ser;
 static double cof[6]={76.18009172947146,-86.50532032941677,
24.01409824083091,-1.231739572450155,
0.1208650973866179e-2,-0.5395239384953e-5};
 int j;

y=x=xx;
tmp=x+5.5;
tmp -= (x+0.5)*log(tmp);
ser=1.000000000190015;
for (j=0;j<=5;j++) 
	ser += cof[j]/++y;
return -tmp+log(2.5066282746310005*ser/x);
}