#include "def.c"

void rkqs(double**, double**, int, int, double* , double, 
	  double, double**, double*, double*, 
	  void (*) (double, double**, double**));
void rkck(double**, double**, int, int, double, double, 
	  double**, double**, 
	  void (*)(double, double**, double**));


void rkck(double **y, double **dydx, int n, int m, double x, double h,
double **yout, double **yerr,
void (*derivs)(double, double**, double**))
{
	int i, k;
	static double a2=0.2, a3=0.3, a4=0.6, a5=1.0, a6=0.875, b21=0.2,
	b31=3.0/40.0, b32=9.0/40.0, b41=0.3, b42=-0.9, b43=1.2,
	b51=-11.0/54.0, b52=2.5, b53=-70.0/27.0, b54=35.0/27.0,
	b61=1631.0/55296.0, b62=175.0/512.0, b63=575.0/13824.0,
	b64=44275.0/110592.0, b65=253.0/4096.0, c1=37.0/378.0,
	c3=250.0/621.0, c4=125.0/594.0, c6=512.0/1771.0,
	dc5=-277.0/14336.0;

	double dc1=c1-2825.0/27648.0, dc3=c3-18575.0/48384.0,
	dc4=c4-13525.0/55296.0, dc6=c6-0.25;
	double  **ak2, **ak3, **ak4, **ak5, **ak6, **ytemp;


	ak2=matrix(0,n-1,0,m-1);
	ak3=matrix(0,n-1,0,m-1);
	ak4=matrix(0,n-1,0,m-1);
	ak5=matrix(0,n-1,0,m-1);
	ak6=matrix(0,n-1,0,m-1);
	ytemp=matrix(0,n-1,0,m-1);


		for (i=0; i<n; i++)
		  for (k=0; k<m; k++)
		  ytemp[i][k]=y[i][k]+b21*h*dydx[i][k];
		(*derivs)(x+a2*h, ytemp, ak2);


		for (i=0; i<n; i++)
		  for (k=0; k<m; k++)
		  ytemp[i][k]=y[i][k]+h*(b31*dydx[i][k]+b32*ak2[i][k]);
		(*derivs)(x+a3*h, ytemp, ak3);

		//printf("in rk %g\n",dtrace(ytemp,X,n/2));

		for (i=0; i<n; i++)
		  for (k=0; k<m; k++)
		  ytemp[i][k]=y[i][k]+h*(b41*dydx[i][k]+b42*ak2[i][k]+b43*ak3[i][k]);
		(*derivs)(x+a4*h, ytemp, ak4);

		//printf("in rk %g\n",dtrace(ytemp,X,n/2));
		for (i=0; i<n; i++)
		  for (k=0; k<m; k++)
		  ytemp[i][k]=y[i][k]+h*(b51*dydx[i][k]+b52*ak2[i][k]+		b53*ak3[i][k]+b54*ak4[i][k]);
		(*derivs)(x+a5*h, ytemp, ak5);

		//printf("in rk %g\n",dtrace(ytemp,X,n/2));

		for (i=0; i<n; i++)
		  for (k=0; k<m; k++)
		  ytemp[i][k]=y[i][k]+h*(b61*dydx[i][k]+b62*ak2[i][k]+
				   b63*ak3[i][k]+b64*ak4[i][k]+b65*ak5[i][k]);
		(*derivs)(x+a6*h, ytemp, ak6);

		//printf("in rk %g\n",dtrace(ytemp,X,n/2));

		for (i=0; i<n; i++)
		  for (k=0; k<m; k++)
		  yout[i][k]=y[i][k]+h*(c1*dydx[i][k]+c3*ak3[i][k]+c4*ak4[i][k]+c6*ak6[i][k]);

		//printf("in rk yout %g\n",dtrace(yout,X,n/2));

		for (i=0; i<n; i++)
		  for (k=0; k<m; k++)
		    { 
		      yerr[i][k]=h*(dc1*dydx[i][k]+dc3*ak3[i][k]+dc4*ak4[i][k]+
				 dc5*ak5[i][k]+dc6*ak6[i][k]);
		    }

		free_matrix(ytemp,0,n-1,0,m-1);
		free_matrix(ak6,0,n-1,0,m-1);
		free_matrix(ak5,0,n-1,0,m-1);
		free_matrix(ak4,0,n-1,0,m-1);
		free_matrix(ak3,0,n-1,0,m-1);
		free_matrix(ak2,0,n-1,0,m-1);
}


void rkqs (double **y, double **dydx, int n, int m, 
double *x, double htry, double eps, 
double **yscal, double *hdid, double *hnext, 
void (*derivs) (double, double**, double**))

{
	int i, k;
	double errmax, h, htemp, xnew;
	double **yerr, **ytemp;

	yerr=matrix(0,n-1,0,m-1);
	ytemp=matrix(0,n-1,0,m-1);
		h=htry;

		for(;;)
			{
			  rkck(y, dydx, n, m, *x, h, ytemp, yerr, derivs);
			  errmax=0.0;
			  for (i=0; i<n; i++)
			    for (k=0; k<m; k++)
			    { 
				errmax=((errmax>cabs(yerr[i][k]/yscal[i][k]))?errmax:fabs(yerr[i][k]/yscal[i][k]));
			    }
			  
			  errmax/=eps;
			  if (errmax<=1.0) break;
			  htemp=SAFETY*h*pow(errmax, PSHRNK);
			  
			  if (h>=0.0)
			    h=((htemp>0.1*h)?htemp:0.1*h);
			else
			  h=((htemp<0.1*h)?htemp:0.1*h);
			  xnew=(*x)+h;
			  if (xnew==*x) nrerror("stepsize undeflow in rkqs \n");
			}
		
		if (errmax>ERRCON) *hnext=SAFETY*h*pow(errmax,PGROW);
		else *hnext=5.0*h;
		*x+=(*hdid=h);
		for (i=0; i<n; i++) 
		  for (k=0; k<m; k++)
		    y[i][k]=ytemp[i][k];

		free_matrix(ytemp,0,n-1,0,m-1);
		free_matrix(yerr,0,n-1,0,m-1);
}
