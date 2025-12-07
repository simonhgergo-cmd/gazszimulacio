# define MAXSTP 1000000000
# define TINY 1.0e-30
#define SAFETY 0.9
#define PGROW -0.2
#define PSHRNK -0.25
#define ERRCON 1.89e-4

#include <complex.h>
#include "def.c"

#define SIGN(a,b) ((b)<0 ? -fabs(a) : fabs(a))

void nrerror(error_text)
char error_text[];
{
	void exit();

	printf("Numerical Recipes run-time error...\n");
	printf("%s\n",error_text);
	printf("...now exiting to system...\n");
	exit(1);
}


/*******************************************************************************
 *	Random number generator						       *
 *******************************************************************************/


#define M1 259200
#define IA1 7141
#define IC1 54773
#define RM1 (1.0/M1)
#define M2 134456
#define IA2 8121
#define IC2 28411
#define RM2 (1.0/M2)
#define M3 243000
#define IA3 4561
#define IC3 51349

float ran1(idum)
int *idum;
{
	static long ix1,ix2,ix3;
	static float r[98];
	float temp;
	static int iff=0;
	int j;
	void nrerror();

	if (*idum < 0 || iff == 0) {
		iff=1;
		ix1=(IC1-(*idum)) % M1;
		ix1=(IA1*ix1+IC1) % M1;
		ix2=ix1 % M2;
		ix1=(IA1*ix1+IC1) % M1;
		ix3=ix1 % M3;
		for (j=1;j<=97;j++) {
			ix1=(IA1*ix1+IC1) % M1;
			ix2=(IA2*ix2+IC2) % M2;
			r[j]=(ix1+ix2*RM2)*RM1;
		}
		*idum=1;
	}
	ix1=(IA1*ix1+IC1) % M1;
	ix2=(IA2*ix2+IC2) % M2;
	ix3=(IA3*ix3+IC3) % M3;
	j=1 + ((97*ix3)/M3);
	if (j > 97 || j < 1) nrerror("RAN1: This cannot happen.");
	temp=r[j];
	r[j]=(ix1+ix2*RM2)*RM1;
	return temp;
}

#undef M1
#undef IA1
#undef IC1
#undef RM1
#undef M2
#undef IA2
#undef IC2
#undef RM2
#undef M3
#undef IA3
#undef IC3



double **matrix(int nrl,int nrh,int ncl,int nch)
{
	int i;
	double **m;

	m=(double **) malloc((unsigned) (nrh-nrl+1)*sizeof(double*));
	if (!m) nrerror("allocation failure 1 in matrix()");
	m -= nrl;

	for(i=nrl;i<=nrh;i++) {
		m[i]=(double *) malloc((unsigned) (nch-ncl+1)*sizeof(double));
		if (!m[i]) nrerror("allocation failure 2 in matrix()");
		m[i] -= ncl;
	}
	return m;
}

double ***vectmatrix(int nrl,int nrh,int ncl,int nch,int nql,int nqh)
{
	int i,k;
	double ***m;

	m=(double ***) malloc((unsigned) (nrh-nrl+1)*sizeof(double**));
	if (!m) nrerror("allocation failure 1 in vectmatrix()");
	m -= nrl;

       for(i=nrl;i<=nrh;i++) {
	      m[i]=(double **) malloc((unsigned) (nch-ncl+1)*sizeof(double*));
	      if (!m[i]) nrerror("allocation failure 2 in vectmatrix()");
	      for(k=ncl;k<=nch;k++) 
		{
m[i][k]=(double *) malloc((unsigned) (nqh-nql+1)*sizeof(double));
	      if (!m[i][k]) nrerror("allocation failure 3 in vectmatrix()");
	      m[i][k] -= nql;
		}
	      m[i] -= ncl;
	}
	return m;
}

int ***ivectmatrix(int nrl,int nrh,int ncl,int nch,int nql,int nqh)
{
	int i,k;
	int ***m;

	m=(int ***) malloc((unsigned) (nrh-nrl+1)*sizeof(int**));
	if (!m) nrerror("allocation failure 1 in vectmatrix()");
	m -= nrl;

       for(i=nrl;i<=nrh;i++) {
	      m[i]=(int **) malloc((unsigned) (nch-ncl+1)*sizeof(int*));
	      if (!m[i]) nrerror("allocation failure 2 in vectmatrix()");
	      for(k=ncl;k<=nch;k++) 
		{
m[i][k]=(int *) malloc((unsigned) (nqh-nql+1)*sizeof(int));
	      if (!m[i][k]) nrerror("allocation failure 3 in vectmatrix()");
	      m[i][k] -= nql;
		}
	      m[i] -= ncl;
	}
	return m;
}

double ****matmatrix(int nrl,int nrh,int ncl,int nch,int nql,int nqh,int nwl,int nwh)
{
	int i,k,l;
	double ****m;

	m=(double ****) malloc((unsigned) (nrh-nrl+1)*sizeof(double***));
	if (!m) nrerror("allocation failure 1 in matmatrix()");
	m -= nrl;

       for(i=nrl;i<=nrh;i++) {
	      m[i]=(double ***) malloc((unsigned) (nch-ncl+1)*sizeof(double**));
	      if (!m[i]) nrerror("allocation failure 2 in matmatrix()");

	      for(k=ncl;k<=nch;k++) 
		{
m[i][k]=(double **) malloc((unsigned) (nqh-nql+1)*sizeof(double*));
	      if (!m[i][k]) nrerror("allocation failure 3 in matmatrix()");

	      for(l=nql;l<=nqh;l++) 
		{
m[i][k][l]=(double *) malloc((unsigned) (nwh-nwl+1)*sizeof(double));
	      if (!m[i][k][l]) nrerror("allocation failure 4 in matmatrix()");
	      m[i][k][l] -= nwl;
		}
	      m[i][k] -= nql;
		}
	      m[i] -= ncl;
	}
	return m;
}

int **imatrix(int nrl,int nrh,int ncl,int nch)
{
	int i;
	int **m;

	m=(int **) malloc((unsigned) (nrh-nrl+1)*sizeof(int*));
	if (!m) nrerror("allocation failure 1 in matrix()");
	m -= nrl;

	for(i=nrl;i<=nrh;i++) {
		m[i]=(int *) malloc((unsigned) (nch-ncl+1)*sizeof(int));
		if (!m[i]) nrerror("allocation failure 2 in matrix()");
		m[i] -= ncl;
	}
	return m;
}

complex double **cmatrix(int nrl,int nrh,int ncl,int nch)
{
	int i;
	complex double **m;

	m=(complex double **) malloc((unsigned) (nrh-nrl+1)*sizeof(complex double*));
	if (!m) nrerror("allocation failure 1 in cmatrix()");
	m -= nrl;

	for(i=nrl;i<=nrh;i++) {
		m[i]=(complex double *) malloc((unsigned) (nch-ncl+1)*sizeof(complex double));
		if (!m[i]) nrerror("allocation failure 2 in cmatrix()");
		m[i] -= ncl;
	}
	return m;
}

double *vector(int nl,int nh)
{
	double *v;

	v=(double *)malloc((unsigned) (nh-nl+1)*sizeof(double));
	if (!v) nrerror("allocation failure in vector()");
	return v-nl;
}

complex double *cvector(int nl,int nh)
{
  
	complex double *v;

	v=(complex double *)malloc((unsigned) (nh-nl+1)*sizeof(complex double));
	if (!v) nrerror("allocation failure in cvector()");
	return v-nl;
}

int *ivector(int nl,int nh)
{
	int *v;

	v=(int *)malloc((unsigned) (nh-nl+1)*sizeof(int));
	if (!v) nrerror("allocation failure in ivector()");
	return v-nl;
}


void free_matrix(double **m,int nrl,int nrh,int ncl,int nch)
{
	int i;

	for(i=nrh;i>=nrl;i--) free((char*) (m[i]+ncl));
	free((char*) (m+nrl));
}

void free_vectmatrix(double ***m,int nrl,int nrh,int ncl,int nch,int nql,int nqh)
{
	int i,k;

	for(i=nrh;i>=nrl;i--)
	  { 
	    for(k=nch;k>=ncl;k--)
	      {
	      free((char*) (m[i][k]+nql));
	      }
	    free((char*) (m[i]+ncl));
	  }
	free((char*) (m+nrl));
}

void free_ivectmatrix(int ***m,int nrl,int nrh,int ncl,int nch,int nql,int nqh)
{
	int i,k;

	for(i=nrh;i>=nrl;i--)
	  { 
	    for(k=nch;k>=ncl;k--)
	      {
	      free((char*) (m[i][k]+nql));
	      }
	    free((char*) (m[i]+ncl));
	  }
	free((char*) (m+nrl));
}

void free_matmatrix(double ****m,int nrl,int nrh,int ncl,int nch,int nql,int nqh,int nwl,int nwh)
{
	int i,k,l;

	for(i=nrh;i>=nrl;i--)
	  { 
	    for(k=nch;k>=ncl;k--) 
	      {
	      for(l=nqh;l>=nql;l--)
		free((char*) (m[i][k][l]+nwl));
	      free((char*) (m[i][k]+nql));
	      }
	    free((char*) (m[i]+ncl));
	  }
	free((char*) (m+nrl));
}

void free_imatrix(int **m,int nrl,int nrh,int ncl,int nch)
{
	int i;

	for(i=nrh;i>=nrl;i--) free((char*) (m[i]+ncl));
	free((char*) (m+nrl));
}

void free_cmatrix(complex **m,int nrl,int nrh,int ncl,int nch)
{
	int i;

	for(i=nrh;i>=nrl;i--) free((char*) (m[i]+ncl));
	free((char*) (m+nrl));
}


void free_vector(double *v,int nl,int nh)
{
	free((char*) (v+nl));
}

void free_ivector(int *v,int nl,int nh)
{
	free((char*) (v+nl));
}

void free_cvector(complex *v,int nl,int nh)
{
	free((char*) (v+nl));
}

double singsimps0 (double ya, double yb,  double *a, int STEPS)
{
int i;
double f,s,h;

h=(yb-ya)/(double)(STEPS+1);

 s=0.0;

 for(i=0; i<=STEPS; i++)
			
			{
			  f=a[i];
			
			if ((i==0)||(i==STEPS))
				s+=f;
			else	 
				s+=(4.0*f*(i % 2)+2.0*f*((i+1) % 2));
			}
			s*=h/3.0;
			
return s;
}

double asysimps (double **a, int S1, int S2, double int1, double int2)
{
int i,k;
double h1,h2,f,si,sk;

h1=int1/(double)S1;
h2=int2/(double)S2;

 for(i=0; i<=S1; i++)
   {
     sk=0.0;
     for(k=0; k<=S2; k++)
       {
	 f=a[i][k];
	 if ((k==0)||(k==S2))
	   sk+=f;
	 else	 
	   sk+=(4.0*f*(k % 2)+2.0*f*((k+1) % 2));
       }
     sk*=h1/3.0;
     if ((i==0)||(i==S1))
       si+=sk;
     else	 
       si+=(4.0*sk*(i % 2)+2.0*sk*((i+1) % 2));
   }
			si*=h2/3.0;
			
			return si;
			
}












