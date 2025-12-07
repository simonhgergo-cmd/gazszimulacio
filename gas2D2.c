#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "additional.c"
#include "int2D.c"
#include "def.c"
#include <string.h>   // for memcpy, memset, memmove

//fordítás gcc-vel: gcc -O3 -march=native -ffast-math -funroll-loops -fno-math-errno -lm -o run gas2D2.c
//futtatás Linux alatt: ./run
//Legyen minden segédfile ugyabban a könyvtárban

FILE *inp;
FILE *out1;


double sigma, epsilon, start, stop, spaceforaball, spfb1D, maxball, dxsav, Dcoeff, g_accel, sim_radius, mass, radius, diameter, eps, h, xlimit, ylimit, maxspeed, minspeed, **distr, **dx, **dy, **pairdist;
int Nballs, inthewall, touching, dim1, dim2, seed, *nok, *nbad;
float wallFup, wallFdown, wallFleft, wallFright;
int vdistribution[10]; // sebesség-eloszlás

void getdiff (double, double**, double**);
void initialize();
void getnumbers();
double energy(double**);
void writetofile(double**, double);

double x;
double molecule_loss, wall_loss, molecule_k, wall_k, molecule_e, wall_e, c;
double e, pi;
int grid[40][40][1000];
int grid_count[40][40];
double box_width;
int selection_x[5];
int selection_y[5];

FILE *fp;

void odeint(FILE*, double**, int, int, double, double, double, double, double, 
	    int*, int*, void(*)(double, double**, double**), 
	     void (*)(double**, double**, int, int, double*, 
		     double, double, double**, double*, double*, 
		     void (*) (double, double**, double**)));


void odeint(FILE *evol, double **ystart, int nvar, int mvar, double x1, double x2, 
double eps, double h1, double hmin, int *nok, int *nbad,
void(*derivs)(double, double**, double**), 
	    void (*rkqs)(double**, double**, int, int, double* , double, 
	     double, double**, double*, double*, 
	     void (*) (double, double**, double**)))
{
  int nstp, i, k, l, jc,jc2,oc=0,tobemodif=0, notyet=1, outnum=0;
  double xsav, hnext, hdid, h;
  double **yscal,**y,**dydx,a,b,sum;

	printf("\n Integration started...\n");
	
	yscal=matrix(0,nvar-1,0,mvar-1);
	y=matrix(0,nvar-1,0,mvar-1);
	dydx=matrix(0,nvar-1,0,mvar-1);
		
		x=x1;
		h=SIGN(h1, x2-x1);
		for(i=0; i<nvar; i++)
		  for(k=0; k<mvar; k++)
		    y[i][k]=ystart[i][k];

    xsav = x1 - dxsav;
		for (nstp=1; nstp<MAXSTP; nstp++)
			{
			 
			(*derivs)(x,y,dydx);

			for (i=0; i<nvar; i++)
			  for(k=0; k<mvar; k++)
			    {
			      yscal[i][k]=0.001+h*dydx[i][k];
			    }

			if (fabs(x-xsav)>=fabs(dxsav))
				{
				  sum=0.0;

          writetofile(y, x);
          

				  fprintf(out1, "%g %g %g %g %g %g %g %g %g %g\n",x, y[0][0], y[1][0], y[2][0], y[3][0], y[0][1], y[1][1], y[2][1], y[3][1],energy(y));

				  printf("\n in odeint t=%g,   h=%g,  wall:%d  touch:%d  energy:%g\n",x,h,inthewall, touching, energy(y));
          printf("up=%g down=%g right=%g left=%g \n", wallFup, wallFdown, wallFright, wallFleft);
          

				 fflush(out1);
				  
				xsav=x;
				}

			
			if ((x+h-x2)>0.0) h=x2-x;
		
			(*rkqs)(y, dydx, nvar, mvar, &x, h, eps, yscal, &hdid, &hnext, derivs);	
			
			if ((x-x2)*(x2-x1)>=0.0)
                                       {
		       
					 for (i=0; i<nvar; i++)
					   for (k=0; k<mvar; k++)
					   {
					     ystart[i][k]=y[i][k];
					    
					   }
					
					
					 fflush(out1);
					
					 
					 free_matrix(dydx,0,nvar-1,0,mvar-1); 
					 free_matrix(y,0,nvar-1,0,mvar-1);
					 free_matrix(yscal,0,nvar-1,0,mvar-1);				

					 return;
				       }
				
			if (fabs(hnext)<=hmin){ printf("h: %f hnext: %f hmin: %f\n", h, hnext, hmin);nrerror("stepsize too small in odeint"); }
			h=hnext;
			}
		nrerror("Too many steps in odeint");
}


void getdiff (double t, double **st, double **ds)
{
  double left, right, up, down, dist, close, F, Fx, Fy;
  int i,k;
  
  inthewall=0;
  touching=0;

   for (k=0; k<dim2; k++)
    {
      ds[0][k]=st[2][k]; //dx/dt=vx
      ds[1][k]=st[3][k]; //dy/dt=vy

      ds[2][k]=0.0;      //ax=0
      ds[3][k]=0.0;      //ay=0
    }

   for (k=0; k<dim2; k++) //falak
    {
      left=st[0][k];
      right=(st[0][k]-xlimit);

      down=st[1][k];
      up=(st[1][k]-ylimit);
      
      if (left<0)
	{
    double Fdmp = (st[2][k] < 0) ? -c*st[2][k] : 0;
	  ds[2][k]+=-Dcoeff*left/mass + Fdmp/mass;
	  inthewall++;
	}
      if (right>0)
	{
    double Fdmp = (st[2][k] > 0) ? -c*st[2][k] : 0;
	  ds[2][k]+=-Dcoeff*right/mass + Fdmp / mass;
	  inthewall++;
	}

      if (down<0)
	{
    double Fdmp = (st[3][k] < 0) ? -c*st[3][k] : 0;
	  ds[3][k]+=-Dcoeff*down/mass + Fdmp/mass;
	  inthewall++;
	}

      if (up>0)
	{
    double Fdmp = (st[3][k] > 0) ? -c*st[3][k] : 0;
	  ds[3][k]+=-Dcoeff*up/mass + Fdmp / mass;
	  inthewall++;
	}
    }

    for(int i = 0; i < (int) (xlimit/box_width); i++){
      for(int j = 0; j < (int) (ylimit/box_width); j++){
        grid_count[i][j] = 0;
      }
    }
    int nx = (int)(xlimit / box_width);
    int ny = (int)(ylimit / box_width);
    int cx, cy;
    for (int p = 0; p < dim2; ++p) {
        int cx = (int)(st[0][p] / box_width);
        int cy = (int)(st[1][p] / box_width);
        if (cx < 0) cx = 0;
        if (cx >= nx) cx = nx - 1;
        if (cy < 0) cy = 0;
        if (cy >= ny) cy = ny - 1;

        int pos = grid_count[cx][cy]++;
        if (pos < dim2) grid[cx][cy][pos] = p;
    }

    for(int l = 0; l < (int) (xlimit / box_width); l++){
      for(int m = 0; m < (int) (ylimit / box_width); m++){
        selection_x[0] = l;     selection_y[0] = m;
        selection_x[1] = l;
        if(m > 0) selection_y[1] = m-1; else selection_y[1] = -2;
        if(l < (int) (xlimit / box_width) - 1) selection_x[2] = l+1; else selection_x[2] = -2;
        selection_y[2] = m;
        if(l > 0) selection_x[3] = l-1; else selection_x[3] = -2;
        if(m > 0) selection_y[3] = m-1; else selection_y[3] = -2;
        if(l < (int) (xlimit / box_width) - 1)selection_x[4] = l+1; else selection_x[4] = -2;
        if(m > 0) selection_y[4] = m-1; else selection_y[4] = -2;
        int count = 0;
        count = 0;
      for(int a = 0; a < 5; a++){
      if(a == 0){
       for (int n=0;n < grid_count[l][m]; n++)
        {
          for (int o=0; o<n; o++)
	    {
        i = grid[l][m][n];
        k = grid[l][m][o];
	      dx[i][k]=st[0][i]-st[0][k];
	      dy[i][k]=st[1][i]-st[1][k];
	      float r=sqrt(dx[i][k]*dx[i][k]+dy[i][k]*dy[i][k]);
        
        double sgm6 = sigma * sigma * sigma * sigma * sigma * sigma;
        double sgm12 = sgm6 * sgm6;

            //Lennard-Jones potenciál
            F = 0;
            if(r < sim_radius){  
              double r6 = r*r*r*r*r*r;
              double r13 = r6 * r6 * r;
              double r7 = r6 * r;
              F = 4*epsilon*((6*sgm6/r7) - (12*sgm12/r13));

              if(r<diameter)
	            {
	              touching++;
              }
            Fx = F * dx[i][k] / r / mass;
            Fy = F * dy[i][k] / r / mass;

	          ds[2][k]+=Fx;
	          ds[2][i]+=-Fx;

	          ds[3][k]+=Fy;
	          ds[3][i]+=-Fy;
            }
          count++;
	    }
    }
  }
  else{
    for (int n=0;n < grid_count[l][m]; n++)
        {
          if(selection_x[a] == -2 || selection_y[a] == -2){
            continue;
          }
          for (int o=0;o < grid_count[selection_x[a]][selection_y[a]]; o++)
	    {
        i = grid[l][m][n];
        k = grid[selection_x[a]][selection_y[a]][o];
	      dx[i][k]=st[0][i]-st[0][k];
	      dy[i][k]=st[1][i]-st[1][k];
	      float r=sqrt(dx[i][k]*dx[i][k]+dy[i][k]*dy[i][k]);
        
        double sgm6 = sigma * sigma * sigma * sigma * sigma * sigma;
        double sgm12 = sgm6 * sgm6;

            //Lennard-Jones potenciál
            F = 0;
            if(r < sim_radius){  
              double r6 = r*r*r*r*r*r;
              double r13 = r6 * r6 * r;
              double r7 = r6 * r;
              F = 4*epsilon*((6*sgm6/r7) - (12*sgm12/r13));

              if(r<diameter)
	            {
	              touching++;
              }
            Fx = F * dx[i][k] / r / mass;
            Fy = F * dy[i][k] / r / mass;

	          ds[2][k]+=Fx;
	          ds[2][i]+=-Fx;

	          ds[3][k]+=Fy;
	          ds[3][i]+=-Fy;
            }
          count++;
      }
  }
}

}
}
        }
        for(int i = 0; i < dim2; i++){
          ds[3][i] += -g_accel;
        }
      }


double energy (double **st)
{
  double kinetic=0.0, potential=0.0, left, right,up,down;
  int i,j,k;
  wallFup = 0.0;
  wallFdown = 0.0;
  wallFleft = 0.0;
  wallFright = 0.0;
  int v;

  for (i=0; i<dim2; i++)
    {
      kinetic+=0.5*mass * 1.6605e-27*(st[2][i]*100*st[2][i]*100+st[3][i]*100*st[3][i]*100);
    }
  
    for(int i = 0; i < (int) (xlimit/box_width); i++){
      for(int j = 0; j < (int) (ylimit/box_width); j++){
        grid_count[i][j] = 0;
      }
    }
    // --- constants for cell geometry ---
    const int nx = (int)(xlimit / box_width);
    const int ny = (int)(ylimit / box_width);
    int cx, cy;
    for (int p = 0; p < dim2; ++p) {
        int cx = (int)(st[0][p] / box_width);
        int cy = (int)(st[1][p] / box_width);
        if (cx < 0) cx = 0; if (cx >= nx) cx = nx - 1;
        if (cy < 0) cy = 0; if (cy >= ny) cy = ny - 1;

        int pos = grid_count[cx][cy]++;
        if (pos < dim2) grid[cx][cy][pos] = p;
    }

    for(int l = 0; l < (int) (xlimit / box_width); l++){
      for(int m = 0; m < (int) (ylimit / box_width); m++){
        selection_x[0] = l;     selection_y[0] = m;
        selection_x[1] = l;
        if(m > 0) selection_y[1] = m-1; else selection_y[1] = -2;
        if(l < (int) (xlimit / box_width) - 1) selection_x[2] = l+1; else selection_x[2] = -2;
        selection_y[2] = m;
        if(l > 0) selection_x[3] = l-1; else selection_x[3] = -2;
        if(m > 0) selection_y[3] = m-1; else selection_y[3] = -2;
        if(l < (int) (xlimit / box_width) - 1)selection_x[4] = l+1; else selection_x[4] = -2;
        if(m > 0) selection_y[4] = m-1; else selection_y[4] = -2;
        int count = 0;
        count = 0;
      for(int a = 0; a < 5; a++){
      if(a == 0){
       for (int n=0;n < grid_count[l][m]; n++)
        {
          for (int o=0; o<n; o++)
	    {
        i = grid[l][m][n];
        k = grid[l][m][o];
        if(i >= k){
          continue;
        }
	      dx[i][k]=st[0][i]-st[0][k];
	      dy[i][k]=st[1][i]-st[1][k];
	      float r=sqrt(dx[i][k]*dx[i][k]+dy[i][k]*dy[i][k]);
        
        double sgm6 = sigma * sigma * sigma * sigma * sigma * sigma;
        double sgm12 = sgm6 * sgm6;
        double r6 = r*r*r*r*r*r;

            //Lennard-Jones potenciál
            if(r < sim_radius){  
              double r6 = r*r*r*r*r*r;
              double r13 = r6 * r6 * r;
              double r7 = r6 * r;
              potential += 4*epsilon*1e-21*((sgm12/ r6 / r6) - (sgm6 / r6));
            }
          count++;
	    }
    }
  }
  else{
    for (int n=0;n < grid_count[l][m]; n++)
        {
          if(selection_x[a] == -2 || selection_y[a] == -2){
            continue;
          }
          for (int o=0;o < grid_count[selection_x[a]][selection_y[a]]; o++)
	    {
        i = grid[l][m][n];
        k = grid[selection_x[a]][selection_y[a]][o];
        if(i >= k){
          continue;
        }
	      dx[i][k]=st[0][i]-st[0][k];
	      dy[i][k]=st[1][i]-st[1][k];
	      float r=sqrt(dx[i][k]*dx[i][k]+dy[i][k]*dy[i][k]);
        
        double sgm6 = sigma * sigma * sigma * sigma * sigma * sigma;
        double sgm12 = sgm6 * sgm6;
        double r6 = r*r*r*r*r*r;

            //Lennard-Jones potenciál
            if(r < sim_radius){  
              double r6 = r*r*r*r*r*r;
              double r13 = r6 * r6 * r;
              double r7 = r6 * r;
              potential += 4*epsilon*1e-21*((sgm12/ r6 / r6) - (sgm6 / r6));
            }
          count++;
      }
  }
}}}}

  for (k=0; k<dim2; k++) //falak
    {
      //idealis gaz

      left=st[0][k] / 1e10;
      right=(st[0][k]-xlimit) / 1e10;
      
      down=st[1][k] / 1e10;
      up=(st[1][k]-ylimit) / 1e10;
      
      //height potential
      potential += mass/1.6605e27 * g_accel * 1e14 * st[1][k];
      
      if (left<0){
	      potential+=0.5*Dcoeff*left*left;
        wallFleft += Dcoeff * left * -1;
      }

      if (right>0){
	      potential+=0.5*Dcoeff*right*right;
        wallFright += Dcoeff * right;
      }

      if (down<0){
	      potential+=0.5*Dcoeff*down*down;
        wallFdown += Dcoeff * down * -1;
      }
       
      if (up>0){
	      potential+=0.5*Dcoeff*up*up;
        wallFup += Dcoeff * up;
      }
    }
  
  return (potential+kinetic);
  
}

void writetofile(double **st, double t){

  for (int k=0; k<dim2; k++) //falak
  {
    double left, right,up,down;
    wallFup = 0.0;
    wallFdown = 0.0;
    wallFleft = 0.0;
    wallFright = 0.0;

    left=st[0][k] / 1e10;
    right=(st[0][k]-xlimit) / 1e10;
    
    down=st[1][k] / 1e10;
    up=(st[1][k]-ylimit) / 1e10;

    if (left<0){
      wallFleft += Dcoeff * left * -1;
    }
    if (right>0){
      wallFright += Dcoeff * right;
    }
    if (down<0){
      wallFdown += Dcoeff * down * -1;
    }
    if (up>0){
      wallFup += Dcoeff * up;
    }
  }

  for(int i= 0; i < Nballs; i++){
    fprintf(fp, "%f %f ", st[0][i], st[1][i]);
  }
  fprintf(fp, "%g %g %g %g %f %d %.6f %g\n", wallFup, wallFdown, wallFleft, wallFright, 0.0, Nballs, t, energy(st));
  printf("%f %f\n", st[0][2], st[0][3]);
}

void initialize ()
{

  int i,j,k, inarow;
  double ran, sign, oneball, x_, y;

  oneball=1*sigma;

  inarow=(int)(xlimit/oneball);
  printf("%d\n", inarow);

  k=0;
  for (j=1; j<inarow*2-1; j++)
    for (i=1; i<inarow/2-1; i++)
      {
	if (k<Nballs)
	  {
	    x_=(i)*oneball;
	    y=(j)*oneball;
	
	    distr[0][k]=x_; 
	    distr[1][k]=y;
	    k++;
	  }
      }
  
  for (k=0; k<dim2; k++)
    {
      /*
      //A lenti két sorban a véletlenszerű elhelyezés gondot okozhat túl közel eső golyók esetén. 
      distr[0][k]=ran1(&seed)*xlimit; //x
      distr[1][k]=ran1(&seed)*ylimit; //y
      */
      distr[2][k]=(ran1(&seed) * 2 - 1)*(maxspeed-minspeed)+minspeed; //vx
      ran=ran1(&seed);
      if (ran!=0.0) distr[2][k]*=ran/fabs(ran);    // kell előjel is a sebességnek
      
      distr[3][k]=(ran1(&seed) * 2 - 1)*(maxspeed-minspeed)+minspeed; //vy

      ran=ran1(&seed);
      if (ran!=0.0) distr[3][k]*=ran/fabs(ran); 

      printf("%f %f \n", distr[2][k], distr[3][k]);
    }
}



void getnumbers()
{
  int i, j, k;
  double attenuate;

  g_accel = 9.81 * 1e-14; // angström/ps^2
  e = 2.71828;
  pi = 3.14159;
  c = (fabs(log(wall_e))/sqrt(pi*pi+log(wall_e)*log(wall_e))*2*sqrt(Dcoeff*mass));
  double log_we = log(wall_e);
  double fabs_log_we = fabs(log_we);
  double denominator = sqrt(pi*pi + log_we*log_we);
  double sqrt_term = sqrt(Dcoeff * mass);

  spfb1D = 2 * sigma;
  spaceforaball = spfb1D * spfb1D;
  printf("Maximal number of balls: %g\n", maxball);
  if (Nballs > maxball)
      nrerror("Too many balls...");

  dim2 = Nballs;

  diameter = sigma;         // Å
}


int main(int argc, char *argv[])
{
  int i,j,k;
  double time, Intensity, prevint, cx, prevcx, cy, prevcy, cintx, cinty, prevdens;

  fp = fopen("output.txt", "w");

  double attenuate;

  Dcoeff = atof(argv[1]);
  sim_radius = atof(argv[2]);
  molecule_e = atof(argv[3]);
  wall_e = atof(argv[4]);
  mass = atof(argv[5]);

  radius = atof(argv[6]);
  dim1 = atof(argv[7]);
  seed = atoi(argv[8]);

  // Simulation box dimensions
  xlimit = atof(argv[9]);
  ylimit = atof(argv[10]);
  box_width = atof(argv[11]);

  spfb1D = 2 * sigma;
  spaceforaball = spfb1D * spfb1D;
  maxball = atof(argv[12]);
  Nballs = atoi(argv[13]);
  
  maxspeed = atof(argv[14]);       // angström/ps
  minspeed = atof(argv[15]);       // angström/ps

  start = atof(argv[16]);
  stop = atof(argv[17]);             // ps total simulation time
  dxsav = atof(argv[18]);
  h = atof(argv[19]);              // integration step (ps)
  eps = atof(argv[20]);            // precision
  dim2 = Nballs;

  // Lennard–Jones potential parameters for water
  epsilon = atof(argv[21]); //*1e-21 J
  sigma   = atof(argv[22]); //angström
  diameter = sigma;         // angström
  
  getnumbers ();

  distr=matrix(0,dim1-1,0,dim2-1); //distr[0][k]: a k-adik atom x koordonátája
                                   //distr[1][k]: a k-adik atom y koordonátája
                                   //distr[2][k]: a k-adik atom vx sebessége
                                   //distr[3][k]: a k-adik atom vy sebessége
  dx=matrix(0,dim2-1,0,dim2-1);
  dy=matrix(0,dim2-1,0,dim2-1);
  pairdist=matrix(0,dim2-1,0,dim2-1);
  
  initialize();

 
  odeint(out1, distr, dim1, dim2,start,stop, eps, h,1e-100, nok, nbad, getdiff, rkqs);
       
  free_matrix(distr,0,dim1-1,0,dim2-1);
  free_matrix(dx,0,dim2-1,0,dim2-1);
  free_matrix(dy,0,dim2-1,0,dim2-1);
  free_matrix(pairdist,0,dim2-1,0,dim2-1);
  
  fflush(out1);
  fclose(out1);
  fclose(fp);
}
