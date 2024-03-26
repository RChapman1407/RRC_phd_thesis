// Ordered Line Integral Method - Midpoint (OLIM-M) for computing the quasi-potential
// for dx = b(x) dt + sqrt(epsilon) Sigma dw where Sigma is a VARIABLE matrix
// Copyright:  Daisy Dahiya and Maria Cameron, March 2018
//
// Edited for AMOC model: Ruth Chapman and Peter Ashwin, November 2023.

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#define PI 3.141592653589793
#define e6 0.166666666666666
#define e9 0.111111111111111
#define mabs(a) ((a) >= 0 ? (a) : -(a))
#define sgn(a) ((a) == 0 ? 0 : ((a) > 0  ? 1 : -1 ))
#define max(a,b) ((a) >= (b) ? (a) : (b))
#define min(a,b) ((a) <= (b) ? (a) : (b))
#define INFTY 1.0e+6
#define TOL 1.0e-12
#define NX 1024
#define NY 1024
#define K 22
#define NLTOL 1.0e-6

//three box model parameters

#include "FAMOUSB2xco2_on.h"


struct myvector {
  double x;  
  double y;
};


struct mymatrix {
  double a11;
  double a12;
  double a21;
  double a22;
};

struct mysol {
  double g;
  char c;
};  

struct sol_info {
	char type; // solution type:	1 = 1ptupdate, 2 = 2ptupdate, 0 = initialization, 'n' = never reached
	int ind0;
	int ind1;
	double s;
};

struct interpdata {
	double fx0;
	double fx1;
	double fx2;
	double fx3;
	double fy0;
	double fy1;
	double fy2;
	double fy3;
};

struct myindex {
	int ix;
	int iy;
};	



int main(void);
struct myvector myfield(struct myvector x); /* B */
void param(void);
void initial_curve(void);
void olim(void);
struct mysol triangle_update(int ind,int ind0,int ind1);
double one_pt_update(int ind,int ind0);
// functions for the binary tree			   
void addtree(int ind); /* adds a node to the binary tree
                                 of the "considered" points */ 
void updatetree(int ind); /* updates the binary tree */
void deltree(void); /* deletes the root of the binary tree */
// linear algebra and misc.
struct mymatrix matrix_inverse(struct mymatrix matr);
struct mymatrix matrix_product(struct mymatrix a,struct mymatrix b);
struct mymatrix matrix_transpose(struct mymatrix matr);
struct mymatrix matrix_sum(struct mymatrix A,struct mymatrix B);
struct mymatrix matrix_difference(struct mymatrix A,struct mymatrix B);
struct mymatrix matrix_lin_comb(struct mymatrix A,struct mymatrix B,double a,double b);
struct myvector matr_vec(struct mymatrix matr,struct myvector vec);
double dot_product(struct myvector a,struct myvector b); 
struct myvector vec_lin_comb(struct myvector v1,struct myvector v2,double a,double b);
struct myvector vec_difference(struct myvector v1,struct myvector v2);
struct myvector vec_sum(struct myvector v1,struct myvector v2);
struct myvector scalar_mult(struct myvector v,double a);
double length_vec(struct myvector x);
double quadform(struct mymatrix matr,struct myvector vec);
double Anorm(struct mymatrix matr,struct myvector vec);
struct mymatrix Qpot_matrix(struct mymatrix M);
struct myvector getpoint(int ind);
struct mymatrix Sigma_matrix(struct myvector x);
struct mymatrix Amatrix(struct myvector x);

// struct mymatrix QpotAnisotropicMatrix(void); // returns the matrix Q s.t. U(x) = x^T Q x for anisotropic diffusion
int getneighbors(int ind,int k,int *nei);
struct myvector myrhs(struct myvector p,struct interpdata gdata,int ind);
void ShootMAP(void);
int get_nearest_meshpoint_index(struct myvector p);
struct interpdata get_interpdata(struct myvector p,int ind);
struct myvector myrhs(struct myvector p,struct interpdata gdata,int ind);
double angle(struct myvector v);

double init(struct myvector x);
struct mysol hybrid_nonlin_solver(double a,double b,double u0,double u1,
			struct myvector x0,struct myvector x1,struct myvector b0,struct myvector b1,
			struct mymatrix A0,struct mymatrix A1,struct myvector x);
double myfun(double s,double u0,double u1,struct myvector x0,struct myvector x1,
						struct myvector b0,struct myvector b1,
						struct mymatrix A0,struct mymatrix A1,
						struct myvector dx,struct myvector db,struct mymatrix dA,
						struct myvector x);
double geometric_action_line(struct myvector x0, struct myvector x1);
double init(struct myvector x);


//double exact_solution(struct myvector x);

/***************************************/

const int nx1=NX-1, ny1=NY-1, nxy=NX*NY;
const int Npathmax = 20*max(NX,NY); //NX*NY;
int Npath; // the number of points representing the instanton
int NCURVE;
int count=0; /* # of considered points */
double h,hx,hy;
int ms[NX*NY]; /* 0 = 'Unknown', 1 = 'Considered', 2 = "in Accepted Front", 3 = "Accepted" */
double g[NX*NY]; /* function to be computed */
double vfx[NX*NY]; /* vector field x */
double vfy[NX*NY]; /* vector field y */

int pos[NX*NY]; /* pos(index of mesh pt) = position in binary tree */
int tree[NX*NY]; /* tree(position in the tree) = index of mesh pt */
const int neii[8]={1, NX+1, NX, NX-1, -1, -NX-1, -NX, -NX+1 }; /* neighbor's indices */
//double Uexact[NX*NY];
struct sol_info solinfo[NX*NY];
struct myvector *path;
struct myvector x_ipoint,x_ShootMAP;
struct mycurve *icurve;

// variables for the potential 
double XMIN,XMAX,YMIN,YMAX; // define the computational domain
struct mymatrix Bi; // Bi = J(ipoint); 
double zeta = 0.5; 
struct mymatrix Qexact; // the exact quasi-potential for this linear field is x^T Qexact x

// FILES
const char *f_qpot_name = "Qpot.txt"; // output file with the quasipotential
const char *f_vfx_name = "VFx.txt"; // VF x compt
const char *f_vfy_name = "VFy.txt"; // VF y compt
const char *f_x_name = "Xaxis.txt"; // output file with the norm of VF
const char *f_y_name = "Yaxis.txt"; // output file with the norm of VF
const char *ifname = "Instanton.txt"; // output file with the MAP (the instanton)


/**************************************/
struct myvector myfield(struct myvector x) {
  struct myvector v;
  //double r2,aux,si,aux1,ffun,hfun;
  double q,c1,c2,SIPv,SNv,STv,aq,ypx,ypy,ymx,ymy;
  
  //r2 = x.y*x.y + x.x*x.x;
  // SNv and STv are only variables for three box model
  SNv = x.x/1000;
  STv = x.y/1000;

  q = (Lambda*(Alpha*(TS - To) + Beta*(SNv - SS0)))/(1+Lambda*Alpha*Mu);
  aq = abs(q);

  c1 = (1 + tanh(Xi*q))/2;
  c2 = 1 - c1;

  SIPv = (SALT-(VN*SNv + VT*STv + VS*SS0 + VB*SB0))/VIP;
  
  // q>0 vector field
  ypx = (q*(STv-SNv)+KN*(STv-SNv)-FN*So)*(Y/(VN));
  ypy = (q*(Gamma*SS0+(1-Gamma)*SIPv-STv)+KS*(SS0-STv)+KN*(SNv-STv)-FT*So)*(Y/(VT));

  // q<0 vector field
  ymx = (aq*(SB0-SNv)+KN*(STv-SNv)-FN*So)*(Y/(VN));
  ymy = (aq*(SNv-STv)+KS*(SS0-STv)+KN*(SNv-STv)-FT*So)*(Y/(VT));

  // mix q>0 and q<0 components
  v.x = c1*ypx + c2*ymx;
  v.y = c1*ypy + c2*ymy;

//   if (q>0)
//   { 
// 	v.x=ypx;
//     v.y=ypy;
//   } else {
// 	v.x=ymx;
//     v.y=ymy;
//   }

  // time units decades
  v.x = v.x*1000*10;
  v.y = v.y*1000*10;
  
  return v;
}

/*************************************/

struct mymatrix Sigma_matrix(struct myvector x) {
  struct mymatrix Sigma;
	//double er;
	
	//er = 1.0/length_vec(x);
	Sigma.a11 = 1; Sigma.a12 = 0;
	Sigma.a21 = 0; Sigma.a22 = 1;
	//Sigma.a11 = 0.2372*0.00001; Sigma.a12 = 0.0021*0.00001;
	//Sigma.a21 = 0; Sigma.a22 = 0.2897*0.00001;
	
	return Sigma;
}


/***********************************/

struct mymatrix Amatrix(struct myvector x) {
	struct mymatrix Sigma;
	
	Sigma = Sigma_matrix(x);

	return matrix_inverse(matrix_product(Sigma,matrix_transpose(Sigma)));
}


/*************************************/

//double exact_solution(struct myvector x) {
//double r2, ue; 

//	r2 = x.x*x.x + x.y*x.y;
//	ue = r2*(-1.0 + r2/18.0) + 4.5 + 2.0*(1.0 - x.x/sqrt(r2));
	
//	return ue;
//}

/*************************************/

void param() {
  int i,j,ind; 
struct mymatrix Sigma,S,Q;

  printf("in param()\n");

  x_ipoint.x= XINIT; 
  x_ipoint.y= YINIT; 
  XMIN = XL; XMAX = XH;
  YMIN = YL; YMAX = YH;
  Bi.a11 = -8.6012e-6; Bi.a12 = 8.6012e-6; Bi.a21 = 3.2521e-11; Bi.a22 = -5.8468e-6;
  //Bi.a11 = 1.; Bi.a12 = 0.; Bi.a21 = 0.; Bi.a22 = 1.;
  x_ShootMAP.x = Xshoot;
  x_ShootMAP.y = Yshoot;

  hx = (XMAX - XMIN)/(NX-1);
  hy = (YMAX - YMIN)/(NY-1);
  h = sqrt(hx*hx + hy*hy);
  for( i=0; i<NX; i++ ) {
	for( j=0; j<NY; j++ ) {
      ind = i + NX*j;
	  ms[ind] = 0;
	  g[ind] = INFTY;
	  solinfo[ind].type = 'n';
	  //Uexact[ind] = exact_solution(getpoint(ind));
    }
  }
}

/************************************/

void ipoint() {
  int i,j,ind,ind0,m,n;
  double gtemp;
  const int isur[4]={0, 1, NX+1, NX};
  struct myvector x,l,b;
  struct mymatrix Q,S,Sigma;
	
  i = floor((x_ipoint.x - XMIN)/hx);
  j = floor((x_ipoint.y - YMIN)/hy);
  ind0 = i+j*NX;
  for( m=0; m<4; m++ ) {
    ind = ind0+isur[m];
	x = vec_difference(getpoint(ind),x_ipoint);	
    Sigma = Sigma_matrix(x_ipoint);
    S = matrix_inverse(Sigma);
    Q = Qpot_matrix(matrix_product(S,matrix_product(Bi,Sigma)));
    Q = matrix_product(matrix_transpose(S),matrix_product(Q,S));
    gtemp = quadform(Q,x);
	g[ind] = gtemp;
	ms[ind] = 1;
	addtree(ind);
	solinfo[ind].type = 0;
  }
}

/*************************************/

struct mymatrix Qpot_matrix(struct mymatrix M) {
	double A,B,C,a,b,c,d,aux1,aux2,aux;
	struct mymatrix Umatr;

	a= M.a11; b = M.a12; c = M.a21; d = M.a22;
	aux1 = c - b;
	aux2 = a + d;
	aux = aux1*aux1 + aux2*aux2;
	aux1 *= aux2/aux;
	aux2 *= aux2/aux; 
	Umatr.a11 = -(a*aux2 + c*aux1);
	Umatr.a12 = -(b*aux2 + d*aux1);
	Umatr.a21 = Umatr.a12;
	Umatr.a22 = -(d*aux2 - b*aux1);

	return Umatr;

}


/**********************************************/

void nvf_compute(void) {
  int i,j,ind;
  struct myvector xx,b;

  ind=0;
  for( j=0; j<NY; j++ ) {
    for( i=0; i<NX; i++ ) {
      ind = i + NX*j;
      xx=getpoint(ind);
      b=myfield(xx);
	  // put vector field into arrays
      vfx[ind]=b.x;
	  vfy[ind]=b.y;
    }
  }
}


/*** ordered line integral method ***/

void olim(void) {
  int i,j,k,m,ind,ind0,ind1,indupdate,imin,n;
  int mycount = 0; /* counter for accepted points */
  double gmin,gtemp,gold;
  int Naf, AFneib[8], Nc, NCneib[8]; /* accepted front neighbors and new considered of the newly accepted point */
  struct mysol sol; 
  struct myvector vec,b,b0,b1,v0,v1,vnew,xm0,xm1; 
  struct mymatrix A0,A1;
  double s, s1;
  char pr = 'n'; // a logic variable: if pr == 'y', solution will be printed out in the terminal window
  int nneib, nnaux, fneib, *fnei, *nei, *naux;
  int fneibmax, nneibmax = 8;
  
  k = 2*K+1;
  fneibmax = k*k - 1;
  
  fnei = (int *)malloc(fneibmax*sizeof(int));	// far neighbors
  nei = (int *)malloc(nneibmax*sizeof(int));	// nearest neighbors
  naux = (int *)malloc(nneibmax*sizeof(int));	// nearest neighbors

  
  printf("in olim()\n");

  while( count > 0 ) {
    ind = tree[1];
    vnew = getpoint(ind);
    j = ind/NX;
    i = ind%NX;
    /* x and y of the newly accepted point */
    ms[ind] = 2;
    deltree();
	mycount++;
    if( i<=0 || i>=nx1 || j<=0 || j>=ny1 || g[ind] >= INFTY-1) {
	  printf("The boundary is reached, %d points are Accepted, Umax = %.4e\n",mycount,g[ind]);
	  break; /* quit if we reach the boundary of the computational domain */
	}
	
	if( pr == 'y' ) {
		//printf("%i: %i is accepted (%i,%i), g=%.4e, uexact = %.4e, stype = %i",
		//		mycount,ind,i,j,g[ind],Uexact[ind],solinfo[ind].type);
		if( solinfo[ind].type == 1 ) printf(", ind0 = %i (%i,%i)\n", 
				solinfo[ind].ind0,solinfo[ind].ind0%NX,solinfo[ind].ind0/NX);
		else if( solinfo[ind].type == 2 ) printf(", ind0 = %i, ind1 = %i, s = %.4e\n", 
				solinfo[ind].ind0,solinfo[ind].ind1,solinfo[ind].s);
		else printf("\n");
 	}

    /* Inspect the neighbors of the new Accepted point */
    Naf = 0;
    Nc = 0;
    nneib = getneighbors(ind,1,nei);

    for( k = 0; k<nneib; k++ ) {
		  ind1 = nei[k]; // neighbor of the new Accepted point
		  // update AcceptedFront
		  if( ms[ind1] == 2 ) {
				m = 0;
				nnaux = getneighbors(ind1,1,naux);
				for( n = 0; n < nnaux; n++ ) {
				   ind0 = naux[n];
				   if( ms[ind0] < 2 ) m++;
				}   
				if( m == 0 ) { /* ind1 has no considered neighbors */
				  ms[ind1] = 3;
				}
				else {
				  AFneib[Naf] = ind1;
				  Naf++;
				}  
		  }
		  else if( ms[ind1] == 0 ) { // the neighbor ind1 will be a new Considered point
		    	vec = getpoint(ind);
				NCneib[Nc] = ind1;
				Nc++;
		  }  
	}
	/* update considered points */
	fneib = getneighbors(ind,K,fnei);

	for( k = 0; k < fneib; k++ ) {
		indupdate = fnei[k];						
		if( ms[indupdate] == 1 ) {
			gold = g[indupdate];
			gtemp  = one_pt_update(indupdate,ind);
			vec = getpoint(indupdate);
			if( gtemp < g[indupdate] ) {
				g[indupdate] = gtemp;
				solinfo[indupdate].type = 1;
				solinfo[indupdate].ind0 = ind;
			}
			xm0 = vec_lin_comb(vnew,vec,0.5,0.5);
			b0 = myfield(xm0); 
			A0 = Amatrix(xm0);
			for( m = 0; m < Naf; m++ ) {
			  ind1 = AFneib[m];
			  v1 = getpoint(ind1);
			  xm1 = vec_lin_comb(v1,vec,0.5,0.5);
			  b1 = myfield(xm1); 
			  A1 = Amatrix(xm1);
			  sol = hybrid_nonlin_solver(0.0,1.0,g[ind],g[ind1],vnew,v1,b0,b1,A0,A1,vec);
			  if( sol.c == 'y' ) {
				  s = sol.g;
				  s1 = 1.0 - s;	
				  gtemp = s1*g[ind] + s*g[ind1] + geometric_action_line(vec_lin_comb(vnew,v1,s1,s),vec);	
				  if( gtemp < g[indupdate] ) {
					g[indupdate] = gtemp;
					solinfo[indupdate].type = 2;
					solinfo[indupdate].ind0 = ind;
					solinfo[indupdate].ind1 = ind1;
					solinfo[indupdate].s = s;
				  }	
			  }	
			}
			if( gold > g[indupdate] ) {
			  updatetree(indupdate);
			}   
		} 
	}		

		
     /* Shift Unknown neighbors of the new Accepted point to Considered and compute values at them */ 			  
	 for( m = 0; m < Nc; m++ ) { /* for all points that switched from unknown to considered */
		   indupdate = NCneib[m];
		   i = indupdate%NX;
		   j = indupdate/NX;
		   vec = getpoint(indupdate);
		   gmin = INFTY;
		   imin = ind;
		   fneib = getneighbors(indupdate,K,fnei);
		   	
		   for( k = 0; k < fneib; k++ ) {
			 ind0 = fnei[k];
			 /* compute tentative values using poins of the accepted front or close accepted poonts */
			 if( ms[ind0] == 2 ) {
				 gtemp = one_pt_update(indupdate,ind0);
				 if( gtemp < gmin ) {
					 gmin = gtemp;
					 imin = ind0;
				 }
			 }
		   }
		   ind0 = imin;	 
		   g[indupdate] = gmin;
		   solinfo[indupdate].type = 1;
		   solinfo[indupdate].ind0 = ind0;
		   v0 = getpoint(ind0);
		   xm0 = vec_lin_comb(v0,vec,0.5,0.5);
		   b0 = myfield(xm0); 
		   A0 = Amatrix(xm0);	
		   
		   nneib = getneighbors(ind0,1,nei);
		   for( k = 0; k < nneib; k++ ) {
				 ind1 = nei[k];
				 if( ms[ind1] == 2 ) {
				     v1 = getpoint(ind1);
				     xm1 = vec_lin_comb(v1,vec,0.5,0.5);
				     b1 = myfield(xm1);
				     A1 = Amatrix(xm1);				     
				     sol = hybrid_nonlin_solver(0.0,1.0,g[ind0],g[ind1],v0,v1,b0,b1,A0,A1,vec);
					  if( sol.c == 'y' ) {
						  s = sol.g;
						  s1 = 1.0 - s;
						  gtemp = s1*g[ind0] + s*g[ind1] + geometric_action_line(vec_lin_comb(v0,v1,s1,s),vec);
						  if( gtemp < g[indupdate] ) {
						  	g[indupdate] = gtemp;
							solinfo[indupdate].type = 2;
							solinfo[indupdate].ind0 = ind0;
							solinfo[indupdate].ind1 = ind1;
							solinfo[indupdate].s = s;
						  }	
					  }				  
				  } /* end if( ms[ind1] == 2 ) */
		  }	      
		  addtree(indupdate);
		  ms[indupdate] = 1;

	} /* end for( m = 0; m < Nc; m++ ) */

  } /* end while ( count > 0 ) */
}

  
/*********************************************/

double one_pt_update(int ind,int ind0) {
  struct myvector x,x0;
  double gtemp;
  
  x = getpoint(ind);
  x0 = getpoint(ind0);
  gtemp = g[ind0] + geometric_action_line(x0,x);

  return gtemp;
}
/*-------------*/

double geometric_action_line(struct myvector x0, struct myvector x1) {
  struct myvector l,b,xm;
  struct mymatrix A;
  
  l = vec_difference(x1,x0);
  xm = vec_lin_comb(x0,x1,0.5,0.5);
  b = myfield(xm);
  A = Amatrix(xm);

  return Anorm(A,b)*Anorm(A,l) - dot_product(b,matr_vec(A,l));
}


/*-------------*/  
  
struct myvector getpoint(int ind) {
  struct myvector x0,x1,l;
  int i,ind0,ind1;
  double s;
  
  if( ind < nxy ) {
	  l.x = hx*(ind%NX) + XMIN;
	  l.y = hy*(ind/NX) + YMIN;
  }
  else {
  	printf("Error in getpoint: ind = %i > NXY\n",ind);
  	exit(1);
  }	
  return l;
}

/*-------------*/  

int getneighbors(int ind,int k,int *nei) {
	int jx,jy,k2 = k*k;
	int i,m,nneib = 0;
	int nex,ney; 
    double mbound;
  	
  	jx = ind%NX;
  	jy = ind/NX;
	for( i = -k; i <= k; i++ ) {
		nex = jx + i;
		mbound = floor(sqrt(k2 - i*i + 1.0) + TOL);
		if( nex >= 0 && nex <= nx1 ) {
			for( m = -mbound; m <= mbound; m++ ) {
				ney = jy + m;
				if( ney >= 0 && ney <= ny1 ) {
			    	if( m != 0 || i != 0 ) {
						nei[nneib] = ney*NX + nex;				
						nneib++;
					}
				}
			}
		}
	}
	return nneib;
}			
/***** N o n l i n e a r   1D    s o l v e r *****/


struct mysol hybrid_nonlin_solver(double a,double b,double u0,double u1,
			struct myvector x0,struct myvector x1,struct myvector b0,struct myvector b1,
			struct mymatrix A0,struct mymatrix A1,struct myvector x) {
	double c, fa, fb, fc, d, fd, dd, df, dm, ds, t;
	struct myvector dx,db;
	struct mymatrix dA;
	struct mysol sol;
	double NONLIN_TOL = NLTOL;
	int iter = 0, itermax = 100;

	dx = vec_difference(x0,x1);
	db = vec_difference(b1,b0);
	dA = matrix_difference(A1,A0);

	c = a;
	fa = myfun(a,u0,u1,x0,x1,b0,b1,A0,A1,dx,db,dA,x); 
	fb = myfun(b,u0,u1,x0,x1,b0,b1,A0,A1,dx,db,dA,x);  
	fc = fa;
	if( (fa > 0 && fb > 0 ) || (fa < 0 && fb < 0) ) {
//	 root is not bracketed 
		sol.c = 'n';
		sol.g = INFTY;
		return sol;
	}
	while( iter < itermax ) {
		if( fabs(fc) < fabs(fb) ) { 
			t = c; c = b; b = t;
			t = fc; fc = fb; fb = t;
			a = c; fa = fc;
		}		
		if( fabs(b - c) < NONLIN_TOL ) break;		
		dm = 0.5*(c - b);
		df = fa - fb;

		if( fabs(df) < NONLIN_TOL ) ds = dm;
		else ds = -fb*(a - b)/df;
		if( (ds > 0 && dm < 0) || (ds < 0 && dm > 0) || (fabs(ds) > fabs(dm)) ) dd = dm;
		else dd = ds;
		
		if( fabs(dd) < NONLIN_TOL ) dd = 0.5*sgn(dm)*NONLIN_TOL;

		d = b + dd;
		fd = myfun(d,u0,u1,x0,x1,b0,b1,A0,A1,dx,db,dA,x); 
		if( fabs(fd) < NONLIN_TOL ) {
			b = d;
			break;
		}
		a = b; b = d; fa = fb; fb = fd;
		if( (fb > 0 && fc > 0) || (fb < 0 && fc < 0)) {
			c = a; fc = fa;
		}
		iter++;
	}
	sol.c = 'y';
	sol.g = b;
	
	return sol;
}


/*--------------*/

double myfun(double s,double u0,double u1,struct myvector x0,struct myvector x1,
						struct myvector b0,struct myvector b1,
						struct mymatrix A0,struct mymatrix A1,
						struct myvector dx,struct myvector db,struct mymatrix dA,
						struct myvector x) {
	double s1 = 1.0 - s,ls,lbs,beta;
	struct myvector xs,xmxs,bs,Ad,Ab;
	struct mymatrix As;
	
	xs = vec_lin_comb(x0,x1,s1,s);
	xmxs = vec_difference(x,xs);
	bs = vec_lin_comb(b0,b1,s1,s);
	As = matrix_lin_comb(A0,A1,s1,s);
	ls = Anorm(As,xmxs);
	lbs = Anorm(As,bs);
	beta = lbs/ls;
	Ad = matr_vec(As,xmxs);
	Ab = matr_vec(As,bs);
	
	if( lbs > TOL ) {	
	return u1 - u0 + (dot_product(dx,Ad) + 0.5*quadform(dA,xmxs))*beta + 
					(dot_product(db,Ab) + 0.5*quadform(dA,bs))/beta 
			         - dot_product(dx,Ab) - dot_product(db,Ad) 
			         - dot_product(xmxs,matr_vec(dA,bs));
	}
	else {
	return u1 - u0 + (dot_product(dx,Ad) + 0.5*quadform(dA,xmxs))*beta + 
			         - dot_product(dx,Ab) - dot_product(db,Ad) 
			         - dot_product(xmxs,matr_vec(dA,bs));
	}
	
}
			



/**********************************************/
/*** linear algebra ***/

double angle(struct myvector v) {
  double ang,lv;
  
  lv = length_vec(v);
  if( lv < 1e-14 ) ang = 0.0;
  else if( v.y >= 0.0 ) ang = acos(v.x/lv);
  else ang = 2.0*PI - acos(v.x/lv);
  return ang;
} 

double length_vec(struct myvector x) {
  return sqrt(x.x*x.x + x.y*x.y);
}

struct mymatrix matrix_inverse(struct mymatrix matr) {
  struct mymatrix mi;
  double rdet;

  rdet=1.0/(matr.a11*matr.a22-matr.a12*matr.a21);
  mi.a11=matr.a22*rdet;
  mi.a12=-matr.a12*rdet;
  mi.a21=-matr.a21*rdet;
  mi.a22=matr.a11*rdet;
  return mi;
}

struct mymatrix matrix_product(struct mymatrix a,struct mymatrix b) {
  struct mymatrix c;

  c.a11=a.a11*b.a11+a.a12*b.a21;
  c.a12=a.a11*b.a12+a.a12*b.a22;
  c.a21=a.a21*b.a11+a.a22*b.a21;
  c.a22=a.a21*b.a12+a.a22*b.a22;
  return c;
}

struct mymatrix matrix_transpose(struct mymatrix matr) {
  struct mymatrix mt;

  mt.a11=matr.a11;
  mt.a22=matr.a22;
  mt.a12=matr.a21;
  mt.a21=matr.a12;
  return mt;
}

struct mymatrix matrix_sum(struct mymatrix A,struct mymatrix B) {
  struct mymatrix C;

  C.a11 = A.a11 + B.a11;
  C.a22 = A.a22 + B.a22;
  C.a12 = A.a12 + B.a12;
  C.a21 = A.a21 + B.a21;
  return C;
}

struct mymatrix matrix_difference(struct mymatrix A,struct mymatrix B) {
  struct mymatrix C;

  C.a11 = A.a11 - B.a11;
  C.a22 = A.a22 - B.a22;
  C.a12 = A.a12 - B.a12;
  C.a21 = A.a21 - B.a21;
  return C;
}

struct mymatrix matrix_lin_comb(struct mymatrix A,struct mymatrix B,double a,double b) {
  struct mymatrix C;

  C.a11 = a*A.a11 + b*B.a11;
  C.a22 = a*A.a22 + b*B.a22;
  C.a12 = a*A.a12 + b*B.a12;
  C.a21 = a*A.a21 + b*B.a21;
  return C;
}


struct myvector matr_vec(struct mymatrix matr,struct myvector vec) {
  struct myvector c;

  c.x=matr.a11*vec.x+matr.a12*vec.y;
  c.y=matr.a21*vec.x+matr.a22*vec.y;
  return c;
}

double dot_product(struct myvector a,struct myvector b) {
   return a.x*b.x+a.y*b.y;
}

struct myvector vec_lin_comb(struct myvector v1,struct myvector v2,double a,double b) {
	struct myvector v;
	
	v.x = a*v1.x + b*v2.x;
	v.y = a*v1.y + b*v2.y;
	
	return v;
}
	
struct myvector vec_difference(struct myvector v1,struct myvector v2) {
	struct myvector v;
	
	v.x = v1.x - v2.x;
	v.y = v1.y - v2.y;
	
	return v;
}
	
struct myvector vec_sum(struct myvector v1,struct myvector v2) {
	struct myvector v;
	
	v.x = v1.x + v2.x;
	v.y = v1.y + v2.y;
	
	return v;
}

struct myvector scalar_mult(struct myvector v,double a) {
	struct myvector w;

	w.x = a*v.x;
	w.y = a*v.y;
	
	return w;
}

double quadform(struct mymatrix matr,struct myvector v) {

	return dot_product(v,matr_vec(matr,v));
}

double Anorm(struct mymatrix matr,struct myvector v) {

	return sqrt(quadform(matr,v)); 
}

/****************************************/
/**************************************************************/
/************ FUNCTIONS RELATED TO THE BINARY TREE ***************/

void addtree(int ind) {
  int loc, ptemp;
  int indp, indc;
  char ch;

  count++;
  tree[count]=ind;
  pos[ind]=count;
  if( count > 1 ) {
    loc=count;
    indc=tree[loc];
    indp=tree[loc/2];
    ch=( g[indc] < g[indp] ) ? 'y' : 'n';
    while( ch == 'y' ) {
      ptemp=pos[indc];
      pos[indc]=pos[indp];
      tree[loc/2]=indc;
      pos[indp]=ptemp;
      tree[loc]=indp;
      loc=loc/2;
      if( loc > 1 ) {
        indc=tree[loc];
        indp=tree[loc/2];
        ch=( g[indc] < g[indp] ) ? 'y' : 'n';
      }
      else ch='n';
    }
  }
}

/*------------------------------------------------------------------*/

void updatetree(int ind) {
  int loc, lcc;
  double g0;

  g0=g[ind];
  loc=pos[ind];
  while( loc > 1 && g0 < g[tree[loc/2]] ) {
    tree[loc]=tree[loc/2];
    pos[tree[loc]]=loc;
    loc=loc/2;
    tree[loc]=ind;
    pos[tree[loc]]=loc;
  }  
  lcc=count;
  while( (loc*2 <= count && g0 > g[tree[loc*2]]) || (loc*2+1 <= count && g0 > g[tree[loc*2+1]]) )  {
    lcc=( loc*2+1 <=count && g[tree[loc*2+1]] < g[tree[loc*2]] ) ? loc*2+1 : loc*2;
    tree[loc]=tree[lcc];
    pos[tree[loc]]=loc;
    loc=lcc;
    tree[loc]=ind; 
    pos[tree[loc]]=loc;
  }
}

/*---------------------------------------------------------------------*/


/* deletes root of the binary tree */
void deltree() {
  int loc, ptemp, ind, lcc, ic, ic1, ic2, mind;
  char chd, ch='n';

  mind=tree[1];
  pos[tree[1]]=0;
  tree[1]=tree[count];
  pos[tree[1]]=1;
  count--;
  loc=1;
  ind=tree[1];
  lcc=2*loc;
  if( lcc < count )  {
    ic1=tree[lcc];
    ic2=tree[lcc+1];
    if( (g[ind]) > (g[ic1]) || (g[ind]) > (g[ic2]) ) {
      if( (g[ic1]) <= (g[ic2]) )  {
        chd='l';
	    ic=ic1;
      }
      else {
        chd='r';
	    ic=ic2;
	    lcc++;
      }
    }
    else chd='n';
  }
  else if( lcc == count ) {
    ic=tree[lcc];
    if( (g[ind]) > (g[ic]) ) {chd='l'; if(ch=='y') printf("left\n");}
    else chd='n';
  }
  else chd='n';
  while( chd != 'n' ) {    
    ptemp=pos[ind];
    pos[ind]=pos[ic];
    tree[loc]=ic;
    pos[ic]=ptemp;
    tree[lcc]=ind;
    loc=lcc;
    lcc=2*loc;
    if( lcc < count )  {
      ic1=tree[lcc];
      ic2=tree[lcc+1];
      if( (g[ind]) > (g[ic1]) || (g[ind]) > (g[ic2]) ) {
        if( (g[ic1]) <= (g[ic2]) )  {
          chd='l';
	      ic=ic1;
        }
        else {
          chd='r';
	      ic=ic2;
	      lcc++;
        }
      }
      else chd='n';
    }
    else if( lcc == count ) {
      ic=tree[lcc];
      if(ch=='y') printf("child: loc(%i)=%i, t1=%.12e\n",ic1,lcc,g[ic1]);
      if( (g[ind]) > (g[ic]) ) { chd='l';if(ch=='y') printf("left\n");}
      else chd='n';
    }
    else chd='n';
  } /* end while( chd != 'n' ) */
}


/********************************************************/	

void ShootMAP() {
// shoots a MAP  
	int j,indprev = -1,ind;
	struct myvector y,k1,k2,k3,k4;
	struct interpdata gdata;
	//double dt = 0.1*h,dt2;
    double dt = 0.05*h,dt2;
	
	path = (struct myvector *)malloc(Npathmax*sizeof(struct myvector));
	
	dt2 = 0.5*dt;
	path[0] = x_ShootMAP;
	path[1].x = path[0].x - 0.5*hx;
	path[1].y = path[0].y + hy;
	j = 1;
	while( length_vec(vec_difference(x_ipoint,path[j])) > 3*h && j < Npathmax-1 ) {
		y = path[j];
		ind = get_nearest_meshpoint_index(y);
		if( ind != indprev ) gdata = get_interpdata(y,ind);
		k1 = myrhs(y,gdata,ind);
		indprev = ind;
		
		y = vec_sum(path[j],scalar_mult(k1,dt2));
		ind = get_nearest_meshpoint_index(y);
		if( ind != indprev ) gdata = get_interpdata(y,ind);
		k2 = myrhs(y,gdata,ind);
		indprev = ind;
		
		y = vec_sum(path[j],scalar_mult(k2,dt2));
		ind = get_nearest_meshpoint_index(y);
		if( ind != indprev ) gdata = get_interpdata(y,ind);
		k3 = myrhs(y,gdata,ind);
		indprev = ind;
		
		y = vec_sum(path[j],scalar_mult(k3,dt));
		ind = get_nearest_meshpoint_index(y);
		if( ind != indprev ) gdata = get_interpdata(y,ind);
		k4 = myrhs(y,gdata,ind);
		indprev = ind;
	
		path[j+1] = vec_sum(path[j],scalar_mult(vec_sum(vec_sum(k1,scalar_mult(vec_sum(k2,k3),2.0)),k4),dt*e6));
		j++;
 		//printf("j = %i, ind %i path: %.4e, %.4e\n",j,ind,path[j].x,path[j].y);

		if( j >= Npathmax ) {
			printf("Failed to compute MAP: the maximal allowed path length is reached: j = %i\n",j);
	//		return;
		}	
	}
	path[j] = x_ipoint;
    printf("j = %i\n, ind %i",j,ind);
    printf(" last but one point: %.4e, %.4e\n",path[j-1].x,path[j-1].y);
    printf(" last point: %.4e, %.4e\n",path[j].x,path[j].y);
	j++;
	Npath = j;
	printf("Npath = %i\n",Npath);
}


/*************************************************/
int get_nearest_meshpoint_index(struct myvector p) {
	int ix,iy;
	
	ix = floor((p.x - XMIN)/hx);
	iy = floor((p.y - YMIN)/hy);
	return ix + iy*NX;
}

/*************************************************/

struct interpdata get_interpdata(struct myvector p,int ind) {
	struct interpdata gdata;
		
	gdata.fx0 = 0.5*(g[ind + 1] - g[ind - 1])/hx; // dg/dx at ind
	gdata.fx1 = 0.5*(g[ind + 2] - g[ind])/hx; // dg/dx at ind + 1
	gdata.fx2 = 0.5*(g[ind + 1 + NX] - g[ind - 1 + NX])/hx; // dg/dx ind + NX
	gdata.fx3 = 0.5*(g[ind + 2 + NX] - g[ind + NX])/hx; // dg/dx ind + NX + 1
	
	gdata.fy0 = 0.5*(g[ind + NX] - g[ind - NX])/hy; // dg/dy at ind
	gdata.fy1 = 0.5*(g[ind + NX + 1] - g[ind - NX + 1])/hy; //dg/dy at ind + 1
	gdata.fy2 = 0.5*(g[ind + 2*NX] - g[ind])/hy; // dg/dy at ind + NX
	gdata.fy3 = 0.5*(g[ind + 2*NX + 1] - g[ind + 1])/hy; // dg/dy at ind + NX + 1

	return gdata;
}	


/*************************************************/
struct myvector myrhs(struct myvector p,struct interpdata gdata,int ind) {
	double ax,ay,ax1,ay1,gnorm;
	struct mymatrix S;
	int i,j;
	
	struct myvector gu,b;
	
	i = ind%NX;
	j = ind/NX;
	if( i > NX - 3 || j > NY - 3) {
		printf("In myrhs: y = (%.4e,%.4e), (i,j) = (%i,%i)\n",p.x,p.y,i,j);
		printf("Too close to the boundary\n");
		exit(1);
	}	

	ax = (p.x - XMIN - i*hx)/hx; ax1 = 1.0 - ax;
	ay = (p.y - YMIN - j*hy)/hy; ay1 = 1.0 - ay;
	
	gu.x = -((gdata.fx0*ax1 + gdata.fx1*ax)*ay1 + (gdata.fx2*ax1 + gdata.fx3*ax)*ay);
	gu.y = -((gdata.fy0*ax1 + gdata.fy1*ax)*ay1 + (gdata.fy2*ax1 + gdata.fy3*ax)*ay);
	b = myfield(p);	
	S = Sigma_matrix(p);
	S = matrix_product(S,matrix_transpose(S));
	gu = vec_difference(matr_vec(S,gu),b);
	
	gnorm = length_vec(gu);
	gu = scalar_mult(gu,1.0/gnorm);

	
	return gu;
}

/********************************************************/

	    
/*** main ***/

 int main() {
  int i,j,ind,k,si; 
  clock_t CPUbegin;
  double cpu,dd,errmax = 0.0,erms = 0.0;
  double xx,yy;
  FILE *fid,*fid2,*fid3;
  
  param();
  ipoint();
  nvf_compute();

  printf("Total salt = %f\n",SALT);

  CPUbegin=clock();
  olim();
  cpu = (clock()-CPUbegin)/((double)CLOCKS_PER_SEC);
  printf("cputime of olim() = %g\n",cpu);  
  
  fid = fopen(f_qpot_name,"w");
  fid2 = fopen(f_vfx_name,"w");
  fid3 = fopen(f_vfy_name,"w");
  ind=0;
  k = 0;
  for( j=0; j<NY; j++ ) {
    yy = hy*j + YMIN;
    for( i=0; i<NX; i++ ) {
        xx = hx*i + XMIN;
    	if( ms[ind] < 2 ) {
    		g[ind] = INFTY;
    		si = -1;
    	}	
        fprintf(fid,"%.12e\t",g[ind]);
        fprintf(fid2,"%.12e\t",vfx[ind]);
        fprintf(fid3,"%.12e\t",vfy[ind]);
        ind++;
    }
    fprintf(fid,"\n");
    fprintf(fid2,"\n");
    fprintf(fid3,"\n");
  }
  fclose(fid);
  fclose(fid2);
  fclose(fid3);

  fid = fopen(f_x_name,"w");
  for( j=0; j<NX; j++ ) {
	  xx = hx*j + XMIN;
          fprintf(fid,"%.12e\t",xx);
  }
  fclose(fid);

  fid = fopen(f_y_name,"w");
  for( j=0; j<NY; j++ ) {
	  yy = hy*j + YMIN;
          fprintf(fid,"%.12e\t",yy);
  }
  fclose(fid);


  printf("NX = %i, NY = %i, K = %i\n",NX,NY,K);
  //printf("errmax = %.4e, erms = %.4e\n",errmax,sqrt(erms/k));
  
// /*   compute the instanton 
  	ShootMAP();
  	fid = fopen(ifname,"w");
  	for( j = 0; j< Npath; j++ ) {
  		fprintf(fid,"%.12e\t%.12e\n",path[j].x,path[j].y);
  	}
  	fclose(fid);

  return 0;
}  
