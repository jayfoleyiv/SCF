#include<math.h>
#include<malloc.h>
#include"blas.h"
#include <stdlib.h>
#include <stdio.h>
#include<string.h>
int DIAG_N(int dim, int number, double *mat, double *en, double *wfn);
void Diagonalize(double*M,long int dim, double*eigval,double*eigvec);
void print_matrix( char* desc, int m, int n, double* a, int lda );
void LoopMM(int dim, double *a, char *transa, double *b, char *transb, double *c);
void BuildDensity(int dim, int occ, double *C, double *D);
void UpdateF(int dim, double *D, double *Hcore, double *EI, double *Fnew);
double E_Total(int dim, double *D, double *HCore, double *F, double Enuc);
int FourDIndx(int i, int j, int k, int l, int dim);
void ReadEI(int dim, FILE *fp, double *EE);
double DensityDiff(int dim, double *D, double *Dnew);
int main() {
  //enuc.dat  eri.dat  s.dat  t.dat  v.dat
  int dim, i, j, ij, kl;
  FILE *Nfp, *Sfp, *Tfp, *NEfp, *EEfp;
  double Enuc, *S, *T, *NE, *EE, *Hcore, *temp, val;
  double *Svals, *Svecs, *SqrtS, *SqrtSvals;
  double *F, *Fnew, *Cp, *C, *D, *Dnew, *eps;
  double sum, ESCF, ESCF_i, deltaE, deltaDD, tolE, tolDD;
  int iter, itermax;

  tolE = 1e-7;
  tolDD = 1e-7;
  dim = 7;
  itermax = 100;

  Nfp  = fopen("enuc.dat","r");
  Sfp  = fopen("s.dat","r");
  Tfp  = fopen("t.dat","r");
  NEfp = fopen("v.dat","r");
  EEfp = fopen("eri.dat","r");

  NE    = (double *)malloc(dim*dim*sizeof(double));
  S     = (double *)malloc(dim*dim*sizeof(double));
  T     = (double *)malloc(dim*dim*sizeof(double));
  Hcore = (double *)malloc(dim*dim*sizeof(double));
  F     = (double *)malloc(dim*dim*sizeof(double));
  Fnew  = (double *)malloc(dim*dim*sizeof(double));
  Cp    = (double *)malloc(dim*dim*sizeof(double));
  C     = (double *)malloc(dim*dim*sizeof(double));
  D     = (double *)malloc(dim*dim*sizeof(double));
  Dnew  = (double *)malloc(dim*dim*sizeof(double));
  eps   = (double *)malloc(dim*sizeof(double));
  Svals = (double *)malloc(dim*sizeof(double));
  SqrtSvals = (double *)malloc(dim*dim*sizeof(double));
  Svecs = (double *)malloc(dim*dim*sizeof(double));
  SqrtS = (double *)malloc(dim*dim*sizeof(double));
  temp  = (double *)malloc(dim*dim*sizeof(double));
  EE    = (double *)malloc(dim*dim*dim*dim*sizeof(double));

  // Read Nuclear Repulsion
  fscanf(Nfp, "%lf",&Enuc);
  
  // Read 2-electron integrals
  ReadEI(dim, EEfp, EE);

  // Read 1-electron matrices
  for (i=0; i<dim; i++) {
    for (j=0; j<=i; j++) {

      fscanf(NEfp,"%i",&ij);
      fscanf(NEfp,"%i",&kl);
      fscanf(NEfp,"%lf",&val);
      NE[i*dim+j] = val;
      NE[j*dim+i] = val;
      
      fscanf(Sfp,"%i",&ij);
      fscanf(Sfp,"%i",&kl);
      fscanf(Sfp,"%lf",&val);
      S[i*dim+j] = val;
      S[j*dim+i] = val;

      fscanf(Tfp,"%i",&ij);
      fscanf(Tfp,"%i",&kl);
      fscanf(Tfp,"%lf",&val);
      T[i*dim+j] = val;
      T[j*dim+i] = val;

      Hcore[i*dim+j] = T[i*dim+j] + NE[i*dim+j];
      Hcore[j*dim+i] = T[j*dim+i] + NE[j*dim+i];

      printf("  %i  %i   %18.14f\n",i+1,j+1,Hcore[i*dim+j]);

    }
  }

  print_matrix("  T  ", dim, dim, T, dim);
  print_matrix("  NE ", dim, dim, NE, dim);
  print_matrix("  S  ", dim, dim, S, dim);
  // Diagonalize overlap 
  DIAG_N(dim, dim, S, Svals, Svecs); 

  for (i=0; i<dim; i++) {
 
    SqrtSvals[i*dim + i] = pow(Svals[i],-1./2);

  }
  // Form S^{-1/2} = L_S s^{-1/2} L_S^t
  LoopMM(dim, SqrtSvals, "n", Svecs, "t", temp);
  LoopMM(dim, Svecs, "n", temp, "n", SqrtS);

  print_matrix( "S^1/2", dim, dim, SqrtS, dim);

  // Form Fock matrix F = S^{-1/2}^t H_core S^{1/2}
  LoopMM(dim, Hcore, "n", SqrtS, "n", temp);
  LoopMM(dim, SqrtS, "t", temp, "n", F);

  print_matrix("  Fock", dim, dim, F, dim);

  //  Get Guess MO matrix from diagonalizing Fock matrix
  // Diag(F) -> MO Coefficients = vecs, MO energies = vals
  DIAG_N(dim, dim, F, eps, Cp);
  print_matrix("  Initial Coefficients", dim, dim, Cp, dim);

  LoopMM(dim, SqrtS, "n", Cp, "n", C);

  BuildDensity(dim, 5, C, D);
  print_matrix("  Coefficients", dim, dim, C, dim);

  print_matrix("  Density Matrix", dim, dim, D, dim);

  ESCF_i = E_Total(dim, D, Hcore, F, Enuc);
  
  printf("  Initial E_SCF is %12.10f\n",ESCF);

  int die = 1;
  iter=0;

  printf("  ITERATION 0:  RHF ENERGY IS %18.14f\n",ESCF_i);
  do {

    // Update Fock matrix
    UpdateF(dim, D, Hcore, EE, Fnew);

    // Form Fock matrix F = S^{-1/2}^t H_core S^{1/2}
    LoopMM(dim, Fnew, "n", SqrtS, "n", temp);
    LoopMM(dim, SqrtS, "t", temp, "n", F);
  
    // Diagonalize new Fock matrix
    DIAG_N(dim, dim, F, eps, Cp);

    // Get new MO coefficients
    LoopMM(dim, SqrtS, "n", Cp, "n", C);

    // Build new Density matrix
    BuildDensity(dim, 5, C, Dnew);
    //print_matrix("  New Density Matrix ", dim, dim, Dnew, dim);

    // Compute new Energy
    ESCF = E_Total(dim, Dnew, Hcore, Fnew, Enuc);

    // Get RMS_D for density matrix, copy new density matrix to D array
    deltaDD = DensityDiff(dim, D, Dnew);

    // get change in energy
    deltaE = ESCF - ESCF_i;
    // call current energy ESCF_i for next iteration
    ESCF_i = ESCF;
    
    if (fabs(deltaE)<tolE && deltaDD<tolDD) die=0;
    else if (iter>itermax) die=0;
    
    iter++;
    printf("  ITERATION %5i:  RHF ENERGY IS %18.14f  DeltaE is %18.14f  DeltaD is %18.14f\n",iter, ESCF, deltaE, deltaDD);

  }while(die);


return 0;
}

double DensityDiff(int dim, double *D, double *Dnew) {
  int m, n;
  double sum;

  sum = 0;

  for (m=0; m<dim; m++) {
    for (n=0; n<dim; n++) {

      sum += (Dnew[m*dim+n]-D[m*dim+n])*(Dnew[m*dim+n]-D[m*dim+n]);
      D[m*dim+n] = Dnew[m*dim+n];

    }
  }   

  return sqrt(sum);

}

void ReadEI(int dim, FILE *fp, double *EE) {
  int i, j, k, l, ij, kl, ijkl;
  double val;

  while(fscanf(fp, "%d %d %d %d %lf",&i,&j,&k,&l,&val) !=EOF) {
    i--;
    j--;
    k--;
    l--;
    // ijkl
    //ij = i*(i+1)/2 + j;
    //kl = k*(k+1)/2 + l;
    //ijkl = ij*(ij+1)/2 + kl;
    EE[FourDIndx(i,j,k,l,dim)] = val;
    EE[FourDIndx(j,i,k,l,dim)] = val;
    EE[FourDIndx(i,j,l,k,dim)] = val;
    EE[FourDIndx(j,i,l,k,dim)] = val;
    EE[FourDIndx(k,l,i,j,dim)] = val;
    EE[FourDIndx(l,k,i,j,dim)] = val;
    EE[FourDIndx(k,l,j,i,dim)] = val;
    EE[FourDIndx(l,k,j,i,dim)] = val;
  }

}

int FourDIndx(int i, int j, int k, int l, int dim) {

  return i*dim*dim*dim+j*dim*dim+k*dim+l;

}


double E_Total(int dim, double *D, double *Hc, double *F, double Enuc) {

  int m, n;  
  double sum;

  /*print_matrix("  D  ", dim, dim, D, dim);
  print_matrix("  Hc ", dim, dim, Hc, dim);
  print_matrix("  F  ", dim, dim, F, dim);
  */
  sum = 0.;
  for (m=0; m<dim; m++) {
    for (n=0; n<dim; n++) {
//      sum += D[m*dim+n]*(Hc[m*dim+n] + F[m*dim+n]);
      sum += D[m*dim+n]*(Hc[m*dim+n] + F[m*dim+n]);
    }
  }
  return sum + Enuc;
}


void UpdateF(int dim, double *D, double *Hcore, double *EI, double *Fnew) {

  int m, n, l, s, mnls, mlns;
  double sum;

  for (m=0; m<dim; m++) {
    for (n=0; n<dim; n++) {

      sum = 0.;
      for (l=0; l<dim; l++) {
        for (s=0; s<dim; s++) {

          mnls = FourDIndx(m, n, l, s, dim);
          mlns = FourDIndx(m, l, n, s, dim); 
          
          sum += D[l*dim+s]*(2*EI[mnls]-EI[mlns]);

        }
      }
      Fnew[m*dim+n] = Hcore[m*dim+n] +  sum;
      //Fnew[m*dim+n] = sum;
    }
  }
}


void BuildDensity(int dim, int occ, double *C, double *D) {
  int i, j, m;
  double sum;

  for (i=0; i<dim; i++) {
    for (j=0; j<dim; j++) {
      sum = 0.;
      for (m=0; m<occ; m++) {
        sum += C[i*dim+m]*C[j*dim+m];
      }
      D[i*dim+j] = sum;
    }
  }
}



void Diagonalize(double*M,long int dim, double*eigval,double*eigvec){
  integer one,info,edim,fdim,*ifail,tdim,i,il,iu,m;
  doublereal zero,vl,vu;
  doublereal tol;
  char N, U;
  doublereal *work;
  integer*iwork;
  edim = 8*dim;
  fdim = 5*dim;
  tdim = 3*dim;
  N    = 'V'; // 'N' for eigenvalues only, 'V' for eigenvectors, too
  U    = 'U';
  one  = dim;   // if N='N', one=1; otherwise, one=dim;
  work  = (doublereal*)malloc(edim*sizeof(doublereal));
  DSPEV(N,U,dim,M,eigval,eigvec,one,work,info);

  //for (i=0; i<dim; i++) printf("  Eig %i  %12.10f\n",i+1,eigval[i]);
  free(work);
}

void print_matrix( char* desc, int m, int n, double* a, int lda ) {
        int i, j;
        printf("\n\n----------------------");
        printf( "\n %s\n", desc );
        for( i = 0; i < m; i++ ) {
                for( j = 0; j < n; j++ ) printf( " %12.9f", a[i+j*lda] );
                printf( "\n" );
        }
        printf("----------------------\n\n");
}

void LoopMM(int dim, double *a, char *transa, double *b, char *transb, double *c) {
  int i, j, k; 
  double sum;

  for (i=0; i<dim; i++) {
    for (j=0; j<dim; j++) {
  
      sum = 0.;
      for (k=0; k<dim; k++) {

        if (strcmp(transa,"t")==0 && strcmp(transb,"t")==0) {
          sum += a[k*dim+i]*b[j*dim+k];  
        }

        else if (strcmp(transa,"n")==0 && strcmp(transb,"t")==0) {
          sum += a[i*dim+k]*b[j*dim+k];
        }
        else if (strcmp(transa,"t")==0 && strcmp(transb,"n")==0) {
          sum += a[k*dim+i]*b[k*dim+j];
        }
        else {
          sum += a[i*dim+k]*b[k*dim+j];
        }
        
      }
      c[i*dim+j] = sum;
    }
  }
}

int DIAG_N(int dim, int number, double *mat, double *en, double *wfn) {
  int i,j,ind, state_max, count;
  double *pacMat, *eigval, *eigvec;

  pacMat = (double *)malloc((dim*(dim+1)/2)*sizeof(double));
  eigval = (double *)malloc(dim*sizeof(double));
  eigvec = (double *)malloc(dim*dim*sizeof(double));

  for (i=0;i<dim;i++) {
    for (j=0;j<dim;j++) {
      if (i<=j) {
        ind =  j*(j+1)/2 + i; // Position(i,j);
        pacMat[ind] = mat[i*dim+j];
      }
    }

  }

  Diagonalize(pacMat,dim,eigval,eigvec);

  count=0;
  for (i=0; i<number; i++) {
    en[i] = eigval[i];
    if (en[i]<=0) count++;
    for (j=0; j<dim; j++) {
      wfn[j*dim+i] = eigvec[i*dim+j];
    }
  }

  return count;
  free(pacMat);
  free(eigval);
  free(eigvec);


}


