#include<math.h>
#include<malloc.h>
#include"memory.h"
#include"blas.h"
#include <stdlib.h>
#include <stdio.h>

int DIAG_N(int dim, int number, double *mat, double *en, double *wfn);
void Diagonalize(double*M,long int dim, double*eigval,double*eigvec);
void print_matrix( char* desc, int m, int n, double* a, int lda );

int main() {

    // integers related to the dimensions of our arrays
    int N, lda; 

    // These will be pointers to our matrices 
    // we will store the matrices as 1-d arrays 
    double *M, *D;

    // this will be a pointer to an array of our eigenvectors
    // we will store this in the same way we store our matrices
    double *eigvec;

    // this will be a 1d array of the eigenvalues for
    // our matrices
    double *eig;

    // scalars for the arguments to the matrix-matrix multiplication function
    double alpha, beta;


    N = 3;
    lda = N;
    alpha = 1.0;
    beta = 0.0;

    //  allocate memory for each array
    //  the matrices have N*N elements
    M      = (double *)malloc(N*N*sizeof(double));
    D      = (double *)malloc(N*N*sizeof(double));
    //  there are N eignevalues
    eig    = (double *)malloc(N*sizeof(double));
    //  there are also N eigenvectors with N elements each, so also N*N
    eigvec = (double *)malloc(N*N*sizeof(double));

    // We will store the following 3x3 matrix as a vector with 9 elements:

    /* 
        1.20000   0.20000  0.11320        [M(0,0)]  [M(0,1)] [M(0,2)]
        0.98100   0.33120  2.20000  <==>  [M(1,0)]  [M(1,1)] [M(1,2)] 
        0.99132   2.20000  3.20000        [M(2,0)]  [M(2,1)] [M(2,2)]  

        will be stored like this in a c program:

         vc(0) vc(1) vc(2) vc(3)   vc(4)  vc(5) vc(6)   vc(7) vc(8)
        [1.2, 0.2, 0.1132, 0.981, 0.3312, 2.2, 0.99132, 2.20, 3.2]

        Note that this is called "Row-major order" and this is
        the convention used in C-programming.
    
        The convention used in Fortran programming is 
        "Column-major ordering", and the same matrix would 
        map to the following vector in Fortran:

        vf[0] 1.2
        vf[1] 0.981
        vf[2] 0.99132
        vf[3] 0.2000
        vf[4] 0.3312
        vf[5] 2.2
        vf[6] 0.1132
        vf[7] 2.2
        vf[8] 3.2

        vf looks like the vectorized version of the transpose of matrix M
        (aka Mt) 

    */
    // M(0,0)
    M[0*N+0] = 1.2;
    // M(0,1) and M(1,0)
    M[0*N+1] = 0.2;
    M[1*N+0] = 0.981;
    /// M(0,2) and M(2,0)
    M[0*N+2] = 0.1132;
    M[2*N+0] = 0.99132;

    // M(1,1)
    M[1*N+1] = 0.3312;
    // M(1,2) and M(2,1) 
    M[1*N+2] = 2.2;
    M[2*N+1] = 2.2;

    // M(2,2)
    M[2*N+2] = 3.2;

    // this will print the matrix in column-major order, reflecting
    // the way the matrix will be interpreted by the
    // matrix library functions, which are written in fortran
    print_matrix("Matrix M",N, N, M, N);

    //  Here we want to take the matrix we formed and
    //  multiply it by its transpose to get a symmetrix product
    //  matrix
    //  this particular F_DGEMM call performs
    //  D = alpha * transpose(M)*M + beta * D, but again,
    //  because of the column-major order of Fortran, a 
    //  c-programmer would interpret this as 
    //  D = alpha * M*transpose(M) + beta * D
    //  (recall also that alpha is a scalar that we have
    //  set equal to 1, beta is a scalar that we have set equal to zero,
    //  'T' in the first argument means take transpose of the first matrix
    //  'N' in the second argument means use the second matrix as is
    F_DGEMM( 'T', 'N', N, N, N, alpha, M, lda, M, lda, beta, D, lda);

    // The product D should be a symmetric matrix... verify that
    // it is actually equal to M*transpose(M) as WE have defined M above
    print_matrix("Matrix M*Mt",N, N, D, N);

    //  Now we will diagonalize D and store its eigenvectors and eigenvalues
    int bound = DIAG_N(N, N, D, eig, eigvec);

    //  Going to print the eigenvectors
    print_matrix("Eigenvectors of (M*Mt)",N,N,eigvec,N);
    printf("\n");

    //  Going to print the eigenvalues
    printf("  eig1 is %f\n",eig[0]);
    printf("  eig2 is %f\n",eig[1]);
    printf("  eig3 is %f\n",eig[2]);

    //  Now perform the matrix multiplication the other way:
    //  D = alpha * M*transpose(M) + beta * D
    //  which in our c-programming convention, will
    //  actually yield D = alpha*transpose(M)*M + beta * D
    F_DGEMM( 'N', 'T', N, N, N, alpha, M, lda, M, lda, beta, D, lda);

    // The product D should be a symmetric matrix again... verify that
    // it is actually equal to transpose(M)*M as WE have defined M above
    print_matrix("Matrix Mt*M",N, N, D, N);

    //  Now we will diagonalize D and store its eigenvectors and eigenvalues
    bound = DIAG_N(N, N, D, eig, eigvec);

    //  Going to print the eigenvectors
    print_matrix("Eigenvectors of (Mt*M)",N,N,eigvec,N);

    //  Going to print the eigenvalues
    printf("  eig1 is %f\n",eig[0]);
    printf("  eig2 is %f\n",eig[1]);
    printf("  eig3 is %f\n",eig[2]);

    return 0;

}

int DIAG_N(int dim, int number, double *mat, double *en, double *wfn) {
  int i,j,ind, state_max, count;
  double *pacMat, *eigval, *eigvec;

  pacMat = VEC_DOUBLE(dim*(dim+1)/2);
  eigval = VEC_DOUBLE(dim);
  eigvec = VEC_DOUBLE(dim*dim);

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
      wfn[i*dim+j] = eigvec[i*dim+j];
    }
  }

  return count;
  free(pacMat);
  free(eigval);
  free(eigvec);


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
                for( j = 0; j < n; j++ ) printf( " %6.2f", a[i+j*lda] );
                printf( "\n" );
        }
        printf("----------------------\n\n");
}
         
