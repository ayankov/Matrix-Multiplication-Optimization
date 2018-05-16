//     *********************************
//     *    Matrix Multiply Project    *
//     *                               *
//     *********************************

//     ** MAIN PROGRAM  **


//     *************************************************
//     ** Any changes you make to this code must      **
//     ** maintain the correctness of the matrix      **
//     ** multiply computed by the original version.
//     **					      **
//     ** An implementation with incorrect results for**
//     ** matrix C earns zero point for this project. **
//     ** 
//     ** To print matrix compile with -DPRINT_MATRIX **
//     ** e.g., g++ -DPRINT_MATRIX                    **
//     **                                             **
//     ** A sample Makefile is provided.
//     ** You may assume m = n = k for your matrices  **
//     *************************************************
//
// Alexei Yankovsky
// Instructor: Dr. Rong Ge
// December 5, 2017
// Project 3 - Speed up Matrix Multiplication



#include <time.h>
#include <stdio.h>
#include <stdlib.h>

#define min(a, b) ((a) < (b)) ? (a) : (b)

double **dmatrix(int nrl,int nrh,int ncl,int nch);
void nerror(char *error_text);

int main(int argc, char** argv)  {
 
      int l,m,n,k;
      int i,j;
      double temp;
      double **A, **B, **C;


     //  ****************************************************
     //  * The following allows matrix parameters to be     *
     //  * entered on the command line to take advantage    *
     //  * of dynamically allocated memory.  You may modify *
     //  * or remove it as you wish.                        *
     //  ****************************************************

     if (argc != 4) {
       nerror("Usage:  <executable> <m-value> <n-value> <k-value>");
     }
      m = atoi(argv[1]);
      n = atoi(argv[2]);
      k = atoi(argv[3]);

      // *********************************************************
      // * Call the dmatrix() subroutine to dynamically allocate *
      // * storage for the matrix sizes specified by m, n, and k *  
      // *********************************************************

      A=dmatrix(0,m-1,0,k-1);
      B=dmatrix(0,k-1,0,n-1);
      C=dmatrix(0,m-1,0,n-1);

      // *********************************************************
      //  * Initialize matrix elements so compiler does not      *
      //  * optimize out                                         *
      // *********************************************************

      for(j=0;j<k;j++) {
        for(i=0;i<m;i++) {
          A[i][j] = i+j+4.0;
        }
      }

      for(j=0;j<n;j++) {
        for(i=0;i<k;i++) {
          B[i][j] = i+j+5.0;
        }
      }

      for(j=0;j<n;j++) {
        for(i=0;i<m;i++) {
          C[i][j] = 0.0;
        }
      }

      // ******************************
      // * Start embedded timing here *
      // ******************************
			clock_t begin = clock();
      // **********************************
      // * Perform simple matrix multiply *
      // **********************************



//loop interchange, loop unroll, indexing optimization, loop fission, blocking

   int b = 48; //blocksize
   int i0,l0,j0;

   for(i0=0;i0<n;i0+=b){   
      for(l0=0;l0<n;l0+=b){
         for(j0=0;j0<n;j0+=b){
            for(i=i0;i<((i0+b)>m?m:(i0+b));i++) {
               double *iRowA;
               iRowA = A[i];
               double *iRowC;
               iRowC = C[i];
               for(l=l0; l<((l0+b)>k?k:(l0+b)); l++) {
                  double *lRowB;
                  lRowB = B[l];
                  double ilA = iRowA[l];
                  int prods[k];
                  for(j=j0; j<((j0+b)>n?n:(j0+b)); j+=4) {
	                  prods[j] = ilA*lRowB[j];
	                  prods[j+1] = ilA*lRowB[j+1];
	                  prods[j+2] = ilA*lRowB[j+2];
	                  prods[j+3] = ilA*lRowB[j+3];
                  }   
                  for(j=j0; j<((j0+b)>n?n:(j0+b)); j+=4) {
	                  iRowC[j] += prods[j];
	                  iRowC[j+1] += prods[j+1];	
	                  iRowC[j+2] += prods[j+2];	
	                  iRowC[j+3] += prods[j+3];		
                  }        
               }
            }
         }               
      }
   }


      // ******************************
      // * Stop embedded timing here  *
      // ******************************
			clock_t end = clock();
			double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
			printf("\nCPU execution time: %f \n", time_spent);


      // **************************************************
      // * Print out a 10 x 10 matrix for testing only    *
      // * Comment out when timing                        *
      // **************************************************

 #ifdef PRINT_MATRIX
      fprintf(stdout, "Here is the matrix A:\n\n");
      for(i=0;i<m;i++) {
        for(j=0;j<k;j++) {
          fprintf(stdout, "%10.2f ",A[i][j]);
        }
        fprintf(stdout, "\n");
      }
      fprintf(stdout, "Here is the matrix B:\n\n");
      for(i=0;i<k;i++) {
        for(j=0;j<n;j++) {
          fprintf(stdout, "%10.2f",B[i][j]);
        }
        fprintf(stdout, "\n");
      }
      fprintf(stdout, "Here is the matrix C:\n\n");
      for(i=0;i<m;i++) {
        for(j=0;j<n;j++) {
          fprintf(stdout, "%10.2f",C[i][j]);
        }
        fprintf(stdout, "\n");
      }
#endif        
        
}
//     **  END MAIN PROGRAM  **

//     ********************************************************
//     *******    BEGIN SUBROUTINES    ************************
//     ********************************************************

double **dmatrix(int nrl,int nrh,int ncl,int nch)
// Allocates a double matrix with range [nrl..nrh][ncl..nch]
{
  int i;
  double **m;

//    Allocate pointers to rows
  m=(double **) malloc((unsigned)(nrh-nrl+1)*sizeof(double *));
  if (!m) nerror("allocation failure in malloc in dmatrix()");
  m -= nrl;
//    Allocate rows and set pointers to them
  for(i=nrl;i<=nrh;i++) {
    m[i]=(double*) malloc((unsigned) (nch-ncl+1)*sizeof(double));
    if (!m[i]) nerror("allocaion failure in malloc in dmatrix()");
    m[i] -= ncl;
  }
  return m;
}

void nerror(char *error_text)
{
  void exit();
  fprintf(stderr, "Run-time error...\n");
  fprintf(stderr,"%s\n",error_text);
  fprintf(stderr,"Exiting...\n");
  exit(1);
}























































