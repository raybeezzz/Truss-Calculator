#include "pch.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

void InverseNxN(double **A, int N, double **Ainv, int *indx);
void ELGS (double **A, int N, int *indx);

void main()
{
   double **A,**Ao,**Ainv,d,*row,*col;
   int *indx;
   int i,j,k,N;

   N=3;

   A    = new double *[N];   
   Ao   = new double *[N];
   Ainv = new double *[N];
   indx = new int [N];
   row  = new double [N];
   col  = new double [N];

   for(i=0;i<N;i++)
   {
      A[i] = new double [N];
      Ao[i] = new double [N];
      Ainv[i] = new double [N];
   }

   A[0][0]=3.0;
   A[0][1]=1.0;
   A[0][2]=0.0;
   A[1][0]=2.0;
   A[1][1]=2.0;
   A[1][2]=1.0;
   A[2][0]=1.0;
   A[2][1]=0.0;
   A[2][2]=2.0;

   printf("Matrix A:\n");
   for(i=0;i<N;i++)
   {
      for(j=0;j<N;j++)
      {
         Ao[i][j]=A[i][j];
         printf("% .4lf%c",A[i][j],(j<N-1 ? '\t' : '\n'));
      }
   }

   printf("\nInverse Matrix Ainv:\n");
   InverseNxN(A,N,Ainv,indx);
   for(i=0;i<N;i++)
   {
      for(j=0;j<N;j++)
      {
         printf("% .4lf%c",Ainv[i][j],(j<N-1 ? '\t' : '\n'));
      }
   }

   printf("\nA x Ainv:\n");
   for(i=0;i<N;i++)
   {
      for(j=0;j<N;j++)
      {
         for(k=0;k<N;k++)
         {
            row[k]=Ao[i][k];
            col[k]=Ainv[k][j];
         }

         d=0.0;
         for(k=0;k<N;k++) d+=row[k]*col[k];
         printf("% .4lf%c",d,(j<N-1 ? '\t' : '\n'));
      }
   }


   for(i=0;i<N;i++)
   {
      delete [] A[i];
      delete [] Ao[i];
      delete [] Ainv[i];
   }
   delete A;
   delete Ao;
   delete Ainv;
   delete indx;
   delete row;
   delete col;

   getchar();
}

/* Function to invert matrix A with the inverse stored in Ainv.*/
void InverseNxN (double **A, int N, double **Ainv, int *indx)
{
   int i,j,k;
   double **b;

   b    = new double *[N];
   for(i=0;i<N;i++)
   {
      b[i] = new double [N];
   }

   
   for(i=0;i<N;i++)
   {
      for(j=0;j<N;j++)
      {
         b[i][j] = 0.0;
      }
   }
   for(i=0;i<N;i++)
   {
      b[i][i] = 1.0;
   }
   
   ELGS (A,N,indx);
   
   for(i=0;i<N-1;i++)
   {
      for(j=i+1;j<N;j++)
      {
         for(k=0;k<N;k++)
         {
            b[indx[j]][k] = b[indx[j]][k]-A[indx[j]][i]*b[indx[i]][k];
         }
      }
   }
   
   for(i=0;i<N;i++)
   {
      Ainv[N-1][i] = b[indx[N-1]][i]/A[indx[N-1]][N-1];
      for (j=N-2;j>= 0;j=j-1)
      {
         Ainv[j][i] = b[indx[j]][i];
         for(k=j+1;k<N;k++)
         {
            Ainv[j][i] = Ainv[j][i]-A[indx[j]][k]*Ainv[k][i];
         }
         Ainv[j][i] = Ainv[j][i]/A[indx[j]][j];
      }
   }

   for(i=0;i<N;i++)
   {
      delete [] b[i];
   }
   delete b;
}

/******************************************************************************************************************* 
Function to perform the partial-pivoting Gaussian elimination. A is the original matrix in the input and transformed
matrix plus the pivoting element ratios below the diagonal in the output.  indx[] records the pivoting order. 
*******************************************************************************************************************/
void ELGS (double **A, int N, int *indx)
{
   int i, j, k, itmp;
   double c1, pi, pi1, pj;
   double *c;
   
   c = new double [N];
   
   // Initialize the index
   for(i=0;i<N;i++)
   {
      indx[i] = i;
   }
   
   // Find the rescaling factors, one from each row
   for(i=0;i<N;i++)
   {
      c1 = 0;
      for(j=0;j<N;j++)
      {
         if(fabs(A[i][j]) > c1) c1=fabs(A[i][j]);
      }
      c[i] = c1;
   }
   
   // Search the pivoting (largest) element from each column 
   for(j=0;j<N-1;j++)
   {
      pi1 = 0;
      k = -1;
      for(i=j;i<N;i++)
      {
         pi = fabs(A[indx[i]][j])/c[indx[i]];
         if(pi>pi1)
         {
            pi1 = pi;
            k = i;
         }
      }
      if(k==-1)
      {
         printf("Pivoting problem\n");
         exit(1);
      }
      
      // Interchange the rows via indx[] to record pivoting order
      itmp = indx[j];
      indx[j] = indx[k];
      indx[k] = itmp;
      for(i=j+1;i<N;i++)
      {
         pj = A[indx[i]][j]/A[indx[j]][j];
         
         // Record pivoting ratios below the diagonal
         A[indx[i]][j] = pj;
         
         // Modify other elements accordingly
         for (k=j+1;k<N;k++)
         {
            A[indx[i]][k] = A[indx[i]][k]-pj*A[indx[j]][k];
         }
      }
   }

   delete [] c;
}
