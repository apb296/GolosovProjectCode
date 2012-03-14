#include "mex.h"
#include <math.h>

/*
MINCOL Returns the minimal value in each column of a matrix
z=mincol(x);
*/

/* returns the minimal values only */
double *mincol(double *A, double *B, int m, int n)
{
  double temp, *Aend;
  int  i;
  Aend=A+m*n;
  for (i=0; i<m; i++) B[i]=*A++;
  while (A<Aend) for (i=0; i<m; i++) {
    temp=*A++; if (B[i]>temp) B[i]=temp;
  }
  return(B);
}


/* returns the minimal values and their column indices */
double *mincol2(double *A, double *B, double *ind, int m, int n)
{
  double temp;
  unsigned int  i, j;
  for (i=0; i<m; i++) {B[i]=*A++; ind[i]=1;}
  for (j=2;j<=n; j++) 
    for (i=0; i<m; i++){
      temp=*A++; if (B[i]>temp) {B[i]=temp; ind[i]=j;}
  }
  return(B);
}

void mexFunction(
   int nlhs, mxArray *plhs[],
   int nrhs, const mxArray *prhs[])
{
  int m, n;
   
  if (nrhs<1) mexErrMsgTxt("No input argument passed.");
  if (!mxIsDouble(prhs[0])) mexErrMsgTxt("Matrix must be double");
  if (mxIsSparse(prhs[0])) mexErrMsgTxt("Matrix must be full (not sparse)");
  if (mxIsComplex(prhs[0])) mexErrMsgTxt("Matrix must be real");
  m=mxGetM(prhs[0]);
  n=mxGetN(prhs[0]);
  plhs[0]=mxCreateDoubleMatrix(m,1,mxREAL);
  if (nlhs<2)
    mincol(mxGetPr(prhs[0]),mxGetPr(plhs[0]),m,n);
  else{
    plhs[1]=mxCreateDoubleMatrix(m,1,mxREAL);
    mincol2(mxGetPr(prhs[0]),mxGetPr(plhs[0]),mxGetPr(plhs[1]),m,n);
  }  
}
