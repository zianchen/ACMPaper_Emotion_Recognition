#include "mex.h"
#include <math.h>

void mexFunction(
int		
nlhs, 	/* number of expected outputs */
mxArray	  *plhs[],	/* mxArray output pointer array */
int		  nrhs, 	/* number of inputs */
const mxArray	  *prhs[]	/* mxArray input pointer array */
)
{
    double *V;
    double *Gx;
    double *Gy;
    double *Dx;
    double *Dy;
    double *Lx;
    double *Ly;
    const int*  sizes;
    int i,j,p,k,m,n;
    double Gpq;
    double dx,dy,lx,ly;
    m=mxGetM(prhs[0]);  //number of row
    n=mxGetN(prhs[0]);  // number colm
//    int m = mxGetM(prhs[0]);
//    int n = mxGetN(prhs[0]);
    V = mxGetPr(prhs[0]);
    Gx = mxGetPr(prhs[1]);
    Gy = mxGetPr(prhs[2]);
 //   mexPrintf("m:%d,n:%d\n",m,n);
    plhs[0]=mxCreateDoubleMatrix(m,n,mxREAL);
    plhs[1]=mxCreateDoubleMatrix(m,n,mxREAL);
    plhs[2]=mxCreateDoubleMatrix(m,n,mxREAL);
    plhs[3]=mxCreateDoubleMatrix(m,n,mxREAL);
    Dx = mxGetPr(plhs[0]);
    Dy = mxGetPr(plhs[1]);
    Lx = mxGetPr(plhs[2]);
    Ly = mxGetPr(plhs[3]);
    for(i=0;i<m;i++){
        for(j=0;j<n;j++){
           for(p=i-1;p<i+1;p++){
               dx = 0;
               dy = 0;
               lx = 0;
               ly = 0;
                for(k = j-1;k<j+1;k++){
                    if(( p>-1 && p<m )&&( k>-1 && k<n )){
                        Gpq = exp(-(pow(i-p,2)+(pow(j-k,2)/8)));
                        dx += Gpq*abs(Gx[j*m+i]);
                        dy += abs(Gpq*Gy[k*m+p]);
                        lx += Gpq*abs(Gx[j*m+i]);
                        ly += abs(Gpq*Gy[k*m+p]);
                    }
                }
            }
        Dx[j*m+i] = dx;
        Dy[j*m+i] = dy;
        Lx[j*m+i] = lx;
        Ly[j*m+i] = ly;
        }
    }
         
 
 //   mexPrintf("hello,world!\n");
}