#include "mex.h"
#include <math.h>
#include <time.h>
#include <sys/time.h>

#define DEBUG 0

/*
 * Quick computation of the similarity structure score
 * (fpereira@princeton.edu)
 
 in:
 k x k x m - similarity matrices (k x k) for voxel 1, voxel 2, etc
             (columnwise)
 k x k - similarity structure matrix
         (columnwise, hence a vector)
 k - #conditions
 m - #voxels

 out:
 1 x m - similarity score
 
 assumes:

%   This file is part of Simitar
%
%   Simitar is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%   Simitar is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%
%   You should have received a copy of the GNU General Public License
%   along with Simitar.  If not, see <http://www.gnu.org/licenses/>.
%
*/

void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
  /* input */
  double *similarity;
  double *simmatrix; 
  int k,m;
  double *ptmp;
  

  double *examples1;
  double *examples2;
  int *labels;             double *plabels; /* original labels */
  int *labelsGroup;        double *plabelsGroup;
  int *voxelsToNeighbours; double *pvoxelsToNeighbours;
  int *numberOfNeighbours; double *pnumberOfNeighbours;
  int nPermutations;       double *pnPermutations;

  /* output */
  double *score;

  /* code */

  int i, v, n;
  double *vptr; double *mptr;
  double accumulated;
  double *sptr;
  

  /* new ones (used to keep the patterns in each neighbourhood) */

  double *tmpe1; double *tmpe2;
  double *mean1; double *mean2;
  double *stdv1; double *stdv2;

  double *eptr1; double *eptr2;
  double *tptr1; double *tptr2;
  double *mptr1; double *mptr2;
  double *sptr1; double *sptr2;
  int e1; int e2;
  int c1; int c2;
  double correlation; 
  int moffset;


  /* rest of variables pertaining to input/output contents */
  int nGroups;
  int itmp;
  int iprev,iseen;
  int group;
  int *groupStarts; double *pgroupStarts; mxArray *groupStarts_mx;
  int *groupEnds;   double *pgroupEnds;   mxArray *groupEnds_mx;
  int *groupSize;   double *pgroupSize;   mxArray *groupSize_mx;
  int *groupNclass; double *pgroupNclass; mxArray *groupNclass_mx;
  int **groupPermutations;
  int p;
  int nClasses;
  int n1,m1,n2,m2;
  int g,c,gidx;
  int *gptr;    /* group pointer */
  double *eptr, *e2ptr; /* example pointer */
  int ngc;
  double *sumX1; mxArray *sumX1_mx; double *x1ptr;
  double *sumX2; mxArray *sumX2_mx; double *x2ptr;
  double *m2ptr;
  double *exampleLargestValue;
  double *examplesCorrect; mxArray *examplesCorrect_mx;
  double *examplesSquared;
  double *means;
  double *meansSquared;
  double *vars;
  double *sumX1train;
  double *sumX2train;
  int gtest, gtrain;
  int *nc;
  int ntest, ntrain, ntrainmo;
  int nn;
  int *nptr, *iptr;
  int neighbour;
  double voxelValue;
  double *cptr;
  int neighbourRadius; double *pneighbourRadius; int neighbourMax;
  
  /* everything else */
  int     e,ig,j;

  double *X,*M,*Ez,*Vz,*sigma;
  double *pEz,*pVz,*pX,*pM,*pL;
  double *pEzM, *pEzM2;
  double **mEz,**mVz,**mX,**mM,**mL;
  double fixed;
  double *log_px;
  double *log_px_pd;
  mxArray *tmp1_mx,*tmp2_mx,*sigma2_mx,*logsigma_mx,*pEzM_mx,*pEzM2_mx;
  double  *tmp1,*tmp2,*sigma2,*logsigma;
  int     status,mrows,ncols;
  int     cX,rX,cM,rM,cEz,rEz,cVz,rVz,cS,rS;
  int     idx;
  char    buffer[10];

  struct timeval tvStart; struct timeval tvEnd; struct timeval tvDuration;
  double  mduration;


  /*  check for proper number of arguments */
  if(nrhs!=5) 
    mexErrMsgTxt("five inputs required");

  /* check to make sure the input arguments are double matrices */
  if( !mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) ) {
    mexErrMsgTxt("inputs must be matrices");
  }
  
  /*  create a pointer to the input matrices */
  similarity           = mxGetPr(prhs[0]);
  simmatrix            = mxGetPr(prhs[1]);
  ptmp                 = mxGetPr(prhs[2]);
  k = (int) round(*ptmp);
  ptmp                 = mxGetPr(prhs[3]);
  m = (int) round(*ptmp);

  /* mexPrintf("k=%d\tm=%d\n",k,m);fflush(stdout); */

  /*  output matrix space is given as an input argument */

  score                = mxGetPr(prhs[4]);

  /*
   * main loop
   * - a pointer will advance from matrix for voxel 1, then 2, then 3
   * - while in one matrix the values encountered are multipled by the corresponding ones in simmatrix
   * - the accumulated value is the score for that voxel
   */

  vptr = similarity;
  sptr = score;
  n = k*k;
  

  for (v=0; v<m; v++)
    {
      mptr = simmatrix;
      accumulated = 0;

      for (i=0; i<n; i++) {
	accumulated = accumulated + (*mptr)*(*vptr);

	vptr++;
	mptr++;
      }

      *sptr = accumulated;
      sptr++;
    }

  return;

}
