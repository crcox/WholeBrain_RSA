#include "mex.h"
#include <math.h>
#include <time.h>
#include <sys/time.h>

#define DEBUG 0

/*
 * Efficient computation of searchlight correlation matrices 
 * (fpereira@princeton.edu)
 
 in:
 n x m  - examples;
 n x 1  - labels;
 n x 1  - labelsGroup;
 m x <radius^3 -1> - voxelsToNeighbours;
 m x 1  - numberOfNeighbours;
 
 out:
 1 x m  - error
 1 x m  - fraction
 
 assumes:
 - examples are ordered by group
 - group numbers are 1:#groups, all of them occurring
 - class labels are 1:#classes, all of them occurring
 
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
  double *dptr;
  double *examples1;
  double *examples2;
  int *labels;             double *plabels; /* original labels */
  int *labelsGroup;        double *plabelsGroup;
  int *voxelsToNeighbours; double *pvoxelsToNeighbours;
  int *numberOfNeighbours; double *pnumberOfNeighbours;
  int nPermutations;       double *pnPermutations;

  /* output */
  double *accuracy;
  double *accuracyCount;
  double *fraction;

  double *similarity;
  double *similarity2;

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
  double correlation; double *sptr;
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
  int n,m;
  int n1,m1,n2,m2;
  int g,c,gidx;
  int *gptr;    /* group pointer */
  double *eptr, *e2ptr; /* example pointer */
  int ngc;
  double *sumX1; mxArray *sumX1_mx; double *x1ptr;
  double *sumX2; mxArray *sumX2_mx; double *x2ptr;
  double *mptr,*m2ptr;
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
  int     e,i,ig,j,k,v;

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
  int     whichMeasure;;

  struct timeval tvStart; struct timeval tvEnd; struct timeval tvDuration;
  double  mduration;


  /*  check for proper number of arguments */
  if(nrhs!=7) 
    mexErrMsgTxt("seven inputs required");

  /* check to make sure the input arguments are double matrices */
  if( !mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) ) {
    mexErrMsgTxt("inputs must be matrices");
  }
  
  /*  create a pointer to the input matrices */
  examples1            = mxGetPr(prhs[0]);
  examples2            = mxGetPr(prhs[1]);
  pneighbourRadius     = mxGetPr(prhs[2]);
  pvoxelsToNeighbours  = mxGetPr(prhs[3]);
  pnumberOfNeighbours  = mxGetPr(prhs[4]);
  similarity           = mxGetPr(prhs[5]);
  dptr                 = mxGetPr(prhs[6]);
  n1  = mxGetN(prhs[0]);  m1 = mxGetM(prhs[0]); /* dimensions */
  n2  = mxGetN(prhs[0]);  m2 = mxGetM(prhs[0]); /* dimensions */

  /* mexPrintf("#examples=(%d,%d) m=(%d,%d)\n",n1,n2,m1,m2);fflush(stdout); */
  m = m1;
  whichMeasure = (int) round(*dptr);

  /*  output matrix space is given as an input argument */

  /** figure out things and create temporary matrices */
     
  gettimeofday(&tvStart,NULL);

  /* some of the inputs are really matrices of integers, convert from double to a local integer copy */

  /* turns neighbour information arguments into arrays of ints, also switch to 0-based indexing */
  neighbourRadius = (int) *pneighbourRadius;
  neighbourMax    = pow( 2*neighbourRadius+1, 3)-1;

  numberOfNeighbours = (int *) mxMalloc( m*sizeof(int) );
  for (j=0; j<m; j++) { numberOfNeighbours[j] = (int) floor(pnumberOfNeighbours[j]); }

  voxelsToNeighbours = (int *) mxMalloc( m*neighbourMax*sizeof(int) );
  for (j=0; j<m; j++) {
    nptr = voxelsToNeighbours  + j*neighbourMax;
    mptr = pvoxelsToNeighbours + j*neighbourMax; 
    for (v=0; v<neighbourMax; v++) { nptr[v] = (int) floor(mptr[v]-1); }
  }

  /* debug examples
     for (i=0; i<10; i++) {
       eptr = examples + i*m;
       for (j=0; j<8; j++) { mexPrintf("%1.2f ",eptr[j]); }; mexPrintf("\n"); fflush(stdout);
     }
  */

  /* mexPrintf("#examples1=%d #examples2=%d #voxels=%d\n",n1,n2,m); fflush(stdout);*/

  /*
  for (j=0; j<m; j++)
    {
      nn   = numberOfNeighbours[j];
      nptr = voxelsToNeighbours + j*neighbourMax;

      mexPrintf("neighbours %d:\t",j);fflush(stdout);
      for (v=0; v<nn; v++)
	{
	  neighbour = nptr[v];
	  mexPrintf("%d ",neighbour);fflush(stdout);
	}
      mexPrintf("\n");fflush(stdout);
    }
  */

  /*
    for (j=0; j<m; j++)
      {
        nptr = voxelsToNeighbours + j*neighbourMax;
        mexPrintf("neighbours %d:\t",j);fflush(stdout);
        for (v=0; v<neighbourMax; v++) { mexPrintf("%d ",nptr[v]);fflush(stdout); }; mexPrintf("\n");fflush(stdout);
      }
  */

  /* debug neighbour info
  for (j=0; j<m; j++)
    {
      nn   = numberOfNeighbours[j];
      nptr = voxelsToNeighbours + j*neighbourMax;

      mexPrintf("neighbours %d:\t",j);fflush(stdout);
      for (v=0; v<nn; v++)
	{
	  neighbour = nptr[v];
	  mexPrintf("%d ",neighbour);fflush(stdout);
	}
      mexPrintf("\n");fflush(stdout);
    }
  */

  /* get all the memory allocation out of the way here (only needs to be done once) */

  /* these will store the bits of each example in the neighbourhood of a voxel, in each dataset */
  tmpe1 = mxMalloc(n1*(neighbourMax+1)*sizeof(double));
  tmpe2 = mxMalloc(n2*(neighbourMax+1)*sizeof(double));

  mean1 = mxMalloc(n1*sizeof(double));
  mean2 = mxMalloc(n2*sizeof(double));
  stdv1 = mxMalloc(n1*sizeof(double));
  stdv2 = mxMalloc(n2*sizeof(double));

  gettimeofday(&tvEnd,NULL);

  timersub(&tvEnd,&tvStart,&tvDuration);
  mduration = (double) tvDuration.tv_usec / 1000000;
  /*
    mexPrintf("finished setting up in %lf second(s)",mduration);mexPrintf("\n");fflush(stdout);
  */

  gettimeofday(&tvStart,NULL);

  /**
   ** Main loop
   **/

  moffset = n1 * n2;

  i = 0;
  sptr  = similarity;

  switch(whichMeasure)
    {
    case 1:
      /* correlation */

      for (j=0; j<m; j++)
	{
	  nn   = numberOfNeighbours[j];
	  nptr = voxelsToNeighbours + j*neighbourMax;

	  /** load the two temporary arrays and compute means/standard deviations:
	      - the mean can be computed as each element is loaded
	      - the standard deviation is computed on a second pass
	  **/

	  /* array 1 */

	  for (e1=0; e1<n1; e1++) {

	    /* transfer the neighbourhood to the temporary array and compute mean */

	    eptr1 = examples1 + e1*m;
	    tptr1 = tmpe1     + e1*(neighbourMax+1);
	    c1 = 0;

	    tptr1[c1] = eptr1[j];
	    mean1[e1] = tptr1[c1];
	    c1++;

	    for (v=0; v<nn; v++)
	      {
		neighbour = nptr[v];
		tptr1[c1] = eptr1[neighbour];
		mean1[e1] = mean1[e1] + tptr1[c1];
		c1++;
	      }

	    mean1[e1] = mean1[e1] / (nn+1);

	    /* now compute standard deviation and subtract means */

	    c1 = 0;
	    tptr1[c1] = tptr1[c1] - mean1[e1];
	    stdv1[e1] = pow(tptr1[c1],2);
	    c1++;

	    for (v=0; v<nn; v++)
	      {
		neighbour = nptr[v];
		tptr1[c1] = tptr1[c1] - mean1[e1];
		stdv1[e1] = stdv1[e1] + pow(tptr1[c1],2);
		c1++;
	      }

	    stdv1[e1] = stdv1[e1] / (nn);
	    stdv1[e1] = sqrt(stdv1[e1]);

	    /* and divide by standard deviation */
	    c1 = 0;
	    tptr1[c1] = tptr1[c1] / stdv1[e1]; c1++;
	    for (v=0; v<nn; v++)
	      {
		tptr1[c1] = tptr1[c1] / stdv1[e1]; c1++;
	      }


	  } /* end of loop over examples */

	  /*mexPrintf("after loading 1\n"); fflush(stdout);*/

	  /* array 2 */


	  for (e2=0; e2<n2; e2++) {

	    /* transfer the neighbourhood to the temporary array and compute mean */

	    eptr2 = examples2 + e2*m;
	    tptr2 = tmpe2     + e2*(neighbourMax+1);
	    c2 = 0;

	    tptr2[c2] = eptr2[j];
	    mean2[e2] = tptr2[c2];
	    c2++;

	    for (v=0; v<nn; v++)
	      {
		neighbour = nptr[v];
		tptr2[c2] = eptr2[neighbour];
		mean2[e2] = mean2[e2] + tptr2[c2];
		c2++;
	      }

	    mean2[e2] = mean2[e2] / (nn+1);

	    /* now compute standard deviation and subtract means */

	    c2 = 0;
	    tptr2[c2] = tptr2[c2] - mean2[e2];
	    stdv2[e2] = pow(tptr2[c2],2);
	    c2++;

	    for (v=0; v<nn; v++)
	      {
		neighbour = nptr[v];
		tptr2[c2] = tptr2[c2] - mean2[e2];
		stdv2[e2] = stdv2[e2] + pow(tptr2[c2],2);
		c2++;
	      }

	    stdv2[e2] = stdv2[e2] / (nn);
	    stdv2[e2] = sqrt(stdv2[e2]);

	    /* and divide by standard deviation */
	    c2 = 0;
	    tptr2[c2] = tptr2[c2] / stdv2[e2]; c2++;
	    for (v=0; v<nn; v++)
	      {
		tptr2[c2] = tptr2[c2] / stdv2[e2]; c2++;
	      }


	  } /* end of loop over examples */
      
	  /** compute correlation **/

	  for (e1=0; e1<n1; e1++) {
	    for (e2=0; e2<n2; e2++){

	      tptr1 = tmpe1 + e1*(neighbourMax+1);
	      tptr2 = tmpe2 + e2*(neighbourMax+1);

	      correlation = (*tptr1)*(*tptr2);

	      for (v=0; v<nn; v++)
		{
		  neighbour = nptr[v];

		  tptr1++;
		  tptr2++;
		  correlation = correlation + (*tptr1)*(*tptr2);
		}

	      *sptr = correlation/nn;
	      sptr++;
	    }
	  }
	} /* end of loop over voxels */
  
      break;


    case 2:

      /* euclidean */

      for (j=0; j<m; j++)
	{
	  nn   = numberOfNeighbours[j];
	  nptr = voxelsToNeighbours + j*neighbourMax;

	  /** 
	      load the two temporary arrays containing examples in neighbourhood
	  **/

	  /* array 1 */

	  for (e1=0; e1<n1; e1++) {

	    /* transfer the neighbourhood to the temporary array */
	    eptr1 = examples1 + e1*m;                /* m voxels per row in examples */
	    tptr1 = tmpe1     + e1*(neighbourMax+1); /* neighbourMax+1 voxels per row */
	    c1 = 0;

	    tptr1[c1] = eptr1[j];
	    c1++;

	    for (v=0; v<nn; v++)
	      {
		neighbour = nptr[v];
		tptr1[c1] = eptr1[neighbour];
		c1++;
	      }
	  } /* end of loop over examples */

	  /* array 2 */

	  for (e2=0; e2<n2; e2++) {

	    /* transfer the neighbourhood to the temporary array */
	    eptr2 = examples2 + e2*m;                /* m voxels per row in examples */
	    tptr2 = tmpe2     + e2*(neighbourMax+1); /* neighbourMax+1 voxels per row */
	    c2 = 0;

	    tptr2[c2] = eptr2[j];
	    c2++;

	    for (v=0; v<nn; v++)
	      {
		neighbour = nptr[v];
		tptr2[c2] = eptr2[neighbour];
		c2++;
	      }
	  } /* end of loop over examples */
      

	  /** compute euclidean distance **/

	  for (e1=0; e1<n1; e1++) {
	    for (e2=0; e2<n2; e2++){

	      tptr1 = tmpe1 + e1*(neighbourMax+1);
	      tptr2 = tmpe2 + e2*(neighbourMax+1);

	      /* correlation = (*tptr1)*(*tptr2); */
	      correlation = pow( (*tptr1)-(*tptr2), 2);

	      for (v=0; v<nn; v++)
		{
		  tptr1++;
		  tptr2++;
		  correlation = correlation + pow( (*tptr1)-(*tptr2), 2);
		}

	      /**sptr = sqrt(correlation);*/
	      *sptr = correlation;
	      sptr++;
	    }
	  }
	} /* end of loop over voxels */
      break;
    }
    
  gettimeofday(&tvEnd,NULL);

  /* give the user an estimate of running time, so that they can abort if necessary */
  timersub(&tvEnd,&tvStart,&tvDuration);
  mduration = (double) tvDuration.tv_usec / 1000000;

  /*
  mexPrintf("finished running one iteration in %lf second(s)",mduration);
  if (nPermutations) { mexPrintf(", estimated running time for the permutation test (%d permutations) is %d seconds.",nPermutations,(int) round(mduration*nPermutations));}
  else { mexPrintf("\tall done!"); }
  mexPrintf("\n");fflush(stdout);
  */
  
  /* get rid of all the matrices obtained with mxMalloc (in the order the mxMallocs appear in the code) */
  /* mxFree(accuracyCount); */

  mxFree(numberOfNeighbours);
  mxFree(voxelsToNeighbours);

  mxFree(tmpe1); mxFree(tmpe2);
  mxFree(mean1); mxFree(stdv1);
  mxFree(mean2); mxFree(stdv2);


  /* no longer used
     mxFree(labels);
     mxFree(labelsGroup);
     mxFree(groupStarts);
     mxFree(groupEnds);
     mxFree(groupSize);  
     mxFree(groupNclass);
     mxFree(sumX1);
     mxFree(sumX2);
     mxFree(examplesSquared);
     mxFree(nc);
     mxFree(means);
     mxFree(vars);
     mxFree(meansSquared);
     mxFree(sumX1train);
     mxFree(sumX2train);
  */
  /* mxFree(exampleLargestValue); */
  /* for (g=0; g<nGroups; g++) { mxFree(groupPermutations[g]); } */
  /* mxFree(groupPermutations); */

  return;


  /*  mexPrintf("here\n");fflush(stdout); */

  /* create other temporary matrices and C pointers to them
  sigma2_mx = mxCreateDoubleMatrix(m,1,mxREAL);
  sigma2 = mxGetPr(sigma2_mx);

  logsigma_mx = mxCreateDoubleMatrix(m,1,mxREAL);
  logsigma = mxGetPr(logsigma_mx);

  tmp1_mx = mxCreateDoubleMatrix(m,1,mxREAL);
  tmp1 = mxGetPr(tmp1_mx);

  pEzM_mx = mxCreateDoubleMatrix(m,1,mxREAL);
  pEzM = mxGetPr(pEzM_mx);
  pEzM2_mx = mxCreateDoubleMatrix(m,1,mxREAL);
  pEzM2 = mxGetPr(pEzM2_mx);
  */


  /*
   * do the computation
   */
  
  /* debug: access  matrices coming in */
  /* NO PROBLEM HERE
  for (e=0; e<n; e++) {
    pEz = Ez + e*l;
    for (i=0; i<l; i++) { mexPrintf("%1.4f ",pEz[i]); }; mexPrintf("\n"); fflush(stdout);
  }
  mexPrintf("\n");

  for (e=0; e<n; e++) {
    pVz = Vz + e*l;
    for (i=0; i<l; i++) { mexPrintf("%1.4f ",pVz[i]); }; mexPrintf("\n"); fflush(stdout);
  }

  mexPrintf("\n");

  for (e=0; e<n; e++) {
    pX = X + e*m;
    for (j=0; j<4; j++) { mexPrintf("%1.4f ",pX[j]); };
    for (j=(m-4); j<m; j++) { mexPrintf("%1.4f ",pX[j]); }; 
    mexPrintf("\n"); fflush(stdout);
  }
  mexPrintf("\n");

  for (i=0; i<l; i++) {
    pM = M + i*m;
    for (j=0; j<4; j++) { mexPrintf("%1.4f ",pM[j]); };
    for (j=(m-4); j<m; j++) { mexPrintf("%1.4f ",pM[j]); }; 
    mexPrintf("\n");
  }
  mexPrintf("\n");

  return;
  */

  /* initialize the result

  *log_px = 0;  
  for (i=0; i<l; i++) {
    pL = log_px_pd + i*m;
    for (j=0; j<m; j++) {
      pL[j] = 0;
    }
  }
  */


  /* precompute things outside the loop

  for (j=0; j<m; j++) { sigma2[j]   = sigma[j]*sigma[j]; }
  for (j=0; j<m; j++) { logsigma[j] = log(sigma[j]); }
  fixed = 0.5 * log(2*PI);
  */

  for (e=0; e<n; e++) {

    /** initialize pointers for this example */
    /*
    pX  = X  + e*m;
    pEz = Ez + e*l;
    pVz = Vz + e*l;
    */

    /** 1) work on the log_px computation for this example */

    /* compute E[zM] and E[(zM)^2]
    for (j=0; j<m; j++) { pEzM[j] = 0;}
    for (j=0; j<m; j++) { pEzM2[j] = 0;}

    for (i=0; i<l; i++) {
      pM = M + i*m;
      for (j=0; j<m; j++) {
	pEzM[j] = pEzM[j] + pEz[i] * pM[j];
      }
    }

    for (i=0; i<l; i++) {
      pM = M + i*m;
      for (j=0; j<m; j++) {
	pEzM2[j] = pEzM2[j] + pVz[i]*(pM[j]*pM[j]);
      }
    }
    for (j=0; j<m; j++) {
      pEzM2[j] = pEzM2[j] + (pEzM[j]*pEzM[j]);
    }
    

    for (j=0; j<m; j++) {
      *log_px = *log_px + (pX[j]*pEzM[j] -0.5*pEzM2[j] -0.5*(pX[j]*pX[j]))/sigma2[j] -logsigma[j] -fixed;
    }
    */


    /*
    for (j=0; j<5; j++) { mexPrintf("%1.3f ",pEzM[j]); };  mexPrintf("\n"); fflush(stdout);
    for (j=0; j<5; j++) { mexPrintf("%1.3f ",pEzM2[j]); }; mexPrintf("\n"); fflush(stdout);
    mexPrintf("%f\n",*log_px);
    mexPrintf("aqui!\n");fflush(stdout); scanf("%s",buffer);
    */

    /* 2) work on the log_px_pd computation for this example */

    /* for (j=0; j<5; j++) { mexPrintf("%1.3f ",pX[j]); }; mexPrintf("\n"); fflush(stdout); */
    /* for (i=0; i<l; i++) { mexPrintf("%1.3f ",pEz[i]); }; mexPrintf("\n"); fflush(stdout); */

    /* precompute Ez*M  for this example */
    
    /*
    for (j=0; j<m; j++) { tmp1[j] = 0; }
    for (i=0; i<l; i++) {
      pM = M + i*m;
      for (j=0; j<m; j++) {	
	tmp1[j] = tmp1[j] + pEz[i] * pM[j];
	idx++;
      }
    }
    */

    /*
    for (j=0; j<5; j++) { mexPrintf("%1.4f ",tmp1[j]); }; mexPrintf("\n"); fflush(stdout);
    mexPrintf("aqui!\n");fflush(stdout); scanf("%s",buffer);
    */
    /* printf("aqui!\n");fflush(stdout); scanf("%s",buffer); */

  } /* end of loop over examples */

  /*
  mxDestroyArray(tmp1_mx);
  mxDestroyArray(sigma2_mx);
  mxDestroyArray(logsigma_mx);
  */
  /* mxFree(tmp1); */
}
