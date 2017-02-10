#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <mkl.h>
#include "gurobi_c.h"
#include "mat.h"
#define MIN(a,b) (a>b?b:a)
#define MAX(a,b) (a>b?a:b)
#define ABS(a) ((a)>0?(a):-1.0*(a))

#define BUFSIZE 256

// Written by Gonghua Li, Ph.D
// Kunming Institute of Zoology, CAS, PR China
// Mailto: ligonghua@mail.kiz.ac.cn
//         ligonghuabio@gmail.com
// 20161016

double norm_vector(double *d, int n)
{
	int i;
	double out;
	out = 0;
	for (i =0;i<n ;i++ ) out = out+ d[i]*d[i]; 
	out = sqrt(out);
	return out;
}

double maxabs_vector(double *v, int n)
{
	int i;
	double out,tmp;
	out = -1000;
	for (i =0; i<n ;i++ )
	{
		tmp = ABS(v[i]);
		out = MAX(out,tmp);
	}
	return out;
}

double** malloc2d(int row,int colum) //malloc 2-dimension double
{
	double **p;
	int i;
	p = (double **)malloc(row*sizeof(double *));
	for (i =0;i<row ;i++ ) p[i]= (double *)malloc(colum*sizeof(double));
	return(p);
}
void free2d(double **p,int row)
{
	int i;
	for (i =0;i<row ;i++ ) free(p[i]);
	free(p);
}

int str2int(const char *cc)
{
	int aa,bb;
	sscanf(cc,"%d",&aa);
	bb=aa;
	return(bb);
}
double str2double(const char *cc)
{
	double aa,bb;
	sscanf(cc,"%lg",&aa);
	bb=aa;
	return(bb);
}

int main(int   argc, char **argv)
{
	GRBenv   *env   = NULL,*modelenv = NULL;
	GRBmodel *model = NULL;
	double  *Svalue = NULL;
	int  *Srowind = NULL, *Scolind = NULL;
	int erro = 0;
	int i,j,k,m,n,a,b,nrxns,nmets,npoints,nsteps,Snnz,len;
	int       error = 0;
    double *sol = NULL;
	double *srow = NULL;
	double **warmupPoints = NULL;
	double **points = NULL;
	//double *S_full = NULL; //vector just for call matlat null function
	double *UB,*LB;
    int       optimstatus;
    double    objval,duration,tmp;
	char infile[200],outfile[200],infile_nullS[200];
	double *N = NULL, *tN = NULL,*tmpvectorN = NULL, *tmpvectorM  = NULL,*tmpvectorFlux;// null_S
	double maxMinTol = 1e-6; // Minimum allowed distance to the closest constraint, default
    double uTol = 1e-6; // Ignore directions where u is really small, default
    double dTol = 1e-7; // Project out of directions that are too close to the boundary,default

	double *data = NULL;
    MATFile *pmat = NULL; //mat file
	mxArray *pa0 = NULL, *pa1 = NULL;  //mat array
	int status; 

	clock_t start, finish;
	start = clock();
	
    
	if (argc<4)
	{
	    fputs ("\nNot enough inputs!\n",stderr);
		fputs ("\nUsage:mcmc_gurobi mpsfile numpoints numsteps outfile_mat!\n",stderr);
		exit (1);
	}
	printf("Number of argc: %i\n",argc);
    strcpy(infile,argv[1]);
	strcpy(infile_nullS,argv[2]);
    npoints=str2int(argv[3]);
	nsteps=str2int(argv[4]);
	strcpy(outfile,argv[5]);
    if (argc>7)
    {
		maxMinTol = str2double(argv[6]);
		uTol = str2double(argv[6]);
		dTol = str2double(argv[7]);
		printf("Use user defined maxMinTol,uTol,dTol:\n");

		//printf("maxMinTol = %e   uTol = %e   dTol = %e \n",maxMinTol,uTol,dTol);
    }else{
		printf("Use Default maxMinTol,uTol,dTol:\n");
	}
	printf("maxMinTol = %e   uTol = %e   dTol = %e \n",maxMinTol,uTol,dTol);
	//set model env
	error = GRBloadenv(&env, "");
	error = GRBsetintparam(env, GRB_INT_PAR_LOGTOCONSOLE, 0);
	error = GRBsetintparam(env, "Method", 0);
	//error = GRBsetintparam(env, "Method", 1);

	//read model
	error = GRBreadmodel(env, infile, &model);
	if (error){printf("Error reading model\n\n");exit(1);}
	//error = GRBsetdblparam(GRBgetenv(model), "LogToConsole", 0);
    	
	error = GRBgetintattr(model,GRB_INT_ATTR_NUMCONSTRS, &nmets);
	error = GRBgetintattr(model,GRB_INT_ATTR_NUMVARS, &nrxns);
	error = GRBgetintattr(model,GRB_INT_ATTR_NUMNZS, &Snnz);

    
	//malloc 
    //S = Malloc_MODELS(nmets,nrxns);
	Svalue = (double *)malloc(Snnz*sizeof(double));
	Srowind = (int *)malloc(Snnz*sizeof(int));
	Scolind = (int *)malloc(Snnz*sizeof(int));

	sol = (double *)malloc(nrxns*sizeof(double));
	srow = (double *)malloc(nrxns*sizeof(double));
	UB = (double *)malloc(nrxns*sizeof(double));
	LB = (double *)malloc(nrxns*sizeof(double));
	//S_full = (double *)malloc(nmets*nrxns*sizeof(double));
	warmupPoints = malloc2d(2*nrxns,nrxns);
	points = malloc2d(npoints,nrxns);

	printf("nmets = %i  nrxns = %i  nnz = %i\n",nmets,nrxns,Snnz);

    // load prepared nullspace of S "mat format"
	pmat = matOpen(infile_nullS, "r");
    if (pmat == NULL) {
        printf("Error reopening file %s\n", infile_nullS);
        return 0;
    }
	//matClose(pmat);

    len = strlen(infile_nullS);
	infile_nullS[len-4] = '\0';
	printf("%s\n",infile_nullS);
	//pa0 = matGetVariable(pmat, infile_nullS);
	//pa0 = mxCreateDoubleMatrix(1,1,mxREAL);
	pa0 = matGetVariable(pmat, infile_nullS);

	m = mxGetM(pa0);
	n = mxGetN(pa0);
    tN = (double *)malloc(nrxns*n*sizeof(double));//rowmajor
	N = (double *)malloc(nrxns*n*sizeof(double));
	tmpvectorM = (double *)malloc(nrxns*sizeof(double));//nrxns
	tmpvectorN = (double *)malloc(n*sizeof(double));//number of 
	tmpvectorFlux = (double *)malloc(nrxns*sizeof(double));


	memcpy(tN,mxGetPr(pa0), m*n*sizeof(double));//matlab using col major
	for (i =0;i<nrxns ;i++ )
	{
		for (j =0;j<n ;j++ )
		{
			N[i*n+j] = tN[j*nrxns+i];
		}
	}

	// get model S
	k = 0; 
	for (i =0;i<nmets ;i++ )
	{
		for (j =0;j<nrxns ;j++ )
		{		
			error = GRBgetcoeff(model, i, j,&tmp);
			//S_full[j*nmets+i] = tmp;
			if (tmp != 0)
			{
				Svalue[k] = tmp;
				Srowind[k] = i;
				Scolind[k] = j;
				k = k+1;
			}
		}
	}



    start = clock();
	printf("Begin to Warmup£ºnumber Warmup = %i\n",2*nrxns);
    //Warmup
	for (i=0;i<nrxns ;i++ )
	{
		error = GRBsetdblattrelement(model, GRB_DBL_ATTR_OBJ, i, 1.0);
		error = GRBsetintattr(model, GRB_INT_ATTR_MODELSENSE, GRB_MAXIMIZE);
		error = GRBoptimize(model);
		//error = GRBupdatemodel(model);
		error = GRBgetdblattrarray(model, GRB_DBL_ATTR_X, 0, nrxns, warmupPoints[i]);
		error = GRBgetdblattr(model, GRB_DBL_ATTR_OBJVAL, UB+i);

	    error = GRBsetintattr(model, GRB_INT_ATTR_MODELSENSE, GRB_MINIMIZE);
		error = GRBoptimize(model);
		//error = GRBupdatemodel(model);
		error = GRBgetdblattrarray(model, GRB_DBL_ATTR_X, 0, nrxns, warmupPoints[i+nrxns]);
		error = GRBgetdblattr(model, GRB_DBL_ATTR_OBJVAL, LB+i);

		error = GRBsetdblattrelement(model, GRB_DBL_ATTR_OBJ, i, 0.0);

	}





	
	finish = clock();
    duration = (double)(finish - start) / CLOCKS_PER_SEC;
	printf( "\nWarmup spend: %f CPU seconds\n", duration );

	//mcmc simulation:ACHRSampler Artificial Centering Hit-and-Run sampler 
	/* method same as matlab */
    //double maxMinTol = 1e-6; // Minimum allowed distance to the closest constraint
    //double uTol = 1e-6; // Ignore directions where u is really small
    //double dTol = 1e-7; // Project out of directions that are too close to the boundary
	double *centerPoint=NULL, *prevPoint = NULL, *curPoint = NULL, *randVector=NULL,*u;
	int stepCount,randPointID,totalStepCount,nWrmup = 2*nrxns;
	double tmprnd,norm_u,maxStep,minStep,stepDist,alpha,beta;
	double *distUB, *distLB;
	char transa = 'N';
	int rndint[2],is_continue;
	VSLStreamStatePtr stream;
	time_t t;

	centerPoint = (double *)malloc(nrxns*sizeof(double));
	prevPoint = (double *)malloc(nrxns*sizeof(double));
	curPoint = (double *)malloc(nrxns*sizeof(double));
	distUB = (double *)malloc(nrxns*sizeof(double));
	distLB = (double *)malloc(nrxns*sizeof(double));
	u = (double *)malloc(nrxns*sizeof(double));

	randVector = (double *)malloc(nsteps*sizeof(double));
	points = malloc2d(npoints,nrxns);

    srand((unsigned) time(&t));

    //centerPoint and pervPoint
	for (i =0;i<nrxns ;i++ ) {centerPoint[i] = 0;}
	for (i = 0;i<nrxns ;i++ )
	{
		for (j = 0;j<2*nrxns ;j++ )
		{
			centerPoint[i] = centerPoint[i] + warmupPoints[j][i];
		}
	}
	for (i =0;i<nrxns ;i++ )
	{
		centerPoint[i] = centerPoint[i]/nrxns/2;
		prevPoint[i] = centerPoint[i];
	}

    // % Move points in same as matlab
	// warmupPts = warmupPts*.33 + .67*centerPoint*ones(1,nPoints);
	for (i =0; i<2*nrxns ;i++ )
	{
		cblas_dscal(nrxns,0.33, warmupPoints[i], 1);
		cblas_daxpy(nrxns, 0.67, centerPoint, 1, warmupPoints[i], 1);	
	}

    printf("Begin to MCMC\n",nmets,nrxns,Snnz);

	vslNewStream(&stream, VSL_BRNG_SFMT19937, 777);

	start = clock();
	//mcmc
	totalStepCount = 0;
	for (i=0;i<npoints ;i++ )
	{
		printf("%i th point\n", i);
		status = vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD_ACCURATE, stream, nsteps, randVector, 0.0, 1.0);
		//for (k =0;k<nsteps ;k++ ) randVector[k] = rand()/(RAND_MAX+1.0);

		stepCount = 0;

		while(stepCount<nsteps)
		{
			
			// Pick a random warmup point
			//tmprnd = rand()/(RAND_MAX+1.0);
			//randPointID = (int)(2*nrxns*tmprnd);
			status = viRngUniform(VSL_RNG_METHOD_UNIFORM_STD, stream, 2, rndint, 0, 2*nrxns);
			randPointID = rndint[0];
			//printf("%i\n",randPointID);
			// Get a direction from the center point to the warmup point
			for (k=0; k<nrxns;k++ ) u[k] = warmupPoints[randPointID][k] - centerPoint[k];
			//norm_u = norm_vector(u,nrxns);
			norm_u = cblas_dnrm2(nrxns, u, 1);
			//for (k=0; k<nrxns;k++ ) u[k] = u[k]/norm_u;
			cblas_dscal(nrxns,1.0/norm_u, u, 1);

			// Figure out the distances to upper and lower bounds
			// Figure out the true max & min step sizes
			maxStep = 1e9;
			minStep = -1e9;
			for (k=0;k<nrxns ;k++ )
			{
				distUB[k] = UB[k] - prevPoint[k];
				distLB[k] = prevPoint[k]-LB[k];
				if (distUB[k] > dTol & distLB[k] > dTol)
				{
					if (u[k]>uTol)
					{
						maxStep = MIN(maxStep,distUB[k]/u[k]);
						minStep = MAX(minStep,-distLB[k]/u[k]);
					}
					if (u[k] < -uTol)
					{
						maxStep = MIN(maxStep,-distLB[k]/u[k]);
						minStep = MAX(minStep,distUB[k]/u[k]);
					}
				}

			}
			
			//if (maxStep >1000 | minStep<-1000)
			//{
			//	continue;
			//}
			if ((ABS(minStep) < maxMinTol & ABS(maxStep) < maxMinTol) | minStep >maxStep)
			{
				continue;

			}

			stepDist = randVector[stepCount]*(maxStep-minStep)+minStep;

			//printf("%i =  %.3f   %.3f %.3f\n",totalStepCount,maxStep,minStep,stepDist);

			
			//for (k=0;k<nrxns ;k++ )
			//{
			//	curPoint[k] = prevPoint[k] + stepDist*u[k];
			//}
			cblas_daxpy(nrxns, stepDist, u, 1, prevPoint, 1);//prevPoint updated

			//printf("111:%f\n", prevPoint[1]);


			// Reproject the current point and go to the next step
			//if ((totalStepCount+1) % 10 == 0)
			//{
			//	MODLES_multiply_vector(S, curPoint, tmpvectorFlux);
			//	if (maxabs_vector(tmpvectorFlux,nrxns)>1e-9)
			//	{
			//		multply_matrix_vector(tN,curPoint,tmpvectorN);
			//		multply_matrix_vector(N,tmpvectorN,curPoint);
			//	}
			//}
			if ((totalStepCount+1) % 100 == 0)
			{
				 is_continue = 0;
				 for (k=0;k<nrxns ;k++ )
			     {
			    	 if (prevPoint[k]>UB[k]) {prevPoint[k] = UB[k]; is_continue = 1;}
			    	 if (prevPoint[k]<LB[k]) {prevPoint[k] = LB[k]; is_continue = 1;}
			     }
				//printf("%i\n",totalStepCount);
				//need sparse multiply  need mXm matrix for multipy so  reult have nrxns value,
				//but just 1:nmets is non zeros, others are zeros.
		//         mkl_dcoogemv (&transa , &nrxns , Svalue , Srowind , Scolind , &Snnz , prevPoint , sol );
		//         if (maxabs_vector(sol,nmets)>1e-9)
		         if (is_continue == 1)
		         {
					cblas_dgemv(CblasRowMajor, CblasNoTrans,n,nrxns,1.0,tN,nrxns,prevPoint,1,0,tmpvectorN,1);
					//printf("222:%f\n", prevPoint[1]);
					cblas_dgemv(CblasRowMajor, CblasNoTrans,nrxns,n,1.0,N,n,tmpvectorN,1,0,prevPoint,1);
					//printf("3333:%f\n", prevPoint[1]);
		         }
			}

			//for (k=0;k<nrxns ;k++ )
			//{
			//	if (prevPoint[k]>UB[k]) prevPoint[k] = UB[k];
			//	if (prevPoint[k]<LB[k]) prevPoint[k] = LB[k];
			//}
			
			

            // for (k=0;k<nrxns ;k++ ) prevPoint[k] = curPoint[k]; prevPoint updated


			stepCount = stepCount + 1;
			totalStepCount = totalStepCount + 1;

			//recalculate the center point
			alpha = 1.0/(nWrmup+totalStepCount+1);
			beta = (nWrmup + totalStepCount)/(nWrmup+totalStepCount+1);
			cblas_daxpby(nrxns, alpha,prevPoint, 1, beta, centerPoint, 1);

			//for (k=0;k<nrxns ;k++ )
			//{
			//	centerPoint[k] = ((nWrmup + totalStepCount)*centerPoint[k] + curPoint[k])/(nWrmup+totalStepCount+1);
			//}
		}

		//memcpy(points[i], curPoint,nrxns*sizeof(double));
		memcpy(points[i], prevPoint,nrxns*sizeof(double));
		//pointCount = pointCount + 1;
	} 


	//write points to mat file example warmup

	//printf("Hello1\n");

	
	data = (double *)malloc(npoints*nrxns*sizeof(double));
	
	for (i = 0;i<npoints ;i++ )
	{
		//for (j = 0;j<nrxns ;j++ )
		//{
		//	data[i*nrxns+j] = points[i][j];
		//}
		//printf("i = %i\n",i);
		memcpy(data+i*nrxns, points[i],nrxns*sizeof(double));
	}
	//printf("Hello2\n");
	pmat = matOpen(outfile, "w");
	pa1 = mxCreateDoubleMatrix(nrxns,npoints,mxREAL);// matlab will transposes
	memcpy((void *)(mxGetPr(pa1)), (void *)data, npoints*nrxns*sizeof(double));
	status = matPutVariable(pmat, "Points", pa1);
    if (status != 0) {
      printf("%s :  Error using matPutVariable on line %d\n", __FILE__, __LINE__);
      return(EXIT_FAILURE);
    } 

   
    /*//write to mat file example warmup
	double *data = NULL;
    MATFile *pmat; //mat file
	mxArray *pa1;  //mat array
	int status; 
	data = (double *)malloc(2*nrxns*nrxns*sizeof(double));
	
	for (i = 0;i<2*nrxns ;i++ )
	{
		memcpy(data+i*nrxns, warmupPoints[i],nrxns*sizeof(double));
	}
	pmat = matOpen(outfile, "w");
	pa1 = mxCreateDoubleMatrix(nrxns,2*nrxns,mxREAL);// matlab will transposes
	memcpy((void *)(mxGetPr(pa1)), (void *)data, 2*nrxns*nrxns*sizeof(double));
	status = matPutVariable(pmat, "warmupPoints", pa1);
    if (status != 0) {
      printf("%s :  Error using matPutVariable on line %d\n", __FILE__, __LINE__);
      return(EXIT_FAILURE);
    } 
	*/



	//printf( "vv  = %.3f", warmupPoints[0][6798] );

	finish = clock();
    duration = (double)(finish - start) / CLOCKS_PER_SEC;
    printf( "\nMCMC spend: %f CPU seconds\n", duration );
	return 0;
}





// subfunctions
