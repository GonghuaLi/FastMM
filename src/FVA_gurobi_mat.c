#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "gurobi_c.h"
#include "mat.h"
//#include <glpk.h>
#include <time.h>
#include <FastKO.h>


int help()
{
	printf("\n---------------------FastMM:FVA-----------------------\n");
	printf("\n---------------------  Version 0.01   -----------------------\n");
    printf("        C version of Fluxvariability analysis   \n");
    printf("Written By Dr Gong-Hua Li.\n");
    printf("State Key Laboratory of Genetic Resources and Evolution,\n");
    printf("  Kunming Institute  of Zoology, Chinese Academy of Sciences,\n");
    printf("  32, Eastern Jiaochang Road, Kunming, Yunnan 650223, China.\n");
    printf("Please mail to: ligonghua@mail.kiz.ac.cn\n");
    printf("This programm is under GNU GPL, you can freely use,modify and\n");
    printf("  distribute. \n");
    printf("--------------------------------------------------\n\n");
	printf("Uasage: FVA\n");
	printf("        -m cobramodel(require)\n");
	printf("        -o outputfile (require)\n");
	printf("        -c constraintfile (option,defent no additional constraint)\n");
	printf("eg: FVA -m cobramodel -o outputfile -c constraintfile\n\n");

	return 0;
}

int main(int argc,char **argv)
{ 
	int i,j,k,m,n,s,num_gen,num_rxns,nnz;
	int error,numRxns, numMets,*coltype;
	double duration,type;
	char buf[2048],objfile[300],outputfile[300],optimizeType[100],constraintfile[300],infile[300],outfile[300];

	int status; 
    clock_t start, finish;
	start = clock();

	GRBenv   *env   = NULL,*modelenv = NULL;
	GRBmodel *model = NULL;
	//FILE *fin, *fout;
	

	error = GRBloadenv(&env, "");
	error = GRBsetintparam(env, GRB_INT_PAR_LOGTOCONSOLE, 0);
	error = GRBsetintparam(env, "Method", 0);


	//read matfile
	MATFile *pmat = NULL,*pmatout = NULL; //mat file
	mxArray *plb = NULL,*pub = NULL,*pb = NULL,*pc = NULL;
	mxArray *pbed = NULL,*pind = NULL,*pval = NULL,*pa1 = NULL,*pa2 = NULL;  //mat array
	
	strcpy(infile,argv[1]);
	strcpy(outfile,argv[2]);
	pmat = matOpen(infile, "r");
    if (pmat == NULL) {
        printf("Error reopening file %s\n", infile);
        return 0;
    }

	/* mat file contains: lb, ub, b, c, (bed,ind,val)*/

	plb = matGetVariable(pmat,"lb");
	pub = matGetVariable(pmat,"ub");
	pb = matGetVariable(pmat,"b");
	pc = matGetVariable(pmat,"c");
	pbed = matGetVariable(pmat,"bed");
	pind = matGetVariable(pmat,"ind");
	pval = matGetVariable(pmat,"val");

	numRxns = mxGetM(plb);
	numMets = mxGetM(pb);
	nnz = mxGetM(pind);

	double *lb,*ub,*c,*b,*bed_double,*ind_double,*val;
	int *bed,*ind;
	char *sense;
	lb = (double *)malloc(numRxns*sizeof(double));
	ub = (double *)malloc(numRxns*sizeof(double));
	c = (double *)malloc(numRxns*sizeof(double));
	b = (double *)malloc(numMets*sizeof(double));
	bed_double = (double *)malloc((numMets+1)*sizeof(double));
	ind_double = (double *)malloc(nnz*sizeof(double));
	val =  (double *)malloc(nnz*sizeof(double));
	bed = (int *)malloc((numMets+1)*sizeof(int));
	ind = (int *)malloc(nnz*sizeof(int));
	sense = (char *)malloc(numMets*sizeof(char));


	memcpy(lb,mxGetPr(plb),numRxns*sizeof(double));
	memcpy(ub,mxGetPr(pub),numRxns*sizeof(double));
	memcpy(c,mxGetPr(pc),numRxns*sizeof(double));
	memcpy(b,mxGetPr(pb),numMets*sizeof(double));
	memcpy(bed_double,mxGetPr(pbed),(numMets+1)*sizeof(double));
	memcpy(ind_double,mxGetPr(pind),nnz*sizeof(double));
	memcpy(val,mxGetPr(pval),nnz*sizeof(double));
	for (i = 0;i<numMets+1 ;i++ )
	{
		bed[i] = floor(bed_double[i] + 0.001);
	}
	for (i = 0;i<nnz ;i++ )
	{
		ind[i] = floor(ind_double[i] + 0.001);
	}
	for (i =0;i<numMets ;i++ )
	{
		sense[i] = GRB_EQUAL;
	}

	error = GRBnewmodel(env, &model, "metModel", numRxns, c, lb, ub, NULL, NULL);
	error = GRBaddconstrs(model,numMets,nnz,bed,ind,val,sense,b,NULL);
	error = GRBupdatemodel(model);

	double *minFlux,*maxFlux;
	minFlux  = (double *)malloc((numRxns+1)*sizeof(double));
	maxFlux  = (double *)malloc((numRxns+1)*sizeof(double));
	for (i=0;i<numRxns ;i++ )
	{
		error = GRBsetdblattrelement(model, GRB_DBL_ATTR_OBJ, i, 1.0);
		error = GRBsetintattr(model, GRB_INT_ATTR_MODELSENSE, GRB_MAXIMIZE);
		error = GRBoptimize(model);
		error = GRBgetdblattr(model, GRB_DBL_ATTR_OBJVAL, maxFlux+i);

	    error = GRBsetintattr(model, GRB_INT_ATTR_MODELSENSE, GRB_MINIMIZE);
		error = GRBoptimize(model);
		error = GRBgetdblattr(model, GRB_DBL_ATTR_OBJVAL, minFlux+i);

		error = GRBsetdblattrelement(model, GRB_DBL_ATTR_OBJ, i, 0.0);
		//glp_set_obj_coef(P,i+1,1.0);
		//glp_simplex(P, &parm);
		//minFlux[i] = glp_get_obj_val(P);

		//glp_set_obj_coef(P,i+1,-1.0);
		//glp_simplex(P, &parm);
		//maxFlux[i] = -1.0*glp_get_obj_val(P);
		//glp_set_obj_coef(P,i+1,0);
	}


	pmatout = matOpen(outfile, "w");
	pa1 = mxCreateDoubleMatrix(numRxns,1,mxREAL);// matlab will transposes
	pa2 = mxCreateDoubleMatrix(numRxns,1,mxREAL);// matlab will transposes
	memcpy((void *)(mxGetPr(pa1)), (void *)minFlux, numRxns*sizeof(double));
	memcpy((void *)(mxGetPr(pa2)), (void *)maxFlux, numRxns*sizeof(double));
	status = matPutVariable(pmatout, "minFlux", pa1);
	status = matPutVariable(pmatout, "maxFlux", pa2);
	matClose(pmat);
	matClose(pmatout);
    
    /*fout = fopen(outputfile,"wb");
	for (i=0;i<numRxns ;i++ )
	{
		fprintf(fout,"%.6f\t%.6f\n",minFlux[i],maxFlux[i]);
	}
	fclose(fout);
	*/

	finish = clock();
    duration = (double)(finish - start) / CLOCKS_PER_SEC;
    printf( "\nFastMM: Runing FVA spend: %f seconds\n", duration );
	return 0;
}


