#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <ilcplex/cplex.h>
#include "matio.h"
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

int get_varLen(matvar_t *matvar){
	int *size;
	size = (int *)malloc(2*sizeof(int));
	
	memcpy(size,matvar->dims,2*sizeof(int));
	if(size[0] > size[1]){
		return(size[0]);
	}else{
		return(size[1]);
	}
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
    
	CPXENVptr     env = NULL;
    CPXLPptr      model = NULL;
	//GRBenv   *env   = NULL,*modelenv = NULL;
	//GRBmodel *model = NULL;
	//FILE *fin, *fout;

	env = CPXopenCPLEX (&status);
	status = CPXsetintparam (env, CPXPARAM_ScreenOutput, CPX_OFF);
	//status = CPXsetintparam (env, CPXPARAM_Read_DataCheck, CPX_ON);

	model = CPXcreateprob (env, &status, "cplex");
	

	//error = GRBloadenv(&env, "");
	//error = GRBsetintparam(env, GRB_INT_PAR_LOGTOCONSOLE, 0);
	//error = GRBsetintparam(env, "Method", 0);


	//read matfile
	mat_t *pmat = NULL,*pmatout = NULL; //mat file
	matvar_t *plb = NULL,*pub = NULL,*pb = NULL,*pc = NULL;
	matvar_t *pbed = NULL,*pind = NULL,*pval = NULL,*pa1 = NULL,*pa2 = NULL;  //mat array
	
	strcpy(infile,argv[1]);
	strcpy(outfile,argv[2]);
	pmat = Mat_Open(infile, MAT_ACC_RDONLY);
    if (pmat == NULL) {
        printf("Error reopening file %s\n", infile);
        return 0;
    }

	/* mat file contains: lb, ub, b, c, (bed,ind,val)*/

	plb = Mat_VarRead(pmat,"lb");
	pub = Mat_VarRead(pmat,"ub");
	pb = Mat_VarRead(pmat,"b");
	pc = Mat_VarRead(pmat,"c");
	pbed = Mat_VarRead(pmat,"bed");
	pind = Mat_VarRead(pmat,"ind");
	pval = Mat_VarRead(pmat,"val");

	numRxns = plb->dims[0];
	numMets = pb->dims[0];
	nnz = pind->dims[0];

	double *lb,*ub,*c,*b,*bed_double,*ind_double,*val;
	int *bed,*ind,*indices;
	char *sense,*lu;
	lb = (double *)malloc(numRxns*sizeof(double));
	ub = (double *)malloc(numRxns*sizeof(double));
	c = (double *)malloc(numRxns*sizeof(double));
	b = (double *)malloc(numMets*sizeof(double));
	bed_double = (double *)malloc((numMets+1)*sizeof(double));
	ind_double = (double *)malloc(nnz*sizeof(double));
	val =  (double *)malloc(nnz*sizeof(double));
	bed = (int *)malloc((numMets+1)*sizeof(int));
	ind = (int *)malloc(nnz*sizeof(int));
	indices = (int *)malloc((numRxns+1)*sizeof(int));
	sense = (char *)malloc(numMets*sizeof(char));
	lu = (char *)malloc(numRxns*sizeof(char));


	memcpy(lb,plb->data,numRxns*sizeof(double));
	memcpy(ub,pub->data,numRxns*sizeof(double));
	memcpy(c,pc->data,numRxns*sizeof(double));
	memcpy(b,pb->data,numMets*sizeof(double));
	memcpy(bed_double,pbed->data,(numMets+1)*sizeof(double));
	memcpy(ind_double,pind->data,nnz*sizeof(double));
	memcpy(val,pval->data,nnz*sizeof(double));
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
		sense[i] = 'E';
	}

    // set coeff and rhs
    status = CPXaddrows (env, model, numRxns, numMets, nnz, b,
                       sense, bed, ind, val, NULL, NULL);
	//set lb ub
	for (i = 0;i<numRxns ;i++ )
	{
		indices[i] = i;
		lu[i] = 'L';
	}
	status = CPXchgbds (env, model, numRxns, indices, lu, lb);
	for (i = 0;i<numRxns ;i++ ) lu[i] = 'U';
	status = CPXchgbds (env, model, numRxns, indices, lu, ub);

	//set method
	status = CPXsetintparam (env, CPXPARAM_LPMethod, CPX_ALG_PRIMAL);

	//CPXaddrows (env, lp, ccnt, rcnt, nzcnt, rhs,
	//      sense, rmatbeg, rmatind, rmatval,
	//      newcolname, newrowname);
	//error = GRBnewmodel(env, &model, "metModel", numRxns, c, lb, ub, NULL, NULL);
	//error = GRBaddconstrs(model,numMets,nnz,bed,ind,val,sense,b,NULL);
	//error = GRBupdatemodel(model);

	double *minFlux,*maxFlux;
	double tmpval_1[2] = {1,1};
	double tmpval_0[2] = {0,0};
	size_t dims1[2],dims2[2];
	minFlux  = (double *)malloc((numRxns+1)*sizeof(double));
	maxFlux  = (double *)malloc((numRxns+1)*sizeof(double));
	for (i=0;i<numRxns ;i++ )
	{
		status = CPXchgobj (env, model, 1, indices+i, tmpval_1);
		CPXchgobjsen (env, model, CPX_MAX);
		status = CPXlpopt (env, model);
		status = CPXgetobjval (env, model, maxFlux+i);

		CPXchgobjsen (env, model, CPX_MIN);
		status = CPXlpopt (env, model);
		status = CPXgetobjval (env, model, minFlux+i);

		status = CPXchgobj (env, model, 1, indices+i, tmpval_0);
		//glp_set_obj_coef(P,i+1,1.0);
		//glp_simplex(P, &parm);
		//minFlux[i] = glp_get_obj_val(P);

		//glp_set_obj_coef(P,i+1,-1.0);
		//glp_simplex(P, &parm);
		//maxFlux[i] = -1.0*glp_get_obj_val(P);
		//glp_set_obj_coef(P,i+1,0);
	}


	//pmatout =Mat_CreateVer(outfile,NULL,MAT_FT_DEFAULT);
	//pa1 = mxCreateDoubleMatrix(numRxns,1,mxREAL);// matlab will transposes
	//pa2 = mxCreateDoubleMatrix(numRxns,1,mxREAL);// matlab will transposes
	dims1[0] = numRxns;
	dims1[1] = 1;

	pmatout =Mat_CreateVer(outfile,NULL,MAT_FT_DEFAULT);
	pa1 = Mat_VarCreate("minFlux",MAT_C_DOUBLE,MAT_T_DOUBLE,2,dims1,minFlux,0);
	pa2 = Mat_VarCreate("maxFlux",MAT_C_DOUBLE,MAT_T_DOUBLE,2,dims1,maxFlux,0);
	Mat_VarWrite(pmatout,pa1,MAT_COMPRESSION_NONE);
	Mat_VarWrite(pmatout,pa2,MAT_COMPRESSION_NONE);



	//memcpy((void *)(mxGetPr(pa1)), (void *)minFlux, numRxns*sizeof(double));
	//memcpy((void *)(mxGetPr(pa2)), (void *)maxFlux, numRxns*sizeof(double));
	//status = matPutVariable(pmatout, "minFlux", pa1);
	//status = matPutVariable(pmatout, "maxFlux", pa2);
	Mat_Close(pmat);
	Mat_Close(pmatout);
    
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


