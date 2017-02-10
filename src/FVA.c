#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <glpk.h>
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
	int i,j,k,m,n,s,num_gen,num_rxns,*is_expr, *changedRxns,contrlputs =0;
	int numRxns, numMets,*coltype,is_constraint=0;
	double duration,type;
	double *lb,*ub;
	char buf[2048],objfile[300],outputfile[300],optimizeType[100],constraintfile[300];
	char **rules;
    clock_t start, finish;
	start = clock();

	glp_prob *P;
    glp_smcp parm;

	CFluxFile *inModelFile,*outModelFile;
	FILE *fin, *fout;

	inModelFile = MallocCFluxFile(300);
	outModelFile = MallocCFluxFile(300);
	is_expr = (int *)malloc(90000*sizeof(int));
	changedRxns = (int *)malloc(90000*sizeof(int));
	
	for (i=1;i<argc-1 ;i++ )
	{
		if (strcmp(argv[i],"-m")==0) {get_model_des(inModelFile,argv[i+1]);contrlputs++;}
		else if (strcmp(argv[i],"-o")==0) {strcpy(outputfile,argv[i+1]);contrlputs++;}
		else if (strcmp(argv[i],"-c")==0) {strcpy(constraintfile,argv[i+1]);is_constraint++;}

		else continue;
    }

    //read gene stats
	P = glp_create_prob();
	glp_read_mps(P, GLP_MPS_DECK, NULL, inModelFile->mps);
	glp_adv_basis(P, 0);
    glp_init_smcp(&parm);
    parm.msg_lev = GLP_MSG_OFF;
	numRxns = glp_get_num_cols(P);
	numMets = glp_get_num_rows(P);

	lb = (double *)malloc((numRxns+1)*sizeof(double));
	ub = (double *)malloc((numRxns+1)*sizeof(double));
	coltype = (int *)malloc((numRxns+1)*sizeof(int));

	for (i = 0;i<numRxns ;i++ ) 
	{
		lb[i] = glp_get_col_lb(P,i+1);
		ub[i] = glp_get_col_ub(P,i+1);
		coltype[i] = glp_get_col_type(P,i+1);
	}

	//read constraintfile
	int *constraint,num_constraint=0;
	char cts[300];
	double *constraint_cutoff;
	constraint = (int *)malloc((numRxns+1)*sizeof(int));
	constraint_cutoff = (double *)malloc((numRxns+1)*sizeof(double));
	if (is_constraint>0)
	{
		fin = fopen(constraintfile,"rb");
		if (fin==NULL){printf("\nError in open constraintfile\n\n");exit(1);}
		i = 0;
		while(fgets(buf, 2048, fin) != NULL)
	    {
		    k = sscanf(buf,"%d%[ ,\t]%lf",constraint+i,cts,constraint_cutoff+i);
			printf("%d\t%f\n",constraint[i],constraint_cutoff[i]);
			if (k<3)
			{
				printf("\nError in read constraintfile\n\n");
				exit(1);
			}
		    i = i+1;
	    }
		num_constraint = i;
	    fclose(fin);
		printf("FastKO: FVA Using canstraint: number of constraint is %d\n",num_constraint);

		for (i = 0;i<num_constraint ;i++ )
		{
			glp_set_obj_coef(P,constraint[i],-1.0);
			glp_simplex(P, &parm);
			constraint_cutoff[i] = -1.0*glp_get_obj_val(P)*constraint_cutoff[i];
			glp_set_obj_coef(P,constraint[i],0);

		}
		for (i = 0;i<num_constraint ;i++ )
		{
			if (constraint_cutoff[i]<ub[constraint[i]-1])
			{
				glp_set_col_bnds(P,constraint[i],coltype[constraint[i]-1],constraint_cutoff[i],ub[constraint[i]-1]);

			}else
			{
				glp_set_col_bnds(P,constraint[i],GLP_FX,constraint_cutoff[i],constraint_cutoff[i]);
			}
		}
	}

	double *minFlux,*maxFlux;
	minFlux  = (double *)malloc((numRxns+1)*sizeof(double));
	maxFlux  = (double *)malloc((numRxns+1)*sizeof(double));
	for (i=0;i<numRxns ;i++ )
	{
		glp_set_obj_coef(P,i+1,1.0);
		glp_simplex(P, &parm);
		minFlux[i] = glp_get_obj_val(P);

		glp_set_obj_coef(P,i+1,-1.0);
		glp_simplex(P, &parm);
		maxFlux[i] = -1.0*glp_get_obj_val(P);
		glp_set_obj_coef(P,i+1,0);
	}
    
    fout = fopen(outputfile,"wb");
	for (i=0;i<numRxns ;i++ )
	{
		fprintf(fout,"%.6f\t%.6f\n",minFlux[i],maxFlux[i]);
	}
	fclose(fout);

	finish = clock();
    duration = (double)(finish - start) / CLOCKS_PER_SEC;
    printf( "\nFastKO: Runing FVA spend: %f seconds\n", duration );
	return 0;
}


