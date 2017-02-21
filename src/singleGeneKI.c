#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <glpk.h>
#include <time.h>
#include <FastKO.h>

int help()
{
	printf("\n---------------------FastMM:singleGeneKI-----------------------\n");
	printf("\n---------------------  Version 0.01   -----------------------\n");
    printf("        C version of singleGene knockout analysis   \n");
    printf("Written By Dr Gong-Hua Li.\n");
    printf("State Key Laboratory of Genetic Resources and Evolution,\n");
    printf("  Kunming Institute  of Zoology, Chinese Academy of Sciences,\n");
    printf("  32, Eastern Jiaochang Road, Kunming, Yunnan 650223, China.\n");
    printf("Please mail to: ligonghua@mail.kiz.ac.cn\n");
    printf("This programm is under GNU GPL, you can freely use,modify and\n");
    printf("  distribute. \n");
    printf("--------------------------------------------------\n\n");
	printf("Uasage: singleGeneKO\n");
	printf("        -m cobramodel(require)\n");
	printf("        -t max or min (require)\n");
	printf("        -o outputfile (require)\n");
	printf("        -f objective function file (option,default is biomass_reaction)\n");
	printf("        -c constraintfile (option)\n");
	printf("eg: singleGeneKI -m cobramodel -t max -o outputfile\n\n");

	return 0;
}


int main(int argc,char **argv)
{ 
	int i,j,k,m,n,s,num_gen,num_rxns,*is_expr, *changedRxns,contrlputs =0;
	int numRxns, numMets,*coltype,is_constraint=0,is_objf=0;
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
	    else if (strcmp(argv[i],"-f")==0) {strcpy(objfile,argv[i+1]);is_objf = 1;}
		else if (strcmp(argv[i],"-t")==0) {strcpy(optimizeType,argv[i+1]);contrlputs++;}
		else if (strcmp(argv[i],"-o")==0) {strcpy(outputfile,argv[i+1]);contrlputs++;}
		else if (strcmp(argv[i],"-c")==0) {strcpy(constraintfile,argv[i+1]);is_constraint++;}

		else continue;
    }
	if (contrlputs<3)
	{
		help();
		exit(1);
	}

	if (strcmp(optimizeType,"min")==0 || strcmp(optimizeType,"MIN")==0) {type = 1.0;}
	else if (strcmp(optimizeType,"max")==0 || strcmp(optimizeType,"MAX")==0) {type = -1.0;}
    else {printf("-t option: the optimizeType must be min or max");help();exit(1);}

    //read gene stats
	fin = fopen(inModelFile->gens,"rb");
	if (fin==NULL){printf("\nError in open modelGeneRulesFile\n\n");exit(1);}
	i = 0;
	while(fgets(buf, 2048, fin) != NULL)
	{
		sscanf(buf+6,"%d",is_expr+i);
		i = i+1;
	}
	fclose(fin);
	num_gen = i;

	P = glp_create_prob();
	glp_read_mps(P, GLP_MPS_DECK, NULL, inModelFile->mps);
	glp_adv_basis(P, 0);
    glp_init_smcp(&parm);
    parm.msg_lev = GLP_MSG_OFF;
	numRxns = glp_get_num_cols(P);
	numMets = glp_get_num_rows(P);
	

	//read objective function file
	int *obj,num_obj;
	obj = (int *)malloc((numRxns+1)*sizeof(int));//maximum number of objective function is numRxns.

	if (is_objf ==1)
	{
		fin = fopen(objfile,"rb");
	}else{
		fin = fopen(inModelFile->objf,"rb");
	}
	if (fin==NULL){printf("\nError in open objectiveFile\n\n");exit(1);}
	i = 0;
	while(fgets(buf, 2048, fin) != NULL)
	{
		sscanf(buf,"%d",obj+i);
		i = i+1;
	}
	fclose(fin);
	num_obj = i;

	printf("FastMM: Number of objective function is %d\n",num_obj);

	lb = (double *)malloc((numRxns+1)*sizeof(double));
	ub = (double *)malloc((numRxns+1)*sizeof(double));
	coltype = (int *)malloc((numRxns+1)*sizeof(int));

	for (i = 0;i<numRxns ;i++ ) 
	{
		lb[i] = glp_get_col_lb(P,i+1);
		ub[i] = glp_get_col_ub(P,i+1);
		coltype[i] = glp_get_col_type(P,i+1);
	}


    rules = malloc2c(numRxns+1,2048);
	
	m = read_rules(inModelFile->ruls,rules);



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
		printf("FastMM: Using canstraint: number of constraint is %d\n",num_constraint);

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

	int *gene_effective,num_effective_genes,**iseffectiveRxns;
	double tmp;
	gene_effective = (int *)malloc((num_gen+1)*sizeof(int));
	iseffectiveRxns = malloc2i(num_obj+1,numRxns+1);
	//iseffectiveRxns = (int *)malloc((numRxns+1)*sizeof(int));
	num_effective_genes = get_effective_genes(rules,numRxns,gene_effective,num_gen);
    printf("FastMM: Number of active genes is %d\n",num_effective_genes);

	fout = fopen(outputfile,"wb");
    //no deletion
	double *val_no_delet;
	val_no_delet = (double *)malloc((num_obj+1)*sizeof(double));
	fprintf(fout,"0\t0");
	for (i =0;i<num_obj ;i++ )
	{
		glp_set_obj_coef(P,obj[i],type);
		glp_simplex(P, &parm);
		val_no_delet[i] = type*glp_get_obj_val(P);
		fprintf(fout,"\t%.6f",val_no_delet[i]);
		for (j =0 ;j<numRxns ;j++ )
	    {
		    tmp = glp_get_col_prim(P,j+1);
		    if (ABS(tmp) <1e-9 )
		    {
			    iseffectiveRxns[i][j] = 0;
		    }else
		    {
			    iseffectiveRxns[i][j] = 1;
		    }
	    }
		glp_set_obj_coef(P,obj[i],0);
	}
	fprintf(fout,"\n");
	

  
	//double **single_delet_vals ;
	int effectiveStat;
	double relaxlb,relaxub;
	//single_delet_vals = malloc2d(num_obj+1,num_gen+1);
	//single_delet_vals = (double *)malloc((num_gen+1)*sizeof(double));

	//single_gene_deletion
	printf("FastMM:  Start single gene knock-in....\n");

	for (i =0;i<num_gen ;i++ )
	{
		fprintf(fout,"%d\t0",i+1);
		if (gene_effective[i]>0)
		{
			is_expr[i] = 0;
			n = get_changedRxns_by_gens(rules,numRxns,is_expr,changedRxns);
			is_expr[i] = 1;
			for (s = 0;s<num_obj ;s++ )
			{
				effectiveStat = 0;
				for (k=0;k<n ;k++ )
				{
					if (iseffectiveRxns[s][changedRxns[k]]==1)
					{
						effectiveStat=1;
						break;
					}
				}

				//if (effectiveStat<1)
			    //{
				    //is_expr[i] = 1;
				//    fprintf(fout,"\t%.6f",val_no_delet[s]);
				    //single_delet_vals[s][i] = val_no_delet[s];
				//    continue;
			    //}

                glp_set_obj_coef(P,obj[s],type);
				for (j = 0;j<n ;j++ ) {
				    //glp_set_col_bnds(P,changedRxns[j]+1,GLP_FX,0,0);
					relaxlb = lb[changedRxns[j]];
					relaxub = ub[changedRxns[j]];
					if (relaxlb<0){
						relaxlb = -1000.0;
					}
					if (relaxub > 0)
					{
						relaxub = 1000.0;
					}
					glp_set_col_bnds(P,changedRxns[j]+1,coltype[changedRxns[j]],relaxlb,relaxub);
			    }
			    glp_simplex(P, &parm);
			    //single_delet_vals[s][i] = type*glp_get_obj_val(P);
			    fprintf(fout,"\t%.6f",type*glp_get_obj_val(P));
			    //is_expr[i] = 1;
			    for (j = 0;j<n ;j++ ) {			
				    glp_set_col_bnds(P,changedRxns[j]+1,coltype[changedRxns[j]],lb[changedRxns[j]],ub[changedRxns[j]]);
			    }
				glp_set_obj_coef(P,obj[s],0);
			}
		}else
		{
			for (s = 0;s<num_obj ;s++ )
			{
				fprintf(fout,"\t%.6f",val_no_delet[s]);
				//single_delet_vals[s][i] = val_no_delet[s];
			}	
		}
		fprintf(fout,"\n");
	}

	printf("FastMM:  Single gene Knock-in: completed!\n\n");
	printf("\n");
	fclose(fout);
	finish = clock();
    duration = (double)(finish - start) / CLOCKS_PER_SEC;
    printf( "\nRuning whole genome gene knock-in spend: %f seconds\n", duration );
	return 0;
}

