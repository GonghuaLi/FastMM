#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <ilcplex/cplex.h>
#include "mat.h"
//#include <glpk.h>
#include <time.h>
#include <FastKO.h>

int help()
{
	printf("\n---------------------FastMM:singleGeneKO-----------------------\n");
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
	printf("        -f objective function file (option,defent is biomass_reaction)\n");
	printf("        -c constraintfile (option)\n");
	printf("        -e effective_rxnsfile (optional,default not use this file)\n");
	printf("eg: singleGeneKO -m cobramodel -t max -o outputfile\n\n");

	return 0;
}


int main(int argc,char **argv)
{ 
	int i,j,k,m,n,s,num_gen,num_rxns,*is_expr, *changedRxns,contrlputs =0,indices[1];
	int numRxns, error, status,numMets,nmets,nrxns,*coltype,is_constraint=0,is_objf=0,is_effective_file = 0;
	double duration,type,tmpval;
	double *lb,*ub;
	char buf[2048],objfile[300],outputfile[300],optimizeType[100],constraintfile[300],effectfile[300];
	char **rules;

	double tmpval_1[2] = {1,1};
	double tmpval_0[2] = {0,0};
	char lu_l[2] = {'L','L'};
	char lu_u[2] = {'U','U'};
	char lu_b[2] = {'B','B'};

    clock_t start, finish;
	start = clock();

	//glp_prob *P;
    //glp_smcp parm;
	//GRBenv   *env   = NULL,*modelenv = NULL;
	//GRBmodel *model = NULL;

	CPXENVptr     env = NULL;
    CPXLPptr      model = NULL;
	env = CPXopenCPLEX (&status);
	status = CPXsetintparam (env, CPXPARAM_ScreenOutput, CPX_OFF);
	//status = CPXsetintparam (env, CPXPARAM_Read_DataCheck, CPX_ON);

	model = CPXcreateprob (env, &status, "cplex");



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
		else if (strcmp(argv[i],"-e")==0) {strcpy(effectfile,argv[i+1]);is_effective_file++;}

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

	//P = glp_create_prob();
	//glp_read_mps(P, GLP_MPS_DECK, NULL, inModelFile->mps);
	//glp_adv_basis(P, 0);
    //glp_init_smcp(&parm);
    //parm.msg_lev = GLP_MSG_OFF;
		//set model env
	//error = GRBloadenv(&env, "");
	//error = GRBsetintparam(env, GRB_INT_PAR_LOGTOCONSOLE, 0);
	//error = GRBsetintparam(env, "Method", 0);
	//error = GRBreadmodel(env, inModelFile->mps, &model);
    printf("%s\n",inModelFile->mps);
	error = CPXreadcopyprob (env, model, inModelFile->mps, NULL);



	if (error){printf("Error reading model\n\n");exit(1);}
	//error = GRBsetdblparam(GRBgetenv(model), "LogToConsole", 0);
    	
	//error = GRBgetintattr(model,GRB_INT_ATTR_NUMCONSTRS, &numMets);
	//error = GRBgetintattr(model,GRB_INT_ATTR_NUMVARS, &numRxns);
	numRxns =  CPXgetnumcols(env, model);
	numMets = CPXgetnumrows(env, model);
	//error = GRBgetintattr(model,GRB_INT_ATTR_NUMNZS, &Snnz);
	for (i=0;i<numRxns ;i++ )
	{
		indices[0] = i;
		status = CPXchgobj (env, model, 1, indices, tmpval_0);
	}

	//numRxns = glp_get_num_cols(P);
	//numMets = glp_get_num_rows(P);
	

	//read objective function file
	int *obj,num_obj,*input_effectRxns;
	obj = (int *)malloc((numRxns+1)*sizeof(int));//maximum number of objective function is numRxns.
	input_effectRxns = (int *)malloc((numRxns+1)*sizeof(int));//maximum number of objective function is numRxns.

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

	if (is_effective_file ==1)
	{
		fin = fopen(effectfile,"rb");
	}
	if (fin==NULL){printf("\nError in open is_effective_file\n\n");exit(1);}
	i = 0;
	while(fgets(buf, 2048, fin) != NULL)
	{
		sscanf(buf,"%d",input_effectRxns+i);
		i = i+1;
	}
	fclose(fin);

	lb = (double *)malloc((numRxns+1)*sizeof(double));
	ub = (double *)malloc((numRxns+1)*sizeof(double));
	coltype = (int *)malloc((numRxns+1)*sizeof(int));

	//error = GRBgetdblattrarray(model, GRB_DBL_ATTR_LB,0,numRxns,lb);
	//error = GRBgetdblattrarray(model, GRB_DBL_ATTR_UB,0,numRxns,ub);
	//error = GRBgetdblattrarray(model, GRB_CHAR_ATTR_VTYPE,0,numRxns,coltype);

	error = CPXgetlb (env, model, lb, 0, numRxns-1);
	error = CPXgetub (env, model, ub, 0, numRxns-1);

	//for (i = 0;i<numRxns ;i++ ) 
	//{
	//	lb[i] = glp_get_col_lb(P,i+1);
	//	ub[i] = glp_get_col_ub(P,i+1);
	//	coltype[i] = glp_get_col_type(P,i+1);
	//}


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
			   indices[0] = constraint[i]-1;
		       error = CPXchgobj(env, model, 1, indices, tmpval_1);
		       CPXchgobjsen(env, model, CPX_MAX);
			   error = CPXlpopt(env, model);
		       error = CPXgetobjval(env, model, &tmpval);
			//error = GRBsetdblattrelement(model, GRB_DBL_ATTR_OBJ, constraint[i]-1, 1.0);
		    //error = GRBsetintattr(model, GRB_INT_ATTR_MODELSENSE, GRB_MAXIMIZE);
			//error = GRBoptimize(model);
			//error = GRBgetdblattr(model, GRB_DBL_ATTR_OBJVAL, &tmpval);
			//glp_set_obj_coef(P,constraint[i],-1.0);
			//glp_simplex(P, &parm);
			constraint_cutoff[i] = 1.0*tmpval*constraint_cutoff[i];
			//error = GRBsetdblattrelement(model, GRB_DBL_ATTR_OBJ, constraint[i]-1, 0.0);
			error = CPXchgobj(env, model, 1, indices, tmpval_0);
			//glp_set_obj_coef(P,constraint[i],0);

		}
		for (i = 0;i<num_constraint ;i++ )
		{
			if (constraint_cutoff[i]<ub[constraint[i]-1])
			{
				indices[0] = constraint[i]-1;
				//GRBsetdblattrelement(model, GRB_DBL_ATTR_UB, constraint[i]-1, constraint_cutoff[i]);
				error = CPXchgbds (env, model, 1, indices, lu_u, constraint_cutoff+i);
				//glp_set_col_bnds(P,constraint[i],coltype[constraint[i]-1],constraint_cutoff[i],ub[constraint[i]-1]);

			}//else
			//{
			//	glp_set_col_bnds(P,constraint[i],GLP_FX,constraint_cutoff[i],constraint_cutoff[i]);
			//}
		}
	}
	

	int *gene_effective,num_effective_genes,**iseffectiveRxns;
	double tmp,*xval,tvalue[1];
	gene_effective = (int *)malloc((num_gen+1)*sizeof(int));
	iseffectiveRxns = malloc2i(num_obj+1,numRxns+1);
	//iseffectiveRxns = (int *)malloc((numRxns+1)*sizeof(int));

	if (num_obj<2 && is_effective_file==1)
    {
		num_effective_genes = get_effective_genes_from_effectrxns(rules,numRxns,gene_effective,num_gen,input_effectRxns);
    }else{
		num_effective_genes = get_effective_genes(rules,numRxns,gene_effective,num_gen);
	}
    printf("FastMM: Number of active genes is %d\n",num_effective_genes);

	fout = fopen(outputfile,"wb");
    //no deletion
	double *val_no_delet;
	val_no_delet = (double *)malloc((num_obj+1)*sizeof(double));
	xval = (double *)malloc((numRxns+1)*sizeof(double));
	fprintf(fout,"0\t0");
	for (i =0;i<num_obj ;i++ )
	{
		//glp_set_obj_coef(P,obj[i],type);
		//error = GRBsetdblattrelement(model, GRB_DBL_ATTR_OBJ, obj[i]-1, 1.0);
		//error = GRBsetintattr(model, GRB_INT_ATTR_MODELSENSE, type);
		//glp_simplex(P, &parm);
		//error = GRBoptimize(model);
		//error = GRBgetdblattr(model, GRB_DBL_ATTR_OBJVAL, &tmpval);
		indices[0] = obj[i]-1;
		error = CPXchgobj (env, model, 1, indices, tmpval_1);
		//CPXchgobjsen (env, model, CPX_MAX);
		CPXchgobjsen (env, model, type);
		error = CPXlpopt (env, model);
		error = CPXgetobjval (env, model, &tmpval);

		printf("%.6f\n",tmpval);
		//val_no_delet[i] = type*tmpval;
		val_no_delet[i] = tmpval;
		fprintf(fout,"\t%.6f",val_no_delet[i]);
		//error = GRBgetdblattrarray(model, GRB_DBL_ATTR_X, 0, numRxns, xval);
		error = CPXgetx (env, model, xval, 0, numRxns-1);
		if (i==0 && is_effective_file==1)
		{
			for (j =0 ;j<numRxns ;j++ )
			{
			    iseffectiveRxns[i][j] = input_effectRxns[j];
			}
		}else{
		    for (j =0 ;j<numRxns ;j++ )
	        {
			//error = GRBsetdblattrelement(model, GRB_DBL_ATTR_X,j,tmp);
			//printf("%i    %.6f\n",j,xval[j]);
		    //tmp = glp_get_col_prim(P,j+1);
		        if (ABS(xval[j]) <1e-9 )
		        {
			        iseffectiveRxns[i][j] = 0;
		        }else
		        {
				    //printf("%j    %.6f\n",j,tmpval);
			        iseffectiveRxns[i][j] = 1;
		        }
	        }
		}
		//glp_set_obj_coef(P,obj[i],0);
		//error = GRBsetdblattrelement(model, GRB_DBL_ATTR_OBJ, obj[i]-1, 0);
		status = CPXchgobj (env, model, 1, indices+i, tmpval_0);

	}
	fprintf(fout,"\n");
	

  
	//double **single_delet_vals ;
	int effectiveStat,thischange;


	//single_delet_vals = malloc2d(num_obj+1,num_gen+1);
	//single_delet_vals = (double *)malloc((num_gen+1)*sizeof(double));

	//single_gene_deletion
	printf("FastMM:  Start single gene deletion....\n");

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

				if (effectiveStat<1)
			    {
				    //is_expr[i] = 1;
				    fprintf(fout,"\t%.6f",val_no_delet[s]);
				    //single_delet_vals[s][i] = val_no_delet[s];
				    continue;
			    }

                //glp_set_obj_coef(P,obj[s],type);
		         indices[0] = obj[s]-1;
		        error = CPXchgobj (env, model, 1, indices, tmpval_1);
		        CPXchgobjsen (env, model, type);
		       //error = CPXlpopt (env, model);
		       //error = CPXgetobjval (env, model, &tmpval);


				//error = GRBsetdblattrelement(model, GRB_DBL_ATTR_OBJ, obj[s]-1, 1.0);
				//error = GRBsetintattr(model, GRB_INT_ATTR_MODELSENSE, type);
				
				for (j = 0;j<n ;j++ ) {
					//thischange = changedRxns[j];

					error = CPXchgbds (env, model, 1, changedRxns+j, lu_b, tmpval_0);
					//error = CPXchgbds (env, model, 1, &thischange, lu_l, tmpval_0);
					//error = CPXchgbds (env, model, 1, &thischange, lu_u, tmpval_0);
					//error = GRBsetdblattrelement(model, GRB_DBL_ATTR_LB,changedRxns[j],0);
					//error = GRBsetdblattrelement(model, GRB_DBL_ATTR_UB,changedRxns[j],0);
				    //glp_set_col_bnds(P,changedRxns[j]+1,GLP_FX,0,0);
			    }
			    //glp_simplex(P, &parm);
				//error = GRBoptimize(model);
				error = CPXlpopt (env, model);
			    //single_delet_vals[s][i] = type*glp_get_obj_val(P);
				//error = GRBgetdblattr(model, GRB_DBL_ATTR_OBJVAL, &tmpval);
				error = CPXgetobjval (env, model, &tmpval);
			    fprintf(fout,"\t%.6f",tmpval);
			    //is_expr[i] = 1;
			    for (j = 0;j<n ;j++ ) {	
					//thischange = changedRxns[j];
					error = CPXchgbds (env, model, 1, changedRxns+j, lu_l, lb+changedRxns[j]);
					error = CPXchgbds (env, model, 1, changedRxns+j, lu_u, ub+changedRxns[j]);
					//error = GRBsetdblattrelement(model, GRB_DBL_ATTR_LB,changedRxns[j],lb[changedRxns[j]]);
					//error = GRBsetdblattrelement(model, GRB_DBL_ATTR_UB,changedRxns[j],ub[changedRxns[j]]);
				    //glp_set_col_bnds(P,changedRxns[j]+1,coltype[changedRxns[j]],lb[changedRxns[j]],ub[changedRxns[j]]);
			    }
				//glp_set_obj_coef(P,obj[s],0);
				//error = GRBsetdblattrelement(model, GRB_DBL_ATTR_OBJ, obj[s]-1, 0);
				error = CPXchgobj (env, model, 1, indices, tmpval_0);
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

	printf("FastMM:  Single gene deletion: completed!\n\n");
	printf("\n");
	fclose(fout);
	finish = clock();
    duration = (double)(finish - start) / CLOCKS_PER_SEC;
    printf( "\nRuning whole genome gene deletion spend: %f seconds\n", duration );
	return 0;
}

