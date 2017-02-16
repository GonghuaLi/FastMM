#include <dirent.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#define MAX(a,b) (a>b?a:b)
#define MIN(a,b) (a>b?b:a)
#define ABS(a) ((a)>0?(a):-1.0*(a))

typedef struct 
{
	char *mps;
	char *rxns;
	char *mets;
	char *gens;
	char *ruls;
	char *objf;
	char *modelName;
	int n;
} CFluxFile;

CFluxFile* MallocCFluxFile(int n)
{
	CFluxFile *p;
	p=(CFluxFile *)malloc(sizeof(CFluxFile));
	p->mps = (char*)malloc(n*sizeof(char));
	p->rxns = (char*)malloc(n*sizeof(char));
	p->mets = (char*)malloc(n*sizeof(char));
	p->gens = (char*)malloc(n*sizeof(char));
	p->ruls = (char*)malloc(n*sizeof(char));
	p->objf = (char*)malloc(n*sizeof(char));
	p->modelName = (char*)malloc(n*sizeof(char));
	return p;
}

void   freeCFluxFile(CFluxFile *p)
{
	free(p->mps); p->mps = NULL;
	free(p->rxns); p->rxns = NULL;
	free(p->mets); p->mets = NULL;
	free(p->gens); p->gens = NULL;
	free(p->ruls); p->ruls = NULL; 
	free(p->objf); p->objf = NULL;
	free(p->modelName); p->modelName = NULL;
	free(p); p=NULL;
}


int get_model_des(CFluxFile *p,char *c)
{
	int i,n,j,k;
	char des='/',name[300];
	n = strlen(c);
	k = 0;
	for (i = n-1;i>=0 ;i-- )
	{
		if (c[i] =='/' || c[i] =='\\')
		{
			k = i+1;
			des = c[i];
			break;
		}
	}

	j=0;

	for (i = k;i<n ;i++ ) name[j++] = c[i];
	name[j] = '\0';
	for (i = 0;i<n ;i++ )
	{
		p->modelName[i] = c[i];
		p->mps[i] = c[i];
		p->gens[i] = c[i];
		p->mets[i] = c[i];
		p->objf[i] = c[i];
		p->rxns[i] = c[i];
		p->ruls[i] = c[i];
	}
	p->modelName[i] = '\0';
    p->mps[i] = des;
	p->gens[i] = des;
	p->mets[i] = des;
	p->objf[i] = des;
	p->rxns[i] = des;
	p->ruls[i] = des;

	i = i+1;

	p->mps[i] = '\0';
	p->gens[i] = '\0';
	p->mets[i] = '\0';
	p->objf[i] = '\0';
	p->rxns[i] = '\0';
	p->ruls[i] = '\0';
	strcat(p->mps, name );
	strcat(p->gens, name);
	strcat(p->mets, name);
	strcat(p->ruls, name);
	strcat(p->rxns, name);
	strcat(p->objf, name);
	strcat(p->mps, ".mps" );
	strcat(p->gens, ".gens" );
	strcat(p->mets, ".mets" );
	strcat(p->ruls, ".ruls" );
	strcat(p->rxns, ".rxns" );
	strcat(p->objf, ".objf" );
	return 0;
}

int copyfile(char *infile,char *outfile)
{
	FILE *fin,*fout;
	char ch;
	fin = fopen(infile,"rb");
	fout = fopen(outfile,"wb");
	if (fin==NULL){printf("Error in open File: %s\n",infile);exit(1);}
	while( ( ch = fgetc(fin) ) != EOF )
      fputc(ch, fout);
	return 0;
}

int  read_rules(char *filename, char **p)
{

	int i=0;
	char buf[2048];
	FILE *fin;
	fin = fopen(filename,"rb");
	if (fin==NULL){printf("\nError in open: %s \n\n",filename);exit(1);}
	while(fgets(buf, 2048, fin) != NULL)
	{
		strcpy(p[i],buf);
		i++;
	}
	fclose(fin);
	return i;
}

int get_changedRxns_by_gens(char **rules,int numRxns,int *is_expr, int *changedRxn)
{

	int i,j,n,k,r,stat,stata,geneid,ctrlBrackets;
	char opera;
	i=0;
	r = 0;
	for (k=0;k<numRxns ;k++ )
	{
		n = strlen(rules[k]);
		if (n>16)
		{
			stat = 1;
			opera = '&';
			for (j =15;j<n ;j++ )
			{
				if (rules[k][j] == '(')
				{
					ctrlBrackets = 1;
					sscanf(rules[k]+j+1,"%d",&geneid);
					stata = is_expr[geneid-1];
					continue;
				}
				if (rules[k][j] == '&' && ctrlBrackets == 1)
				{
					sscanf(rules[k]+j+2,"%d",&geneid);
					stata = MIN(stata,is_expr[geneid-1]);
					continue;
				}
				if (rules[k][j] == '|' && ctrlBrackets == 1)
				{
					sscanf(rules[k]+j+2,"%d",&geneid);
					stata = MAX(stata,is_expr[geneid-1]);
					continue;
				}
				if (rules[k][j] == '|' && ctrlBrackets == 1)
				{
					sscanf(rules[k]+j+2,"%d",&geneid);
					stata = MAX(stata,is_expr[geneid-1]);
					continue;
				}
				if (rules[k][j] == ')')
				{
					if (opera == '&') stat = MIN(stat,stata);
					if (opera == '|') stat = MAX(stat,stata);
					ctrlBrackets = 0;
					continue;
				}

				if (rules[k][j] == '|' && ctrlBrackets == 0) opera = '|';
				if (rules[k][j] == '&' && ctrlBrackets == 0) opera = '&';
			}

			if (stat == 0) 
			{
				changedRxn[r] = i;
				r = r+1;
			}
			//is_rxns[i] = stat;
		}
		i = i+1;
	}
	return r;
}

int get_changed_effectiveRxns_by_gens(char **rules,int numRxns,int *is_expr, int *changedRxn,int *iseffetiveRxns,int *stat_effective)
{

	int i,j,n,k,r,stat,stata,geneid,ctrlBrackets;
	char opera;
	stat_effective[0] = 0;
	i=0;
	r = 0;
	for (k=0;k<numRxns ;k++ )
	{
		n = strlen(rules[k]);
		if (n>16)
		{
			stat = 1;
			opera = '&';
			for (j =15;j<n ;j++ )
			{
				if (rules[k][j] == '(')
				{
					ctrlBrackets = 1;
					sscanf(rules[k]+j+1,"%d",&geneid);
					stata = is_expr[geneid-1];
					continue;
				}
				if (rules[k][j] == '&' && ctrlBrackets == 1)
				{
					sscanf(rules[k]+j+2,"%d",&geneid);
					stata = MIN(stata,is_expr[geneid-1]);
					continue;
				}
				if (rules[k][j] == '|' && ctrlBrackets == 1)
				{
					sscanf(rules[k]+j+2,"%d",&geneid);
					stata = MAX(stata,is_expr[geneid-1]);
					continue;
				}
				if (rules[k][j] == '|' && ctrlBrackets == 1)
				{
					sscanf(rules[k]+j+2,"%d",&geneid);
					stata = MAX(stata,is_expr[geneid-1]);
					continue;
				}
				if (rules[k][j] == ')')
				{
					if (opera == '&') stat = MIN(stat,stata);
					if (opera == '|') stat = MAX(stat,stata);
					ctrlBrackets = 0;
					continue;
				}

				if (rules[k][j] == '|' && ctrlBrackets == 0) opera = '|';
				if (rules[k][j] == '&' && ctrlBrackets == 0) opera = '&';
			}

			if (stat == 0) 
			{
				changedRxn[r] = i;
				if (iseffetiveRxns[i] == 1)
				{
					stat_effective[0] = stat_effective[0] + 1;
				}
				r = r+1;
			}
			//is_rxns[i] = stat;
		}
		i = i+1;
	}
	return r;
}

int get_effective_genes(char **rules,int numRxns,int *gene_effective,int num_genes)
{

	int i,j,n,k,r,stat,stata,geneid,ctrlBrackets;
	char opera;
    for (i = 0;i<num_genes ;i++ )
		gene_effective[i] = 0;
	i=0;
	for (k=0;k<numRxns ;k++ )
	{
		n = strlen(rules[k]);
		if (n>16)
		{
			stat = 1;
			opera = '&';
			for (j =15;j<n ;j++ )
			{
				if (rules[k][j] == '(')
				{
					ctrlBrackets = 1;
					sscanf(rules[k]+j+1,"%d",&geneid);
					gene_effective[geneid-1]=1;
				}
			}
		}
	}
	r = 0;
	for (i=0;i<num_genes ;i++ )
	{
		if(gene_effective[i]>0) r++;
	}
	return r;
}

int get_changedRxns_by_one_met(int *IA,double *VA, int lenA, double *lb,double *ub, int *tmprun,int *changedRxn)
{
	int i, j,k;
	k = 0;
	for (i =1 ;i<lenA+1 ;i++ )
	{
		if ((lb[IA[i]-1]>-1e-6 && VA[i] < -1e-6) || (ub[IA[i]-1]<1e-6 && VA[i] > 1e-6))
		{
			if (tmprun[IA[i]-1]==0)
			{
				changedRxn[k] = IA[i]-1;
				k=k+1;
				tmprun[IA[i]-1] = 1;
			}
		}
	}
	for (i = 0;i<k ;i++ )
	{
		tmprun[changedRxn[i]] = 0;
	}
	return k;
}

int get_changedRxns_by_two_mets(int *IA,double *VA, int lenA, int *IB,double *VB, int lenB,double *lb,double *ub, int *tmprun,int *changedRxn)
{
	int i, j,k;
	k = 0;
	for (i =1 ;i<lenA+1 ;i++ )
	{
		if ((lb[IA[i]-1]>-1e-6 && VA[i] < -1e-6) || (ub[IA[i]-1]<1e-6 && VA[i] > 1e-6))
		{
			if (tmprun[IA[i]-1]==0)
			{
				changedRxn[k] = IA[i]-1;
				k=k+1;
				tmprun[IA[i]-1] = 1;
			}
		}
	}
	for (i =1 ;i<lenB+1 ;i++ )
	{
		if ((lb[IB[i]-1]>-1e-6 && VB[i] < -1e-6) || (ub[IB[i]-1]<1e-6 && VB[i] > 1e-6))
		{
			if (tmprun[IB[i]-1]==0)
			{
				changedRxn[k] = IB[i]-1;
				k=k+1;
				tmprun[IB[i]-1] = 1;
			}
		}
	}
	for (i = 0;i<k ;i++ )
	{
		tmprun[changedRxn[i]] = 0;
	}
	return k;
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


int** malloc2i(int row,int colum) //malloc 2-dimension int
{
	int **p;
	int i;
	p = (int **)malloc(row*sizeof(int *));
	for (i =0;i<row ;i++ ) p[i]= (int *)malloc(colum*sizeof(int));
	return(p);
}
void free2i(int **p,int row)
{
	int i;
	for (i =0;i<row ;i++ ) free(p[i]);
	free(p);
}

char** malloc2c(int row,int colum) //malloc 2-dimension int
{
	char **p;
	int i;
	p = (char **)malloc(row*sizeof(char *));
	for (i =0;i<row ;i++ ) p[i]= (char *)malloc(colum*sizeof(char));
	return(p);
}

void free2c(char **p,int row)
{
	int i;
	for (i =0;i<row ;i++ ) free(p[i]);
	free(p);
}

