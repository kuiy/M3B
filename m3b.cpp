

///////////////////////////////////////////////////////////////////////////////
////////////////////////////////////includes///////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <math.h>
#include <time.h>

#include "chisq.h"

///////////////////////////////////////////////////////////////////////////////
//////////////////////////////////defines//////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

#define MAX_LINE_LENGTH 1000000
#define MAX_N_STATES 15
#define MAX_STATE_LENGTH 100

#define N_CASES_PER_DF 5

///////////////////////////////////////////////////////////////////////////////
////////////////////////////////global definitions/////////////////////////////
///////////////////////////////////////////////////////////////////////////////

int n_cases, n_vars, max_max_cond_size;

double alpha;

clock_t clock_init;

int *n_states;

char *data_file;

int **data;



int ci_times=0;//record times of independency tests

int answerList[5000]={0};


///////////////////////////////////////////////////////////////////////////////
///////////////////////////////functions' prototypes///////////////////////////
///////////////////////////////////////////////////////////////////////////////

void parse_data(void);//read data
int state_index(int var, char *state, char ***states);

int next_cond_index(int n_pc, int cond_size, int *cond_index);
double compute_dep(int var, int target, int *cond);


FILE *fp_mb;
FILE *fp_pc;
FILE *fp_times;

void adjv_superset(int target, int *pc,int **sep);
double new_subset_dep(int var, int target, int *pc, int *order,int *sep2, int subset_size);


void m3b(int targetNode);
void RecSearch(int oldNode, int targetNode, int *answerList, int deep,
                     int *AdjForOldNode, int** SepForOldNode);
void printAnswerList(int *answerList, int length);
bool inAnswerList(int *answerList, int target, int length);
bool isVStruct (int oldNode, int targetNode, int AdjNode, int* AdjForOldNode,
                int* SepForOldNodeAndAdjNode);
bool inSepSet(int target, int *SepForOldNodeAndAdjNode);



///////////////////////////////////////////////////////////////////////////////
////////////////////////////////functions' body////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

int main(int argc, char* argv[])

{


int i, target;


if(argc!=6)
	{
	 printf("\n\n Usage: data_file n_cases n_vars target(-1=all) alpha. \n\n");
	 exit(0);
	}



srand((unsigned)time(NULL));

clock_init=clock();




data_file=argv[1];
n_cases=atoi(argv[2]);
n_vars=atoi(argv[3]);
target=atoi(argv[4]);
alpha=atof(argv[5]);


n_states=new int[n_vars];

data=new int*[n_cases];
for(i=0;i<n_cases;i++)
    data[i]=new int[n_vars];




printf("starting  to read data.....\n");


parse_data();//read data sets

printf("read data is over...............\n");
printf("starting  to preform algoirhtms.....\n");


         m3b(target);

        printf("\n\n ci_times=%d",ci_times);


fp_times=fopen("ci_time.txt","a+");

fprintf(fp_times,"%d ",ci_times);

fprintf(fp_times,"\n");
fclose(fp_times);

printf(" time:%f\n",(double)(clock()-clock_init)/CLOCKS_PER_SEC);

return 0;

}



////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////


void parse_data(void)

{int i, j, min_n_states;
 char buffer[(int)MAX_LINE_LENGTH];
 char *token;
 const char *separators=" ,\t\n\r";
 FILE *f_in;
 char ***states;


states=new char**[n_vars];
for(i=0;i<n_vars;i++)
	{states[i]=new char*[(int)MAX_N_STATES];
	 for(j=0;j<(int)MAX_N_STATES;j++)
		states[i][j]=new char[(int)MAX_STATE_LENGTH];
	}

for(i=0;i<n_vars;i++)
	n_states[i]=0;

f_in=fopen(data_file,"r");

for(i=0;i<n_cases;i++)
	{fgets(buffer,(int)MAX_LINE_LENGTH,f_in);
	 if((int)strlen(buffer)==(int)MAX_LINE_LENGTH-1)
		{printf("\n\n Exceeding MAX_LINE_LENGTH. \n\n");
		 exit(0);
		}

	 token=strtok(buffer,separators);
	 for(j=0;j<n_vars;j++)
		{if((int)strlen(token)>(int)MAX_STATE_LENGTH-1)
			{printf("\n\n Exceeding MAX_STATE_LENGTH. \n\n");
			 exit(0);
			}

		 data[i][j]=state_index(j,token,states);

		 token=strtok(0,separators);
		}
	}

fclose(f_in);

min_n_states=0;
for(i=0;i<n_vars;i++)
	if(n_states[i]>1 && (min_n_states==0 || n_states[i]<min_n_states))
		min_n_states=n_states[i];

if(min_n_states==0)
	{printf("\n\n Only one state for all the variables. \n\n");
	 exit(0);
	}

max_max_cond_size=(int)(log((double)n_cases/(double)((int)N_CASES_PER_DF*(min_n_states-1)*(min_n_states-1)))/log((double)min_n_states)); //min_n_states may not be representative.

for(i=0;i<n_vars;i++)
	{for(j=0;j<(int)MAX_N_STATES;j++)
		delete [] states[i][j];
	 delete [] states[i];
	}
delete [] states;
}


////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////

int state_index(int var, char *state, char ***states)

{int i;

for(i=0;i<n_states[var];i++)
	if(strcmp(states[var][i],state)==0)
		return i;

if(i>(int)MAX_N_STATES-1)
	{printf("\n\n Exceeding MAX_N_STATES. \n\n");
	 exit(0);
	}

strcpy(states[var][i],state);
n_states[var]++;

return i;
}




////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////

void m3b(int targetNode)
{   /* A pc_spouse function written by Liang */
    // be careful no Loop in graph
    printf("\n\n\n-----------------------------------------------\n");
    printf("start m3b:targetNode is %d\n",targetNode);
	answerList[0] = targetNode;

	
	int deep = 1, i;
	int *Adj = new int[n_vars];
	int **Sep = new int*[n_vars];
	for(i=0; i<n_vars; i++)  // init the sep array;
	    Sep[i]=new int[max_max_cond_size];


	adjv_superset(targetNode, Adj, Sep); // get Adj and Sep for targetNode;

    printf("\npc:");
	fp_mb=fopen("result_mb.txt","a+");
	fp_pc=fopen("result_pc.txt","a+");

	for (i = 0 ; i<n_vars; i++) 
		if (Adj[i]==1)
		{
			printf("%d,",i);
			fprintf(fp_mb,"%d ",i);
			fprintf(fp_pc,"%d ",i);

		}
	printf("\n");

	//fprintf(fp_mb,"\n");
	  fprintf(fp_pc,"\n");
    
    fclose(fp_mb);
	fclose(fp_pc);

    


	//for (every AdjNode in target's Adj)
    for (i = 0 ; i<n_vars ; i++) {
        if (Adj[i]==1) {
            int AdjNode = i;
            RecSearch(targetNode, AdjNode, answerList, deep, Adj, Sep);
        }
	}



		fp_mb=fopen("result_mb.txt","a+");
	    fprintf(fp_mb,"\n");
        fclose(fp_mb);

    // remember free the memory;
	delete [] Adj;
    for(i=0;i<n_vars;i++)
        delete [] Sep[i];
    delete [] Sep;
}


void RecSearch(int oldNode, int targetNode, int* answerList, int deep,
                     int *AdjForOldNode, int** SepForOldNode) {
    //printf("\nnow in dfs function, oldNode=%d, targetNode=%d, deep=%d\n",oldNode, targetNode, deep);

    
	
    fp_mb=fopen("result_mb.txt","a+");

	answerList[deep] = targetNode; // add node to answer list;

    int *Adj = new int[n_vars],i;
	int **Sep = new int*[n_vars];
	for(i=0; i<n_vars; i++) // init the sep array;
	    Sep[i]=new int[max_max_cond_size];

    adjv_superset(targetNode, Adj, Sep); // get Adj and Sep for targetNode;

    for (i = 0 ; i<n_vars ; i++) {
        if (Adj[i]==1 && !inAnswerList(answerList, i, deep)) { //
            int AdjNode = i; // get Adj node;
            if ( isVStruct(oldNode, targetNode, AdjNode,  AdjForOldNode, SepForOldNode[AdjNode]) ) {
                 RecSearch(targetNode, AdjNode, answerList, deep+1, Adj, Sep);
            }
        }
    }
    if (deep >= 2)
    {   printAnswerList(answerList, deep);}


	fclose(fp_mb);
    // remember free the memory
	delete [] Adj;
    for(i=0;i<n_vars;i++)
        delete [] Sep[i];
    delete [] Sep;

}


bool inAnswerList(int *answerList, int target, int length) // 0 ~ length
{
    for (int i = 0 ; i<=length; i++) {
        if (target == answerList[i])
        {    return true;    }
    }
    return false;
}


bool isVStruct (int oldNode, int targetNode, int AdjNode, int* AdjForOldNode, int* SepForOldNodeAndAdjNode)
{
    if (oldNode==targetNode || targetNode==AdjNode || oldNode==AdjNode) return false;
    // 3 nodes must be different

    if (AdjForOldNode[AdjNode] == 1) return false;
    // oldNode and AdjNode can't be connect

    if (inSepSet(targetNode, SepForOldNodeAndAdjNode)) return false;
    // targetNode can't in the sep of oldNode and AdjNode

    int *cond=new int[max_max_cond_size];
    int n_conds = 0;
    for (int i = 0 ; i<max_max_cond_size; i++){
        cond[i] = SepForOldNodeAndAdjNode[i];
        if(cond[i]!=-1)   n_conds++;
    }
    if( n_conds < max_max_cond_size ){
        cond[n_conds]=targetNode;
        if(compute_dep(AdjNode, oldNode, cond)>=(double)1.0){
            delete [] cond;
            return true;
        }
    }
    delete [] cond;
    return false;
}

bool inSepSet(int target, int *SepForOldNodeAndAdjNode)
{
    for (int i = 0 ; i<max_max_cond_size; i++) {
        if (target == SepForOldNodeAndAdjNode[i])  return true;
    }
    return false;
}
void printAnswerList(int *answerList, int length)
{
   
    fp_mb=fopen("result_mb.txt","a+");
	printf("\nstart to print the answer:\n");
    printf("%d->",answerList[0]);
	fprintf(fp_mb,"%d ",answerList[0]);

    for (int i = 1 ; i<=length-2; i++) {
        printf("%d<->",answerList[i]);
		fprintf(fp_mb,"%d ",answerList[i]);
    }

    printf("%d",answerList[length-1]);
	fprintf(fp_mb,"%d ",answerList[length-1]);
    printf("<-%d\n",answerList[length]);
    fprintf(fp_mb,"%d ",answerList[length]);
	//fprintf(fp_out,"\n");
	fclose(fp_mb);
}






///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
void adjv_superset(int target, int *pc,int **sep)

{int i,j,t2,t3,condsize,n_pc;
 int *pc1,*cond;
 double *dep,cond_dep,t1;
 int *order,*sep2;

sep2=new int[max_max_cond_size];
pc1=new int[n_vars];
dep=new double[n_vars];
cond=new int[max_max_cond_size];


order=new int[n_vars];
dep=new double[n_vars];
for(i=0;i<n_vars;i++)
{
    order[i]=i;
	dep[i]=(double)0.0;
}


for(i=0;i<n_vars;i++)
	if(n_states[i]>1 && n_states[target]>1)
         pc[i]=1;
	else
	     pc[i]=0;



pc[target]=0;


for(i=0;i<max_max_cond_size;i++)

	cond[i]=-1;


//printf("\n");
for(i=0;i<n_vars;i++)
{
    if(pc[i]==1)
         dep[i]=compute_dep(i,target,cond);

    if(dep[i]<=(double)(-1.0))
	   {
		 pc[i]=0;
		 for(j=0;j<max_max_cond_size;j++)
				 sep[i][j]=cond[j];
	   }
}

//printf("\n");

for(i=0;i<n_vars;i++)
      pc1[i]=pc[i];





//sort dep by descending order


  for(j=0; j<n_vars; ++j)
    {
        for(i=0; i<n_vars-j-1; ++i)
        {
            if(dep[i]<dep[i+1])
            {

				t1=dep[i];
                dep[i]=dep[i+1];
                dep[i+1]=t1;

                t3=pc1[i];
				pc1[i]=pc1[i+1];
				pc1[i+1]=t3;


                t2= order[i];
                order[i]= order[i+1];
                order[i+1]=t2;
            }
        }
    }




condsize=1;


n_pc=0;
for(i=0;i<n_vars;i++)
	if(pc1[i]==1)
		n_pc++;




while (condsize<=4)
//while (condsize<=n_pc)

{
	  for (i=n_vars-1;i>=0;--i)
	  {
         if (pc1[i]==1)
		 {  pc1[i]=0;

              cond_dep=new_subset_dep(order[i],target,pc1,order,sep2, condsize);
              

		    pc1[i]=1;

			if (cond_dep==88888)
				   break;

			if(cond_dep<=(double)(-1.0))
			 {
			    pc1[i]=0;

          	    pc[order[i]]=0;

			    for(j=0;j<max_max_cond_size;j++)
				     sep[order[i]][j]=sep2[j];
			}
		 }
	  }


n_pc=0;
for(i=0;i<n_vars;i++)
	if(pc1[i]==1)
		n_pc++;

condsize++;


}




delete [] pc1;
delete [] cond;
delete [] order;
delete [] dep;

delete [] sep2;
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////


int next_cond_index(int n_pc, int cond_size, int *cond_index)


{int i, j, stop;

stop=1;
for(i=cond_size-1;i>=0;i--)
	{if(cond_index[i]<n_pc+i-cond_size)
		{cond_index[i]++;

		 if(i<cond_size-1)
			for(j=i+1;j<cond_size;j++)
				cond_index[j]=cond_index[j-1]+1;

		 stop=0;
		 i=-1;
		}
	}

return stop;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////


double compute_dep(int var, int target, int *cond)



{int i, j, k, n_cond_states, cond_state, df, df_target, df_var;
 double statistic, dep, aux;
 int *ss_cond;
 int **ss_target, **ss_var;
 int ***ss;




n_cond_states=1;
for(i=0;i<max_max_cond_size;i++)
	if(cond[i]!=-1)
		n_cond_states=n_cond_states*n_states[cond[i]];

if(n_cond_states*(n_states[target]-1)*(n_states[var]-1)*(int)N_CASES_PER_DF>n_cases)
{
  
   return (double)0.0;
 

}

ci_times++;

 



ss_cond=new int[n_cond_states];
for(i=0;i<n_cond_states;i++)
	ss_cond[i]=0;

ss_target=new int*[n_cond_states];
for(i=0;i<n_cond_states;i++)
	{ss_target[i]=new int[n_states[target]];
	 for(j=0;j<n_states[target];j++)
		ss_target[i][j]=0;
	}

ss_var=new int*[n_cond_states];
for(i=0;i<n_cond_states;i++)
	{ss_var[i]=new int[n_states[var]];
	 for(j=0;j<n_states[var];j++)
		ss_var[i][j]=0;
	}

ss=new int**[n_cond_states];
for(i=0;i<n_cond_states;i++)
	{ss[i]=new int*[n_states[target]];
	 for(j=0;j<n_states[target];j++)
		{ss[i][j]=new int[n_states[var]];
		 for(k=0;k<n_states[var];k++)
			 ss[i][j][k]=0;
		}
	}

for(i=0;i<n_cases;i++)
	{cond_state=0;
	 for(j=0;j<max_max_cond_size;j++)
		if(cond[j]!=-1)
			cond_state=cond_state*n_states[cond[j]]+data[i][cond[j]];

	 ss_cond[cond_state]++;
	 ss_target[cond_state][data[i][target]]++;
	 ss_var[cond_state][data[i][var]]++;
	 ss[cond_state][data[i][target]][data[i][var]]++;
	}



//G^2 statistic based on observed and expected frequencies.

statistic=(double)0.0;
for(i=0;i<n_cond_states;i++)
	if(ss_cond[i]>0)
		for(j=0;j<n_states[target];j++)
			if(ss_target[i][j]>0)
				for(k=0;k<n_states[var];k++)
					if(ss[i][j][k]>0)
						{aux=(double)(ss_target[i][j]*ss_var[i][k])/(double)ss_cond[i];
						 statistic=statistic+(double)ss[i][j][k]*(log((double)ss[i][j][k])-log(aux));
						}
statistic=statistic*(double)2.0;


//Reduced df due to zero marginals.

df=0;
for(i=0;i<n_cond_states;i++)
	if(ss_cond[i]>0)
		{df_target=0;
		 for(j=0;j<n_states[target];j++)
			if(ss_target[i][j]>0)
				df_target++;

		 df_var=0;
		 for(k=0;k<n_states[var];k++)
			 if(ss_var[i][k]>0)
				df_var++;

		 df=df+(df_target-1)*(df_var-1);
		}


delete [] ss_cond;
for(i=0;i<n_cond_states;i++)
	{delete [] ss_target[i];
	 delete [] ss_var[i];
	 for(j=0;j<n_states[target];j++)
		delete [] ss[i][j];
	 delete [] ss[i];
	}
delete [] ss_target;
delete [] ss_var;
delete [] ss;

if(statistic<(double)0.0) //statistic sometimes takes values -0.0. Due to loss of precision ?
	statistic=(double)0.0;

if(df<=0) //Naive ?
	df=1;

dep=gammq((double)0.5*(double)df,(double)0.5*statistic);

if(dep<(double)0.0) //Just in case.
	dep=(double)0.0;
else
	if(dep>(double)1.0)
		dep=(double)1.0;

if(dep<=alpha)
	{dep=(double)2.0-dep;

	 if(dep==(double)2.0)
		dep=(double)2.0+statistic/(double)df;

	 return dep;
	}
else
	{dep=(double)(-1.0)-dep;

	 if(dep==(double)(-2.0))
		dep=(double)(-2.0)-statistic/(double)df;

	 return dep;
	}
}




////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double new_subset_dep(int var, int target, int *pc, int *order,int *sep2, int subset_size)



{int i, j, n_pc, stop, cond_size;//max_cond_size;
 double dep1;
 int *cond, *code, *cond_index;

cond=new int[max_max_cond_size];
for(i=0;i<max_max_cond_size;i++)
	{cond[i]=-1;
	 sep2[i]=-1;
	}

n_pc=0;
for(i=0;i<n_vars;i++)
	if(pc[i]==1)
		n_pc++;


if (subset_size>n_pc)
{
	dep1=88888;
	return dep1;
}





code=new int[n_pc];
j=0;
for(i=0;i<n_vars;i++)
	if(pc[i]==1)
		{code[j]=order[i];

	      
		 j++;
		}



dep1=(double)(0.0);

cond_size=subset_size;

cond_index=new int[cond_size];
for(i=0;i<cond_size;i++)
		cond_index[i]=i;

 do
	{stop=0;

		for(i=0;i<cond_size;i++)
	    	cond[i]=code[cond_index[i]];

		dep1=compute_dep(var,target,cond);


	     if(dep1<=(double)(-1.0))
				{for(i=0;i<max_max_cond_size;i++)
					sep2[i]=cond[i];
   				    stop=1;
				
				   
				}


		 if(stop==0)
			stop=next_cond_index(n_pc,cond_size,cond_index);
		}
	 while(stop==0);


delete [] cond_index;
delete [] code;
delete [] cond;

return dep1;
}
