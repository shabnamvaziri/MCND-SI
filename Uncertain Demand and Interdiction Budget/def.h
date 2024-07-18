#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <limits.h>
#include <string.h>
#include <sys/types.h>
//#include <ilcplex/cplex.h>
#include <../../../../../encs/pkg/cplex-22.1.0/root/cplex/include/ilcplex/cplex.h>
#include <sys/time.h>
struct timeval start, stop;
struct timeval startTotal, stopTotal, startImprove, endImprove;

#define getrandom( min, max ) ((rand() % (int) (((max)+1)-(min)))+(min))
#define true 1
#define false 0
#define NOT_AVAIL -1
#define NONE -1


//char outfile[20];
//FILE  *Outy = NULL;
time_t	ttt;
struct tm	*tm;

typedef struct AVAL
{
	int index;
	double value;
	int index1;
	int index2;
} AVAL;
AVAL* Sorted;


void read_instance(const char *name);
double solve_MasterProblem(void);
double solve_ParetoCut(void);
double warm_start(void);
double warm_start2(void);
double warm_start3(void);
double solve_thirdlevel(void);
double solve_SubProblem(void);
double solve_branch_benders_cut(void);
double SPViolated(double *temp_x);
void Print_solution(void);
FILE *Open_File(const char *, const char *);
void free_memory(void);
void Initialize_memory(void);
void i_vector(int **vector, int n, char *s);
void d_vector(double **vector, int n, char *s);
void c_vector(char **vector, int n, char *s);
int **create_int_matrix(int, int);
double **create_double_matrix(int, int);
double ***create_double_matrix3D(int, int, int);
int *create_int_vector(int);
double *create_double_vector(int);
void createScenario(void);
void rec2(int, int, int);
void rec(int, int);
int Comparevalue(const void* , const void* );
char* create_char_vector(int);
struct cutinfo {
	CPXLPptr lp;
	int    numcols;
	double bestlb;
	////////////Ibrahim add by me////////
	long nodeid;
	double nodeobjval;
	int objsen;
	////////////Ibrahim add by me////////
};

typedef struct cutinfo CUTINFO, *CUTINFOptr;  //Declaring the cutinfo and the cutinfoPtr 


static int makeusercuts(CPXENVptr env, CPXLPptr lp, CUTINFOptr cutinfo);
static int CPXPUBLIC user_cut_callback1(CPXCENVptr env, void *cbdata, int wherefrom, void *info, int *useraction_p);
static int makelazyconstraint(CPXENVptr env, CPXLPptr lp, CUTINFOptr cutinfo);
static int CPXPUBLIC lazy_feas_callback1(CPXCENVptr env, void *cbdata, int wherefrom, void *info, int *useraction_p);


int			vioptcuts;//Number of violated optimality cuts found

/*******Global Reporting*******/
double		Upper_bound;
double		best_lower_bound;    //Value of the best lower bound so far
int		nodecount;
double		cputimeBD;