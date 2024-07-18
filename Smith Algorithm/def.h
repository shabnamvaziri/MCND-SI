#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <limits.h>
#include <string.h>
#include <sys/types.h>
//#include <ilcplex/cplex.h>
#include <../../../../../encs/pkg/cplex-12.10.0/root/cplex/include/ilcplex/cplex.h>
#include <sys/time.h>
struct timeval start, stop;
struct timeval startTotal, stopTotal;

#define getrandom( min, max ) ((rand() % (int) (((max)+1)-(min)))+(min))
#define true 1
#define false 0
#define NOT_AVAIL -1
#define NONE -1


//char outfile[20];
//FILE  *Outy = NULL;
time_t	ttt;
struct tm	*tm;


void read_instance(const char *name);
double solve_MasterProblem(void);
double solve_ParetoCut(void);
double solve_SubProblem(void);
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