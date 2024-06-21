#include "def.h"

extern double     **pi, *alpha, *varphi, *y;
extern double        *x1;
extern double     *w;
extern int        S, iter, K, N, A, *from_node_arc, *to_node_arc, *source, *destination, sglobal, max_iter, *b;
extern int        *B;
extern double     *r, *f, *p, *d, *c, M, obj, *vc, rho;
extern double     *theta, theta0, Q;
extern double     MAX_DOUBLE;
extern double     MAX_DOUBLE2;
extern double	  MasterTime, SubTime;
extern int		  *flag;
extern double     *piVio, *alphaVio;
extern double	violation;
extern int		numCutRelax;
extern int		numCutBranchLazy;
extern int		numCutBranchuser;
clock_t	   startTemp, endTemp;
extern FILE *Out;
extern char OutName[100];
extern char outfile[20];


double SPViolated(double *temp_x)
{
	int i, k, s, h, l;
	int index, index1;  // auxiliar indices to fill in the constraint matrix
	double best_upper_bound, best_lower_bound;
	int nodecount;     //Variables to call cplex
	CPXLPptr  lpVio;      // data strucutre to store a problem in cplex ...................
	CPXENVptr envVio;     // cplex environment.............................................
	int       numcols; // number of variables ..........................................
	int       numrows; // number of constraints.........................................
	int       numnz;   // number of non-zero elements in the matrix ....................
	int       objsen;  // optimization sense (min:1, max:-1 ) ..........................
	double    *obj;    // objective function coefficients ..............................
	double    *rhs;    // right and side of constraints ................................
	char      *sense;  // constraints sense (<=: 'L', =:'E', >=:'G') ...................
	int       *matbeg; // index of first non-zero element in each row...................
	int       *matind; // associated column of each non-zelo element ...................
	double    *matval; // coefficient values fo the non-zero elements of constraints....
	double    *lb;     // lower bounds of variables.....................................
	double    *ub;     // upper bounds of variables.....................................
	int       status;  // optimization status......................... .................
	double    *x;      // solution vector (double, even if the problem is integer) .....
	char		probname[16]; // problem name for cplex .......................................
	char      *ctype;  // variable type ('C', 'I', 'B') only if integer.................
	double    value;   // objevtive value of solution ..................................
	double		temp1 = 0, temp2 = 0, temp3 = 0, *temp4;
	double      num_alfa_var, num_pi_var, num_varphi_var;
	double      *pos_alfa, *pos_pi, *pos_varphi, *pos_x;
	int          num_x_var;
	int auxilary = 0;
	int auxilary1 = 0;
	pos_x = create_double_vector(A + K);
	pos_alfa = create_double_vector(A + K);
	pos_pi = create_double_vector(N*K);
	pos_varphi = create_double_vector(A + K);
	auxilary = 0;
	auxilary1 = 0;

	//Initialize CPLEX environment
	envVio = CPXopenCPLEX(&status);
	if (envVio == NULL) {
		char  errmsg[1024];
		printf("Could not open CPLEX. \n");
		CPXgeterrorstring(envVio, status, errmsg);
		printf("%s", errmsg);
	}

	// Create the problem in CPLEX 
	strcpy(probname, "UFLP");
	lpVio = CPXcreateprob(envVio, &status, probname);
	if (envVio == NULL) {
		char  errmsg[1024];
		printf("Could not create LP. \n");
		CPXgeterrorstring(envVio, status, errmsg);
		printf("%s", errmsg);
	}

	CPXchgobjsen(envVio, lpVio, CPX_MAX);

	//Define x variables
	index1 = 0;  // index of columns
	numcols = A + K;
	d_vector(&obj, numcols, "open_cplex:1");
	d_vector(&lb, numcols, "open_cplex:8");
	d_vector(&ub, numcols, "open_cplex:9");
	c_vector(&ctype, numcols, "open_cplex:01");

	for (h = 0; h < A + K; h++){
		pos_x[h] = index1;
		obj[index1] = 0;
		ctype[index1] = 'B';
		lb[index1] = 0;
		ub[index1] = 1;
		index1++;
	}

	status = CPXnewcols(envVio, lpVio, index1, obj, lb, ub, ctype, NULL);
	if (status)
		fprintf(stderr, "CPXnewcols failed.\n");
	free(obj);
	free(lb);
	free(ub);
	free(ctype);
	num_x_var = index1;

	//Define alfa variables
	index1 = 0;  // index of columns
	numcols = A + K;
	d_vector(&obj, numcols, "open_cplex:1");
	d_vector(&lb, numcols, "open_cplex:8");
	d_vector(&ub, numcols, "open_cplex:9");
	c_vector(&ctype, numcols, "open_cplex:01");

	for (h = 0; h < A + K; h++){
		pos_alfa[h] = index1 + num_x_var;
		obj[index1] = -1 * (c[h] * w[h]) - 0.0000001;
		ctype[index1] = 'C';
		lb[index1] = 0;
		ub[index1] = CPX_INFBOUND;
		index1++;
	}

	status = CPXnewcols(envVio, lpVio, index1, obj, lb, ub, ctype, NULL);
	if (status)
		fprintf(stderr, "CPXnewcols failed.\n");
	free(obj);
	free(lb);
	free(ub);
	free(ctype);
	num_alfa_var = index1;

	//Define pi variables
	index1 = 0;  // index of columns
	numcols = N*K;
	d_vector(&obj, numcols, "open_cplex:1");
	d_vector(&lb, numcols, "open_cplex:8");
	d_vector(&ub, numcols, "open_cplex:9");
	c_vector(&ctype, numcols, "open_cplex:01");

	for (l = 0; l < N; l++){
		for (k = 0; k < K; k++){
			pos_pi[l*K + k] = index1 + num_x_var + num_alfa_var;
			if (destination[k] == l){
				obj[index1] = -1 * d[k];
			}
			else if (source[k] == l){
				obj[index1] = d[k];
			}
			else
				obj[index1] = 0;
			ctype[index1] = 'C';
			lb[index1] = -CPX_INFBOUND;
			ub[index1] = CPX_INFBOUND;
			index1++;
		}
	}

	status = CPXnewcols(envVio, lpVio, index1, obj, lb, ub, ctype, NULL);
	if (status)
		fprintf(stderr, "CPXnewcols failed.\n");
	free(obj);
	free(lb);
	free(ub);
	free(ctype);
	num_pi_var = index1;

	//Add constraint pi-pi-alfa <= vc  
	numrows = (A + K)*K;
	numnz = (A + K) + 2 * K*(A + K)*N;
	d_vector(&rhs, numrows, "open_cplex:2");
	c_vector(&sense, numrows, "open_cplex:3");
	i_vector(&matbeg, numrows, "open_cplex:4");
	i_vector(&matind, numnz, "open_cplex:6");
	d_vector(&matval, numnz, "open_cplex:7");

	index = 0;
	index1 = 0;
	for (h = 0; h < A + K; h++){
		for (k = 0; k < K; k++){
			sense[index1] = 'L';
			rhs[index1] = vc[h];
			matbeg[index1++] = index;
			matind[index] = pos_alfa[h];
			matval[index++] = -1;
			matind[index] = pos_pi[(int)from_node_arc[h] * K + k];
			matval[index++] = 1;
			matind[index] = pos_pi[to_node_arc[h] * K + k];
			matval[index++] = -1;
			matind[index] = pos_x[h];
			matval[index++] = -M;
		}
	}

	status = CPXaddrows(envVio, lpVio, 0, index1, index, rhs, sense, matbeg, matind, matval, NULL, NULL);
	if (status)
		fprintf(stderr, "CPXaddrows failed.\n");
	free(matbeg);
	free(matind);
	free(matval);
	free(sense);
	free(rhs);

	//Add constraint budget
	numrows = 1;
	numnz = A + K;
	d_vector(&rhs, numrows, "open_cplex:2");
	c_vector(&sense, numrows, "open_cplex:3");
	i_vector(&matbeg, numrows, "open_cplex:4");
	i_vector(&matind, numnz, "open_cplex:6");
	d_vector(&matval, numnz, "open_cplex:7");

	index = 0;
	index1 = 0;
	sense[index1] = 'L';
	rhs[index1] = B[sglobal];
	matbeg[index1++] = index;
	for (h = 0; h < A + K; h++){
		matind[index] = pos_x[h];
		matval[index++] = b[h];
	}

	status = CPXaddrows(envVio, lpVio, 0, index1, index, rhs, sense, matbeg, matind, matval, NULL, NULL);
	if (status)
		fprintf(stderr, "CPXaddrows failed.\n");
	free(matbeg);
	free(matind);
	free(matval);
	free(sense);
	free(rhs);

		//Add constraint pi origin=0
	numrows = N * K;
	numnz = N * K;
	d_vector(&rhs, numrows, "open_cplex:2");
	c_vector(&sense, numrows, "open_cplex:3");
	i_vector(&matbeg, numrows, "open_cplex:4");
	i_vector(&matind, numnz, "open_cplex:6");
	d_vector(&matval, numnz, "open_cplex:7");

	index = 0;
	index1 = 0;
	for (l = 0; l < N; l++) {
		for (k = 0; k < K; k++) {
			if (source[k] == l) {
				sense[index1] = 'E';
				rhs[index1] = 0;
				matbeg[index1++] = index;
				matind[index] = pos_pi[l * K + k];
				matval[index++] = 1;
			}
		}
	}

	status = CPXaddrows(envVio, lpVio, 0, index1, index, rhs, sense, matbeg, matind, matval, NULL, NULL);
	if (status)
		fprintf(stderr, "CPXaddrows failed.\n");
	free(matbeg);
	free(matind);
	free(matval);
	free(sense);
	free(rhs);



	//CPXwriteprob(envVio, lpVio, "model3.lp", NULL);                          //write the model in .lp format if needed (to debug)
	CPXsetintparam(envVio, CPX_PARAM_THREADS, 1);
	CPXsetintparam(envVio, CPX_PARAM_SCRIND, CPX_OFF); //output display
	//CPXsetintparam(env,CPX_PARAM_INTSOLLIM,1);    //stops after finding first integer sol.
	CPXsetintparam(envVio, CPX_PARAM_MIPDISPLAY, 3); //different levels of output display
	//CPXsetintparam(env,CPX_PARAM_MIPEMPHASIS,1);//0:balanced; 1:feasibility; 2:optimality,3:bestbound, 4:hiddenfeas
	CPXsetdblparam(envVio, CPX_PARAM_TILIM, 7200); // time limit
	CPXsetdblparam(envVio, CPX_PARAM_EPGAP, 0.0001); // e-optimal solution (%gap)
	//CPXsetdblparam(envVio, CPX_PARAM_EPGAP, 0.000001); // e-optimal solution (%gap)
	CPXsetintparam(envVio, CPX_PARAM_NODEFILEIND, 0);
	status = CPXsetintparam(envVio, CPX_PARAM_MEMORYEMPHASIS, 1);	//conserve memory where possible
	//CPXsetdblparam(env,CPX_PARAM_TRELIM, 14000); // B&B memory limit
	//CPXsetdblparam(env,CPX_PARAM_EPGAP, 0.0000000001); // e-optimal solution (%gap)
	//CPXsetdblparam(env,CPX_PARAM_EPAGAP, 0.0000000001); // e-optimal solution (absolute value)
	//CPXsetdblparam(env,CPX_PARAM_EPINT, 0.0000000001); // integer precision
	//CPXsetintparam(env,CPX_PARAM_THREADS, 1); // Number of threads to use
	//CPXsetdblparam(env,CPX_PARAM_EPRHS, 0.0000001);
	//CPXsetintparam(env,CPX_PARAM_REDUCE, 0);  // only needed when adding lazy constraints
	//CPXsetintparam(env,CPX_PARAM_HEURFREQ, -1); //heuristic frequency and intensisty 
	//CPXsetdblparam(env, CPX_PARAM_CUTSFACTOR, 1.0);  //limit the number of cuts added by cplex 1.0002
	//CPXsetdblparam(env,CPX_PARAM_CUTUP,UpperBound+.01); // provide an initial upper bound
	//CPXsetintparam(envVio, CPX_PARAM_MIPEMPHASIS, CPX_MIPEMPHASIS_OPTIMALITY);  // MIP emphasis: optimality, feasibility, moving best bound
	CPXsetintparam(envVio, CPX_PARAM_MIPEMPHASIS, CPX_MIPEMPHASIS_OPTIMALITY);  // MIP emphasis: optimality, feasibility, moving best bound
	//CPXsetintparam(env,CPX_PARAM_PARALLELMODE, 1); 
	//CPXsetintparam(env,CPX_PARAM_PREIND,0);
	//CPXsetintparam(env,CPX_PARAM_MIPORDIND,CPX_ON); // Turn on or off the use of priorities on bracnhing variables
	//CPXsetintparam(env,CPX_PARAM_MIPEMPHASIS,CPX_MIPEMPHASIS_BESTBOUND);  // MIP emphasis: optimality, feasibility, moving best bound

	//gettimeofday(&start, NULL);
	//startTemp = clock();
	CPXmipopt(envVio, lpVio);  //solve the integer program
	//gettimeofday(&stop, NULL);
	//endTemp = clock();
	//SubTime += (double)(endTemp - startTemp) / (double)(CLOCKS_PER_SEC);
	//SubTime += ((double)(stop.tv_sec - start.tv_sec) * 1000 + (double)(stop.tv_usec - start.tv_usec) / 1000) / 1000;

	i = CPXgetstat(envVio, lpVio);
	if (i == 101)
		printf("Optimal solution found\n");
	else if (i == 102)
		printf("e-optimal solution found\n");
	else if (i == 103)
		printf(" infeasible solution\n");
	else if (i == 107)
		printf("Time limit reached\n");
	else
		printf("Unknown stopping criterion (%d)\n", i);
	// retrive solution values
	CPXgetmipobjval(envVio, lpVio, &value);
	//printf("Upper bound: %.2f   ", value);
	best_upper_bound = value;
	// If CPLEX was able to find the optimal solution, the previous function provides the optimal solution value
	//if not, it provides the best upper bound
	CPXgetbestobjval(envVio, lpVio, &value);  //best lower bound in case the problem was not solved to optimality
	best_lower_bound = value;
	//printf("Lower bound: %.2f  \n", value);

	nodecount = CPXgetnodecnt(envVio, lpVio);
	//printf("Number of BB nodes : %ld  \n", nodecount);

	numcols = CPXgetnumcols(envVio, lpVio);
	d_vector(&x, numcols, "open_cplex:0");
	CPXgetmipx(envVio, lpVio, x, 0, numcols - 1);  // obtain the values of the decision variables

	if (lpVio != NULL) {
		status = CPXfreeprob(envVio, &lpVio);
		if (status) {
			fprintf(stderr, "CPXfreeprob failed, error code %d.\n", status);
		}
	}
	if (envVio != NULL) {
		status = CPXcloseCPLEX(&envVio);
		if (status) {
			char  errmsg[1024];
			fprintf(stderr, "Could not close CPLEX environment.\n");
			CPXgeterrorstring(envVio, status, errmsg);
			fprintf(stderr, "%s", errmsg);
		}
	}

	index = 0;
	for (h = 0; h < A + K; h++){
		x1[h*S + sglobal] = x[index];
		index++;
	}
	double sum_alfa = 0;
	for (h = 0; h < A + K; h++){
		alphaVio[h] = x[index];
		if (alphaVio[h] > 0){
			sum_alfa += alphaVio[h];
		}
		index++;
	}

	for (l = 0; l < N; l++){
		for (k = 0; k < K; k++){
			piVio[(l*K + k)] = x[index];
			index++;
		}
	}

	//printf("best_upper_bound = %f\n", best_upper_bound);
	free(x);


	index = 0;
	for (h = 0; h < A + K; h++){
		w[h] = temp_x[index];
		index++;
	}
	for (h = 0; h < A + K; h++){
		for (k = 0; k < K; k++){
			y[h* K + k] = temp_x[index];
			index++;
		}
	}
	for (s = 0; s < S; s++){
		theta[s] = temp_x[index];
		index++;
	}

	if (best_upper_bound + 0.0000001*sum_alfa - theta[sglobal] >= violation)//?????????
		flag[sglobal] = 1;
	else 
		flag[sglobal] = 0;


	return best_upper_bound + 0.0000001*sum_alfa;
}



/* This simple routine frees up the pointer *ptr, and sets *ptr
to NULL */

static void free_and_null(char **ptr)
{
	if (*ptr != NULL) {
		free(*ptr);
		*ptr = NULL;
	}
} /* END free_and_null */

