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
extern double	violation;

clock_t	   startTemp, endTemp;

double solve_SubProblem(void)
{
	int i, k, s, h, l;
	int index, index1;  // auxiliar indices to fill in the constraint matrix
	double best_upper_bound, best_lower_bound;
	int nodecount;     //Variables to call cplex
	CPXLPptr  lp2;      // data strucutre to store a problem in cplex ...................
	CPXENVptr env2;     // cplex environment.............................................
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
	env2 = CPXopenCPLEX(&status);
	if (env2 == NULL) {
		char  errmsg[1024];
		printf("Could not open CPLEX. \n");
		CPXgeterrorstring(env2, status, errmsg);
		printf("%s", errmsg);
	}

	// Create the problem in CPLEX 
	strcpy(probname, "UFLP");
	lp2 = CPXcreateprob(env2, &status, probname);
	if (env2 == NULL) {
		char  errmsg[1024];
		printf("Could not create LP. \n");
		CPXgeterrorstring(env2, status, errmsg);
		printf("%s", errmsg);
	}

	CPXchgobjsen(env2, lp2, CPX_MAX);

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

	status = CPXnewcols(env2, lp2, index1, obj, lb, ub, ctype, NULL);
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

	status = CPXnewcols(env2, lp2, index1, obj, lb, ub, ctype, NULL);
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

	status = CPXnewcols(env2, lp2, index1, obj, lb, ub, ctype, NULL);
	if (status)
		fprintf(stderr, "CPXnewcols failed.\n");
	free(obj);
	free(lb);
	free(ub);
	free(ctype);
	num_pi_var = index1;

	//Define varphi variables
	index1 = 0;  // index of columns
	numcols = A + K;
	d_vector(&obj, numcols, "open_cplex:1");
	d_vector(&lb, numcols, "open_cplex:8");
	d_vector(&ub, numcols, "open_cplex:9");
	c_vector(&ctype, numcols, "open_cplex:01");

	for (h = 0; h < A + K; h++){
		pos_varphi[h] = index1 + num_x_var + num_pi_var + num_alfa_var;
		obj[index1] = (c[h] * w[h]);
		ctype[index1] = 'C';
		lb[index1] = 0;
		ub[index1] = CPX_INFBOUND;
		index1++;
	}

	status = CPXnewcols(env2, lp2, index1, obj, lb, ub, ctype, NULL);
	if (status)
		fprintf(stderr, "CPXnewcols failed.\n");
	free(obj);
	free(lb);
	free(ub);
	free(ctype);
	num_varphi_var = index1;


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
		}
	}

	status = CPXaddrows(env2, lp2, 0, index1, index, rhs, sense, matbeg, matind, matval, NULL, NULL);
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

	status = CPXaddrows(env2, lp2, 0, index1, index, rhs, sense, matbeg, matind, matval, NULL, NULL);
	if (status)
		fprintf(stderr, "CPXaddrows failed.\n");
	free(matbeg);
	free(matind);
	free(matval);
	free(sense);
	free(rhs);

	//Add constraint linearization varphi-alfa <= 0  
	numrows = A + K;
	numnz = 2 * (A + K);
	d_vector(&rhs, numrows, "open_cplex:2");
	c_vector(&sense, numrows, "open_cplex:3");
	i_vector(&matbeg, numrows, "open_cplex:4");
	i_vector(&matind, numnz, "open_cplex:6");
	d_vector(&matval, numnz, "open_cplex:7");

	index = 0;
	index1 = 0;
	for (h = 0; h < A + K; h++){
		sense[index1] = 'L';
		rhs[index1] = 0;
		matbeg[index1++] = index;
		matind[index] = pos_varphi[h];
		matval[index++] = 1;
		matind[index] = pos_alfa[h];
		matval[index++] = -1;
	}

	status = CPXaddrows(env2, lp2, 0, index1, index, rhs, sense, matbeg, matind, matval, NULL, NULL);
	if (status)
		fprintf(stderr, "CPXaddrows failed.\n");
	free(matbeg);
	free(matind);
	free(matval);
	free(sense);
	free(rhs);

	//Add constraint linearization varphi - M*x <= 0 
	numrows = A + K;
	numnz = 2 * (A + K);
	d_vector(&rhs, numrows, "open_cplex:2");
	c_vector(&sense, numrows, "open_cplex:3");
	i_vector(&matbeg, numrows, "open_cplex:4");
	i_vector(&matind, numnz, "open_cplex:6");
	d_vector(&matval, numnz, "open_cplex:7");

	index = 0;
	index1 = 0;
	for (h = 0; h < A + K; h++){
		sense[index1] = 'L';
		rhs[index1] = 0;
		matbeg[index1++] = index;
		matind[index] = pos_varphi[h];
		matval[index++] = 1;
		matind[index] = pos_x[h];
		matval[index++] = -M;
	}

	status = CPXaddrows(env2, lp2, 0, index1, index, rhs, sense, matbeg, matind, matval, NULL, NULL);
	if (status)
		fprintf(stderr, "CPXaddrows failed.\n");
	free(matbeg);
	free(matind);
	free(matval);
	free(sense);
	free(rhs);



	//Add constraint linearization varphi - alfa - M*x >= -M 
	/*numrows = A + K;
	numnz = 3 * (A + K);
	d_vector(&rhs, numrows, "open_cplex:2");
	c_vector(&sense, numrows, "open_cplex:3");
	i_vector(&matbeg, numrows, "open_cplex:4");
	i_vector(&matind, numnz, "open_cplex:6");
	d_vector(&matval, numnz, "open_cplex:7");

	index = 0;
	index1 = 0;
	for (h = 0; h < A + K; h++){
	sense[index1] = 'G';
	rhs[index1] = -M;
	matbeg[index1++] = index;
	matind[index] = pos_varphi[h];
	matval[index++] = 1;
	matind[index] = pos_x[h];
	matval[index++] = -M;
	matind[index] = pos_alfa[h];
	matval[index++] = -1;
	}

	status = CPXaddrows(env2, lp2, 0, index1, index, rhs, sense, matbeg, matind, matval, NULL, NULL);
	if (status)
	fprintf(stderr, "CPXaddrows failed.\n");
	free(matbeg);
	free(matind);
	free(matval);
	free(sense);
	free(rhs);
	*/

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

	status = CPXaddrows(env2, lp2, 0, index1, index, rhs, sense, matbeg, matind, matval, NULL, NULL);
	if (status)
		fprintf(stderr, "CPXaddrows failed.\n");
	free(matbeg);
	free(matind);
	free(matval);
	free(sense);
	free(rhs);

	


	//CPXwriteprob(env2, lp2, "model3.lp", NULL);                          //write the model in .lp format if needed (to debug)
	CPXsetintparam(env2, CPX_PARAM_THREADS, 1);
	CPXsetintparam(env2, CPX_PARAM_SCRIND, CPX_OFF); //output display
	//CPXsetintparam(env,CPX_PARAM_INTSOLLIM,1);    //stops after finding first integer sol.
	CPXsetintparam(env2, CPX_PARAM_MIPDISPLAY, 3); //different levels of output display
	//CPXsetintparam(env,CPX_PARAM_MIPEMPHASIS,1);//0:balanced; 1:feasibility; 2:optimality,3:bestbound, 4:hiddenfeas
	CPXsetdblparam(env2, CPX_PARAM_TILIM, 3600); // time limit
	CPXsetdblparam(env2, CPX_PARAM_EPGAP, 0.005); // e-optimal solution (%gap)
	//CPXsetdblparam(env2, CPX_PARAM_EPGAP, 0.000001); // e-optimal solution (%gap)
	CPXsetintparam(env2, CPX_PARAM_NODEFILEIND, 0);
	status = CPXsetintparam(env2, CPX_PARAM_MEMORYEMPHASIS, 1);	//conserve memory where possible
	CPXsetintparam(env2, CPX_PARAM_MIPEMPHASIS, CPX_MIPEMPHASIS_OPTIMALITY);  // MIP emphasis: optimality, feasibility, moving best bound
	
	gettimeofday(&start, NULL);
	//startTemp = clock();
	CPXmipopt(env2, lp2);  //solve the integer program
	gettimeofday(&stop, NULL);
	//endTemp = clock();
	//SubTime += (double)(endTemp - startTemp) / (double)(CLOCKS_PER_SEC);
	SubTime += ((double)(stop.tv_sec - start.tv_sec) * 1000 + (double)(stop.tv_usec - start.tv_usec) / 1000) / 1000;

	i = CPXgetstat(env2, lp2);
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
	CPXgetmipobjval(env2, lp2, &value);
	//printf("Upper bound: %.2f   ", value);
	best_upper_bound = value;
	// If CPLEX was able to find the optimal solution, the previous function provides the optimal solution value
	//if not, it provides the best upper bound
	CPXgetbestobjval(env2, lp2, &value);  //best lower bound in case the problem was not solved to optimality
	best_lower_bound = value;
	//printf("Lower bound: %.2f  \n", value);

	nodecount = CPXgetnodecnt(env2, lp2);
	//printf("Number of BB nodes : %ld  \n", nodecount);

	numcols = CPXgetnumcols(env2, lp2);
	d_vector(&x, numcols, "open_cplex:0");
	CPXgetmipx(env2, lp2, x, 0, numcols - 1);  // obtain the values of the decision variables

	if (lp2 != NULL) {
		status = CPXfreeprob(env2, &lp2);
		if (status) {
			fprintf(stderr, "CPXfreeprob failed, error code %d.\n", status);
		}
	}
	if (env2 != NULL) {
		status = CPXcloseCPLEX(&env2);
		if (status) {
			char  errmsg[1024];
			fprintf(stderr, "Could not close CPLEX environment.\n");
			CPXgeterrorstring(env2, status, errmsg);
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
		alpha[(h*S + sglobal) + (iter*(A + K)*S)] = x[index];
		if (alpha[(h*S + sglobal) + (iter*(A + K)*S)] > 0){
			sum_alfa += alpha[(h*S + sglobal) + (iter*(A + K)*S)];
		}
		index++;
	}

	for (l = 0; l < N; l++){
		for (k = 0; k < K; k++){
			pi[(l*K + k) + (sglobal*N*K)][iter] = x[index];
			index++;
		}
	}

	for (h = 0; h < A + K; h++){
		varphi[(h*S + sglobal) + (iter*(A + K)*S)] = x[index];
		index++;
	}	

	//printf("best_upper_bound = %f\n", best_upper_bound);
	free(x);

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

