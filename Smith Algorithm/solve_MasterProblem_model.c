#include "def.h"

extern double     **pi, *alpha, *varphi, *y;
extern int        *x1;
extern int        *w;
extern int        S, iter, K, N, A, *from_node_arc, *to_node_arc, *source, *destination, sglobal, max_iter, *b;
extern int        *B;
extern double     *r, *f, *p, *d, *c, M, obj, *vc, rho;
extern double     *theta, theta0, Q;
extern double     MAX_DOUBLE;
extern double     MAX_DOUBLE2;
extern double	  MasterTime, SubTime;
clock_t	   startTemp, endTemp;
extern double installation, pre;

double solve_MasterProblem(void)
{
	int i, j, k, s, l, h;
	int index, index1;  // auxiliar indices to fill in the constraint matrix
	double best_upper_bound, best_lower_bound;
	int nodecount;     //Variables to call cplex
	CPXLPptr  lp;      // data strucutre to store a problem in cplex ...................
	CPXENVptr env;     // cplex environment.............................................
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
	char      probname[16]; // problem name for cplex .......................................
	char      *ctype;  // variable type ('C', 'I', 'B') only if integer.................
	double    value;   // objevtive value of solution ..................................
	double    num_w_var, num_y_var, num_theta_var;
	double    *pos_theta;
	double    *pos_y;
	int       *pos_w;
	double tempV = 0;
	int jcounter = 0;
	int ii = 0;
	double *temp1;
	double *temp2, *temp3;
	pos_theta = create_double_vector(S);
	pos_y = create_double_vector((A + K)* K);
	pos_w = create_int_vector(A + K);
	temp1 = create_double_vector(S);
	temp2 = create_double_vector((A + K)*S);

	//Initialize CPLEX environment
	env = CPXopenCPLEX(&status);
	if (env == NULL) {
		char  errmsg[1024];
		printf("Could not open CPLEX. \n");
		CPXgeterrorstring(env, status, errmsg);
		printf("%s", errmsg);
	}

	// Create the problem in CPLEX 
	strcpy(probname, "UFLP");
	lp = CPXcreateprob(env, &status, probname);
	if (env == NULL) {
		char  errmsg[1024];
		printf("Could not create LP. \n");
		CPXgeterrorstring(env, status, errmsg);
		printf("%s", errmsg);
	}
	CPXchgobjsen(env, lp, CPX_MIN);

	//Define w_h variables
	index1 = 0;  // index of columns
	numcols = A + K;
	d_vector(&obj, numcols, "open_cplex:1");
	d_vector(&lb, numcols, "open_cplex:8");
	d_vector(&ub, numcols, "open_cplex:9");
	c_vector(&ctype, numcols, "open_cplex:01");

	for (h = 0; h < A + K; h++){
		pos_w[h] = index1;
		obj[index1] = f[h];
		ctype[index1] = 'B';
		lb[index1] = 0;
		ub[index1] = 1;
		index1++;
	}

	status = CPXnewcols(env, lp, index1, obj, lb, ub, ctype, NULL);
	if (status)
		fprintf(stderr, "CPXnewcols failed.\n");
	free(obj);
	free(lb);
	free(ub);
	free(ctype);
	num_w_var = index1;

	//Define y_hk variables
	index1 = 0;  // index of columns
	numcols = (A + K) * K;
	d_vector(&obj, numcols, "open_cplex:1");
	d_vector(&lb, numcols, "open_cplex:8");
	d_vector(&ub, numcols, "open_cplex:9");
	c_vector(&ctype, numcols, "open_cplex:01");

	for (h = 0; h < A + K; h++){
		for (k = 0; k < K; k++){
			pos_y[h*K + k] = index1 + num_w_var;
			obj[index1] = vc[h] * rho;
			ctype[index1] = 'C';
			lb[index1] = 0;
			ub[index1] = CPX_INFBOUND;
			index1++;
		}
	}
	status = CPXnewcols(env, lp, index1, obj, lb, ub, ctype, NULL);
	if (status)
		fprintf(stderr, "CPXnewcols failed.\n");
	free(obj);
	free(lb);
	free(ub);
	free(ctype);
	num_y_var = index1;

	//Define theta variables
	index1 = 0;  // index of columns
	numcols = S;
	d_vector(&obj, numcols, "open_cplex:1");
	d_vector(&lb, numcols, "open_cplex:8");
	d_vector(&ub, numcols, "open_cplex:9");
	c_vector(&ctype, numcols, "open_cplex:01");

	for (s = 0; s < S; s++){
		pos_theta[s] = index1 + num_w_var + num_y_var;
		obj[index1] = (1 - rho)*p[s];
		ctype[index1] = 'C';
		lb[index1] = -CPX_INFBOUND;
		ub[index1] = CPX_INFBOUND;
		index1++;
	}
	status = CPXnewcols(env, lp, index1, obj, lb, ub, ctype, NULL);
	if (status)
		fprintf(stderr, "CPXnewcols failed.\n");
	free(obj);
	free(lb);
	free(ub);
	free(ctype);
	num_theta_var = index1;

	//Add constraint 1 // outward flow - inward fow = 0, d, -d
	numrows = K*N;
	numnz = 2 * (A + K)*K*N;
	d_vector(&rhs, numrows, "open_cplex:2");
	c_vector(&sense, numrows, "open_cplex:3");
	i_vector(&matbeg, numrows, "open_cplex:4");
	i_vector(&matind, numnz, "open_cplex:6");
	d_vector(&matval, numnz, "open_cplex:7");

	index = 0;
	index1 = 0;
	for (k = 0; k < K; k++){
		for (l = 0; l < N; l++){
			sense[index1] = 'G';
			if (destination[k]==l)
				rhs[index1] = -d[k];
			else if (source[k] == l)
				rhs[index1] = d[k];
			else
				rhs[index1] = 0;
			matbeg[index1++] = index;
			for (h = 0; h < A + K; h++){
				matind[index] = pos_y[h*K + k];
				if (from_node_arc[h] == l){
					matval[index++] = 1;
				}
				else if (to_node_arc[h] == l){
					matval[index++] = -1;
				}
			}
		}
	}

	status = CPXaddrows(env, lp, 0, index1, index, rhs, sense, matbeg, matind, matval, NULL, NULL);
	if (status)
		fprintf(stderr, "CPXaddrows failed.\n");
	free(matbeg);
	free(matind);
	free(matval);
	free(sense);
	free(rhs);

	//Add constraint 2 // capacity
	numrows = A + K;
	numnz = K * (A + K) + A + K;
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
		for (k = 0; k < K; k++){
			matind[index] = pos_y[h*K + k];
			matval[index++] = 1;
		}
			matind[index] = pos_w[h];
			matval[index++] = -1 * c[h];
		}

	status = CPXaddrows(env, lp, 0, index1, index, rhs, sense, matbeg, matind, matval, NULL, NULL);
	if (status)
		fprintf(stderr, "CPXaddrows failed.\n");
	free(matbeg);
	free(matind);
	free(matval);
	free(sense);
	free(rhs);


	//Add constraint theta >= theta 0
	numrows = S;
	numnz = S;
	d_vector(&rhs, numrows, "open_cplex:2");
	c_vector(&sense, numrows, "open_cplex:3");
	i_vector(&matbeg, numrows, "open_cplex:4");
	i_vector(&matind, numnz, "open_cplex:6");
	d_vector(&matval, numnz, "open_cplex:7");

	index = 0;
	index1 = 0;
	for (s = 0; s < S; s++){
		sense[index1] = 'G';
		rhs[index1] = theta0;
		matbeg[index1++] = index;
		matind[index] = pos_theta[s];
		matval[index++] = 1;
	}

	status = CPXaddrows(env, lp, 0, index1, index, rhs, sense, matbeg, matind, matval, NULL, NULL);
	if (status)
		fprintf(stderr, "CPXaddrows failed.\n");
	free(matbeg);
	free(matind);
	free(matval);
	free(sense);
	free(rhs);

	if (iter > 0){
		for (ii = 0; ii < iter; ii++){
			for (s = 0; s < S; s++){
				temp1[s] = 0;
			}
			for (s = 0; s < S; s++){
				for (k = 0; k < K; k++){
					for (l = 0; l < N; l++){
						if (destination[k] == l){
							temp1[s] -= d[k] * pi[(l*K + k) + (s*N*K)][ii];
						}
						else if (source[k] == l){
							temp1[s] += d[k] * pi[(l*K + k) + (s*N*K)][ii];
						}
						else
							temp1[s] += 0;
					}
				}
			}

			//// add opimality cut, single cut version////
			numrows = S*iter;
			numnz = S + S*(A + K);
			d_vector(&rhs, numrows, "open_cplex:2");
			c_vector(&sense, numrows, "open_cplex:3");
			i_vector(&matbeg, numrows, "open_cplex:4");
			i_vector(&matind, numnz, "open_cplex:6");
			d_vector(&matval, numnz, "open_cplex:7");

			index = 0;
			index1 = 0;
			for (s = 0; s < S; s++){
				sense[index1] = 'G';
				rhs[index1] = temp1[s];
				matbeg[index1++] = index;
				matind[index] = pos_theta[s];
				matval[index++] = 1;
				for (h = 0; h < A + K; h++){
					matind[index] = pos_w[h];
					matval[index++] = c[h]*alpha[(h*S + s) + (ii*(A + K)*S)];
				}
			}
			status = CPXaddrows(env, lp, 0, index1, index, rhs, sense, matbeg, matind, matval, NULL, NULL);
			if (status)
				fprintf(stderr, "CPXaddrows failed.\n");
			free(matbeg);
			free(matind);
			free(matval);
			free(sense);
			free(rhs);

		}
	}


	CPXsetintparam(env, CPX_PARAM_SCRIND, CPX_ON); //output display
	CPXsetintparam(env, CPX_PARAM_THREADS, 4);
	CPXsetintparam(env, CPX_PARAM_MIPDISPLAY, 3); //different levels of output display
	CPXsetdblparam(env, CPX_PARAM_TILIM, 86400); // time limit
	status = CPXsetintparam(env, CPX_PARAM_MEMORYEMPHASIS, 1);	//conserve memory where possible
	CPXsetintparam(env, CPX_PARAM_NODEFILEIND, 0);
	CPXsetdblparam(env, CPX_PARAM_EPGAP, 0.000001); // e-optimal solution (%gap)

	gettimeofday(&start, NULL);
	//startTemp = clock();
	CPXmipopt(env, lp);  //solve the integer program
	gettimeofday(&stop, NULL);
	//endTemp = clock();
	//MasterTime = (double)(endTemp - startTemp) / (double)(CLOCKS_PER_SEC);
	MasterTime = ((double)(stop.tv_sec - start.tv_sec) * 1000 + (double)(stop.tv_usec - start.tv_usec) / 1000) / 1000;

	i = CPXgetstat(env, lp);
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
	CPXgetmipobjval(env, lp, &value);
	printf("Upper bound: %.2f   ", value);
	best_upper_bound = value;
	
	CPXgetbestobjval(env, lp, &value);  //best lower bound in case the problem was not solved to optimality
	best_lower_bound = value;
	printf("Lower bound: %.2f  \n", value);

	nodecount = CPXgetnodecnt(env, lp);
	printf("Number of BB nodes : %ld  \n", nodecount);

	numcols = CPXgetnumcols(env, lp);
	d_vector(&x, numcols, "open_cplex:0");
	CPXgetmipx(env, lp, x, 0, numcols - 1);  // obtain the values of the decision variables

	if (lp != NULL) {
		status = CPXfreeprob(env, &lp);
		if (status) {
			fprintf(stderr, "CPXfreeprob failed, error code %d.\n", status);
		}
	}
	if (env != NULL) {
		status = CPXcloseCPLEX(&env);
		if (status) {
			char  errmsg[1024];
			fprintf(stderr, "Could not close CPLEX environment.\n");
			CPXgeterrorstring(env, status, errmsg);
			fprintf(stderr, "%s", errmsg);
		}
	}

	int number = 0;
	index = 0;
	for (h = 0; h < A + K; h++){
		if (x[index]>0.5)
			w[h] = 1;
		else
			w[h] = 0;
		//printf("w[%d]= %d\n", h, w[h]);
		index++;
	}

	for (h = 0; h < A + K; h++){
		for (k = 0; k < K; k++){
			y[h* K + k] = x[index];
			index++;
		}
	}

for (h = 0; h < A + K; h++){
installation += w[h]*f[h];
}
for (h = 0; h < A + K; h++){
		for (k = 0; k < K; k++){
pre += vc[h]*y[h* K + k] ;
}
}
	for (s = 0; s < S; s++){
		theta[s] = x[index];
		index++;
	}

	free(x);
	return best_upper_bound;
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

