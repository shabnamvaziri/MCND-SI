#include "def.h"

extern double     **pi, *alpha, *varphi, *y, **pi_pareto, *alpha_pareto, *piVio_pareto, *alphaVio_pareto;
extern double        *x1;
extern double        *w, *w_pareto;
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
extern double bestUB, bestLB, *dj;

double solve_MasterProblem(void)
{
	int i, j, k, s, l, h, pareto_iter;
	int index, index1;  // auxiliar indices to fill in the constraint matrix
	double best_upper_bound, best_lower_bound;
	int nodecount;     //Variables to call cplex
	CPXLPptr  lpPre;      // data strucutre to store a problem in cplex ...................
	CPXENVptr envPre;     // cplex environment.............................................
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
	int    num_w_var, num_y_var, num_theta_var;
	double    *pos_theta;
	double    *pos_y;
	double       *pos_w;
	double tempV = 0;
	int jcounter = 0;
	int ii = 0;
	double *temp1;
	double *temp2, *temp3;
	pos_theta = create_double_vector(S);
	pos_y = create_double_vector((A + K)* K);
	pos_w = create_double_vector(A + K);
	temp1 = create_double_vector(S);
	temp2 = create_double_vector((A + K)*S);

	//Initialize CPLEX envPreironment
	envPre = CPXopenCPLEX(&status);
	if (envPre == NULL) {
		char  errmsg[1024];
		printf("Could not open CPLEX. \n");
		CPXgeterrorstring(envPre, status, errmsg);
		printf("%s", errmsg);
	}

	// Create the problem in CPLEX 
	strcpy(probname, "UFlpPre");
	lpPre = CPXcreateprob(envPre, &status, probname);
	if (envPre == NULL) {
		char  errmsg[1024];
		printf("Could not create lpPre. \n");
		CPXgeterrorstring(envPre, status, errmsg);
		printf("%s", errmsg);
	}
	CPXchgobjsen(envPre, lpPre, CPX_MIN);

	//Define w_h variables
	index1 = 0;  // index of columns
	numcols = A + K;
	d_vector(&obj, numcols, "open_cplex:1");
	d_vector(&lb, numcols, "open_cplex:8");
	d_vector(&ub, numcols, "open_cplex:9");
	//c_vector(&ctype, numcols, "open_cplex:01");
	ctype = create_char_vector((A + K) + ((A + K) * K) + S);


	for (h = 0; h < A + K; h++){
		pos_w[h] = index1;
		obj[index1] = f[h];
		ctype[index1] = 'C';
		lb[index1] = 0;
		ub[index1] = 1;
		index1++;
	}

	status = CPXnewcols(envPre, lpPre, index1, obj, lb, ub, NULL, NULL);
	if (status)
		fprintf(stderr, "CPXnewcols failed.\n");
	free(obj);
	free(lb);
	free(ub);
	//free(ctype);
	num_w_var = index1;

	//Define y_hk variables
	index1 = 0;  // index of columns
	numcols = (A + K) * K;
	d_vector(&obj, numcols, "open_cplex:1");
	d_vector(&lb, numcols, "open_cplex:8");
	d_vector(&ub, numcols, "open_cplex:9");
	//c_vector(&ctype, numcols, "open_cplex:01");

	for (h = 0; h < A + K; h++){
		for (k = 0; k < K; k++){
			pos_y[h*K + k] = index1 + num_w_var;
			obj[index1] = vc[h] * rho;
			ctype[index1 + num_w_var] = 'C';
			lb[index1] = 0;
			ub[index1] = CPX_INFBOUND;
			index1++;
		}
	}
	status = CPXnewcols(envPre, lpPre, index1, obj, lb, ub, NULL, NULL);
	if (status)
		fprintf(stderr, "CPXnewcols failed.\n");
	free(obj);
	free(lb);
	free(ub);
	//free(ctype);
	num_y_var = index1;

	//Define theta variables
	index1 = 0;  // index of columns
	numcols = S;
	d_vector(&obj, numcols, "open_cplex:1");
	d_vector(&lb, numcols, "open_cplex:8");
	d_vector(&ub, numcols, "open_cplex:9");
	//c_vector(&ctype, numcols, "open_cplex:01");

	for (s = 0; s < S; s++){
		pos_theta[s] = index1 + num_w_var + num_y_var;
		obj[index1] = (1 - rho)*p[s];
		ctype[index1 + num_w_var + num_y_var] = 'C';
		lb[index1] = -CPX_INFBOUND;
		ub[index1] = CPX_INFBOUND;
		index1++;
	}
	status = CPXnewcols(envPre, lpPre, index1, obj, lb, ub, NULL, NULL);
	if (status)
		fprintf(stderr, "CPXnewcols failed.\n");
	free(obj);
	free(lb);
	free(ub);
	//free(ctype);
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
			if (destination[k] == l)
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

	status = CPXaddrows(envPre, lpPre, 0, index1, index, rhs, sense, matbeg, matind, matval, NULL, NULL);
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

	status = CPXaddrows(envPre, lpPre, 0, index1, index, rhs, sense, matbeg, matind, matval, NULL, NULL);
	if (status)
		fprintf(stderr, "CPXaddrows failed.\n");
	free(matbeg);
	free(matind);
	free(matval);
	free(sense);
	free(rhs);

	//Add constraint number of open arcs > A/2
	numrows = 1;
	numnz = A;
	d_vector(&rhs, numrows, "open_cplex:2");
	c_vector(&sense, numrows, "open_cplex:3");
	i_vector(&matbeg, numrows, "open_cplex:4");
	i_vector(&matind, numnz, "open_cplex:6");
	d_vector(&matval, numnz, "open_cplex:7");

	index = 0;
	index1 = 0;
	sense[index1] = 'G';
	rhs[index1] = A * 0.5;
	matbeg[index1++] = index;
	for (h = 0; h < A; h++){
		matind[index] = pos_w[h];
		matval[index++] = 1;
	}

	status = CPXaddrows(envPre, lpPre, 0, index1, index, rhs, sense, matbeg, matind, matval, NULL, NULL);
	if (status)
		fprintf(stderr, "CPXaddrows failed.\n");
	free(matbeg);
	free(matind);
	free(matval);
	free(sense);
	free(rhs);

	//Add constraint w dummy =1
	numrows = K;
	numnz = A + K;
	d_vector(&rhs, numrows, "open_cplex:2");
	c_vector(&sense, numrows, "open_cplex:3");
	i_vector(&matbeg, numrows, "open_cplex:4");
	i_vector(&matind, numnz, "open_cplex:6");
	d_vector(&matval, numnz, "open_cplex:7");

	index = 0;
	index1 = 0;
	for (h = A; h < A + K; h++){
		sense[index1] = 'E';
		rhs[index1] = 1;
		matbeg[index1++] = index;
		matind[index] = pos_w[h];
		matval[index++] = 1;
	}

	status = CPXaddrows(envPre, lpPre, 0, index1, index, rhs, sense, matbeg, matind, matval, NULL, NULL);
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

	status = CPXaddrows(envPre, lpPre, 0, index1, index, rhs, sense, matbeg, matind, matval, NULL, NULL);
	if (status)
		fprintf(stderr, "CPXaddrows failed.\n");
	free(matbeg);
	free(matind);
	free(matval);
	free(sense);
	free(rhs);

	//////pareto cut//////
	for (pareto_iter = 3; pareto_iter < (iter + 1); pareto_iter++){
		for (s = 0; s < S; s++){
			temp1[s] = 0;
		}
		for (k = 0; k < K; k++){
			for (l = 0; l < N; l++){
				for (s = 0; s < S; s++){
					if (destination[k] == l){
						temp1[s] -= d[k] * pi_pareto[(l*K + k) + (s*N*K)][pareto_iter];
					}
					else if (source[k] == l){
						temp1[s] += d[k] * pi_pareto[(l*K + k) + (s*N*K)][pareto_iter];
					}
					else
						temp1[s] += 0;
				}
			}
		}

		//// add opimality cut, single cut version////
		numrows = S*(iter + 1);
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
			for (h = 0; h < (A + K); h++){
				matind[index] = pos_w[h];
				matval[index++] = c[h] * (alpha_pareto[(h*S + s) + (pareto_iter*(A + K)*S)]);
			}
		}
		status = CPXaddrows(envPre, lpPre, 0, index1, index, rhs, sense, matbeg, matind, matval, NULL, NULL);
		if (status)
			fprintf(stderr, "CPXaddrows failed.\n");
		free(matbeg);
		free(matind);
		free(matval);
		free(sense);
		free(rhs);
	}
	//////single optimality cut
	if (iter > 0){
		for (ii = 0; ii < iter; ii++){
			for (s = 0; s < S; s++){
				temp1[s] = 0;
			}
			for (k = 0; k < K; k++){
				for (l = 0; l < N; l++){
					for (s = 0; s < S; s++){
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
				for (h = 0; h < (A + K); h++){
					matind[index] = pos_w[h];
					matval[index++] = c[h] * (alpha[(h*S + s) + (ii*(A + K)*S)]);
				}
			}
			status = CPXaddrows(envPre, lpPre, 0, index1, index, rhs, sense, matbeg, matind, matval, NULL, NULL);
			if (status)
				fprintf(stderr, "CPXaddrows failed.\n");
			free(matbeg);
			free(matind);
			free(matval);
			free(sense);
			free(rhs);
		}
	}


	//CPXwriteprob(envPre, lpPre, "model2.lp", NULL);                          //write the model in .lpPre format if needed (to debug)

	CPXsetintparam(envPre, CPX_PARAM_SCRIND, CPX_OFF); //output display
	CPXsetintparam(envPre, CPX_PARAM_THREADS, 1);
	//CPXsetintparam(envPre,CPX_PARAM_INTSOLLIM,1);    //stops after finding first integer sol.
	CPXsetintparam(envPre, CPX_PARAM_MIPDISPLAY, 3); //different levels of output display
	//CPXsetintparam(envPre,CPX_PARAM_MIPEMPHASIS,CPX_MIPEMPHASIS_OPTIMALITY);//0:balanced; 1:feasibility; 2:optimality,3:bestbound, 4:hiddenfeas
	CPXsetdblparam(envPre, CPX_PARAM_TILIM, 3600); // time limit
	//CPXsetdblpPrearam(envPre,CPX_PARAM_TRELIM, 14000); // B&B memory limit
	status = CPXsetintparam(envPre, CPX_PARAM_MEMORYEMPHASIS, 1);	//conserve memory where possible
	CPXsetintparam(envPre, CPX_PARAM_NODEFILEIND, 0);
	CPXsetdblparam(envPre, CPX_PARAM_EPGAP, 0.0001); // e-optimal solution (%gap)
	
	gettimeofday(&start, NULL);
	//startTemp = clock();
	CPXlpopt(envPre, lpPre);
	gettimeofday(&stop, NULL);
	//endTemp = clock();
	//MasterTime = (double)(endTemp - startTemp) / (double)(CLOCKS_PER_SEC);
	MasterTime = ((double)(stop.tv_sec - start.tv_sec) * 1000 + (double)(stop.tv_usec - start.tv_usec) / 1000) / 1000;


	numcols = CPXgetnumcols(envPre, lpPre);
	dj = create_double_vector(numcols);
	d_vector(&x, numcols, "open_cplex:0");
	for (i = 0; i < numcols; i++) {
		dj[i] = 0;
	}


	CPXsolution(envPre, lpPre, &status, &value, x, NULL, NULL, dj);

	best_upper_bound = value;

	if (lpPre != NULL) {
		status = CPXfreeprob(envPre, &lpPre);
		if (status) {
			fprintf(stderr, "CPXfreeprob failed, error code %d.\n", status);
		}
	}
	if (envPre != NULL) {
		status = CPXcloseCPLEX(&envPre);
		if (status) {
			char  errmsg[1024];
			fprintf(stderr, "Could not close CPLEX envPreironment.\n");
			CPXgeterrorstring(envPre, status, errmsg);
			fprintf(stderr, "%s", errmsg);
		}
	}

	int number = 0;
	index = 0;
	for (h = 0; h < A + K; h++){
		w[h] = x[index];
		//printf("w[%d]= %d\n", h, w[h]);
		index++;
	}

	for (h = 0; h < A + K; h++){
		for (k = 0; k < K; k++){
			y[h* K + k] = x[index];
			index++;
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

