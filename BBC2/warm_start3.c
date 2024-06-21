#include "def.h"

extern double** pi, * alpha, * varphi, * y;
extern double* x1;
extern double* w;
extern int  S, iter, K, N, A, * from_node_arc, * to_node_arc, * source, * destination, sglobal, max_iter, * b;
extern int* B;
extern double* r, * f, * p, * d, * c, M, obj, * vc, rho;
extern double* theta, theta0, Q;
extern double     MAX_DOUBLE;
extern double     MAX_DOUBLE2;
extern double	  MasterTime, SubTime;
extern int* flag;
extern double	violation;
clock_t	   startTemp, endTemp;
extern double* dh;


double warm_start3(void)
{
	int i, j, k, s, l, h;
	int index, index1;  // auxiliar indices to fill in the constraint matrix
	double best_upper_bound, best_lower_bound;
	int nodecount;     //Variables to call cplex
	CPXLPptr  lpWarm;      // data strucutre to store a problem in cplex ...................
	CPXENVptr envWarm;     // cplex environment.............................................
	int       numcols; // number of variables ..........................................
	int       numrows; // number of constraints.........................................
	int       numnz;   // number of non-zero elements in the matrix ....................
	int       objsen;  // optimization sense (min:1, max:-1 ) ..........................
	double* obj;    // objective function coefficients ..............................
	double* rhs;    // right and side of constraints ................................
	char* sense;  // constraints sense (<=: 'L', =:'E', >=:'G') ...................
	int* matbeg; // index of first non-zero element in each row...................
	int* matind; // associated column of each non-zelo element ...................
	double* matval; // coefficient values fo the non-zero elements of constraints....
	double* lb;     // lower bounds of variables.....................................
	double* ub;     // upper bounds of variables.....................................
	int       status;  // optimization status......................... .................
	double* x;      // solution vector (double, even if the problem is integer) .....
	char      probname[16]; // problem name for cplex .......................................
	char* ctype;  // variable type ('C', 'I', 'B') only if integer.................
	double    value;   // objevtive value of solution ..................................
	double    num_w_var, num_y_var, num_theta_var;
	double* pos_theta;
	double* pos_y;
	int* pos_w;
	double tempV = 0;
	int jcounter = 0;
	int ii = 0;
	double* temp1;
	double* temp2, * temp3;


	int counter = 0;
	char* btype;
	int* bind;
	double* realub;


	pos_theta = create_double_vector(S);
	pos_y = create_double_vector((A + K) * K);
	pos_w = create_int_vector(A + K);
	temp1 = create_double_vector(S);
	temp2 = create_double_vector((A + K) * S);

	//Initialize CPLEX envPreironment
	envWarm = CPXopenCPLEX(&status);
	if (envWarm == NULL) {
		char  errmsg[1024];
		printf("Could not open CPLEX. \n");
		CPXgeterrorstring(envWarm, status, errmsg);
		printf("%s", errmsg);
	}
	// Create the problem in CPLEX 
	strcpy(probname, "UFlpPre");
	lpWarm = CPXcreateprob(envWarm, &status, probname);
	if (envWarm == NULL) {
		char  errmsg[1024];
		printf("Could not create lpWarm. \n");
		CPXgeterrorstring(envWarm, status, errmsg);
		printf("%s", errmsg);
	}
	CPXchgobjsen(envWarm, lpWarm, CPX_MIN);



	//Define w_h variables
	index1 = 0;  // index of columns
	numcols = A + K;
	d_vector(&obj, numcols, "open_cplex:1");
	d_vector(&lb, numcols, "open_cplex:8");
	d_vector(&ub, numcols, "open_cplex:9");
	c_vector(&ctype, numcols, "open_cplex:01");

	for (h = 0; h < A + K; h++) {
		pos_w[h] = index1;
		obj[index1] = f[h];
		ctype[index1] = 'B';
		lb[index1] = 0;
		ub[index1] = 1;
		index1++;
	}

	status = CPXnewcols(envWarm, lpWarm, index1, obj, lb, ub, ctype, NULL);
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

	for (h = 0; h < A + K; h++) {
		for (k = 0; k < K; k++) {
			pos_y[h * K + k] = index1 + num_w_var;
			obj[index1] = vc[h] * rho;
			ctype[index1] = 'C';
			lb[index1] = 0;
			ub[index1] = CPX_INFBOUND;
			index1++;
		}
	}
	status = CPXnewcols(envWarm, lpWarm, index1, obj, lb, ub, ctype, NULL);
	if (status)
		fprintf(stderr, "CPXnewcols failed.\n");
	free(obj);
	free(lb);
	free(ub);
	free(ctype);
	num_y_var = index1;

	//Add constraint 1 // outward flow - inward fow = 0, d, -d
	numrows = K * N;
	numnz = 2 * (A + K) * K * N;
	d_vector(&rhs, numrows, "open_cplex:2");
	c_vector(&sense, numrows, "open_cplex:3");
	i_vector(&matbeg, numrows, "open_cplex:4");
	i_vector(&matind, numnz, "open_cplex:6");
	d_vector(&matval, numnz, "open_cplex:7");

	index = 0;
	index1 = 0;
	for (k = 0; k < K; k++) {
		for (l = 0; l < N; l++) {
			sense[index1] = 'G';
			if (destination[k] == l)
				rhs[index1] = -d[k];
			else if (source[k] == l)
				rhs[index1] = d[k];
			else
				rhs[index1] = 0;
			matbeg[index1++] = index;
			for (h = 0; h < A + K; h++) {
				matind[index] = pos_y[h * K + k];
				if (from_node_arc[h] == l) {
					matval[index++] = 1;
				}
				else if (to_node_arc[h] == l) {
					matval[index++] = -1;
				}
			}
		}
	}

	status = CPXaddrows(envWarm, lpWarm, 0, index1, index, rhs, sense, matbeg, matind, matval, NULL, NULL);
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
	for (h = 0; h < A + K; h++) {
		sense[index1] = 'L';
		rhs[index1] = 0;
		matbeg[index1++] = index;
		for (k = 0; k < K; k++) {
			matind[index] = pos_y[h * K + k];
			matval[index++] = 1;
		}
		matind[index] = pos_w[h];
		matval[index++] = -1 * c[h];
	}

	status = CPXaddrows(envWarm, lpWarm, 0, index1, index, rhs, sense, matbeg, matind, matval, NULL, NULL);
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
	for (h = A; h < A + K; h++) {
		sense[index1] = 'E';
		rhs[index1] = 1;
		matbeg[index1++] = index;
		matind[index] = pos_w[h];
		matval[index++] = 1;
	}

	status = CPXaddrows(envWarm, lpWarm, 0, index1, index, rhs, sense, matbeg, matind, matval, NULL, NULL);
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
	for (h = 0; h < A; h++) {
		matind[index] = pos_w[h];
		matval[index++] = 1;
	}

	status = CPXaddrows(envWarm, lpWarm, 0, index1, index, rhs, sense, matbeg, matind, matval, NULL, NULL);
	if (status)
		fprintf(stderr, "CPXaddrows failed.\n");
	free(matbeg);
	free(matind);
	free(matval);
	free(sense);
	free(rhs);

	//__________________________________________________________________________________________________________________________________________________
	//CPXwriteprob(envWarm, lpWarm, "model2_fixed.lp", NULL);

	Sorted = (AVAL*)calloc(A, sizeof(AVAL));
	
	counter = 0;
	for (h = 0; h < A; h++) {
		//w_hat[h] = w[h];
		Sorted[counter].index = 0;
		for (k = 0; k < K; k++) {
			Sorted[counter].value += y[h * K + k] + w[h];
			if (y[h * K + k] > 0) {
				Sorted[counter].index++;
			}
		}
		Sorted[counter].index2 = h;
		counter++;
	}
	qsort((AVAL*)Sorted, counter, sizeof(Sorted[0]), Comparevalue);
	for (h = 0; h < (int)(A * 0.2); h++) {
		i_vector(&bind, 1, "open_cplex:4");
		c_vector(&btype, 1, "open_cplex:3");
		d_vector(&realub, 1, "open_cplex:3");
		realub[0] = 0;
		btype[0] = 'U';
		bind[0] = pos_w[Sorted[h].index2];
		status = CPXchgbds(envWarm, lpWarm, 1, bind, btype, realub);//fix variable to maximum zero
		if (status) fprintf(stderr, "CPXchange bounds failed.\n");

	}



	//CPXwriteprob(envWarm, lpWarm, "model2_fixed.lp", NULL);                          //write the model in .lpWarm format if needed (to debug)

	CPXsetintparam(envWarm, CPX_PARAM_SCRIND, CPX_OFF); //output display
	CPXsetintparam(envWarm, CPX_PARAM_THREADS, 1);
	//CPXsetintparam(envWarm,CPX_PARAM_INTSOLLIM,1);    //stops after finding first integer sol.
	CPXsetintparam(envWarm, CPX_PARAM_MIPDISPLAY, 3); //different levels of output display
	//CPXsetintparam(envWarm,CPX_PARAM_MIPEMPHASIS,1);//0:balanced; 1:feasibility; 2:optimality,3:bestbound, 4:hiddenfeas
	CPXsetdblparam(envWarm, CPX_PARAM_TILIM, 3600); // time limit
	//CPXsetdblpPrearam(envPre,CPX_PARAM_TRELIM, 14000); // B&B memory limit
	//status = CPXsetintparam(envWarm, CPX_PARAM_MEMORYEMPHASIS, 1);	//conserve memory where possible
	CPXsetintparam(envWarm, CPX_PARAM_NODEFILEIND, 0);
	CPXsetintparam(envWarm, CPX_PARAM_MIPEMPHASIS, CPX_MIPEMPHASIS_OPTIMALITY);  // MIP emphasis: optimality, feasibility, moving best bound
	CPXsetdblparam(envWarm, CPX_PARAM_EPGAP, 0.005); // e-optimal solution (%gap)
	CPXsetintparam(envWarm, CPX_PARAM_PRELINEAR, 0);

	status = CPXsetdblparam(envWarm, CPX_PARAM_WORKMEM, 12288);			//trigger tree size in MB to create node file (Default 128 MB)


	gettimeofday(&start, NULL);
	//startTemp = clock();
	CPXmipopt(envWarm, lpWarm);  //solve the integer program
	//CPXlpopt(envWarm, lpWarm);
	gettimeofday(&stop, NULL);
	//endTemp = clock();
	//MasterTime = (double)(endTemp - startTemp) / (double)(CLOCKS_PER_SEC);
	MasterTime = ((double)(stop.tv_sec - start.tv_sec) * 1000 + (double)(stop.tv_usec - start.tv_usec) / 1000) / 1000;

	i = CPXgetstat(envWarm, lpWarm);
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
	CPXgetmipobjval(envWarm, lpWarm, &value);
	//printf("Upper bound: %.2f   ", value);
	best_upper_bound = value;

	//best_upper_bound += 2 * x1F + x2F;
	// If CPLEX was able to find the optimal solution, the previous function provides the optimal solution value
	//if not, it provides the best upper bound
	CPXgetbestobjval(envWarm, lpWarm, &value);  //best lower bound in case the problem was not solved to optimality
	best_lower_bound = value;
	//printf("Lower bound: %.2f  \n", value);

	nodecount = CPXgetnodecnt(envWarm, lpWarm);
	//printf("Number of BB nodes : %ld  \n", nodecount);

	numcols = CPXgetnumcols(envWarm, lpWarm);
	d_vector(&x, numcols, "open_cplex:0");
	CPXgetmipx(envWarm, lpWarm, x, 0, numcols - 1);  // obtain the values of the decision variables

	if (lpWarm != NULL) {
		status = CPXfreeprob(envWarm, &lpWarm);
		if (status) {
			fprintf(stderr, "CPXfreeprob failed, error code %d.\n", status);
		}
	}
	if (envWarm != NULL) {
		status = CPXcloseCPLEX(&envWarm);
		if (status) {
			char  errmsg[1024];
			fprintf(stderr, "Could not close CPLEX envWarmironment.\n");
			CPXgeterrorstring(envWarm, status, errmsg);
			fprintf(stderr, "%s", errmsg);
		}
	}

	int number = 0;
	index = 0;
	for (h = 0; h < A + K; h++) {
		w[h] = x[index];
		index++;
	}

	for (h = 0; h < A + K; h++) {
		for (k = 0; k < K; k++) {
			y[h * K + k] = x[index];
			/*if (y[h * K + k] > 0) {
				printf("y[%d][%d]= %lf\n", h, k, y[h * K + k]);
			}*/
			index++;
		}
	}
	
	free(Sorted);
	free(x);
	return best_upper_bound;
}

