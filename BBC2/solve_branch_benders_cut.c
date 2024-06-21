#include "def.h"
int	  			cur_numcols;
double			best_lb;
int				globaldepth = -1;
int				maxdepth = -1;
//double			*lbfix;  //lower bounds of variables at a given node
//double			*ubfix;// upper bounds of variables at a given node
int				last_node;
int				usercutsnum;
int				rep_node;
int				flag_upd_av;
int				terminate;
int				size_av;																					//How many of the last solutions are to be evaluated at the root node to calculate the improvement 
double			tol;
int				nodemult, depmult;
double			*average;
int				flag_incumb = 0;
double			*pos_theta;
double			*pos_w;
int				chgobj, rLimmtCBobj, rLimmtCB, tLimmtCBobj, tLimmtCB; 
int countbackLim;
int countbackObj;
//int countbackLimL;
//int countbackObjL;
int doLazy;

extern double     **pi, *alpha, *varphi, *y, *varphiVio;
extern double        *x1;
extern double        *w;
extern int        S, iter, K, N, A, *from_node_arc, *to_node_arc, *source, *destination, sglobal, max_iter, *b;
extern int        *B;
extern double     *r, *f, *p, *d, *c, M, obj, *vc, rho;
extern double     *theta, theta0, Q;
extern double     MAX_DOUBLE;
extern double     MAX_DOUBLE2;
extern double	  MasterTime, SubTime;
extern int		  *flag;
extern double     *piVio, *alphaVio;
extern double	  violation;
extern int		numCutRelax;
extern int		numCutBranchLazy;
extern int		numCutBranchuser;
clock_t			  startTemp, endTemp;
extern double gap;
extern 	int nodecount;     //Variables to call cplex
extern double checkTime;

extern double bestUB, bestLB, *dj;
extern int NumberFixedW;

char* btype;
int* bind;
double* realub;


extern FILE *Out;
extern char OutName[100];
extern char outfile[20];


int		startOptimalityCuts;
int		finishOptimalityCuts;
int		row;
double* slack;
char* sense2;


double solve_branch_benders_cut(void)
{
	int i, j, k, s, l, h;
	int index, index1;  // auxiliar indices to fill in the constraint matrix
	double best_upper_bound, best_lower_bound;
	//int nodecount;     //Variables to call cplex
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
	int    num_w_var, num_y_var, num_theta_var;
	gap=0;
	double    *pos_y;

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

	int TimeLimitation2 = 2 * 60 * 60;
	startOptimalityCuts = 0;
	finishOptimalityCuts = 0;
	depmult = 1;
	nodemult = 100;
	size_av = 2;
	tol = 0.005;
	rLimmtCBobj = 10;
	rLimmtCB = 10;
	tLimmtCBobj = 1;
	tLimmtCB = 1;
	countbackLim = 0;
	countbackObj = 0;
	//countbackLimL = 0;
	//countbackObjL = 0;
	doLazy = 1; 
	/*Declare the structures for the cut callbacks*/
	/************************************************/
	CUTINFO lazyconinfo;// These are the ones that are necessary to the model for the solution to actually be feasible. These will be verified at each integer node only.
	CUTINFO usercutinfo;     // These are the ones that are not necessary to the model for the solution to actually be feasible. These appear at every node in order to improve the lower bound
	/*user cuts info initialized*/

	average = create_double_vector(size_av); //added to calculate the average improvement
	for (i = 0; i<size_av; i++) average[i] = 0;


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
	//c_vector(&ctype, numcols, "open_cplex:01");
	ctype = create_char_vector((A + K)+ ((A + K) * K) + S);

	for (h = 0; h < A + K; h++){
		pos_w[h] = index1;
		obj[index1] = f[h];
		ctype[index1] = 'B';
		lb[index1] = 0;
		ub[index1] = 1;
		index1++;
	}

	status = CPXnewcols(env, lp, index1, obj, lb, ub, NULL, NULL);
	//status = CPXnewcols(env, lp, index1, obj, lb, ub, ctype, NULL);
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
	status = CPXnewcols(env, lp, index1, obj, lb, ub, NULL, NULL);
	//status = CPXnewcols(env, lp, index1, obj, lb, ub, ctype, NULL);
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
	
	status = CPXnewcols(env, lp, index1, obj, lb, ub, NULL, NULL);
	//status = CPXnewcols(env, lp, index1, obj, lb, ub, ctype, NULL);
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
	startOptimalityCuts = startOptimalityCuts + index1;
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
	startOptimalityCuts = startOptimalityCuts + index1;
	status = CPXaddrows(env, lp, 0, index1, index, rhs, sense, matbeg, matind, matval, NULL, NULL);
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
	startOptimalityCuts = startOptimalityCuts + index1;
	status = CPXaddrows(env, lp, 0, index1, index, rhs, sense, matbeg, matind, matval, NULL, NULL);
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
	startOptimalityCuts = startOptimalityCuts + index1;
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
	startOptimalityCuts = startOptimalityCuts + index1;
	status = CPXaddrows(env, lp, 0, index1, index, rhs, sense, matbeg, matind, matval, NULL, NULL);
	if (status)
		fprintf(stderr, "CPXaddrows failed.\n");
	free(matbeg);
	free(matind);
	free(matval);
	free(sense);
	free(rhs);

	finishOptimalityCuts = startOptimalityCuts;
	for (ii = 0; ii < iter+1; ii++){
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
		//// add opimality cut, multi cut version////
		numrows = S*(iter+1);
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
				matval[index++] = c[h] * (alpha[(h*S + s) + (ii*(A + K)*S)]);
			}
		}
		finishOptimalityCuts = finishOptimalityCuts + index1;
		status = CPXaddrows(env, lp, 0, index1, index, rhs, sense, matbeg, matind, matval, NULL, NULL);
		numCutRelax += index1;
		if (status)
			fprintf(stderr, "CPXaddrows failed.\n");
		free(matbeg);
		free(matind);
		free(matval);
		free(sense);
		free(rhs);
	}


	for (h = 0; h < A + K; h++) {
		if (bestLB + dj[h] > bestUB + (bestUB - bestLB) + 0.1) {
			i_vector(&bind, 1, "open_cplex:4");
			c_vector(&btype, 1, "open_cplex:3");
			d_vector(&realub, 1, "open_cplex:3");
			realub[0] = 0;
			btype[0] = 'U';
			bind[0] = pos_w[h];
			status = CPXchgbds(env, lp, 1, bind, btype, realub);//fix variable to maximum zero
			if (status) fprintf(stderr, "CPXchange bounds failed.\n");
			NumberFixedW++;
		}
	}


	CPXsetintparam(env, CPX_PARAM_SCRIND, CPX_OFF); //output display
	CPXsetintparam(env, CPX_PARAM_THREADS, 1);
	//CPXsetintparam(env,CPX_PARAM_INTSOLLIM,1);    //stops after finding first integer sol.
	CPXsetintparam(env, CPX_PARAM_MIPDISPLAY, 3); //different levels of output display
	CPXsetintparam(env, CPX_PARAM_MIPEMPHASIS, 2);//0:balanced; 1:feasibility; 2:optimality,3:bestbound, 4:hiddenfeas
	CPXsetdblparam(env, CPX_PARAM_TILIM, 86400 - checkTime); // time limit
	CPXsetdblparam(env, CPX_PARAM_EPGAP, 0.0001); // e-optimal solution (%gap)
	
	status = CPXsetintparam(env, CPX_PARAM_PRELINEAR, 0);											/* Do not use presolve */
	//status = CPXsetintparam (env,CPXPARAM_Preprocessing_Reduce,1);									/*No dual reductions*/
	status = CPXsetintparam(env, CPX_PARAM_MIPSEARCH, CPX_MIPSEARCH_TRADITIONAL); 					/* Turn on traditional search for use with control callbacks */
	status = CPXsetintparam(env, CPX_PARAM_MIPCBREDLP, CPX_OFF);
	status = CPXsetintparam(env, CPXPARAM_MIP_Strategy_HeuristicFreq, -1);
	status = CPXsetintparam(env, CPXPARAM_MIP_Strategy_CallbackReducedLP, CPX_OFF);


	double SumSlacks = 0;
	int iterLP = 0;
	int whileIter = 1;
	double UBlp = 999999999;
	double checkTimeImprove = 0;
	if (whileIter == 1) {
		
		j = 0;
		double OBJlp = 0;
		UBlp = 999999999;
		int iterLP = 0;
		double checkTime;
		int		NOimprovement = 0;
		double oldOBJlp = 0;
		int imp_it = 0;
		double tempRelax = 0;
		

		gettimeofday(&startImprove, NULL);
		//checkTimeImprove = ((double)(endImprove .tv_sec - startImprove.tv_sec) * 1000 + (double)(endImprove.tv_usec - startImprove.tv_usec) / 1000) / 1000;

		while (iterLP < 1 && (fabs((UBlp - OBJlp) / UBlp) >= 0.005) && checkTimeImprove < TimeLimitation2 && (NOimprovement < 1)) {// || imp_it < 1)) {

			numcols = CPXgetnumcols(env, lp);

			status = CPXlpopt(env, lp);
			i = CPXgetstat(env, lp);
			oldOBJlp = OBJlp;
			status = CPXgetobjval(env, lp, &OBJlp);
			if (OBJlp - oldOBJlp > 0.1) { imp_it++; }
			if (OBJlp - oldOBJlp < 0.1 && imp_it >= 1) { NOimprovement++; }

			d_vector(&x, numcols, "open_cplex:0");
			status = CPXgetx(env, lp, x, 0, numcols - 1);
			index = 0;
			for (h = 0; h < A + K; h++) {
				w[h] = x[index];
				index++;
			}
			for (h = 0; h < A + K; h++) {
				for (k = 0; k < K; k++) {
					y[h * K + k] = x[index];
					index++;
				}
			}
			for (s = 0; s < S; s++) {
				theta[s] = x[index];
				index++;
			}

			for (sglobal = 0; sglobal < S; sglobal++){
		SPViolated(x);
		if (flag[sglobal] == 1){
			tempRelax = 0;
			for (k = 0; k < K; k++){
				for (l = 0; l < N; l++){
					if (destination[k] == l){
						tempRelax -= d[k] * piVio[(l*K + k)];
					}
					else if (source[k] == l){
						tempRelax += d[k] * piVio[(l*K + k)];
					}
				}
			}
			numrows = 1;
			numnz = 1 + (A + K);
			d_vector(&rhs, numrows, "open_cplex:2");
			c_vector(&sense, numrows, "open_cplex:3");
			i_vector(&matbeg, numrows, "open_cplex:4");
			i_vector(&matind, numnz, "open_cplex:6");
			d_vector(&matval, numnz, "open_cplex:7");

			index = 0;
			index1 = 0;
			sense[index1] = 'G';
			rhs[index1]= tempRelax ;
			matbeg[index1++] = index;
			matind[index] = pos_theta[sglobal];
			matval[index++] = 1;
			for (h = 0; h < A + K; h++){
				matind[index] = pos_w[h];
				matval[index++] =c[h] * (alphaVio[(h)]);
			}

			finishOptimalityCuts = finishOptimalityCuts + index1;
			status = CPXaddrows(env, lp, 0, index1, index, rhs, sense, matbeg, matind, matval, NULL, NULL);
			numCutBranchuser += index1;
			if (status)
				fprintf(stderr, "CPXaddrows failed.\n");
			free(matbeg);
			free(matind);
			free(matval);
			free(sense);
			free(rhs);
		}
	}
			for (h = 0; h < A + K; h++) {
				if (w[h] > 0.01)
					w[h] = 1;
				else
					w[h] = 0;
			}
			UBlp = 0;
			for (h = 0; h < A + K; h++) {
				UBlp += f[h] * w[h];
				for (k = 0; k < K; k++) {
					UBlp += y[h * K + k] * vc[h] * rho;
				}
			}
			for (s = 0; s < S; s++) {
				UBlp += theta[s] * (1 - rho) * p[s];
			}

			gettimeofday(&endImprove, NULL);
			checkTimeImprove = ((double)(endImprove.tv_sec - startImprove.tv_sec) * 1000 + (double)(endImprove.tv_usec - startImprove.tv_usec) / 1000) / 1000;
			iterLP++;
		}


		printf("# of iteration = %d\n", iterLP);
	}
	
	status = CPXlpopt(env, lp);
	i = CPXgetstat(env, lp);


	numrows = CPXgetnumrows(env, lp);
	slack = create_double_vector(numrows);
	sense2 = create_char_vector(numrows);

	//// Get Constraint Senses:
	status = CPXgetsense(env, lp, sense2, 0, numrows - 1);


	//// Get Constraint Slacks:
	status = CPXgetslack(env, lp, slack, 0, numrows - 1);


	//// Amend slack vector so a positive slack means a non-binding constraint
	for (row = 0; row < numrows; row++) { if (sense2[row] == 'G') { slack[row] *= -1.0; } }

	
	//CPXwriteprob(env, lp, "modelBB.lp", NULL);
	j = 0;
	for (row = numrows - 1; row >= startOptimalityCuts ; row--) {
		SumSlacks += slack[row];
		j++;
	}
	iterLP = 0; 
	for (row = numrows - 1; row >= startOptimalityCuts; row--) {		// We must erase backwards to avoid loosing the indices!
		if (slack[row] > (SumSlacks/(j))) {		// We only erase non-binding constraints !
			status = CPXdelrows(env, lp, row, row);	// One at a time.
			iterLP++;
		}
	}

	Out = Open_File(outfile, "a+");
	fprintf(Out, "deleted cuts = %d \t", iterLP);
	fclose(Out);
	//printf("# of deleted cuts = %d\n\n\n", iterLP);
	
	status = CPXcopyctype(env, lp, ctype);

	/* Code to use usercut constraints*/
	/************************************************/
	status = makeusercuts(env, lp, &usercutinfo); 												//Linking the LP defined as the LP in the usercutinfo structure and get sizes
	status = CPXsetusercutcallbackfunc(env, user_cut_callback1, &usercutinfo); 							//Now create the user cuts using our cutcallback which solves a separation problem

	/* Code to use lazy constraints*/
	/************************************************/
	status = makelazyconstraint(env, lp, &lazyconinfo); 											//Linking the LP defined as the LP in the lazyconinfo structure and get sizes
	status = CPXsetlazyconstraintcallbackfunc(env, lazy_feas_callback1, &lazyconinfo);

	//CPXwriteprob(env, lp, "modelBB.lp", NULL);                          //write the model in .lp format if needed (to debug)
	gettimeofday(&start, NULL);
	//startTemp = clock();
	CPXmipopt(env, lp);  //solve the integer program
	gettimeofday(&stop, NULL);
	//endTemp = clock();
	//MasterTime = (double)(endTemp - startTemp) / (double)(CLOCKS_PER_SEC);
	MasterTime = ((double)(stop.tv_sec - start.tv_sec) * 1000 + (double)(stop.tv_usec - start.tv_usec) / 1000) / 1000;

	i = CPXgetstat(env, lp);

	if (i == 101 || i == 102 || i == 106 || i == 105 || i == 107 || i == 108)
	{
		//endTemp = clock();
		CPXgetmipobjval(env, lp, &Upper_bound);
		CPXgetbestobjval(env, lp, &best_lower_bound);
		gap = (Upper_bound - best_lower_bound) / Upper_bound;
		nodecount = CPXgetnodecnt(env, lp);
		printf("Number of BB nodes : %ld  \n", nodecount);
		//Printing to the screen
		if (i == 101 || i == 102){
			printf("Optimal Solution found with value: %.2f\n", Upper_bound);
		}
		else{
			printf("Time limit reached best solution with value:%.2f\n", Upper_bound);
		}
	}
	else{
		//endTemp = clock();
		Upper_bound = -1;
		best_lower_bound = -1;
		nodecount = -1;
		printf("Unknown stopping criterion (%d) \n", status);
	}

	// retrive solution values
	//CPXgetmipobjval(env, lp, &value);
	//printf("Upper bound: %.2f   ", value);
	//best_upper_bound = value;
	//best_upper_bound += 2 * x1F + x2F;
	// If CPLEX was able to find the optimal solution, the previous function provides the optimal solution value
	//if not, it provides the best upper bound
	//CPXgetbestobjval(env, lp, &value);  //best lower bound in case the problem was not solved to optimality
	//best_lower_bound = value;
	//printf("Lower bound: %.2f  \n", value);

	//nodecount = CPXgetnodecnt(env, lp);

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

	index = 0;
	for (h = 0; h < A + K; h++){
		w[h] = x[index];
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
	return Upper_bound;
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

/*Defining make lazycuts*/
static int makelazyconstraint(CPXENVptr  env, CPXLPptr   lp, CUTINFOptr lazyconinfo)
{
	int status = 0;
	lazyconinfo->nodeid = -1;
	cur_numcols = CPXgetnumcols(env, lp);
	lazyconinfo->lp = lp;
	lazyconinfo->numcols = cur_numcols;
	lazyconinfo->nodeobjval = 0.0;
	status = CPXgetbestobjval(env, lp, &best_lb);
	return (status);
}

/*Defining make usercuts*/
/************************************************/
static int makeusercuts(CPXENVptr  env, CPXLPptr   lp, CUTINFOptr usercutinfo)
{
	usercutinfo->nodeid = -1;
	int status = 0;
	cur_numcols = CPXgetnumcols(env, lp);
	usercutinfo->lp = lp;
	usercutinfo->nodeobjval = 0.0;
	usercutinfo->numcols = cur_numcols;
	status = CPXgetbestobjval(env, lp, &best_lb);
	return (status);
}

/*Defining the usercuts that we are going to use for fractional solutions*/

/*********************************************/
static int CPXPUBLIC user_cut_callback1(CPXCENVptr env, void *cbdata, int wherefrom, void *cbhandle, int *useraction_p)
{
	int status = 0;
	/*Declaring the local cut structure*/
	/************************************************/
	CUTINFOptr cutinfo = (CUTINFOptr)cbhandle;
	int      numcols = cutinfo->numcols;
	/************************************************/
	/*Variables used for creating the cut*/
	/************************************************/
	int      *cutind;
	double   *cutval;
	int      *cutind2;
	double   *cutval2;
	//double   rhs;
	double   cutvio;
	int      addcuts = 0;
	int		num_cuts = 0;
	int		cutnz;
	int		depth;
	int      optvio = 0;
	/************************************************/
	int      i, j, k, index = 0, SEQNUM;
	double      *temp_x;
	double	*objval;
	double   obval;
	double	change;
	int		Flag_feas = 1;
	int		Found_Vio = 0;
	int		com_check;
	double value;
	/************************************************/
	double    rhs;    // right and side of constraints ................................
	char      *sense;  // constraints sense (<=: 'L', =:'E', >=:'G') ...................
	int       *matbeg; // index of first non-zero element in each row...................
	int       *matind; // associated column of each non-zelo element ...................
	double    *matval; // coefficient values fo the non-zero elements of constraints....
	double		temp1;
	int			l, h;
	//int       numcols; // number of variables ..........................................
	int       numrows; // number of constraints.........................................
	int       numnz;   // number of non-zero elements in the matrix ....................
	int		 index1;  // auxiliar indices to fill in the constraint matrix
	int		flagVio = 0;
	
	objval = &obval;
	
	d_vector(&temp_x, cur_numcols, "open_cplex:0");											//Create the variable that will be used to store the solution at the node
	objval = &obval;																		//Assign the memory slot to the pointer so that we may use it
	numcols = cur_numcols;																	//Transfer the information on the number of variables which we obtained at the make lazy constraint callback
	*useraction_p = CPX_CALLBACK_DEFAULT;													//At this point, don't add any cuts. We are only starting to solve the separation problem

        int oldnodeid = SEQNUM;
	double oldnodeobjval = obval;
	status = CPXgetcallbacknodex(env, cbdata, wherefrom, temp_x, 0, numcols - 1);			//Obtain the solution at the node
	status = CPXgetcallbacknodeobjval(env, cbdata, wherefrom, objval);	
	status = CPXgetcallbacknodeobjval(env, cbdata, wherefrom, &cutinfo->nodeobjval);						//Obtain the objective value of the solution at the curent node
	status = CPXgetcallbacknodeinfo(env, cbdata, wherefrom, 0, CPX_CALLBACK_INFO_NODE_SEQNUM, &SEQNUM);
	
	status = CPXgetcallbacknodeinfo(env, cbdata, wherefrom, 0, CPX_CALLBACK_INFO_NODE_SEQNUM, &cutinfo->nodeid);
	if (status) {
		fprintf(stderr, "Failed to get node id.\n");
		goto Terminate;
	}
	
	status = CPXgetcallbacknodeinfo(env, cbdata, wherefrom, 0, CPX_CALLBACK_INFO_NODE_DEPTH, &depth);
	if (status) {
		fprintf(stderr, "Failed to get depth.\n");
		goto Terminate;
	}
	
	//printf("Best lb is %lf\n", best_lb);getchar();
	globaldepth = depth;
	if (globaldepth>maxdepth) maxdepth = globaldepth;
	
	if (oldnodeid == SEQNUM) {
		countbackLim++;
		double objchg = ((obval - oldnodeobjval) / oldnodeobjval);
		if (objchg < tol) {
			countbackObj++;
		}
		else {
			countbackObj = 0;
		}
	}
	else {
		countbackLim = 0;
		countbackObj = 0;
	}
	if (SEQNUM == 0) {
		if (countbackObj > rLimmtCBobj || countbackLim > rLimmtCB) {
			*useraction_p = CPX_CALLBACK_ABORT_CUT_LOOP;
			goto Terminate;
		}
	}
	else {
		if (countbackObj > tLimmtCBobj || countbackLim > tLimmtCB) {
			*useraction_p = CPX_CALLBACK_ABORT_CUT_LOOP;
			goto Terminate;
		}
	}
	if (depth > depmult) {
		//if (SEQNUM % nodemult > 0){
			*useraction_p = CPX_CALLBACK_ABORT_CUT_LOOP;
			goto Terminate;
		//}
	}
	violation = 0.0001;
	
	/*else{ violation=LPval*.1;}*/
	//cutsettol = 0.9999; //Cutset inequality

	//printf("Opened the following arcs\n");
	/*for (j = 0; j<M; j++){                  //fill up the values of ybar_frac
	ybar_frac[edges[j].i][edges[j].j] = temp_x[pos_y[edges[j].i][edges[j].j]];
	if (depth == 0) LPsol[j] = temp_x[pos_y[edges[j].i][edges[j].j]];//Keeping the LP solution value
	}*/
	//printf("\n");
	//	if (1 == 1){

	/*if (num_cuts > 0 && Flag_feas == 1) {
	*useraction_p = CPX_CALLBACK_SET;
	//printf("Exited Callback\n");
	goto Terminate;
	}
	else{*/
	//printf("\nEntered to check for optimality\n"); getchar();
	rep_node++;
	//				flag_upd_av = 1;
	/******Updating the stabilizer only if there has not been an incumbent yet********/
	/*if (depth>0 && flag_incumb == 0){
	for (j = 0; j<M; j++){
	stabilizer[j] = ini_core[j];
	core_point[edges[j].i][edges[j].j] = ini_core[j];
	}
	flag_incumb = 1;
	}*/

	///////////SHAAAAAAAAAAAAAAAAAAAAAAAABIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII//////
	for (sglobal = 0; sglobal < S; sglobal++){
		cutvio = SPViolated(temp_x);
		if (flag[sglobal] == 1){
			temp1 = 0;
			for (k = 0; k < K; k++){
				for (l = 0; l < N; l++){
					if (destination[k] == l){
						temp1 -= d[k] * piVio[(l*K + k)];
					}
					else if (source[k] == l){
						temp1 += d[k] * piVio[(l*K + k)];
					}
				}
			}
			numrows = 1;
			numnz = 1 + (A + K);
			//d_vector(&rhs, numrows, "open_cplex:2");
			//c_vector(&sense, numrows, "open_cplex:3");
			i_vector(&matbeg, numrows, "open_cplex:4");
			i_vector(&matind, numnz, "open_cplex:6");
			d_vector(&matval, numnz, "open_cplex:7");

			index = 0;
			index1 = 0;
			//sense[index1] = 'G';
			rhs = temp1;
			matbeg[index1++] = index;
			matind[index] = pos_theta[sglobal];
			matval[index++] = 1;
			for (h = 0; h < A + K; h++){
				matind[index] = pos_w[h];
				matval[index++] =c[h] * (alphaVio[(h)]);
			}

			//status = CPXaddrows(env, lp, 0, index1, index, rhs, sense, matbeg, matind, matval, NULL, NULL);
			status = CPXcutcallbackadd(env, cbdata, wherefrom, index, rhs, 'G', matind, matval, CPX_USECUT_FORCE);
			numCutBranchuser++;
			//status = CPXcutcallbackadd(env, cbdata, wherefrom, optcuts[k].numnz, optcuts[k].RHS, 'G', optcuts[k].ind, optcuts[k].coeff, CPX_USECUT_PURGE);

			if (status)
				fprintf(stderr, "CPXaddrows failed.\n");
			free(matbeg);
			free(matind);
			free(matval);
			//free(sense);
			//free(rhs);
			flagVio = 1;
		}
	}
	/////////////SHAAAAAAAAAAAAAAAAAAAAAAAABIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII///////////////

	/*
	if (flagVio>0){
	for (k = optcount - vioptcuts; k<optcount; k++){					//Add the violated cuts
	status = CPXcutcallbackadd(env, cbdata, wherefrom, optcuts[k].numnz, optcuts[k].RHS, 'G', optcuts[k].ind, optcuts[k].coeff, CPX_USECUT_PURGE);
	//status = CPXcutcallbackaddlocal (env, cbdata, wherefrom,  optcuts[k].numnz, optcuts[k].RHS, 'G',  optcuts[k].ind, optcuts[k].coeff);
	cuts_SC[3]++;
	}
	}*/

	if (flagVio>0){
		*useraction_p = CPX_CALLBACK_SET;
		usercutsnum++;
		//printf("added %d optimality cuts at fractional solution\n", vioptcuts);
		goto Terminate;
	}
	goto Terminate;
	//}
	//}
Terminate:

//Out = Open_File(outfile, "a+");
//fprintf(Out, "check2");
//fclose(Out);
	/*printf("Exited fractional separation\n");
	if(obval>=3088089.25)getchar();*/
	if (*useraction_p != CPX_CALLBACK_SET){  /*printf("Exited without any violation \n");*/ /*getchar();*/ terminate = 1; }
	//if (Flag_Print == 1){ printf("\n"); /*getchar();*/ }
	free(temp_x);
	//free(cutval);
	//free(cutind);
	return (status);
}

static int CPXPUBLIC lazy_feas_callback1(CPXCENVptr env, void *cbdata, int wherefrom, void *cbhandle, int *useraction_p)
{
	int status = 0;
	/*Declaring the local cut structure*/
	/************************************************/
	CUTINFOptr cutinfo = (CUTINFOptr)cbhandle;
	int      numcols = cutinfo->numcols;
	int cutnz;
	/************************************************/
	/*Variables used for creating the cut*/
	/************************************************/
	int      *cutind;
	double   *cutval;
	int      *cutind2;
	double   *cutval2;
	double   cutvio;
	int      addcuts = 0;
	int		num_cuts = 0;
	int		addcuts2 = 0;
	int      optvio = 0;
	/************************************************/
	int      i, j, k, index = 0, SEQNUM;
	double      *temp_x;
	double	*objval;
	double   obval;
	int			depth;
	int			com_check; //Indicator as to what commodity will be brought to trial
	/*Obtaining the solution information necessary at the current node*/
	/************************************************/
	double    rhs;    // right and side of constraints ................................
	char      *sense;  // constraints sense (<=: 'L', =:'E', >=:'G') ...................
	int       *matbeg; // index of first non-zero element in each row...................
	int       *matind; // associated column of each non-zelo element ...................
	double    *matval; // coefficient values fo the non-zero elements of constraints....
	double		temp1;
	int			l, h;
	//int       numcols; // number of variables ..........................................
	int       numrows; // number of constraints.........................................
	int       numnz;   // number of non-zero elements in the matrix ....................
	int		 index1;  // auxiliar indices to fill in the constraint matrix
	int		flagVio = 0;
	/************************************************/

	//cutind = create_int_vector(M + K);
	//cutval = create_double_vector(M + K);
	d_vector(&temp_x, cur_numcols, "open_cplex:0");
	//One percent test
	/********************************************************/
	//if (Flag_Print == 1){ printf("Entered lazy cutcallback\n"); }
	//printf("Entered lazy cutcallback\n"); getchar();
	//flag_frac = 0;
	//last_feas = feascount;
	/*flag_MWeval=0;*/															//Create the variable that will be used to store the solution at the node
	objval = &obval;																							//Assign the memory slot to the pointer so that we may use it
	numcols = cur_numcols;	//Transfer the information on the number of variables which we obtained at the make lazy constraint callback
	//lazycalled++;
	*useraction_p = CPX_CALLBACK_DEFAULT;
	
	int    oldnodeid = cutinfo->nodeid;
	double oldnodeobjval = cutinfo->nodeobjval;
	

	//At this point, don't add any cuts. We are only starting to solve the separation problem
	status = CPXgetcallbacknodex(env, cbdata, wherefrom, temp_x, 0, numcols - 1);									//Obtain the solution at the node
	status = CPXgetcallbacknodeobjval(env, cbdata, wherefrom, objval);
	status = CPXgetcallbacknodeobjval(env, cbdata, wherefrom, &cutinfo->nodeobjval);		//Obtain the objective value of the solution at the curent node
	status = CPXgetcallbacknodeinfo(env, cbdata, wherefrom, 0, CPX_CALLBACK_INFO_NODE_SEQNUM, &SEQNUM);
	
	status = CPXgetcallbacknodeinfo(env, cbdata, wherefrom, 0, CPX_CALLBACK_INFO_NODE_SEQNUM, &cutinfo->nodeid);
	if (status) {
		fprintf(stderr, "Failed to get node id.\n");
		goto Terminate;
	}
	status = CPXgetcallbacknodeinfo(env, cbdata, wherefrom, 0, CPX_CALLBACK_INFO_NODE_DEPTH, &depth);
	if (status) {
		fprintf(stderr, "Failed to get depth.\n");
		goto Terminate;
	}
		for (sglobal = 0; sglobal < S; sglobal++){
			cutvio = SPViolated(temp_x);
			if (flag[sglobal] == 1){
				temp1 = 0;
				for (k = 0; k < K; k++){
					for (l = 0; l < N; l++){
						if (destination[k] == l){
							temp1 -= d[k] * piVio[(l*K + k)];
						}
						else if (source[k] == l){
							temp1 += d[k] * piVio[(l*K + k)];
						}
					}
				}
				numrows = 1;
				numnz = 1 + (A + K);
				i_vector(&matbeg, numrows, "open_cplex:4");
				i_vector(&matind, numnz, "open_cplex:6");
				d_vector(&matval, numnz, "open_cplex:7");
				index = 0;
				index1 = 0;
				rhs = temp1;
				matbeg[index1++] = index;
				matind[index] = pos_theta[sglobal];
				matval[index++] = 1;
				for (h = 0; h < A + K; h++){
					matind[index] = pos_w[h];
					matval[index++] = c[h] * (alphaVio[(h)]);
				}
				status = CPXcutcallbackadd(env, cbdata, wherefrom, index, rhs, 'G', matind, matval, CPX_USECUT_FORCE);
				numCutBranchLazy++;
				if (status)
					fprintf(stderr, "CPXaddrows failed.\n");
				free(matbeg);
				free(matind);
				free(matval);
				flagVio = 1;
			}
		}
		if (flagVio > 0){
			*useraction_p = CPX_CALLBACK_SET;
			//printf("added %d optimality cuts at integer solution\n", vioptcuts);
			goto Terminate;
		}
		goto Terminate;
	//}
Terminate:
	/******Updating the stabilizer only if there we found a new integer solution********/
	/*if(*useraction_p!=CPX_CALLBACK_SET){
	for(j=0;j<M;j++){
	stabilizer[j]=temp_x[j];
	core_point[edges[j].i][edges[j].j]=temp_x[j];
	}
	}
	flag_incumb=1;*/
	/****************************************/
	//if (base_fact - decrement*(disc) == 0)disc = 0;
	free(temp_x);
	//free(cutval);
	//free(cutind);
	return (status);
}
/* End of lazycut callback */
