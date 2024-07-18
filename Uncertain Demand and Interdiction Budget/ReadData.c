#include "def.h"


extern double** pi, * alpha, * varphi, * y;
extern double* x1;
extern double* w;
extern int        S, Sperim, iter, K, N, A, * from_node_arc, * to_node_arc, * source, * destination, sglobal, max_iter, * b;
extern int* B;
extern double* r, * f, * p, * d, * c, M, obj, * vc, rho, *demand, * pPerim;
extern double* theta, theta0, Q;
extern double     MAX_DOUBLE;
extern double     MAX_DOUBLE2;
extern double	  MasterTime, SubTime;
extern int* flag;
clock_t	   startTemp, endTemp;
extern double* piVio, * alphaVio;
extern double	violation;

// Function used to read the input parameters of the instance being solved

void read_instance(const char* name)
{
	int i, j, k, s, l, h, ss;
	FILE* in;

	in = Open_File(name, "r");

	if (fscanf(in, "%d %d %d %d %d", &N, &A, &K, &S, &Sperim) != 5) {
		fprintf(stderr, "ERROR: Cannot read instance size \n");
		exit(1);
	}
	//printf("Instance size: N:%d M:%d \n", N, M);

	Initialize_memory();

	for (h = 0; h < A; h++) {
		fscanf(in, "%d", &from_node_arc[h]);
	}
	for (h = 0; h < A; h++) {
		fscanf(in, "%d", &to_node_arc[h]);
	}
	for (h = 0; h < A; h++) {
		fscanf(in, "%lf", &vc[h]);
	}
	for (h = 0; h < A; h++) {//capacity
		fscanf(in, "%lf", &c[h]);
	}
	for (h = 0; h < A; h++) {
		fscanf(in, "%lf", &f[h]);
	}
	for (k = 0; k < K; k++) {
		fscanf(in, "%d", &source[k]);
	}
	for (k = 0; k < K; k++) {
		fscanf(in, "%d", &destination[k]);
	}
	for (ss = 0; ss < Sperim; ss++) {
		fscanf(in, "%lf", &pPerim[ss]);
		for (k = 0; k < K; k++) {
			fscanf(in, "%lf", &d[k * Sperim + ss]);
		}
	}
	for (h = 0; h < A; h++) {
		for (k = 0; k < K; k++) {
			fscanf(in, "%lf", &r[h * K + k]);
		}
	}
	for (s = 0; s < S; s++) {
		fscanf(in, "%lf", &p[s]);
	}
	for (s = 0; s < S; s++) {
		fscanf(in, "%d", &B[s]);
	}
	fscanf(in, "%lf", &rho);
	for (h = 0; h < A; h++) {
		//fscanf(in, "%d", &b[h]);
		b[h] = 1;
	}
	fclose(in);
}

//Function used to allocate memory to arrays based on the size of the input parameter of the specific instance being solved

void Initialize_memory(void)
{
	max_iter = 1000;
	d = create_double_vector(K * Sperim);
	c = create_double_vector(A + K);
	vc = create_double_vector(A + K);
	y = create_double_vector((A + K) * K * Sperim);
	f = create_double_vector(A + K);
	p = create_double_vector(S);
	r = create_double_vector((A + K) * K);
	destination = create_int_vector(K);
	source = create_int_vector(K);
	from_node_arc = create_int_vector(A + K);
	to_node_arc = create_int_vector(A + K);
	x1 = create_double_vector((A + K) * S);
	w = create_double_vector(A + K);
	alpha = create_double_vector((A + K) * S * max_iter);
	pi = create_double_matrix(N * K * S, max_iter);
	varphi = create_double_vector((A + K) * S * max_iter);
	B = create_int_vector(S);
	theta = create_double_vector(S);
	b = create_int_vector(A + K);
	flag = create_int_vector(S);
	demand = create_double_vector(K);
	piVio = create_double_vector(N * K);
	alphaVio = create_double_vector(A + K);
	pPerim = create_double_vector(Sperim);
	//piVio = create_double_vector(N * K * S);
	//alphaVio = create_double_vector((A + K) * S);
}


//Function used to release memory of the arrays based on the size of the input parameter of the specific instance being solved


void free_memory(void)
{
	int i;
}


// Function used to open data file

FILE* Open_File(const char* name, const char* mode)  // This function opens the file "name" under mode "mode" (reading, writing, etc)
{
	FILE* file;

	if ((file = fopen(name, mode)) == NULL) {
		printf("\nError: File cannot be opened \n");
		//OK=1;
		exit(8);
	}
	return file;
}

// Functions to allocate memory to one and two dimensional arrays


int** create_int_matrix(int rows, int Columns)
{
	int i;
	int** ptr;

	if ((ptr = (int**)calloc(rows, sizeof(int*))) == NULL) {
		printf("\nError: Insuficient memory \n");
		exit(8);
	}
	for (i = 0; i < rows; i++)
		ptr[i] = create_int_vector(Columns);
	return ptr;
}

double** create_double_matrix(int rows, int Columns)
{
	int i;
	double** ptr;

	if ((ptr = (double**)calloc(rows, sizeof(double*))) == NULL) {
		printf("\nError: Insuficient memory \n");
		exit(8);
	}
	for (i = 0; i < rows; i++) {
		ptr[i] = create_double_vector(Columns);
	}
	return ptr;
}

double*** create_double_matrix3D(int rows, int Columns, int Columns2)
{
	int i;
	double*** ptr;

	if ((ptr = (double***)calloc(rows, sizeof(double*))) == NULL) {
		printf("\nError: Insuficient memory \n");
		exit(8);
	}
	for (i = 0; i < rows; i++) {
		ptr[i] = create_double_matrix(Columns, Columns2);
	}
	return ptr;
}


int* create_int_vector(int dim)
{
	int* ptr;

	if ((ptr = (int*)calloc(dim, sizeof(int))) == NULL) {
		printf("\nError: Insuficient memory \n");
		exit(8);
	}
	return ptr;
}

double* create_double_vector(int dim)
{
	double* ptr;

	if ((ptr = (double*)calloc(dim, sizeof(double))) == NULL) {
		printf("\nError: Insuficient memory \n");
		exit(8);
	}
	return ptr;
}


// CPLEX functions to allocate memeory to arrays

void i_vector(int** vector, int n, char* s)
{
	if ((*vector = (int*)calloc(n, sizeof(int))) == NULL)
		//error(s);
		printf("Error: Insuficient memory \n");
	return;
}

void d_vector(double** vector, int n, char* s)
{
	if ((*vector = (double*)calloc(n, sizeof(double))) == NULL)
		// error(s);
		printf("Error: Insuficient memory \n");
	return;
}

void c_vector(char** vector, int n, char* s)
{
	if ((*vector = (char*)calloc(n, sizeof(char))) == NULL)
		//error(s);
		printf("Error: Insuficient memory \n");
	return;
}

char* create_char_vector(int cells) {
	char* ptr = (char*)calloc(cells, sizeof(char));
	if (ptr == NULL) {
		fprintf(stderr, "ERROR: Unable to allocate memory for char_vector\n");
	}
	return ptr;
}

