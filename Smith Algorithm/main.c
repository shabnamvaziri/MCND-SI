#include "def.h"

double     **pi, *alpha, *varphi, *y;
int        *x1;
int     *w;
int        S, iter, K, N, A, *from_node_arc, *to_node_arc, *source, *destination, sglobal, max_iter, *b;
int        *B;
double     *r, *f, *p, *d, *c, M, obj, *vc, rho;
double     *theta, theta0, Q;
double     MAX_DOUBLE = 10000000000;
double     MAX_DOUBLE2 = 10000000;
double     timeTotal = 0;
double	   MasterTime, SubTime;
clock_t	   startTemp, endTemp;
double TimeLimitation = 48 * 60 * 60;
double checkTime;
double installation, pre;
char OutName[100];
FILE *Out = NULL;
char outfile[20];
double tempgap = 0;

void main(int argc, char *argv[])
{
	char	instance[20];
	char	path[50];
	FILE		*ini;
	clock_t  start, end;
	double obj_value_ch, obj_value_ls, opt_value;
	double cputime;
	int i, j, k, s, h, l, MaxNumInst, numSample;
	start = clock();
	M = 0;
	if (argc == 1) {
		printf("Error: Input file not specified \n");
		exit(8);
	}
	ini = Open_File(argv[1], "r");
	fscanf(ini, "%d", &MaxNumInst);
	fscanf(ini, "%s", &outfile);

	ttt = time(NULL);
	tm = localtime(&ttt);
	Out = Open_File(outfile, "a+");
	fprintf(Out, "\n %s\n", asctime(tm));
	fclose(Out);

	for (numSample = 1; numSample < MaxNumInst; numSample++){
		fscanf(ini, "%s", &instance);
		sprintf(path, "./Data/");
		strcat(path, instance);
		
		read_instance(path);
		//read_instance("Input1.txt");
		
	
		Out = Open_File(outfile, "a+");					//Writing what instance we are solving
		fprintf(Out, "\n%s | %s ;\t", argv[1], instance);
		fprintf(Out, "%d\t%d\t%d\t%d\t%lf\t|\t", N, A, K,S,rho);
		fclose(Out);
		double tempM = 0;
		double TotalDemand = 0;
		for (k = 0; k < K; k++){
			TotalDemand += d[k];
		}
		//subproblem is lowerbound and master problem is upperbound. The gap = (UB-LB)/LB

		for (k = 0; k < K; k++){
			from_node_arc[A + k] = source[k];
			to_node_arc[A + k] = destination[k];
			c[A + k] = TotalDemand;
			vc[A + k] = 10000;
			f[A + k] = 0;
			r[(A + k)*K + k] = 0;
			b[A + k] = 1000;
		}
		gettimeofday(&startTotal, NULL);
		for (h = 0; h < A + K; h++){
			M += vc[h];
		}

		double Total_theta = 0;
		//read_instance("Input1.txt");
		for (s = 0; s < S; s++){
			theta[s] = 0;
			Total_theta += theta[s];
		}
		theta0 = -1 * MAX_DOUBLE;
		float xx = 0;
		Q = 10;
		iter = 0;
		checkTime = 0;
		while (fabs((Q - Total_theta) / Q) >= 0.000001 && checkTime < TimeLimitation){
			Total_theta = 0;
			installation = 0;
			 pre = 0;
			obj = solve_MasterProblem();
			for (s = 0; s < S; s++){
				Total_theta += theta[s]*p[s];
			}
			Q = 0;
			for (sglobal = 0; sglobal < S; sglobal++){
				Q += p[sglobal] * solve_SubProblem();
			}
			Out = Open_File(outfile, "a+");
			fprintf(Out, "%d\t%f\t%f\t%f\t%f\t%f\t|\t", iter + 1, obj, Total_theta, Q, MasterTime, SubTime);
			SubTime=0;
			fclose(Out);
			iter++;
			gettimeofday(&stop, NULL);
			checkTime = ((double)(stop.tv_sec - startTotal.tv_sec) * 1000 + (double)(stop.tv_usec - startTotal.tv_usec) / 1000) / 1000;

			//tempgap = fabs(theta - Q);
		}
		//end = clock();
		gettimeofday(&stopTotal, NULL);
		//timeTotal = (double)(end - start) / (double)(CLOCKS_PER_SEC);
		timeTotal = ((double)(stopTotal.tv_sec - startTotal.tv_sec) * 1000 + (double)(stopTotal.tv_usec - startTotal.tv_usec) / 1000) / 1000;

		Out = Open_File(outfile, "a+");
		fprintf(Out, "%f\t", timeTotal);
		fprintf(Out, "\t%lf\t%lf\t", installation, pre);

		fprintf(Out, "\t|\t w: \t|\t");
		for (h = 0; h < A; h++){
				if (w[h]>0.5)
					fprintf(Out, "%d \t", h);
			}

		fprintf(Out, "\t|\t x: \t|\t");
		for (h = 0; h < A; h++){
			for (s = 0; s<S; s++){
				if (x1[h*S + s]>0.5)
					fprintf(Out, "%d\t%d\t|\t", h,s);
			}
		}

		fclose(Out);



		//Print_solution();
		free_memory();

	}
	fclose(ini);
}


void Print_solution(void)
{
	/*fprintf(Out, "%d\n", N);
	fprintf(Out, "%d\n", Nperim);
	fprintf(Out, "%f\n", LB);
	for (int j = 0; j < M; j++){
	fprintf(Out, "%f\n", UB[j]);
	}
	for (int j = 0; j < M; j++){
	fprintf(Out, "%f\n", gap[j]);
	}
	fprintf(Out, "%f\n", timeTotal);
	fprintf(Out, "___________________\n\n");
	fprintf(Out, "%f\n", CLL);
	fprintf(Out, "%f\n", CLU);
	for (int j = 0; j < M; j++){
	fprintf(Out, "%f\n", CUL[j]);
	fprintf(Out, "%f\n", CUU[j]);
	}*/

}
