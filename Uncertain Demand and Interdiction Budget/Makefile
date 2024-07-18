LIB = ../../../../../../encs/pkg/cplex-22.1.0/root/cplex
BBCPenaltymulti : main.o ReadData.o SPViolated.o solve_MasterProblem_model.o solve_SubProblem_model.o solve_branch_benders_cut.o 
	gcc -O2 -I$(LIB)/include main.o ReadData.o SPViolated.o solve_MasterProblem_model.o solve_SubProblem_model.o solve_branch_benders_cut.o -L$(LIB)/lib/x86-64_linux/static_pic  -lcplex -lm -lpthread -ldl -o BBCPenaltymulti

main.o: main.c def.h
	gcc -O2 -c main.c
ReadData.o: ReadData.c def.h
	gcc -O2 -c ReadData.c
SPViolated.o: SPViolated.c def.h
	gcc -O2 -c SPViolated.c
solve_MasterProblem_model.o: solve_MasterProblem_model.c def.h
	gcc -O2 -c solve_MasterProblem_model.c
solve_SubProblem_model.o: solve_SubProblem_model.c def.h
	gcc -O2 -c solve_SubProblem_model.c
solve_branch_benders_cut.o: solve_branch_benders_cut.c def.h
	gcc -O2 -c solve_branch_benders_cut.c
clean :
	rm BBCPenaltymulti main.o ReadData.o SPViolated.o solve_MasterProblem_model.o solve_SubProblem_model.o solve_branch_benders_cut.o 