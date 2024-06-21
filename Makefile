LIB = ../../../../../../encs/pkg/cplex-22.1.0/root/cplex
BBCParetomulti : main.o ReadData.o SPViolated.o solve_MasterProblem_model.o solve_SubProblem_model.o solve_branch_benders_cut.o warm_start.o warm_start2.o warm_start3.o solve_ParetoCut.o
	gcc -O2 -I$(LIB)/include main.o ReadData.o SPViolated.o solve_MasterProblem_model.o solve_SubProblem_model.o solve_branch_benders_cut.o warm_start.o warm_start2.o warm_start3.o solve_ParetoCut.o -L$(LIB)/lib/x86-64_linux/static_pic  -lcplex -lm -lpthread -ldl -o BBCParetomulti

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
warm_start.o: warm_start.c def.h
	gcc -O2 -c warm_start.c
warm_start2.o: warm_start2.c def.h
	gcc -O2 -c warm_start2.c
warm_start3.o: warm_start3.c def.h
	gcc -O2 -c warm_start3.c
solve_ParetoCut.o: solve_ParetoCut.c def.h
	gcc -O2 -c solve_ParetoCut.c
clean :
	rm BBCParetomulti main.o ReadData.o SPViolated.o solve_MasterProblem_model.o solve_SubProblem_model.o warm_start.o warm_start2.o warm_start3.o solve_branch_benders_cut.o solve_ParetoCut.o