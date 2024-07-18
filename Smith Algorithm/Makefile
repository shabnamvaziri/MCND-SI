LIB = ../../../../../../encs/pkg/cplex-12.10.0/root/cplex
Penaltymulti : main.o ReadData.o solve_MasterProblem_model.o solve_SubProblem_model.o 
	gcc -O2 -I$(LIB)/include main.o ReadData.o solve_MasterProblem_model.o solve_SubProblem_model.o -L$(LIB)/lib/x86-64_linux/static_pic  -lcplex -lm -lpthread -ldl -o Penaltymulti

main.o: main.c def.h
	gcc -O2 -c main.c
ReadData.o: ReadData.c def.h
	gcc -O2 -c ReadData.c
solve_MasterProblem_model.o: solve_MasterProblem_model.c def.h
	gcc -O2 -c solve_MasterProblem_model.c
solve_SubProblem_model.o: solve_SubProblem_model.c def.h
	gcc -O2 -c solve_SubProblem_model.c
clean :
	rm Penaltymulti main.o ReadData.o solve_MasterProblem_model.o solve_SubProblem_model.o 
