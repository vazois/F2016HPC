CC=mpicc
EXEC=sieve
EXEC2=sorg

EXP=10
P=8

all: cc

cc: 
	$(CC) main.c -lm -o $(EXEC)
	$(CC) sieve_original.c -lm -o $(EXEC2)

clean:
	rm -rf $(EXEC)
	rm -rf $(EXEC2)
	
run: 
	mpiexec -np=$(P) ./$(EXEC) $(EXP)
	#mpiexec -np=$(P) ./$(EXEC2) $(EXP)
	
