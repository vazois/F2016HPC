CC=mpicc
EXEC=sieve
EXEC2=sorg

EXP=6
P=8

all: cc

cc: 
	$(CC) main.c -lm -o $(EXEC)

clean:
	rm -rf $(EXEC)
	
run: 
	mpiexec -np=$(P) ./$(EXEC) $(EXP)
	
