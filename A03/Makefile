CC=mpicc
EXEC=sieve
EXEC2=sorg

EXP=2
P=1

all: cc

cc: 
	$(CC) main.c -lm -o $(EXEC)

ccO3:
	$(CC) -O3 main.c -lm -o $(EXEC)

clean:
	rm -rf $(EXEC)
	
run: 
	mpiexec -np=$(P) ./$(EXEC) $(EXP)
	
