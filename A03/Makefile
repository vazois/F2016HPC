CC=mpicc
EXEC=sieve
EXEC2=sorg

EXP=10
P=256
C=10000

all: cc

cc: 
	$(CC) main.c -lm -o $(EXEC)

ccO2:
	$(CC) -O2 main.c -lm -o $(EXEC)
	
ccO3:
	$(CC) -O3 main.c -lm -o $(EXEC)

clean:
	rm -rf $(EXEC)
	
run: 
	mpiexec -np=$(P) ./$(EXEC) $(EXP) $(C)
	
