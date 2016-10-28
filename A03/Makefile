CC=mpicc
EXEC=sieve
EXEC2=sorg

EXP=10
P=32
C=1000000

all: cc

cc: 
	$(CC) -O3 main.c -lm -o $(EXEC)

ccO2:
	$(CC) -O2 main.c -lm -o $(EXEC)
	
ccO3:
	$(CC) -O3 main.c -lm -o $(EXEC)

clean:
	rm -rf $(EXEC)
	
run: 
	mpiexec -np=$(P) ./$(EXEC) $(EXP) $(C)
	
