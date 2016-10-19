CC=mpicc
EXEC=sieve


all: cc

cc: 
	$(CC) main.c -lm -o $(EXEC)

clean:
	rm -rf $(EXEC)
	
run: 
	mpiexec -np=128 ./$(EXEC) 8


