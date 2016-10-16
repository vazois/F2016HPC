CC=g++ 
FLAGS= -O3

EXEC=dgemm


all: c99

cc: 
	$(CC) -std=c++11 main.cpp -o $(EXEC)
	
c99:
	$(CC) main.cpp -o $(EXEC)

c11o3:
	$(CC) -O3 -std=c++11 main.cpp -o $(EXEC)

c99o0:
	$(CC) -O0 main.cpp -o $(EXEC)
	
c99o1:
	$(CC) -O1 main.cpp -o $(EXEC)

c99o2:
	$(CC) -O2 main.cpp -o $(EXEC)

c99o3:
	$(CC) -O3 main.cpp -o $(EXEC)
	

clean:
	rm -rf $(EXEC)

