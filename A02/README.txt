#############################
###COMPILATION INSTRUCTION###
#############################
To compile the code just type make. If you want to use optimizations you need to compile by typing 
make c99oX where X is the optimization level that you want. The default compilation mode calls make 
c99 which does not have any optimization flags enabled.

############################
###EXECUTION INSTRUCTIONS###
############################
To execute the program type ./dgemm to get the execution instructions.
For the part where we use only registers you need to provide only the size of the matrices (i.e. ./dgemm -md=0 -n=2048).
For the blocking solution, we need to provide an additional argument indicating the block size (i.e. ./dgemm -md=1 -n=2048 -b=64).
The hybrid solution requires similar arguments, with md only changing which specifies what code to run (i.e. ./dgemm -md=2 -n=2048 -b=64).


