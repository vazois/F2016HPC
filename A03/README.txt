+++++++++++++++++++++++++
+Compilation Instruction+
+++++++++++++++++++++++++
Just type make

++++++++++++++++++++++++
+Execution Instructions+
++++++++++++++++++++++++

"make run" gives you an indication how to run
the program

mpiexec ./sieve 10 1000000

the first value is the exponent for the power of
10 the user wishes to count to

the second value is the block size used only for
the cache aware version. Any given value will be 
divided by 2, so k/2 elements are considered each
time.
 
