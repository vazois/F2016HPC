rm -rf sieve
mpicc main.c -lm -o sieve
Executing original sieve
There are {455052511} primes less than or equal to 10000000000
Elapsed time of <Original Sieve> for (256) processes   7.286523
Executing odd numbers sieve
There are {455052511} primes less than or equal to 10000000000
Elapsed time of <Sieve on odd numbers> for (256) processes   3.550466
Executing local odd numbers sieve
There are {455052511} primes less than or equal to 10000000000
Elapsed time of <Sieve with local sieve> for (256) processes   3.384981
Executing odd numbers sieve with cache awareness 
Cached Array Size:10000
There are {455052623} primes less than or equal to 10000000000
Elapsed time of <Sieve with cache awareness> for (256) processes   2.813120
Comparison between original to odd version,(PASS)
Comparison between original to local odd version,(PASS)
Comparison between original to local cache aware version,(FAIL)
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


Executing original sieve
There are {455052511} primes less than or equal to 10000000000
Elapsed time of <Original Sieve> for (256) processes   7.311476
Executing odd numbers sieve
There are {455052511} primes less than or equal to 10000000000
Elapsed time of <Sieve on odd numbers> for (256) processes   3.559590
Executing local odd numbers sieve
There are {455052511} primes less than or equal to 10000000000
Elapsed time of <Sieve with local sieve> for (256) processes   3.392860
Executing odd numbers sieve with cache awareness 
Cached Array Size:100000
There are {455052511} primes less than or equal to 10000000000
Elapsed time of <Sieve with cache awareness> for (256) processes   1.026791
Comparison between original to odd version,(PASS)
Comparison between original to local odd version,(PASS)
Comparison between original to local cache aware version,(PASS)
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


Executing original sieve
There are {455052511} primes less than or equal to 10000000000
Elapsed time of <Original Sieve> for (256) processes   7.280681
Executing odd numbers sieve
There are {455052511} primes less than or equal to 10000000000
Elapsed time of <Sieve on odd numbers> for (256) processes   3.558429
Executing local odd numbers sieve
There are {455052511} primes less than or equal to 10000000000
Elapsed time of <Sieve with local sieve> for (256) processes   3.396591
Executing odd numbers sieve with cache awareness 
Cached Array Size:1000000
There are {455052511} primes less than or equal to 10000000000
Elapsed time of <Sieve with cache awareness> for (256) processes   0.871949
Comparison between original to odd version,(PASS)
Comparison between original to local odd version,(PASS)
Comparison between original to local cache aware version,(PASS)
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


Executing original sieve
There are {455052511} primes less than or equal to 10000000000
Elapsed time of <Original Sieve> for (256) processes   7.287101
Executing odd numbers sieve
There are {455052511} primes less than or equal to 10000000000
Elapsed time of <Sieve on odd numbers> for (256) processes   3.566859
Executing local odd numbers sieve
There are {455052511} primes less than or equal to 10000000000
Elapsed time of <Sieve with local sieve> for (256) processes   3.392363
Executing odd numbers sieve with cache awareness 
Cached Array Size:10000000
There are {455052511} primes less than or equal to 10000000000
Elapsed time of <Sieve with cache awareness> for (256) processes   2.793496
Comparison between original to odd version,(PASS)
Comparison between original to local odd version,(PASS)
Comparison between original to local cache aware version,(PASS)
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

