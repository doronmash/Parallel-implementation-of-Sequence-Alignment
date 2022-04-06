                                           Parallel implementation of Sequence Alignment

MPI+OpenMP+CUDA Integration 
At first we create two processes. 
The input file initially is known for one machine only, this machine is related to the first process (with rank 0). 
The same machine at the end of the run writes the results to the output file. 
The first process reads the data from the input file and then sends the weights, the first sequence (Seq1), and half of the query sequences (Seq2) to the second process (with rank 1).
The communication between the two processes is done using MPI. 
Both processes perform the exact same computations on their half of the sequences - partially with OpenMP, partially with CUDA. 
The results are sent from the second process to the first process, which after writes the results to the output file.

Sequence Alignment Algorithm Implementation
My implementation for finding the offset and n, k with the best alignment score for a given sequence is composed of the following steps:
1.	Calculate the alphabet matrix according to our weights (compare in char in alphabet to all other).
2.	Calculate all n, k according to the number of mutants.
3.	For each mutant, we iterate through all offsets and calculate its score according to the formula (CUDA).
4.	We save the best score and offset for each mutant.
5.	After finding all mutants scores, we iterate through all mutants and find the best score (OMP).

