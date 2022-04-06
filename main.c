#include "mpi.h"
#include "omp.h"
#include "header.h"
#include <stddef.h>
#define RESULT_ATTRS_NUM 5
#define BUNDLE_ATTRS_NUM 4

int main(int argc, char* argv[])
{
	int rank, num_procs;
	Bundle bundle;	
	Result* results;
	
	MPI_Datatype MPI_RESULT;
	MPI_Datatype MPI_BUNDLE;
	MPI_Status status;
	
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank); 
	MPI_Comm_size(MPI_COMM_WORLD, &num_procs); 

	
///////////////////////////////////////////////////////////////////////////////////////////////////////////// create result type
	
	int blockSize1[RESULT_ATTRS_NUM] = { 1, 1, 1, 1, 1 };
	MPI_Aint disp1[RESULT_ATTRS_NUM];
	MPI_Datatype types1[RESULT_ATTRS_NUM] = { MPI_FLOAT, MPI_INT, MPI_INT, MPI_INT, MPI_INT };

	disp1[0] = offsetof(Result, score);
	disp1[1] = offsetof(Result, offset);
	disp1[2] = offsetof(Result, n);
	disp1[3] = offsetof(Result, k);
	disp1[4] = offsetof(Result, mutant_num);
	
	MPI_Type_create_struct(RESULT_ATTRS_NUM, blockSize1, disp1, types1, &MPI_RESULT);
	MPI_Type_commit(&MPI_RESULT);
	
	
///////////////////////////////////////////////////////////////////////////////////////////////////////////// create bundle type
	
	int blockSize2[BUNDLE_ATTRS_NUM] = { SEQ1_LEN, 1, MAX_SEQS2 * SEQ2_LEN, WEIGHTS_NUM };
	MPI_Aint disp2[BUNDLE_ATTRS_NUM];
	MPI_Datatype types2[BUNDLE_ATTRS_NUM] = { MPI_CHAR, MPI_INT, MPI_CHAR, MPI_FLOAT };

	disp2[0] = offsetof(Bundle, seq1);
	disp2[1] = offsetof(Bundle, seq2Numbers);
	disp2[2] = offsetof(Bundle, seq2s);
	disp2[3] = offsetof(Bundle, weights);
	
	MPI_Type_create_struct(BUNDLE_ATTRS_NUM, blockSize2, disp2, types2, &MPI_BUNDLE);
	MPI_Type_commit(&MPI_BUNDLE);
	
	if (num_procs != NUM_PROCS)
	{
		printf("please use this program with %d processes\n", NUM_PROCS);
		MPI_Abort(MPI_COMM_WORLD, __LINE__);
	}
	int start_time = omp_get_wtime();
/////////////////////////////////////////////////////////////////////////////////////////////////////////////master work
	if (rank == MASTER) 
	{
		int numSeq2s = 0;
		char seq1[SEQ1_LEN];
		float weights[WEIGHTS_NUM];
		
		//getint from file input data
		char** seq2s = readFromFile(READ_FILE_NAME, weights, seq1, &numSeq2s); 	

		//allocate memory for results
		results = (Result*) malloc(numSeq2s * sizeof(Result)); 

		populateBundleData(seq1, seq2s + numSeq2s/NUM_PROCS, numSeq2s/NUM_PROCS, weights, &bundle);

		MPI_Send(&bundle, 1, MPI_BUNDLE, SLAVE, 0, MPI_COMM_WORLD); 	
		
		#pragma omp parallel for
		for (int i = 0; i < numSeq2s/NUM_PROCS; i++) // calculate best results
			results[i] = findBestMutant(seq1, seq2s[i], weights);

		//recieve results from slave
		for (int i = numSeq2s/NUM_PROCS; i < numSeq2s; i++) 
			MPI_Recv(&results[i], 1, MPI_RESULT, SLAVE, 0, MPI_COMM_WORLD, &status);
			
		//write results to output file
		writeToFile(WRITE_FILE_NAME, results, numSeq2s); 	
	}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////slave work
	else 
	{	
		//get data from master
		MPI_Recv(&bundle, 1, MPI_BUNDLE, MASTER, 0, MPI_COMM_WORLD, &status); 
		
		//allocate memory for all results
		results = (Result*) malloc(bundle.seq2Numbers * sizeof(Result)); 

		#pragma omp parallel for
		for (int i = 0; i < bundle.seq2Numbers; i++) // get best results
			results[i] = findBestMutant(bundle.seq1, bundle.seq2s[i], bundle.weights);
		
		for (int i = 0; i < bundle.seq2Numbers; i++) // send to master the best results
			MPI_Send(&results[i], 1, MPI_RESULT, MASTER, 0, MPI_COMM_WORLD);
	}
	
	printf("process %d finished work after - %.2f seconds\n", rank, omp_get_wtime() - start_time);
		
	MPI_Finalize();
	return 0;
}




