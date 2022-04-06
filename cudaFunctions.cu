#include <cuda_runtime.h>
#include <helper_cuda.h>
#include "header.h"


#define CHECK_ERR(err,msg) (\
		{if (err != cudaSuccess) { \
			fprintf(stderr, msg " - %s\n", cudaGetErrorString(err)); \
			exit(EXIT_FAILURE); \
		} \
	})

//////////////////////////////////////////////////////////////////////////////////////////////////calculate mutant score

__device__ float calcMutantScore(char* seq1, char* seq2, float* weights, float* lettersGrid, int len2, int n, int k)
{
	float score = 0;
	int i = 0, j = 0;
	for (i = 0; i < len2; i++, j++)
	{
		if (j == n || j == k) 
			j++;
		int fc_idx = seq1[i] - 65;
		int sc_idx = seq2[j] - 65;
		score += lettersGrid[fc_idx*ABC_NUM + sc_idx];
	}	

	return score;	
}

//////////////////////////////////////////////////////////////////////////////////////////////// find mutant best score

__global__ void getMutantBestScore(char* d_seq1, char* d_seq2, float* d_weights, 
						float* d_bestScores, int* d_bestOffsets, int* d_nkArr, float* d_lettersGrid,
						int num_mutants, int maxOffset, int len2)
{
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	float bestScore = -10000;
	int offset = 0;
	
	if (i < num_mutants)
	{
		for (int j = 0; j < maxOffset; j++)
		{
			float score = calcMutantScore(d_seq1 + j,d_seq2, d_weights, d_lettersGrid, len2, d_nkArr[i], d_nkArr[i+num_mutants]);
			if (score > bestScore)
			{
				bestScore = score;
				offset = j;
			}
		}
		d_bestScores[i] = bestScore;
		d_bestOffsets[i] = offset;
	}		
} 

void calcBestScoreCUDA(char* seq1, char* seq2, float weights[], float* bestScores, int* bestOffsets, int num_mutants, int len2, int* nkArr, float* lettersGrid)
{
	int len1 = strlen(seq1);
	int maxOffset = len1 - (len2-2) + 1;

////////////////////////////////////////////////////////////////////////////////////////////// allocate seq1 memory

	char* d_seq1 = NULL; 
	
	cudaError_t err = cudaSuccess;
	size_t  arrSize = len1 * sizeof(char);
	err = cudaMalloc((void**)&d_seq1, arrSize);
	CHECK_ERR(err, "Failed to allocate device memory");
	err = cudaMemcpy(d_seq1, seq1, arrSize, cudaMemcpyHostToDevice);
	CHECK_ERR(err, "Failed to copy data from host to device"); 
	
////////////////////////////////////////////////////////////////////////////////////////////// allocate seq2 memory	

	char* d_seq2 = NULL; 
	
	err = cudaSuccess;
	arrSize = len2 * sizeof(char);
	cudaMalloc((void**)&d_seq2, arrSize);
	CHECK_ERR(err, "Failed to allocate device memory");
	cudaMemcpy(d_seq2, seq2, arrSize, cudaMemcpyHostToDevice);
	CHECK_ERR(err, "Failed to copy data from host to device"); 

////////////////////////////////////////////////////////////////////////////////////////////// allocate weights memory

	float* d_weights = NULL; 
	
	err = cudaSuccess;
	arrSize = WEIGHTS_NUM * sizeof(float);
	cudaMalloc((void**)&d_weights, arrSize);
	CHECK_ERR(err, "Failed to allocate device memory");
	cudaMemcpy(d_weights, weights, arrSize, cudaMemcpyHostToDevice);
	CHECK_ERR(err, "Failed to copy data from host to device");
	
////////////////////////////////////////////////////////////////////////////////////////////// allocate weights memory

	int* d_nkArr = NULL; 

	err = cudaSuccess;
	arrSize = num_mutants*2 * sizeof(int);
	cudaMalloc((void**)&d_nkArr, arrSize);
	CHECK_ERR(err, "Failed to allocate device memory");
	cudaMemcpy(d_nkArr, nkArr, arrSize, cudaMemcpyHostToDevice);
	CHECK_ERR(err, "Failed to copy data from host to device"); 

////////////////////////////////////////////////////////////////////////////////////////////// allocate weights memory

	float* d_lettersGrid = NULL; 
	
	err = cudaSuccess;
	arrSize = ABC_NUM*ABC_NUM * sizeof(float);
	cudaMalloc((void**)&d_lettersGrid, arrSize);
	CHECK_ERR(err, "Failed to allocate device memory");
	cudaMemcpy(d_lettersGrid, lettersGrid, arrSize, cudaMemcpyHostToDevice);
	CHECK_ERR(err, "Failed to copy data from host to device"); 

//////////////////////////////////////////////////////////////////////////////////////////// allocate best scores memory

	float* d_bestScores = NULL;  
	
	err = cudaSuccess;
	arrSize = num_mutants * sizeof(float);
	cudaMalloc((void**)&d_bestScores, arrSize);
	CHECK_ERR(err, "Failed to allocate device memory");
	cudaMemcpy(d_bestScores, bestScores, arrSize, cudaMemcpyHostToDevice);
	CHECK_ERR(err, "Failed to copy data from host to device"); 

//////////////////////////////////////////////////////////////////////////////////////////// allocate best offsets memory

	int* d_bestOffsets = NULL;
	
	err = cudaSuccess;
	arrSize = num_mutants * sizeof(int);
	cudaMalloc((void**)&d_bestOffsets, arrSize);
	CHECK_ERR(err, "Failed to allocate device memory");
	cudaMemcpy(d_bestOffsets, bestOffsets, arrSize, cudaMemcpyHostToDevice);
	CHECK_ERR(err, "Failed to copy data from host to device");
	
	int threads = 256;
	int blocks = (num_mutants + threads-1) / threads;
	
	getMutantBestScore<<<blocks, threads>>>(d_seq1, d_seq2, d_weights, d_bestScores, d_bestOffsets, d_nkArr, d_lettersGrid, num_mutants, maxOffset, len2);
	
	cudaMemcpy(bestScores, d_bestScores, num_mutants * sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(bestOffsets, d_bestOffsets, num_mutants * sizeof(int), cudaMemcpyDeviceToHost);
	
	cudaFree(d_seq1);
	cudaFree(d_seq2);
	cudaFree(d_weights);
	cudaFree(d_nkArr);
	cudaFree(d_bestScores);
	cudaFree(d_bestOffsets);
}

// ========================

