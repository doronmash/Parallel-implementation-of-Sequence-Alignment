#include "header.h"
#include "omp.h"
#include <float.h>

///////////////////////////////////////////////////////////////////////////////////////////////// find best mutant from all seq2 mutants

Result findBestMutantAndOffsetOMP(float* mutantsBestScores, int* mutantsBestOffsets, int num_mutants, int seq2_len, int* nkArr)
{
	int bestMutant = -1;
	Result result;
	result.score = -FLT_MAX;
	#pragma omp parallel for
	for(int i = 0; i < num_mutants; i++) // find best score of all mutants
	{
		if (mutantsBestScores[i] > result.score)
		{
			result.score = mutantsBestScores[i]; 
			result.offset = mutantsBestOffsets[i];
			result.n = nkArr[i];
			result.k = nkArr[i+num_mutants];
		}
	}
	return result;
}


///////////////////////////////////////////////////////////////////////////////////////////////// main work function

Result findBestMutant(char* seq1, char* seq2, float weights[])
{
	int len2 = strlen(seq2);
	int num_mutants = len2 * (len2 - 1) / 2;
	int* nkArr = createNK(num_mutants);
	float* lettersGrid = createLettersGrid(weights);
	float* mutantsBestScores = (float*) malloc(num_mutants * sizeof(float));	
	int* mutantsBestOffsets = (int*) malloc(num_mutants * sizeof(int));	
	calcBestScoreCUDA(seq1, seq2, weights, mutantsBestScores, mutantsBestOffsets, num_mutants, len2, nkArr, lettersGrid);
	Result result = findBestMutantAndOffsetOMP(mutantsBestScores, mutantsBestOffsets, num_mutants, len2, nkArr);
	
	return result;
}






