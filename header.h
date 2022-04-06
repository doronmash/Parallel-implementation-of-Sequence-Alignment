#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <string.h>


#define READ_FILE_NAME "input.txt"
#define WRITE_FILE_NAME "output.txt"
#define MASTER 0
#define SLAVE 1
#define NUM_PROCS 2
#define SEQ1_LEN 5000
#define SEQ2_LEN 3000
#define CONSERVATIVES_LEN 9
#define SEMI_CONSERVATIVES_LEN 11
#define WEIGHTS_NUM 4
#define MAX_SEQS2 30
#define ABC_NUM 26

typedef struct
{
	float score;
	int offset, n, k;
	int mutant_num;
} Result;

typedef struct
{
	char seq1[SEQ1_LEN];
	int seq2Numbers;
	char seq2s[MAX_SEQS2][SEQ2_LEN];
	float weights[WEIGHTS_NUM];
	
} Bundle;




/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////auxiliaryfunctions

void populateBundleData(char* seq1, char* seqs2[], int numSeqs2, float weights[], Bundle* bundle);
int* createNK(int num_mutants);
float* createLettersGrid(float weights[]);
float compareChars(char first, char second, float weights[]);
int CheckIf2CharInSameGroup(char first, char second, const char** group, int len);

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////cudaFunctions

void calcBestScoreCUDA(char* seq1, char* seq2, float weights[], float* mutantsBestScores, int* mutantsBestOffsets, int num_mutants, int len2, int* nkArr, float* ABCgrid);

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////fileUtil

char** readFromFile(const char* fileName, float weights[], char* seq1, int* numOfSeq2);
void writeToFile(const char* fileName, Result* results, int num_seqs);

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////cFunctions

Result findBestMutantAndOffsetOMP(float* mutantsBestScores, int* mutantsBestOffsets, int num_mutants, int seq2_len, int* nkArr);
Result findBestMutant(char* seq1, char* seq2, float weights[]);




