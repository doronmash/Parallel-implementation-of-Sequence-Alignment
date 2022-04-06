#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "header.h"

/////////////////////////////////////////////////////////////////////////////read data from file

char** readFromFile(const char* fileName, float weights[], char* seq1, int* numOfSeq2)
{
	FILE* fp;
	char* line = NULL;
	
	if ((fp = fopen(fileName, "r")) == 0) 
	{
		printf("cannot open file  %sfor reading\n", fileName);
		exit(0);
	}
	
	//read weights
	for (int i = 0; i < WEIGHTS_NUM; i++) 
		fscanf(fp, "%f", &weights[i]);	
		
	//read seq1
	fscanf(fp, "%s", seq1); 
	//read number of seq2
	fscanf(fp, "%d",numOfSeq2); 
	
	//allocate memory for seq1
	char** seq2s = (char**)malloc(sizeof(char*) * *numOfSeq2); 
	
	//allocate memory for all seq2
	for (int i = 0; i < *numOfSeq2; i++) 
		seq2s[i] = (char*)malloc(sizeof(char) * SEQ2_LEN);
	//read seq2
	for (int i = 0; i < *numOfSeq2; i++) 
		fscanf(fp, "%s", seq2s[i]);

	fclose(fp);
	
	return seq2s;
}

///////////////////////////////////////////////////////////////////////////write results to file

void writeToFile(const char* fileName, Result* results, int num_seqs)
{
	FILE* fp;
	
	
	if ((fp = fopen(fileName, "w")) == 0) 
	{
		printf("cannot open file  %sfor writing\n", fileName);
		exit(0);
	}
	
	//write best result of seq2 to output file
	for (int i = 0; i < num_seqs; i++) 
		fprintf(fp, "seq %d\toffset - %d\t score - %.2f \tMS(%d, %d)\n", i, results[i].offset, results[i].score, results[i].n, results[i].k);
	
	fclose(fp);
}

// ========================
