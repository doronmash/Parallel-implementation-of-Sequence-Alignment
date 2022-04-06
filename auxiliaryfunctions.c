#include "header.h"

const char* conservatives[CONSERVATIVES_LEN] = {
	"NDEQ", "NEQK", "STA", "MILV", "QHRK", "NHQK", "FYW", "HY", "MILF"
};

const char* semi_conservatives[SEMI_CONSERVATIVES_LEN] = {
	"SAG", "ATV", "CSA", "SGND", "STPA", "STNK", "NEQHRK", "NDEQHK", "SNDEQK", "HFY", "FVLIM"
};




///////////////////////////////////////////////////////////////////////////////////////// populate bundle data for MPI

void populateBundleData(char* seq1, char* seq2s[], int seq2Numbers, float weights[], Bundle* bundle)
{
	int len = strlen(seq1);
	for (int i = 0; i < len; i++)
		bundle->seq1[i] = seq1[i];

	bundle->seq2Numbers = seq2Numbers;
	for (int i = 0; i < seq2Numbers; i++)
	{
		len = strlen(seq2s[i]);
		for (int j = 0; j < len; j++)
			bundle->seq2s[i][j] = seq2s[i][j];
	}

	for (int i = 0; i < WEIGHTS_NUM; i++)
		bundle->weights[i] = weights[i];
}

/////////////////////////////////////////////////////////////////////////////////////////assign weight for chars 

float compareChars(char first, char second, float weights[])
{
	if (first == second)
		return weights[0];
	else if (CheckIf2CharInSameGroup(first, second, conservatives, CONSERVATIVES_LEN))
		return -weights[1];
	else if (CheckIf2CharInSameGroup(first, second, semi_conservatives, SEMI_CONSERVATIVES_LEN))
		return -weights[2];

	return -weights[3];
}


int* createNK(int num_mutants)
{
	int* nkArr = (int*) malloc(num_mutants * 2 * sizeof(int));
	int n = 0, k = 1;
	for (int i = 0; i < num_mutants; i++) 
	{
		nkArr[i] = n;
		nkArr[i+num_mutants] = k;
		n++;
		if (n == k)
		{
			n = 0;
			k++;
		} 
	}
	return nkArr;
}


/////////////////////////////////////////////////////////////////////////////////////////checks if both chars in the group

int CheckIf2CharInSameGroup(char first, char second, const char** group, int len)
{
	for (int i = 0; i < len; i++)
	{
		if (strchr(group[i], first) != NULL && strchr(group[i], second) != NULL)
			return 1;
	}
	return 0;
}


///////////////////////////////////////////////////////////////////////////////////////creates letters grid for char comparison

float* createLettersGrid(float weights[])
{
	float* lettersGrid = (float*)malloc(ABC_NUM * ABC_NUM * sizeof(float));
	const char* letters = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
	
	for (int i = 0; i < ABC_NUM; i++)
	{
		for (int j = 0; j < ABC_NUM; j++)
			lettersGrid[i*ABC_NUM + j] = compareChars(letters[i],letters[j],weights);
	}
	
	return lettersGrid;
}





