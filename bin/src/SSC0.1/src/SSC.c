/*
 *  ssc.c
 *  Spacer Scoring for CRISPR
 *
 *  Created by Han Xu on 10/02/14.
 *  Copyright 2014 Dana Farber Cancer Institute. All rights reserved.
 *
 */

#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "performance.h"
#include "words.h"

#define MAX_SEQ_LEN 100
#define MAX_COLUMN  20
#define NUCLEOTIDE_NUM 4

//Read matrix file
int ReadMatrix(char *fileName, double *matrix, int col, double *intercept);

//Process the spacer sequence file and output
int ProcessSeqFile(char *input, char *output, double *matrix, double intercept, int col, int row, int mode);

//Compute score for a spacer sequence
double ComputeSeqScore(char *src, double *matrix, double intercept, int col, int row, int mode);

//Convert sequence character to index
int SeqToIndex(char c);

//print the usage of Command
void PrintCommandUsage(const char *command);

double matrix[NUCLEOTIDE_NUM*MAX_SEQ_LEN];

int main (int argc, const char * argv[]) 
{
	int seqLen=-1;
	char inputFile[1000], outputFile[1000], matrixFile[1000];
	int i, seqNum = 0;
	int mode;
	double intercept;
	
	inputFile[0]=0;
	outputFile[0]=0;
	matrixFile[0]=0;
    mode = 0;
	
	for (i=2;i<argc;i++)
	{
		if (strcmp(argv[i-1], "-l")==0)
		{
			seqLen = atoi(argv[i]);
		}
		if (strcmp(argv[i-1], "-m")==0)
		{
			strcpy(matrixFile, argv[i]);
		}
		if (strcmp(argv[i-1], "-i")==0)
		{
			strcpy(inputFile, argv[i]);
		}
		if (strcmp(argv[i-1], "-o")==0)
		{
			strcpy(outputFile, argv[i]);
		}
		if (strcmp(argv[i-1], "-b")==0)
		{
			mode = atoi(argv[i]);
		}
	}
	
	if ((inputFile[0]==0)||(outputFile[0]==0)||(matrixFile[0]==0)||(seqLen<0))
	{
		printf("Command error!\n");
		PrintCommandUsage(argv[0]);
		return -1;
	}
	
	if (seqLen>MAX_SEQ_LEN)
	{
		printf("length of sequence should be less than %d.\n", MAX_SEQ_LEN);
		return -1;
	}
	
	SetTimer();
	
	if (ReadMatrix(matrixFile, matrix, NUCLEOTIDE_NUM,&intercept)!=seqLen)
	{
		printf("Error: matrix does not match sequence length. \n");
		return -1;
	}
	
	seqNum = ProcessSeqFile(inputFile,outputFile,matrix, intercept, NUCLEOTIDE_NUM,seqLen, mode);
	
	if (seqNum>0)
	{
		printf("%d sequence processed.\nElapsedSeconds = %f\n",seqNum, ElapsedSeconds());
		return 1;
	}
	else 
	{
		printf("Processing failed.\n");
		return -1;	
	}	
}

//Read matrix file
int ReadMatrix(char *fileName, double *matrix, int col, double *intercept)
{
	FILE *fh;
	int rowNum = 0;
	char tmpS[100];
	int i;
	int flag = 0;
	
	fh = (FILE *)fopen(fileName, "r");
	
	if (!fh)
	{
		printf("Cannot open %s\n", fileName);
		return -1;
	}
	
	fscanf(fh, "%s", tmpS);
	fscanf(fh, "%s", tmpS);
	*intercept = atof(tmpS);

	for (i=0;i<col;i++)
	{
		fscanf(fh, "%s", tmpS);
	}
	
	while (!feof(fh))
	{
		flag = 1;
		
		for (i=0;i<col;i++)
		{
			if (feof(fh))
			{
				flag = 0;
				break;
			}
			
			fscanf(fh,"%s",tmpS);
			
			matrix[rowNum*col+i] = atof(tmpS);
		}
		
		if (!flag)
		{
			break;
		}
		
		rowNum++;
	}
	
	fclose(fh);
	
	return rowNum;
}

//Process the spacer sequence file and output
int ProcessSeqFile(char *input, char *output, double *matrix, double intercept, int col, int row, int mode)
{
	FILE *inputFh, *outputFh;
    char **words, tmpLine[MAX_COLUMN*(MAX_SEQ_LEN+1)];
	int count = 0, itemNum, i;
	
	inputFh = (FILE *)fopen(input, "r");
	
	if (!inputFh)
	{
		printf("Cannot open %s\n", input);
		return -1;
	}
	
	outputFh = (FILE *)fopen(output, "w");
	
	if (!outputFh)
	{
		printf("Cannot open %s\n", output);
		return -1;
	}
	
    words = AllocWords(MAX_COLUMN, MAX_SEQ_LEN);
	
	while (fgets(tmpLine, MAX_COLUMN*(MAX_SEQ_LEN+1), inputFh)!=NULL)
	{
        itemNum = StringToWords(words, tmpLine, MAX_SEQ_LEN, MAX_COLUMN, " \t\r\n\v\f");
		if ((itemNum>0)&&(strlen(words[0])==row))
		{
            for (i=0;i<itemNum;i++)
            {
                fprintf(outputFh,"%s\t",words[i]);
            }
			fprintf(outputFh,"%f\n", ComputeSeqScore(words[0], matrix, intercept, col, row, mode));
			count++;
		}
	}
	
    FreeWords(words, MAX_COLUMN);
    
	fclose(inputFh);
	fclose(outputFh);
	
	return count;
}

//Compute score for a spacer sequence
double ComputeSeqScore(char *src, double *matrix, double intercept, int col, int row, int mode)
{
	int len = strlen(src);
	int i,index;
	double score = intercept;
	
	for (i=0;i<len;i++)
	{
		index = SeqToIndex(src[i]);
		
		if ((index>=0)&&(index<col))
		{
			score += matrix[i*col+index];
		}
	}
	
	if (mode)
	{
		return 1/(1.0+exp(-score));
	}
	else
	{
		return score;
	}
}

//Convert sequence character to index
int SeqToIndex(char c)
{
	if ((c=='A')||(c=='a'))
	{
		return 0;
	}
	if ((c=='C')||(c=='c'))
	{
		return 1;
	}
	if ((c=='G')||(c=='g'))
	{
		return 2;
	}
	if ((c=='T')||(c=='t'))
	{
		return 3;
	}
	return -1;
}

//print the usage of Command
void PrintCommandUsage(const char *command)
{
	//prdouble the options of the command
	printf("%s - Spacer Scoring for CRISPR. This program scans spacer sequence of CRISPR to predict effectiveness of guide RNA. \n", command);
	printf("Usage:\n");
	printf("-l <length of input sequence>. Maximum sequence length is %d.\n", MAX_SEQ_LEN);
	printf("-m <matrix file>. The number of rows in the matrix should match the length of sequence. \n");
	printf("-b <mode>. 0:linear; 1:logistic. Default:0\n");
	printf("-i <input sequence file>. One sequence per line. \n");
	printf("-o <output scoring file>. \n");
}
