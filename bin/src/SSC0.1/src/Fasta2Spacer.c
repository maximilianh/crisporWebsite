/*
 *  fasta2spacer.c
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

#define MAX_SEQ_LEN 10000
#define MAX_SPACER_LEN 100
#define MAX_COLUMN 10
#define MAX_SEQ_NUM 1000000

//Extract spacers from sequence and save to file
int ExtractSpacers2File(char *seq, FILE *fh, int len5Prime, int len3Prime, char strand, char *info);

// compute reverse complement of a sequence
int ReverseSeq(char *src, char *dest, int maxLen);

//compute reverse complement of a nucleotide
char ReverseChar(char c);

//print the usage of Command
void PrintCommandUsage(const char *command);

//Procee a fasta file and output to a sequence file
int ProcessFastaFile(char *input, char *output, int len5Prime, int len3Prime);

int main (int argc, const char * argv[]) 
{
   	 char inputFile[1000], outputFile[1000];
   	 int i, len5Prime, len3Prime;
	
   	 inputFile[0]=0;
   	 outputFile[0]=0;
   	 len5Prime = 20;
   	len3Prime = 10;
	
	for (i=2;i<argc;i++)
	{
		if (strcmp(argv[i-1], "-5")==0)
		{
			len5Prime = atoi(argv[i]);
		}
        if (strcmp(argv[i-1], "-3")==0)
        {
            len3Prime = atoi(argv[i]);
        }
		if (strcmp(argv[i-1], "-i")==0)
		{
			strcpy(inputFile, argv[i]);
		}
		if (strcmp(argv[i-1], "-o")==0)
		{
			strcpy(outputFile, argv[i]);
		}
	}
	
	if ((inputFile[0]==0)||(outputFile[0]==0)||(len5Prime<0)||(len3Prime<0)||(len5Prime+len3Prime>MAX_SPACER_LEN))
	{
		printf("Command error!\n");
		PrintCommandUsage(argv[0]);
		return -1;
	}
	
	SetTimer();
    
    ProcessFastaFile(inputFile,outputFile, len5Prime, len3Prime);

    printf("ElapsedSeconds = %f\n",ElapsedSeconds());
    return 1;
}


//Procee a fasta file and output to a sequence file
int ProcessFastaFile(char *input, char *output, int len5Prime, int len3Prime)
{
	FILE *inputFh, *outputFh;
    char **words, tmpLine[MAX_SEQ_LEN+1], tmpSeq[MAX_SEQ_LEN+1], tmpRevSeq[MAX_SEQ_LEN+1], info[MAX_SEQ_LEN];
    int seqCount = 0, spacerCount = 0, itemNum = 0;
	
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
    
    tmpSeq[0] = 0;
    strcpy(info, "NULL");
	
	while (fgets(tmpLine, MAX_COLUMN*(MAX_SEQ_LEN+1), inputFh)!=NULL)
	{
        itemNum = StringToWords(words, tmpLine, MAX_SEQ_LEN, MAX_COLUMN, " \t\r\n\v\f");
      
	if (itemNum<=0)
	{
		continue;
	} 
	if (words[0][0]=='#')
	{
		continue;
	}
 
        if (words[0][0]!='>')
        {
            strcat(tmpSeq, words[0]);
            continue;
        }
        
        if ((strlen(tmpSeq)>0)&&(ReverseSeq(tmpSeq, tmpRevSeq, MAX_SEQ_LEN)==strlen(tmpSeq)))
        {
            spacerCount += ExtractSpacers2File(tmpSeq, outputFh, len5Prime, len3Prime, '+', info);
            spacerCount += ExtractSpacers2File(tmpRevSeq, outputFh, len5Prime, len3Prime, '-', info);
            tmpSeq[0]=0;
        }
                                 
        strcpy(info, words[0]+1);
        seqCount ++;
	}
    
    if ((strlen(tmpSeq)>0)&&(ReverseSeq(tmpSeq, tmpRevSeq, MAX_SEQ_LEN)==strlen(tmpSeq)))
    {
        spacerCount += ExtractSpacers2File(tmpSeq, outputFh, len5Prime, len3Prime, '+', info);
        spacerCount += ExtractSpacers2File(tmpRevSeq, outputFh, len5Prime, len3Prime, '-', info);
        tmpSeq[0]=0;
    }
	
    FreeWords(words, MAX_COLUMN);
    
	fclose(inputFh);
	fclose(outputFh);
	
    printf("%d sequences processed, %d spacers identified.\n", seqCount, spacerCount);
    
    return spacerCount;
}

//Extract spacers from sequence and save to file
int ExtractSpacers2File(char *seq, FILE *fh, int len5Prime, int len3Prime, char strand, char *info)
{
    char *p, spacer[MAX_SPACER_LEN];
    int seqLen, spacerLen;
    int count=0, i, start, end;
    
    seqLen = (int)strlen(seq);
    spacerLen = len5Prime+len3Prime;
    
    if ((spacerLen<=0)||(seqLen<spacerLen))
    {
        return 0;
    }

    for (i=0;i<seqLen-(len3Prime<3?len5Prime+3:spacerLen)+1;i++)
    {
        p = seq+i;
        
        if (((p[len5Prime+1]=='G')||(p[len5Prime+1]=='g'))
            &&((p[len5Prime+2]=='G')||(p[len5Prime+2]=='g')))
        {
            memcpy(spacer, p, spacerLen*sizeof(char));
            spacer[spacerLen]=0;
            
            if (strand=='+')
            {
                start = i;
                end = i+spacerLen-1;
            }
            else
            {
                start = seqLen-spacerLen-i;
                end = seqLen-i-1;
            }
            
            fprintf(fh, "%s\t%d\t%d\t%c\t%s\n", spacer, start, end, strand, info);
            count ++;
        }
    }
    
    return count;
}

// compute reverse complement of a sequence
int ReverseSeq(char *src, char *dest, int maxLen)
{
    int i, len;
    char c;
    
    len = (int)strlen(src);
    
    for (i=0;i<len;i++)
    {
        c = ReverseChar(src[len-1-i]);
        
        if (c<=0)
        {
            return -1;
        }
        
        dest[i]=c;
    }
    dest[len]=0;
    
    return len;
}

//compute reverse complement of a nucleotide
char ReverseChar(char c)
{
    if (c=='A')
    {
        return 'T';
    }
    if (c=='a')
    {
        return 't';
    }
    if (c=='C')
    {
        return 'G';
    }
    if (c=='c')
    {
        return 'g';
    }
    if (c=='G')
    {
        return 'C';
    }
    if (c=='g')
    {
        return 'c';
    }
    if (c=='T')
    {
        return 'A';
    }
    if (c=='t')
    {
        return 'a';
    }
    if (c=='N')
    {
        return 'N';
    }
    if (c=='n')
    {
        return 'n';
    }
    
    return 0;
}

//print the usage of Command
void PrintCommandUsage(const char *command)
{
	//the options of the command
	printf("%s - extract spacers containing NGG PAM motif from the fasta file \n", command);
	printf("\tmaximum length of a single sequence: 10000\n");
	printf("\tmaximum number of sequences in a fasta file: 1000000\n");
	printf("Usage:\n");
	printf("-5 <length of sequence to the 5' of the PAM>, default: 20\n");
	printf("-3 <length of sequence to the 3' of the PAM, including PAM>, default: 10\n");
	printf("-i <input fasta file>. \n");
	printf("-o <output spacer file>. \n");
}
