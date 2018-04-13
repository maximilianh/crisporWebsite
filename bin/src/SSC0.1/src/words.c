/*
 *  words.c
 *	Word manipulations
 *
 *  Created by Han Xu on 10/12/12.
 *  Copyright 2012 Dana Farber Cancer Institute. All rights reserved.
 *
 */

#include <string.h>
#include <dirent.h>
#include <assert.h>
#include "words.h"
#include "stdlib.h"

//allocate 2d array of characters. Return the pointer to the array, and NULL if failure
char **AllocWords(int wordNum, int wordLen)
{
	char **ptr;
	int i;
	
	ptr = (char **)malloc(wordNum*sizeof(char *));
	
	if (!ptr)
	{
		return NULL;
	}
	
	for (i=0;i<wordNum;i++)
	{
		ptr[i] = malloc(wordLen*sizeof(char));
	}
	
	return ptr;
}

//Free 2d array of characters
void FreeWords(char **ptr, int wordNum)
{
	int i;
	
	for (i=0;i<wordNum;i++)
	{
		free(ptr[i]);
	}
	
	free(ptr);
}

//Extract word from a string. Words seperated by the deliminators. return number of words extracted. return -1 if failure.
//Example: wordNum = StringToWords(words, string, maxWordLen, maxWordNum, " \t\r\n\v\f");
int StringToWords(char **words, char *str, int maxWordLen, int maxWordNum, const char *delim)
{
	char *pch, *tmpStr;
	int wordNum = 0;
	int strlength = strlen(str);
	
	if (strlength<=0)
	{
		return -1;
	}
	
	tmpStr = (char *)malloc((strlength+1)*sizeof(char));
	strcpy(tmpStr, str);
	
	pch = strtok(tmpStr, delim);
	
	while (pch!=NULL) 
	{
		if (strlen(pch)>=maxWordLen)
		{
			free(tmpStr);
			return -1;
		}
		
		strcpy(words[wordNum], pch);
		wordNum++;
		
		if (wordNum >=maxWordNum)
		{
			break;
		}
		pch = strtok(NULL, delim);
	}
	
	free(tmpStr);
	
	return wordNum;
}

//Read a list of file names from a directory to word structure. Return the number of files read. Return -1 if failure.
//Read files end with ext. If ext=NULL, read all files.
int DirToWords(char **words, char *dirName, int maxWordLen, int maxWordNum, const char *ext)
{
	DIR *dp;
	struct dirent *ep;
	int count = 0;
	
	dp = opendir(dirName);
	
	if (dp==NULL)
	{
		return -1;
	}
	
	while ((ep = readdir(dp))!=NULL)
	{
		if ((ext != NULL)&&(strcmp(ext, ep->d_name+strlen(ep->d_name)-strlen(ext))))
		{
			continue;
		}
			
		assert(strlen(ep->d_name)<=maxWordLen);
		assert(count<maxWordNum);
		
		if (strlen(ep->d_name)>maxWordLen)
		{
			continue;
		}
		
		if (count>=maxWordNum)
		{
			break;
		}
		
		strcpy(words[count], ep->d_name);
		count++;
	}
	
	closedir(dp);
	
	return count;
}
