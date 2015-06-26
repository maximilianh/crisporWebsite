/*
 *  words.h
 *	Word manipulations
 *
 *  Created by Han Xu on 10/12/12.
 *  Copyright 2012 Dana Farber Cancer Institute. All rights reserved.
 *
 */

//allocate 2d array of characters. Return the pointer to the array, and NULL if failure
char **AllocWords(int wordNum, int wordLen);

//Free 2d array of characters
void FreeWords(char **ptr, int wordNum);

//Extract word from a string. Words seperated by the deliminators. return number of words extracted. return -1 if failure.
//Example: wordNum = StringToWords(words, string, maxWordLen, maxWordNum, " \t\r\n\v\f");
int StringToWords(char **words, char *str, int maxWordLen, int maxWordNum, const char *delim);

//Read a list of file names from a directory to word structure. Return the number of files read. Return -1 if failure.
//Read files end with ext. If ext=NULL, read all files.
int DirToWords(char **words, char *dirName, int maxWordLen, int maxWordNum, const char *ext);