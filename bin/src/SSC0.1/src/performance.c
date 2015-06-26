/*
 *  performance.c
 *  measurement of computing cost in an application
 *
 *  Created by Han Xu on 10/17/12.
 *  Copyright 2012 Dana Farber Cancer Institute. All rights reserved.
 *
 */

#include "performance.h"

struct timeval t0;

//Set the timer
void SetTimer()
{
	gettimeofday(&t0, NULL);
}

//Compute the time elapsed from timer, in seconds
double ElapsedSeconds()
{
	struct timeval t;
	double second1, second2;
	
	gettimeofday(&t, NULL);
	
	second1 = (double)t0.tv_sec+(double)t0.tv_usec/1000000;
	second2 = (double)t.tv_sec+(double)t.tv_usec/1000000;
	
	return second2-second1;
}