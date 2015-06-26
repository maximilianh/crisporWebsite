/*
 *  performance.h
 *  measurement of computing cost in an application
 *
 *  Created by Han Xu on 10/17/12.
 *  Copyright 2012 Dana Farber Cancer Institute. All rights reserved.
 *
 */

#include <time.h>
#include <sys/time.h>

//Set the timer
void SetTimer();

//Compute the time elapsed from timer, in seconds
double ElapsedSeconds();