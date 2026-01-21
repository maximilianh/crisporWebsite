#ifndef DRAW2D_H
#define DRAW2D_H

#include "geometry2D.h"

namespace draw2D {

//coord is a structure which keeps track of each nucleotide's coordinate
struct coordinates{
	coordinates(int size);
    ~coordinates();
	int *x,*y;
	int **num; //[maxbases/10+1][2];
    short int bases;
	void set(int index, int x, int y);
	void set(int index, double x, double y);
};

void place(int number,structure *ct,coordinates *out,int height,int width);
void place(int structnum,structure *ct,coordinates *coord, int diameter, float bondLength, float helixSpacing, float loopSpacing);
void placepk(structure *ct,coordinates *out,int height,int width);

void sortxy(coordinates *coord,bool counter, int height, int width);

} //namespace draw2D
#endif//DRAW2D_H
