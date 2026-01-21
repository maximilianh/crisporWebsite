#ifndef DRAW_H
#define DRAW_H


//coord is a structure which keeps track of each nucleotide's coordinate
struct coordinates{
	coordinates(int size);
   ~coordinates();
	int *x,*y;
	int **num;//[maxbases/10+1][2];
   short int bases;
};



void place(int number,structure *ct,coordinates *out,int height,int width);
void placepk(structure *ct,coordinates *out,int height,int width);

void sortxy(coordinates *coord,bool counter, int height, int width);


#endif//DRAW_H
