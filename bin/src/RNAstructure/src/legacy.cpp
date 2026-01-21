#define amax 400 //this is a maximum line length for void linout (below)
#define col 80  //this is the number of columns in an output file
#include "algorithm.h"

#include "structure.h"


void linout(structure *ct,char *file) {
//Function to output a CT file to a structure form
//	Writes to file *file

char dash[2],bl[2],number[numlen];
int counter,xx,ip,jp,half,k,j,i,ll,kk,jj,l,stz,countr,ii,pos,p;
stackstruct stack;
arraystruct table;

strcpy (dash,"-");
strcpy (bl," ");

//Open output file
ofstream out(file);


for (counter = 1;counter<=ct->numofstructures;counter++) {

	out << "Structure# " <<counter<<"\n";
	out << ct->ctlabel[counter]<<"\n";
   gcvt((float (ct->energy[counter]))/10,6,number);
	if (ct->energy[counter]!=0)
		out << "energy = "<<number<<"\n";


	
	stz = 0;
	i = 1;
	j = ct->numofbases;
	countr = 0;
	xx = 0;
   stack.sp = 0;
   for (k=1;k<=amax-2;k++) { //empty the array
		for (p=1;p<=6;p++) {
			strcpy(table.array[p][k],bl);
		}
	}

   loop:
		//look for dangling ends
	ip = i;
	jp = j;
	while (ct->basepr[counter][ip]==0) {

	ip++;
		if (ip>=j) {

			if (i>=j) goto outputportion;

			half = ((j-i+2)/(2));//this is a hairpin loop
			ii = i - 1;
			jj = j + 1;
			for (k=1;k<=half;k++) {

				ii++;
				jj--;
				countr++;
				if ((ii%10)==0) digit(1,ii,countr,& table);
				if ((jj%10)==0) digit(6,jj,countr,& table);
				if (k!=half) {
					strcpy(table.array[2][countr],tobase(ct->numseq[ii]));
					strcpy(table.array[5][countr],tobase(ct->numseq[jj]));
				}
				else {
					if (ii<jj) {
					strcpy(table.array[3][countr],tobase(ct->numseq[ii]));
               }
					strcpy(table.array[4][countr],tobase(ct->numseq[jj]));

				}
			}

				goto outputportion;
		}
	}
	while (ct->basepr[counter][jp] ==0) jp--;

	k = ip - i;
	if ((j-jp)>k) k = j-jp;
	if (k!=0) {

		ii = ip;
		jj = jp;
		pos = countr + k + 1;
		for (kk=1;kk<=k;kk++) {
			pos--;
			ii--;
			jj++;
			if (ii>=i) {
				strcpy(table.array[2][pos],tobase(ct->numseq[ii]));
				if ((ii%10)==0) digit(1,ii,pos,& table);
			}
			else strcpy(table.array[2][pos],dash);
			if (jj<=j) {
				strcpy(table.array[5][pos],tobase(ct->numseq[jj]));
				if ((jj%10)==0) digit(6,jj,pos,& table);
			}
			else strcpy(table.array[5][pos],dash);
		}
		countr = countr + k;
	}
	//Stacking or bifuraction

	i = ip;
	j = jp;
	if (ct->basepr[counter][i]!=j) { //bifurcation has occured

		countr = countr+2;
		push(&stack,(ct->basepr[counter][i]+1),j,countr,0);

		j = ct->basepr[counter][i];
	}
	while (ct->basepr[counter][i]==j) {

		countr++;
		strcpy(table.array[3][countr],tobase(ct->numseq[i]));
		strcpy(table.array[4][countr],tobase(ct->numseq[j]));
		if ((i%10)==0) digit(1,i,countr,& table);
		if ((j%10)==0) digit(6,j,countr,& table);
//?		if (i==iret&&j==jret) {
//?		strcpy(table.array[2][countr],"|");
//?		strcpy(table.array[5][countr],"^");
//?		}
		i++;
		j--;
	}

	goto loop;

	outputportion: while (1) {


		out << "\n";

		ll = countr;
		if (col<ll) ll=col;
		for (k=1;k<=6;k++) {
	for (l=1;l<=ll;l++) {

				out << table.array[k][l];

				//cout << table.array[k][l];

}
		out <<"\n";
		}
		if (countr<=col) break;
		for (k=1;k<=5;k++) {
			for (j=1;j<=6;j++) {
				strcpy (table.array[j][k],bl);
				strcpy (table.array[j][k+5],bl);
			}
			strcpy(table.array[2][k+5],".");
			strcpy(table.array[5][k+5],".");
		}
		k = 10;
		ll = col+1;
		for (i=ll;i<=countr;i++) {
			k++;
			for (j=1;j<=6;j++) strcpy(table.array[j][k],table.array[j][i]);
		}
		countr = k;
	}

   pull (&stack,&i,&j,&countr,&xx,&stz);
   

	if (stz==0) {

		for (k=1;k<=amax-2;k++) { //empty the array
			for (p=1;p<=6;p++) {
				strcpy(table.array[p][k],bl);
			}
		}


   goto loop;
   }
}
return;
}


