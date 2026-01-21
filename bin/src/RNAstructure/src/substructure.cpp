#include "substructure.h"

substructure::substructure() {
	num=0;
}

//	Put a 5' nuc and its x and y coordinates onto the stack
void substructure::putsub(int putnuc5) {
	//cout << "putting "<<putnuc5<<"\n";
   nuc5[num]=putnuc5;
	num++;
   //cout << "num = "<<num<<"\n";
}

//	Return the next 5' nuc and its x and y coordinates
bool substructure::getsub(int *getnuc5){
  	if(num==0) return false;
   else{
   	num--;
      *getnuc5=nuc5[num];
      //cout << "sending "<<*getnuc5<<"\n";
      //cout << "num = "<<num<<"\n";
      return true;
   }
}
