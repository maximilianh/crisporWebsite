#ifndef DRAW_SUBSTRUCTURE_H
#define DRAW_SUBSTRUCTURE_H

#define maxdomain 100 //the maximum number of domains followed in class
							//substructure
//! substructure

//! this is a stack that keeps track of domains yet
//! to be assigned coordinates in void coordinates
class substructure {
   private:
      int nuc5[maxdomain],num;
      //nuc5 is the number of the 5' base leading into a domain
      //num is the number of domains being stored

   public:
      substructure();
      bool getsub(int *getnuc5);
      void putsub(int putnuc5);
};

#endif //DRAW_SUBSTRUCTURE_H