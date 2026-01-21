#if !defined(GLOBALS_H)
#define GLOBALS_H


//get the names of the free energy data files
void getdat(char *loop, char *stackf, char *tstackh, char *tstacki,
		char *tloop, char *miscloop, char *danglef, char *int22,
      char *int21,char *coax, char *tstackcoax,
      char *coaxstack, char *tstack, char *tstackm, char *triloop,
      char *int11, char *hexaloop, char *tstacki23, char *tstacki1n, char *datapath, bool isRNA);


//get the names of the enathlpy data files
void getenthalpydat(char *loop, char *stackf, char *tstackh, char *tstacki,
	char *tloop, char *miscloop, char *danglef, char *int22,
	char *int21,char *coax, char *tstackcoax,
    char *coaxstack, char *tstack, char *tstackm, char *triloop,
    char *int11, char *hexaloop, char *tstacki23, char *tstacki1n,
	char *datapath, bool isRNA);



#endif