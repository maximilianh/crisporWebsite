
#ifdef _WINDOWS_GUI
	#include "../Windows_interface_deprecated/platform.h"
#else
	#include "platform.h"
#endif //_WINDOWS

#include "rna_library.h"
#include <cmath>
#include <cstdlib>
#include <functional>
#include <iomanip>
#include <cstring>
#include <stdlib.h> 
#include <numeric>

using namespace std;

#define DEBUG_OPENDAT 0
#define outprecision 0

#ifdef COUNTING
	#define NO_COUNT(...) __VA_ARGS__._value
#else
	#define NO_COUNT(...) __VA_ARGS__
#endif //COUNTING


//***********************************code for Structures:

double operator*(const std::vector<double> &a, std::vector<double> &x){
  int i;
  
  double prod = 0;

//  if (a.size() != x.size())
//    cout << "Mismatched Vectors" << a.size() << " x " << x.size() << endl;

//  prod = std::inner_product(a.begin(), a.end(), x.begin(), 0);
  for(i = 0; i < a.size(); i++)
//      prod = fma(a[i], x[i], prod);
      prod += a[i]*x[i];
  return prod;
}

//constuctor
datatable::datatable() : 
    poppen(5, 0) , eparam(11, 0) 
{
	temperature=TEMP_37C; // set default temperature. this can be changed by calling ScaleToTemperature()
    RT=(float)RT_37C;
	alphabetName="";         // the name of the alphabet loaded by a call to opendat (e.g. "rna" "dna" etc)
	loadedAlphabet=false;    // whether the alphabet was loaded (by a call to  
	loadedTables=false;      // whether the thermodynamic parameters were loaded (by a call to opendat with the loadThermo parameter set to true) opendat)
	allowUnknownBases=false; // whether basetonum should return 0 instead of -1 for unknown bases.

#ifdef COUNTING
    loadedCov=false;
#endif
};

//void datatable::setAllowUnknownBases(const bool allow) {
//	allowUnknownBases = allow;
//}
//bool datatable::getAllowUnknownBases() {
//	return allowUnknownBases;
//}
string datatable::GetAlphabetName() {
	return alphabetName;
}
/****fixes log*****
in dev1.4.3 version: 
Datatype: data variables used short int in old version keep using short int type,
Only exceptions are tloop, triloop and hexaloop, the vectors are built of int type.

Matrix dimensional order: 4D vectors use old stack reading order, 
6D, 7D, 8D vectors use old int11, int21, int22 reading order.
read_dangle() uses old dang reading oder
*/


//check if a line contains empty space  or "=", for specification.dat reading
//e.g: X = x = N = n
inline bool is_space_or_equal(char c) {return c==' ' || c=='=' || c=='\r';}

//function floor a float data to a int number by a conversionfactor=10.
// short int floor_to_short(float inum) {
//     short int onum = (short) (floor(inum * conversionfactor + .5));
//     return onum;
// }  
//function floor a char data to a short int number by a conversionfactor=10.
inline int floor_entry_to_short(const char * inum) {
    short int onum;
    if (strcmp(inum,".")){
        onum = (short) (floor(atof(inum) * conversionfactor + .5));
    }
    else{
        onum = INFINITE_ENERGY ;
    }
    return onum;
}  

bool datatable::can_pair(int i, int j, short* sequence)
{   
	return pairing[sequence[i]][sequence[j]];
}



//function searches alphabet table vector for a input base, and return a number which is the index of it.
//The index is based on specification.dat Bases section Alphabet order.
int datatable::basetonum(char base) {
    
    for (int i=0; i<alphabet.size(); i++) {
       //statement that check if a element is in a vector.
        if (std::find(alphabet[i].begin(), alphabet[i].end(), base) != alphabet[i].end()) {
            return i;
        }
    }
    return allowUnknownBases ? 0 : -1;
    
   /*
    unordered_map<char, int>::iterator it; 
    it = alphabet_map.find(base);
    if (it != alphabet_map.end())
        return alphabet_map[base];
    else
        //If the base is not found, return -1
        return allowUnknownBases ? 0 : -1;
    */
}

//function returns the first alphabet entry for a given integer
//This returns '?' for items not found, but this could, of course, be a valid nucleotide in some alphabets 
char datatable::numtobase(int i) {
	if (i>0&&i<alphabet.size()) {
		return alphabet[i][0];
	}
	else return '?';
}

//read in specification.dat to get alphabet vector info, paring bool vector info and not_paring char 1D vector.
bool datatable::read_spec_file(char* ifname) {
    // if the alphabet was read previously, reset it.
    if (alphabet.size()!=0) {
        alphabet.clear();
        pairing.clear();   
        not_pairing.clear();
        non_interacting.clear();
        linker.clear();
		LinkerInts.clear();
    }
    
    std::ifstream spec;
    spec.open(ifname);
    if (!spec) {
        //std::cerr << "\nMissing Data table file!\n";
        return false;
    }

    //declare line and process the lines in the specification.dat file
    std::string line;
    int stage = -1;
    //int count = 0;
    // c0 and c1 are the two indeces for the point where base paring was checked in bool paring matrix.
    int c0;
    int c1;
    //while loop processing the lines of "specification.dat" file
    while(getline(spec, line)) {
        //erase the space and '=' from a line. (also removes '\r' which could be present on Windows)
        // e.g. "X = N = x = n" becomes "XNxn"
        line.erase(std::remove_if(line.begin(), line.end(), is_space_or_equal), line.end());

        //skip empty lines and comments (lines starting with '#')
        if (line.empty() || line[0] == '#') continue;

        //read line is "Bases" then stage=0,     
        if (line=="Bases") {
            stage = 0;
        }
        //read line is "Paring" then stage=1,  
        else if (line=="Pairing") {
            stage = 1;
			//set the pairing array size, it needs to be square with size of the alphabet
			pairing.resize(alphabet.size());
			for (int i=0; i<alphabet.size(); ++i)
				pairing[i].resize(alphabet.size(),false);
        }
        //read line is "Single" then stage=2, 
        else if (line=="Single") {  
            stage = 2;
        }
        //read line is "Non-interacting" then stage=3,
        else if (line=="Non-interacting") {
            stage = 3;
        }
        //read line is "Linker" then stage=4;
        else if (line=="Linker") {
            stage = 4;
        }
        else {
            switch(stage) {
                //in stage=0 case, read in alphabet info.
                case 0:{
                    //count++;
                    //std::vector<char> v;
                   // for(int i=0; i<line.size();i++) {
                    //    v.push_back(line[i]);
                    //}
                    //alphabet.push_back(v);
					alphabet.resize(alphabet.size()+1);
					LinkerInts.resize(alphabet.size()+1);
					alphabet[alphabet.size()-1].resize(line.size());
					for(int i=0; i<line.size();i++) {
                        alphabet[alphabet.size()-1][i]=line[i];

						//identify the special cases of U and A, to help identify nucleotides in pairs that require end penalties.
                        if (line[i] == 'U')
                            baseU=alphabet.size()-1;
						if (line[i]=='A')
							baseA = alphabet.size() - 1;
                    }
                }
                break;

                //in stage=1 case, read in paring info vector.
                case 1:{
                    //same as basetonum(line[0]): get the Paring section's each line's 
                    //first base's corresponding indexing number
                    for (int i=0; i<alphabet.size(); i++) {
                        if (std::find(alphabet[i].begin(), alphabet[i].end(), line[0])!= alphabet[i].end()) {
                            c0 = i;
                            break;
                        }
                    }
                    //same as basetonum(line[1]):get the Pairing section's each line's 
                    //second base's corresponding indexing number
                    for (int i=0; i<alphabet.size(); i++) {
                        if (std::find(alphabet[i].begin(), alphabet[i].end(), line[1]) != alphabet[i].end()) {
                            c1 = i;
                            break;
                        }
                    }
                    //set base paring bool to be true for defined bases pairs in specification.dat file.
                    pairing[c0][c1] = true;
                    pairing[c1][c0] = true;
                    }
                break;

                //in stage=2 case, read in not_paring (Single) info vetor.
                case 2:
                {
                    not_pairing.push_back(line[0]);
                }
                break;

                //in stage=3 case, read in non_interacting info vetor.
                case 3:
                {
                    non_interacting.push_back(line[0]);
                }
                break;       

                //in stage=4 case, read in linker info vector.
                case 4:
                {
                    linker.push_back(line[0]);		
                }
                break;            
            }
        }
    }
    spec.close();

    // construct unordered map from char to index i
    for (int i = 0; i < alphabet.size(); i++) {
        for (int j = 0; j < alphabet[i].size(); j++) {
            alphabet_map.insert(pair<char, int>(alphabet[i][j], i));
        }
    }

    //Post process to put linker information into LinkerInts array:
	for (int i=0;i<LinkerInts.size();++i) {
		LinkerInts[i]=false;
	}
	
	for (int i=0;i<linker.size();++i) {

		LinkerInts[basetonum(linker[i])]=true;

	}

    return true;
}

//function converts a sequence to a number with an alphabet.size()-
//determined factor. 
//This is currently used to encode hairpin loops of 3,4,and 6 unpaired nucleotides (plus closing pair) as an int.
int datatable::seqtonum(std::string seq) {
    
    int n = seq.size();
    int digit_value, seq_value = 0, digit_factor = 1;
    for (int i = 0; i < n; i++) {
        digit_value = basetonum(seq[i]) * digit_factor;
        seq_value += digit_value;
        digit_factor = digit_factor * alphabet.size();
    }
    return seq_value;
}

//function that reads lines from datatable file and store them in to a string vector.
//lines starts with '#' or empty line will be skipped.
//read in data table file to a string vector.
bool datatable::read_datatable_lines(const char* const ifname, std::vector<std::string> &lines){
    std::ifstream ifs(ifname);
    // verify that the file was opened successfully.
    if (!ifs) {
        std::cerr << "\nCritical Error -- Missing Data File: " << ifname << endl;
        return false;
    }
    
    //loop over lines in datatable file.
    std::string line;
    while(getline(ifs, line)) {
        //trim the left empty space of the line, for next step's '#' sign check.
        trimLeft(line);
        //skip the line with '#' and empty.
        if (line.empty() || line[0] == '#') 
            continue;
        //append line to the end of vector v_datatable_lines.
        lines.push_back(line);
    }
    return true;
}

//read in loop.dat type file, get inter, bulge and hairpin info
//bool datatable::read_loop(char* ifname, std::vector<short int> &inter, std::vector<short int> &bulge,
//               std::vector<short int> &hairpin) {
bool datatable::read_loop(char* ifname, std::vector<COUNT_CLASS(short int) > &inter,
                   std::vector<COUNT_CLASS(short int) > &bulge, std::vector<COUNT_CLASS(short int) > &hairpin){
    std::string temp0, temp1, temp2, temp3;
    //declare a 1D int vector with (size of lines in datatable file) to store 
    //internal bulge hairpin data. The data stored in vector indexing from "1"
    //so occupy 0 position with INFINITE_ENERGY
    inter.push_back(INFINITE_ENERGY);
    bulge.push_back(INFINITE_ENERGY);
    hairpin.push_back(INFINITE_ENERGY);
    //declare a vector that holds the lines from the datatable file.
    std::vector<std::string> v_datatable_lines;
    //read in lines from datatable file, return false if file does not exist
    if(!read_datatable_lines(ifname, v_datatable_lines)){
        return false;
    }

    //loop over the vector v_datatable_lines which holds the datatable lines.
    for (int i=0; i<v_datatable_lines.size(); i++) {
        std::istringstream line(v_datatable_lines[i]);
        //get thre values from each line, and send to temp1, temp2, temp3
        line >> temp0 >> temp1 >> temp2 >> temp3;
        //read in data to internal bulge hairpin vector.
        inter.push_back(floor_entry_to_short(temp1.c_str()));
        bulge.push_back(floor_entry_to_short(temp2.c_str()));
        hairpin.push_back(floor_entry_to_short(temp3.c_str()));
    }
    return true;
}

//read in tloop.dat type files, get a 2D matrix info. e.g: triloop.dat, tloop.dat, hexaloop.dat
bool datatable::read_xloop(char* ifname, intVector2D &matrix) {

    std::string seq;
    std::string temp;
    //declare a vector that holds the lines from the datatable file.
    std::vector<std::string> v_datatable_lines;
    //read in lines from datatable file, return false if file does not exist
    if(!read_datatable_lines(ifname, v_datatable_lines)){
        return false;
    }

    //declare a 2D int vector with (size of lines in datatable file) * 2 to store 
    //triloop, tloop, or hexaloop data.
    //matrix = std::vector<std::vector<int> >(v_datatable_lines.size(), std::vector<int>(2,0));
    matrix = std::vector<std::vector<COUNT_CLASS(int) > >(v_datatable_lines.size(), std::vector<COUNT_CLASS(int) >(2,0));

    //loop over the vector v_datatable_lines which holds the datatable lines.
    for (int i=0; i<v_datatable_lines.size(); i++) {
        std::istringstream line(v_datatable_lines[i]);
        //
        line >> seq >> temp;
        //convert sequence to a number by function seqtonum(). 
        //read in data to matrix
        matrix[i][0] = seqtonum(seq);
        //floor_entry_to_short(const char* inum), convert a entry char number to a short int by 
        //conversionfactor 10.
        matrix[i][1] = floor_entry_to_short(temp.c_str());
    }
    return true;
}

//read in dangle.dat type files, get a 4D dangle matrix info
bool datatable::read_dangle(char* ifname, shortVector4D &matrix) {
    //declare a vector that holds the lines from the datatable file.
    std::vector<std::string> v_datatable_lines;
    //read in lines from datatable file, return false if file does not exist
    if(!read_datatable_lines(ifname, v_datatable_lines)){
        return false;
    }

    //initialize dangle matrix vector, index d starts from 0, and old matrix d starts from 1.
    /*
    matrix = std::vector<std::vector<std::vector<std::vector<short int> > > > 
                (alphabet.size(),std::vector<std::vector<std::vector<short int> > > 
                (alphabet.size(),std::vector<std::vector<short int> > 
                (alphabet.size(),std::vector<short int> 
                (3,INFINITE_ENERGY))));
*/

    matrix = std::vector<std::vector<std::vector<std::vector<COUNT_CLASS(short int) > > > > 
                (alphabet.size(),std::vector<std::vector<std::vector<COUNT_CLASS(short int) > > >
                (alphabet.size(),std::vector<std::vector<COUNT_CLASS(short int) > >
                (alphabet.size(),std::vector<COUNT_CLASS(short int) >
                (3,INFINITE_ENERGY))));                
    
    //from a to d, "a" changes slowest, and "d" changes fastest.
    //5'a- c3'  a indicate the base X. 
    //     d    d indicates the location of X, "1" is at 3' end, "2" is at 5' end.
    //3'b- c5'
    //loop over the vector v_datatable_lines which holds the datatable lines.
    int a, b, c, d;

    for (int i=0; i<v_datatable_lines.size(); ) {
        //read the first headline(e.g.:"AX") and get b.
        a = basetonum(v_datatable_lines[i][0]);
        //check length of the first headline, and get d.
        //if it is "AX" format size=2: d=1; else: d=2., old d index starts from 1
        if (v_datatable_lines[i].size() == 2) {
           

			//also make sure the second character is X (not a /r or other whitespace)
			if (v_datatable_lines[i][1]=='X') {
				 d = 1;
			}
			else d = 2;

        }
        else {
            d = 2;
        }
        i++;
        //read the second headline(e.g.:"A"), and get a.
        b = basetonum(v_datatable_lines[i][0]);      
        i++;
        //declare a vector storing the bases (conveted to num) at X position, 
        //using basetonum().
        std::vector <int> v;
        //read the third head_line (e.g.: x line "A   C   G   U ").
        std::istringstream xline(v_datatable_lines[i]);
        i++;
        //loop over the xline items(bases), push_back the number of base to vecor v.
        char base_x;
        while (xline >> base_x) {
            v.push_back(basetonum(base_x));
        }
        //read the data line (e.g.: data line "-0.8 -0.5 -0.8 -0.6").
        std::istringstream data_line(v_datatable_lines[i]);
        i++;
        for (int x=0; x<v.size(); ++x) {
            //get c, the number of base X in vector v. 
            c = v[x];
            std::string temp;
            data_line >> temp;
            //read in data to matrix
            //floor_entry_to_short(const char* inum), convert a entry char number to a short int by 
            //conversionfactor 10.
            matrix[a][b][c][d] = floor_entry_to_short(temp.c_str());
            // const char * entry = temp.c_str();
            // if (strcmp(entry,".")) {
                //convert the float data to int data by CONVERSION_FACTOR
                // matrix[a][b][c][d] = floor_to_short(atof(entry));
            // }    
            // else {
            //     matrix[a][b][c][d] = INFINITE_ENERGY;
            // }
        }
    }
    return true;    
}

//read in stack.dat type files, get a 4D matrix info: stack, tstackh, tstacki, 
//coaxial, tstackcoax, coaxstack, tstack, tstackm, tstacki23, stacki1n
bool datatable::read_4D_table(char* ifname, shortVector4D &matrix) {
    //declare a vector that holds the lines from the datatable file.
    std::vector<std::string> v_datatable_lines;
    //read in lines from datatable file, return false if file does not exist
    if(!read_datatable_lines(ifname, v_datatable_lines)){
        return false;
    } 

    //initialize vector 4D matrix 
    /*
    matrix = std::vector<std::vector<std::vector<std::vector<short int> > > > 
                (alphabet.size(),std::vector<std::vector<std::vector<short int> > > 
                (alphabet.size(),std::vector<std::vector<short int> > 
                (alphabet.size(),std::vector<short int> 
                (alphabet.size(),INFINITE_ENERGY))));
*/

    matrix = std::vector<std::vector<std::vector<std::vector<COUNT_CLASS(short int) > > > >
                (alphabet.size(),std::vector<std::vector<std::vector<COUNT_CLASS(short int) > > >
                (alphabet.size(),std::vector<std::vector<COUNT_CLASS(short int) > >
                (alphabet.size(),std::vector<COUNT_CLASS(short int) >
                (alphabet.size(),INFINITE_ENERGY))));                
    
    //use old stack reading order
    //5'a-c3'
    //3'b-d5'
    //loop over the vector v_datatable_lines which holds the datatable lines.
    int a, b, c, d;

    for (int i=0; i<v_datatable_lines.size(); ) {
        //read the first headline(e.g.:"AX"), and get a.
        a = basetonum(v_datatable_lines[i][0]);
        i++;
        //read the second headline(e.g.:"AY"), and get c.
        b = basetonum(v_datatable_lines[i][0]);      
        i++;

        //declare a vector storing the bases (conveted to num) at Y position, 
        //using basetonum().
        std::vector <int> v;
        //read the third headline (e.g.:"Y" title line "A   C   G   U").
        std::istringstream yline(v_datatable_lines[i]);
        //loop over the yline items(bases), push_back the number of base to 
        //vector v. the v.size determines the next loop size for reading the data.
        char base_y;
        while (yline >> base_y) {
            v.push_back(basetonum(base_y));
        }
        i++;
        for (int x=0; x<v.size(); ++x) {
            char base_x;
            std::istringstream xline(v_datatable_lines[i]);           
            //read each xline's first entry which is the base_x, and get b.
            xline >> base_x; 
            c = basetonum(base_x);
            i++;             
            //read the rest of xline (data) into matrix of v.size()*v.size()
            for (int y=0; y<v.size(); ++y) {
                //get d
                d = v[y];
                std::string temp;
                xline >> temp;
                //read data into matrix
                //floor_entry_to_short(const char* inum), convert a entry char number to a short int by 
                //conversionfactor 10.
                matrix[a][b][c][d] = floor_entry_to_short(temp.c_str());
                // const char * entry = temp.c_str();
                // if (strcmp(entry,".")) {
                //     //convert the float data to int data by CONVERSION_FACTOR
                //     matrix[a][b][c][d] = int (floor((conversionfactor * (atof(entry))) + .5));
                // }    
                // else {
                //     matrix[a][b][c][d] = INFINITE_ENERGY;
                // }
            }
        }
    }
    return true;
}

//read in int11.dat type file, get a 6D matrix info
bool datatable::read_6D_table(char* ifname, shortVector6D &matrix) {
    //declare a vector that holds the lines from the datatable file.
    std::vector<std::string> v_datatable_lines;
    //read in lines from datatable file, return false if file does not exist
    if(!read_datatable_lines(ifname, v_datatable_lines)){
        return false;
    }

    //initialize 6D vector marix.
    /*
    matrix = std::vector<std::vector<std::vector<std::vector<std::vector<std::vector<short int> > > > > >
                (alphabet.size(),std::vector<std::vector<std::vector<std::vector<std::vector<short int> > > > >
                (alphabet.size(),std::vector<std::vector<std::vector<std::vector<short int> > > > 
                (alphabet.size(),std::vector<std::vector<std::vector<short int> > > 
                (alphabet.size(),std::vector<std::vector<short int> > 
                (alphabet.size(),std::vector<short int> 
                (alphabet.size(),INFINITE_ENERGY))))));

*/

    matrix = std::vector<std::vector<std::vector<std::vector<std::vector<std::vector<COUNT_CLASS(short int) > > > > > >
                (alphabet.size(),std::vector<std::vector<std::vector<std::vector<std::vector<COUNT_CLASS(short int) > > > > >
                (alphabet.size(),std::vector<std::vector<std::vector<std::vector<COUNT_CLASS(short int) > > > >
                (alphabet.size(),std::vector<std::vector<std::vector<COUNT_CLASS(short int) > > >
                (alphabet.size(),std::vector<std::vector<COUNT_CLASS(short int) > >
                (alphabet.size(),std::vector<COUNT_CLASS(short int) >
                (alphabet.size(),INFINITE_ENERGY))))));
    //use old int11 reading order.
    //5'abc3'
    //3'def5'
    //loop over the vector v_datatable_lines which holds the datatable lines.
    int a, b, c, d, e, f;

    for (int i=0; i<v_datatable_lines.size(); ) {
        //skip the headline1 as 'X'
        i++;
        //read the first headline2(e.g.:"A A") of the datatable.
        std::istringstream headline2(v_datatable_lines[i]);
        char temp1; 
        char temp2;
        headline2 >> temp1 >> temp2;
        //get a and c from headline2
        a = basetonum(temp1);
        c = basetonum(temp2);
        i++;
        //read the second headline2(e.g.:"U U") of the datatable.
        std::istringstream headline3(v_datatable_lines[i]);
        headline3 >> temp1 >> temp2;
        ///get d and f from headline3
        d = basetonum(temp1); 
        f = basetonum(temp2);      
        //skip the headline4 as 'Y'
        i++;
        //declare a vector storing the number of Y position bases, 
        //with basetonum function.
        i++;
        std::vector <int> v;
        //read the third line (e.g.:"Y" title line) from v_datatable_lines
        std::istringstream yline(v_datatable_lines[i]);
        //loop over the yline items(bases), push_back the number of base to 
        //vector v. the v.size determines the next loop size for reading the data.
        char base_y;
        while (yline >> base_y) {
            v.push_back(basetonum(base_y));
        }
        i++;
        for (int x=0; x<v.size(); ++x) {
            char base_x;
            std::istringstream xline(v_datatable_lines[i]);           
            //read each xline's first entry which is the base_x
            xline >> base_x; 
            //get b from xline fist char, base_x
            b = basetonum(base_x);
            i++;             
            //read the data into matrix of v.size()*v.size()
            for (int y=0; y<v.size(); ++y) {
                //get e from base_y position.
                e = v[y];
                std::string temp;
                xline >> temp;
                //read data into matrix
                //floor_entry_to_short(const char* inum), convert a entry char number to a short int by 
                //conversionfactor 10.
                matrix[a][b][c][d][e][f] = floor_entry_to_short(temp.c_str());
                // const char * entry = temp.c_str();
                // if (strcmp(entry,".")) {
                //     //convert the float data to int data by CONVERSION_FACTOR
                //     matrix[a][b][c][d][e][f] = floor_to_short(atof(entry));
                // }    
                // else {
                //     matrix[a][b][c][d][e][f] = INFINITE_ENERGY;
                // }
            }
        }
    }
    return true;
}

//read in int21.dat type file, get a 7D matrix info
bool datatable::read_7D_table(char* ifname, shortVector7D &matrix) {
    //declare a vector that holds the lines from the datatable file.
    std::vector<std::string> v_datatable_lines;
    //read in lines from datatable file, return false if file does not exist
    if(!read_datatable_lines(ifname, v_datatable_lines)){
        return false;
    }

    //initialize 7D vector marix.
    /*
    matrix = std::vector<std::vector<std::vector<std::vector<std::vector<std::vector<std::vector<short int> > > > > > >
            (alphabet.size(),std::vector<std::vector<std::vector<std::vector<std::vector<std::vector<short int> > > > > >
            (alphabet.size(),std::vector<std::vector<std::vector<std::vector<std::vector<short int> > > > >
            (alphabet.size(),std::vector<std::vector<std::vector<std::vector<short int> > > > 
            (alphabet.size(),std::vector<std::vector<std::vector<short int> > > 
            (alphabet.size(),std::vector<std::vector<short int> > 
            (alphabet.size(),std::vector<short int> 
            (alphabet.size(),INFINITE_ENERGY)))))));
    */

    matrix = std::vector<std::vector<std::vector<std::vector<std::vector<std::vector<std::vector<COUNT_CLASS(short int) > > > > > > >
            (alphabet.size(),std::vector<std::vector<std::vector<std::vector<std::vector<std::vector<COUNT_CLASS(short int) > > > > > >
            (alphabet.size(),std::vector<std::vector<std::vector<std::vector<std::vector<COUNT_CLASS(short int) > > > > >
            (alphabet.size(),std::vector<std::vector<std::vector<std::vector<COUNT_CLASS(short int) > > > >
            (alphabet.size(),std::vector<std::vector<std::vector<COUNT_CLASS(short int) > > >
            (alphabet.size(),std::vector<std::vector<COUNT_CLASS(short int) > >
            (alphabet.size(),std::vector<COUNT_CLASS(short int) >
            (alphabet.size(),INFINITE_ENERGY)))))));

    //use old int21 reading order.
    //5'ac f3'
    //3'bdeg5'
    //   
    //loop over the vector v_datatable_lines which holds the datatable lines.
    int a, b, c, d, e, f, g;

    for (int i=0; i<v_datatable_lines.size(); ) {
        //skip the healdline1 as 'X'
        i++;
        //read the first headline2(e.g.:"A A") of the datatable.
        std::istringstream headline2(v_datatable_lines[i]);
        char temp1; 
        char temp2;
        headline2 >> temp1 >> temp2;
        //get a and f from headline2
        a = basetonum(temp1);
        f = basetonum(temp2);
        i++;
        //read the second headline3(e.g.:"U U") of the datatable.
        std::istringstream headline3(v_datatable_lines[i]);
        headline3 >> temp1 >> temp2;
        //get b and g from headline3
        b = basetonum(temp1); 
        g = basetonum(temp2);      
        //read the headline4 as 'YA'
        i++;
        // get e from headline4
        e = basetonum(v_datatable_lines[i][1]);
        //declare a vector storing the number of Y position bases, 
        //with basetonum function.
        i++;
        std::vector <int> v;
        //read the third line (e.g.:"Y" title line) from v_datatable_lines
        std::istringstream yline(v_datatable_lines[i]);
        //loop over the yline items(bases), push_back the number of base to 
        //vector v. the v.size determines the next loop size for reading the data.
        char base_y;
        while (yline >> base_y) {
            //  << base << endl;
            v.push_back(basetonum(base_y));
        }
        i++;
        for (int x=0; x<v.size(); ++x) {
            char base_x;
            std::istringstream xline(v_datatable_lines[i]);           
            //read each xline's first entry which is the base_x
            xline >> base_x; 
            //get c from xline's first char, base_x
            c = basetonum(base_x);
            i++;  
                       
            //read the data into matrix of v.size()*v.size()
            for (int y=0; y<v.size(); ++y) {
                //get d from base_y location
                d = v[y];
                std::string temp;
                xline >> temp;              
                //read data into matrix
                //floor_entry_to_short(const char* inum), convert a entry char number to a short int by 
                //conversionfactor 10.
                matrix[a][b][c][d][e][f][g] = floor_entry_to_short(temp.c_str());
                // const char * entry = temp.c_str();
                // if (strcmp(entry,".")) {
                //     //convert the float data to int data by CONVERSION_FACTOR
                //     matrix[a][b][c][d][e][f][g] = int (floor((conversionfactor * (atof(entry))) + .5));
                // }
                // else {
                //     matrix[a][b][c][d][e][f][g] = INFINITE_ENERGY;
                // }
            }
        }
    }
    return true;
}

//read in int22.dat type file, get a 8D matrix info
bool datatable::read_8D_table(char* ifname, shortVector8D &matrix) {
    //declare a vector that holds the lines from the datatable file.
    std::vector<std::string> v_datatable_lines;
    //read in lines from datatable file, return false if file does not exist
    if(!read_datatable_lines(ifname, v_datatable_lines)){
        return false;
    }

    //initialize 8D vector marix.
    /*
    matrix = std::vector<std::vector<std::vector<std::vector<std::vector<std::vector<std::vector<std::vector<short int> > > > > > > >
            (alphabet.size(),std::vector<std::vector<std::vector<std::vector<std::vector<std::vector<std::vector<short int> > > > > > >
            (alphabet.size(),std::vector<std::vector<std::vector<std::vector<std::vector<std::vector<short int> > > > > >
            (alphabet.size(),std::vector<std::vector<std::vector<std::vector<std::vector<short int> > > > >
            (alphabet.size(),std::vector<std::vector<std::vector<std::vector<short int> > > > 
            (alphabet.size(),std::vector<std::vector<std::vector<short int> > > 
            (alphabet.size(),std::vector<std::vector<short int> > 
            (alphabet.size(),std::vector<short int> 
            (alphabet.size(),INFINITE_ENERGY))))))));
    */

    matrix = std::vector<std::vector<std::vector<std::vector<std::vector<std::vector<std::vector<std::vector<COUNT_CLASS(short int) > > > > > > > >
            (alphabet.size(),std::vector<std::vector<std::vector<std::vector<std::vector<std::vector<std::vector<COUNT_CLASS(short int) > > > > > > >
            (alphabet.size(),std::vector<std::vector<std::vector<std::vector<std::vector<std::vector<COUNT_CLASS(short int) > > > > > >
            (alphabet.size(),std::vector<std::vector<std::vector<std::vector<std::vector<COUNT_CLASS(short int) > > > > >
            (alphabet.size(),std::vector<std::vector<std::vector<std::vector<COUNT_CLASS(short int) > > > >
            (alphabet.size(),std::vector<std::vector<std::vector<COUNT_CLASS(short int) > > >
            (alphabet.size(),std::vector<std::vector<COUNT_CLASS(short int) > >
            (alphabet.size(),std::vector<COUNT_CLASS(short int) >
            (alphabet.size(),INFINITE_ENERGY))))))));
            
    //use olf int22 reading order
    //5'a e f b3'
    //3'c g h d5'
    //loop over the vector v_datatable_lines which holds the datatable lines.
    int a, b, c, d, e, f, g, h;

    for (int i=0; i<v_datatable_lines.size(); ) {
        //read the first line(e.g.:"A X1 Y1 A") of the datatable.
        std::istringstream headline1(v_datatable_lines[i]);
        std::vector <std::string> v_header1;
        std::string temp1;
        while (headline1 >> temp1) {
             v_header1.push_back(temp1);
        }
        //get a and b based on headline1's first base and last base.
        char base_a = v_header1[0][0]; //get fist char from a string.
        char base_b = v_header1[3][0];
        //get a and b from 
        a = basetonum(base_a);
        b = basetonum(base_b);
        i++;
        //read the second line(e.g.:"AY") of the datatable.
        std::istringstream headline2(v_datatable_lines[i]);
        std::vector <std::string> v_header2;
        std::string temp2;
        while (headline2 >> temp2) {
            v_header2.push_back(temp2);
        }
        //get c and d based on headline2's first base and last base.
        char base_c = v_header2[0][0];
        char base_d = v_header2[3][0];
        c = basetonum(base_c);
        d = basetonum(base_d);    

        i++;
        //declare a vector storing the number of Y1 position bases, 
        //with basetonum function.
        std::vector <int> v_y1;
        //read the third line (e.g.:"Y1" line) from v_datatable_lines
        std::istringstream yline1(v_datatable_lines[i]);
        //loop over the yline items(bases), push_back the number of base to 
        //vector v. the v.size determines the next loop size for reading the data.
        char base_y1;
        while (yline1 >> base_y1) {
            v_y1.push_back(basetonum(base_y1));
        }
        i++;
        //declare a vector storing the number of Y2 position bases, 
        //with basetonum function.
        std::vector <int> v_y2;
        //read the fourth line (e.g.:"Y2" line) from v_datatable_lines
        std::istringstream yline2(v_datatable_lines[i]);
        //loop over the yline items(bases), push_back the number of base to 
        //vector v. the v.size determines the next loop size for reading the data.
        char base_y2;
        while (yline2 >> base_y2) {
            v_y2.push_back(basetonum(base_y2));
        }

        i++;
        for (int x=0; x<v_y1.size(); ++x) {
            std::string base_x;
            char base_x1;
            char base_x2;
            std::istringstream xline(v_datatable_lines[i]);           
            //read each xline's first entry which is the base_x
            xline >> base_x; 
            //get e and g from xline first and scond bases
            base_x1 = base_x[0];
            base_x2 = base_x[1];
            e = basetonum(base_x1);
            g = basetonum(base_x2);
            i++;             
            //read the data into matrix of v.size()*v.size()
            for (int y=0; y<v_y1.size(); ++y) {
                //get f and h from bases of first and second yline 
                f = v_y1[y];
                h = v_y2[y];
                std::string temp;
                xline >> temp;
                //read data into matrix
                //floor_entry_to_short(const char* inum), convert a entry char number to a short int by 
                //conversionfactor 10.
                matrix[a][b][c][d][e][f][g][h] = floor_entry_to_short(temp.c_str());
                // const char * entry = temp.c_str();
                // if (strcmp(entry,".")) {
                //     //convert the float data to int data by CONVERSION_FACTOR
                //     matrix[a][b][c][d][e][f][g][h] =\
                //      int (floor((conversionfactor * (atof(entry))) + .5));
                // }    
                // else {
                //     matrix[a][b][c][d][e][f][g][h] = INFINITE_ENERGY;
                // }
            }
        }
    }
    return true;
}

/*
//read in miscloop.dat
bool datatable::read_miscloop(char* ifname, float &prelog, short int &maxpen, short int &efn2a, 
                   short int &efn2b, short int &efn2c, short int &strain, short int &auend,
                   short int &gubonus, short int &cslope, short int &cint, short int &c3, short int &init,
                   short int &singlecbulge, std::vector<short int> &poppen,
                   std::vector<short int> &eparam) {
*/
bool datatable::read_miscloop(char* ifname, COUNT_CLASS(float)& prelog, COUNT_CLASS(short int)& maxpen,
	COUNT_CLASS(short int)& efn2a, COUNT_CLASS(short int)& efn2b, COUNT_CLASS(short int)& efn2c,
	COUNT_CLASS(short int)& strain, COUNT_CLASS(short int)& auend, bool& AUappliestoGU,  COUNT_CLASS(short int)& gubonus,
                   COUNT_CLASS(short int) &cslope, COUNT_CLASS(short int) &cint, COUNT_CLASS(short int) &c3, COUNT_CLASS(short int) &init,
                   COUNT_CLASS(short int) &singlecbulge, std::vector<COUNT_CLASS(short int) > &poppen,
                   std::vector<COUNT_CLASS(short int) > &eparam){                       

    //


    // float temp;  
    std::string temp;

    //pass a pointer to the vector to the function
    //declare a vector that holds the lines from the datatable file.
    std::vector<std::string> v_datatable_lines;
    //read in lines from datatable file, return false if file does not exist
    if(!read_datatable_lines(ifname, v_datatable_lines)){
        return false;
    }

    //loop over the vector v_datatable_lines which holds the datatable lines.
    for (int i=0; i<v_datatable_lines.size(); ) {

        //read the prelog.
        // temp = atof(v_datatable_lines[i].c_str());
        // prelog = temp * conversionfactor;
        temp = v_datatable_lines[i];
        prelog = (atof(temp.c_str())) * conversionfactor;
        i++;
    
        //read maxpen
        temp = v_datatable_lines[i];
        // maxpen = (short) (temp * CONVERSION_FACTOR + .5);
        // maxpen = floor_to_short(temp);
        maxpen = floor_entry_to_short(temp.c_str());
        i++;

        //read poppen
        std::istringstream poppen_line(v_datatable_lines[i]);
        i++;
        for (int k=1; k<=4; k++){
            poppen_line >> temp;
            // poppen[k] = (short) (temp * CONVERSION_FACTOR + .5);
            poppen[k] = floor_entry_to_short(temp.c_str());
        }

        //read eparam
        std::istringstream eparam_line(v_datatable_lines[i]);
        i++;
        eparam[1] = 0;            // assign some variables that are
        eparam[2] = 0;            // "hard-wired" into code
        eparam[3] = 0;
        eparam[4] = 0;
        eparam_line >> temp;
        // eparam[5] = (short) (floor(temp * CONVERSION_FACTOR + .5));  //constant multi-loop penalty
        eparam[5] = floor_entry_to_short(temp.c_str());
        eparam_line >> temp;
        // eparam[6] = (short) (floor(temp * CONVERSION_FACTOR + .5));  //constant multi-loop penalty
        eparam[6] = floor_entry_to_short(temp.c_str());
        eparam[7] = 30;
        eparam[8] = 30;//This is now deprecated...  This was a limit on internal loop size.
        eparam[9] = -500;
        eparam_line >> temp;
        // eparam[10] = (short) (floor(temp * CONVERSION_FACTOR + .5));
        eparam[10] = floor_entry_to_short(temp.c_str());
        
        //read efn
        std::istringstream efn_line(v_datatable_lines[i]);
        i++;
        efn_line >> temp;
        //(short) (floor(temp * CONVERSION_FACTOR + .5));  //constant multi-loop penalty for efn2
        efn2a = floor_entry_to_short(temp.c_str());
        efn_line >> temp;
        // efn2b = (short) (floor(temp * CONVERSION_FACTOR + .5)); 
        efn2b = floor_entry_to_short(temp.c_str());
        efn_line >> temp;
        // efn2c = (short) (floor(temp * CONVERSION_FACTOR + .5)); 
        efn2c = floor_entry_to_short(temp.c_str());

        //read mlasym
        temp = v_datatable_lines[i];
        i++;
        mlasym = floor_entry_to_short(temp.c_str());
        // (short) (floor(temp * CONVERSION_FACTOR + .5));

        //read strain
        temp = v_datatable_lines[i];
        i++;
        strain = floor_entry_to_short(temp.c_str());
        // (short) (floor(temp * CONVERSION_FACTOR + .5));

        //read the terminal AU penalty
        temp = v_datatable_lines[i];
        i++;
        auend = floor_entry_to_short(temp.c_str());
        // (short) (floor(temp * CONVERSION_FACTOR + .5)); 

		//read whether the terminal AU end penalty also applies to GU pairs
		temp = v_datatable_lines[i];
		i++;
		temp = temp.substr(0, 4);
		if (temp == "true") {
			AUappliestoGU = true;
		}
		else if (temp == "fals") {
			AUappliestoGU = false;
		}
		else {
			//There was a problem here
			return false;
		}

        //read the GGG hairpin bonus
        temp = v_datatable_lines[i];
        i++;
        gubonus = floor_entry_to_short(temp.c_str());
        //(short) (floor(temp * CONVERSION_FACTOR + .5));

        //read the poly c hairpin penalty slope
        temp = v_datatable_lines[i];
        i++;
        cslope = floor_entry_to_short(temp.c_str());
        // (short) (floor(temp * CONVERSION_FACTOR + .5));

        //read the poly c hairpin penalty intercept
        temp = v_datatable_lines[i];
        i++;
        cint = floor_entry_to_short(temp.c_str());
        // (short) (floor(temp * CONVERSION_FACTOR + .5));

        //read the poly c penalty for a loop of 3
        temp = v_datatable_lines[i];
        i++;
        c3 = floor_entry_to_short(temp.c_str());
        // (short) (floor(temp * CONVERSION_FACTOR + .5));

        //read int
        temp = v_datatable_lines[i];
        i++;
        init = floor_entry_to_short(temp.c_str());
        // (short) (floor(temp * CONVERSION_FACTOR + .5));

        //read the single C bulge bonus
        temp = v_datatable_lines[i];
        i++;
        singlecbulge = floor_entry_to_short(temp.c_str());

    }
    return true;
}

//Allocate the data tables
void datatable::allocate_data_tables() {
//	int i,j,k,l,m,n,o,p;
	int i, j, k, l, m, n, o;

	//First allocate the loop tables:
	inter.resize(31);
	bulge.resize(31);
	hairpin.resize(31);

	//now the dangle table:
	dangle.resize(alphabet.size());
	for (i=0;i<alphabet.size();++i) {
		dangle[i].resize(alphabet.size());
		for (j=0;j<alphabet.size();++j) {
			dangle[i][j].resize(alphabet.size());
			for (k=0;k<alphabet.size();++k) dangle[i][j][k].resize(3);
		}
	}

	

	//Now all the 4x4 tables:
	size4D(stack);
	size4D(tstkh);
	size4D(tstki);
	size4D(coax);
	size4D(tstackcoax);
	size4D(coaxstack);
	size4D(tstack);
	size4D(tstkm);
	size4D(tstki23);
	size4D(tstki1n);

	//Now the 1x1, 2x1, and 2x2 internal loops
	iloop11.resize(alphabet.size());
	iloop21.resize(alphabet.size());
	iloop22.resize(alphabet.size());
	for (i=0;i<alphabet.size();++i) {
		iloop11[i].resize(alphabet.size());
		iloop21[i].resize(alphabet.size());
		iloop22[i].resize(alphabet.size());

		for (j=0;j<alphabet.size();++j) {
			iloop11[i][j].resize(alphabet.size());
			iloop21[i][j].resize(alphabet.size());
			iloop22[i][j].resize(alphabet.size());
			
			for (k=0;k<alphabet.size();++k) {
				iloop11[i][j][k].resize(alphabet.size());
				iloop21[i][j][k].resize(alphabet.size());
				iloop22[i][j][k].resize(alphabet.size());

				for (l=0;l<alphabet.size();++l) {
					iloop11[i][j][k][l].resize(alphabet.size());
					iloop21[i][j][k][l].resize(alphabet.size());
					iloop22[i][j][k][l].resize(alphabet.size());

					for (m=0;m<alphabet.size();m++) {
						iloop11[i][j][k][l][m].resize(alphabet.size());
						iloop21[i][j][k][l][m].resize(alphabet.size());
						iloop22[i][j][k][l][m].resize(alphabet.size());

						for (n=0;n<alphabet.size();n++) {
							iloop21[i][j][k][l][m][n].resize(alphabet.size());
							iloop22[i][j][k][l][m][n].resize(alphabet.size());
							
							for (o=0;o<alphabet.size();o++) {
								iloop22[i][j][k][l][m][n][o].resize(alphabet.size());
								
								
							}

						}	

					}
				}
			}
		}
	}

}

//size the 4D tables
void datatable::size4D(shortVector4D &matrix) {
	int i,j,k;
	
	matrix.resize(alphabet.size());
	for (i=0;i<alphabet.size();++i) {
		matrix[i].resize(alphabet.size());
		for (j=0;j<alphabet.size();++j) {
			matrix[i][j].resize(alphabet.size());
			for (k=0;k<alphabet.size();++k) matrix[i][j][k].resize(alphabet.size());
		}
	}

}

// Reset properties in preparation for (re)loading data tables.
void datatable::reset() {
    RT = (float)RT_37C; // we are reading the data-tables from disk so reset RT to its value at 310.15
    temperature=TEMP_37C; // reset temperature (in case it was previously set to some other value)
	loadedAlphabet = loadedTables = false; // reset in case of failure
}

// utilty class to construct data file names. More efficient than string
struct datfile {
    datfile(const char*const directory, const char*const alphabet, const bool useEnthalpies) {
        const int MAX_FILENAME_SIZE = 30; // includes extension
        size_t baseLen = strlen(directory)+strlen(alphabet)+2;
        buffer = new char[baseLen+MAX_FILENAME_SIZE+1];
        // Full file name = directory/alphabet.name.ext (ext includes final dot)
        ext=useEnthalpies?EXT_ENTHALPY:EXT_ENERGY;
        strcpy(buffer, directory);strcat(buffer, "/");strcat(buffer, alphabet);strcat(buffer, ".");
        start = buffer+baseLen;
    }
    char* get(const char*const dataFileName, const char*const fileExtension=NULL) { 
        strcpy(start,dataFileName); strcat(start,fileExtension==NULL?ext:fileExtension); 
        return buffer; 
    }
    //data_directory+"/"+alphabet+"."+
    ~datfile() { delete[] buffer; }
    private:
        char *buffer, *start; const char *ext;
};

// Function opens data files to read thermodynamic data
int datatable::opendat(const char *const directory, const char *const alphabet,  const bool useEnthalpies, const bool skipThermoTables) {
    reset(); // Reset properties in preparation for (re)loading data tables.
    data_directory=is_blank(directory)?getDataPath(alphabet):directory; // store this for later use by ScaleToTemperature etc.
	alphabetName = alphabet;

    // convenience class for generating data file names
    datfile dat(data_directory.c_str(),alphabet,useEnthalpies);

    #if DEBUG_OPENDAT
      cout << "Loading Thermo Tables. opendat(" << data_directory << "," << alphabetName << ")" << (useEnthalpies?" (Enthalpy)":"") << " spec=" << dat.get(F_spec, EXT_SPEC); cout << "  fstack=" << dat.get(F_stack) << endl;
    #endif
	
    //read fspec
    if (!read_spec_file(dat.get(F_spec, EXT_SPEC))) return 0;

	if (!skipThermoTables) {
		//read info from floop for internal loops, hairpin loops, and bulge loops
		if (!read_loop(dat.get(F_loop), inter, bulge, hairpin)) return 0;

		 //read info from fdangle
		if (!read_dangle(dat.get(F_dangle), dangle)) return 0;

		//Read info from fstack
		if (!read_4D_table(dat.get(F_stack), stack)) return 0;

		//Read info from ftstkh
		if (!read_4D_table(dat.get(F_tstackh), tstkh)) return 0;

		//Read info from ftstki
		if (!read_4D_table(dat.get(F_tstacki), tstki)) return 0;

		//Read info from ftstki23
		if (!read_4D_table(dat.get(F_tstacki23), tstki23)) return 0;

		//Read info from ftstki1n
		if (!read_4D_table(dat.get(F_tstacki1n), tstki1n)) return 0;

		//Read info from fcoax 
		if (!read_4D_table(dat.get(F_coax), coax)) return 0;

		//Read info from ftstackcoax 
		if (!read_4D_table(dat.get(F_tstackcoax), tstackcoax)) return 0;

		//Read info from fcoaxstack
		if (!read_4D_table(dat.get(F_coaxstack), coaxstack)) return 0;

		// Read info from ftstack 
		if (!read_4D_table(dat.get(F_tstack), tstack)) return 0; //dangle cases

		//Read info from ftstackm 
		if (!read_4D_table(dat.get(F_tstackm), tstkm)) return 0;

		//Read the 1x1 internal filoop11 
		if (!read_6D_table(dat.get(F_int11), iloop11)) return 0;
		
		//Read the 2x1 internal filoop21 
		if (!read_7D_table(dat.get(F_int21), iloop21)) return 0;
		
		//Read the 2x2 internal filoop22
		if (!read_8D_table(dat.get(F_int22), iloop22)) return 0;
		
		//Read info from ftloops
		if (!read_xloop(dat.get(F_tloop), tloop)) return 0;
		
		//Read info from ftriloops
		if (!read_xloop(dat.get(F_triloop), triloop)) return 0;
		
		//Read info from fhexaloops
		if (!read_xloop(dat.get(F_hexaloop), hexaloop)) return 0;
		
		//Read information from fmiscloop
		if (!read_miscloop(dat.get(F_miscloop), prelog, maxpen, efn2a, efn2b, efn2c, strain,
						auend, AUappliestoGU, gubonus, cslope, cint, c3, init, singlecbulge, poppen,
						eparam)) return 0;
	} else  { // loadThermoTables
		// allocate just what is needed for processDat
		//the dangle table:
		const int sz = this->alphabet.size();
		dangle.resize(sz);
		for (int i=0;i<sz;++i) {
			dangle[i].resize(sz);
			for (int j=0;j<sz;++j) {
				dangle[i][j].resize(sz);
				for (int k=0;k<sz;++k) dangle[i][j][k].resize(3);
			}
		}
		//Now all the 4x4 tables:
		size4D(tstack);  size4D(tstkh);
		size4D(tstki);   size4D(tstki1n);
		size4D(tstki23); size4D(tstkm);
		//size4D(stack);size4D(tstkh);size4D(tstki);size4D(coax);size4D(tstackcoax);
		//size4D(coaxstack);size4D(tstack);size4D(tstkm);size4D(tstki23);size4D(tstki1n);
	}
	processDat();
	loadedAlphabet = true;
	loadedTables = !skipThermoTables;
	return 1;
}

void datatable::processDat() {
    // //Fill out some variables used for backwards compatability:
    numoftloops = (short) tloop.size();
    numoftriloops = (short) triloop.size();
    numofhexaloops = (short) hexaloop.size();

	
	//handle X, a nucleotide that cannot pair or stack; these are listed in non_interacting:
	//First, build a list of ints to which this applies:
	std::vector<int> non_interacting_ints;
	non_interacting_ints.resize(non_interacting.size());
	for (int i=0; i<non_interacting.size();++i) {
		non_interacting_ints[i]=basetonum(non_interacting[i]);
	}
	//also, make a bool array for non-interacters:
	//can_interact[i] is false for nucleotide types i that cannot intercat
	std::vector<bool> can_interact;
	can_interact.resize(alphabet.size());
	for (int i=0;i<can_interact.size();++i) can_interact[i]=true;
	for (int i=0;i<non_interacting.size();++i) {
		can_interact[basetonum(non_interacting[i])]=false;
	}
	//also, note that linker nucs also cannot interact:
	for (int i=0;i<linker.size();++i) {
		can_interact[basetonum(linker[i])]=false;
	}

	//now handle individual tables:

	//First handle dangle:
	for (int i=0;i<alphabet.size();++i) {
		for (int j=0;j<alphabet.size();++j) {
			for (int k=0; k<non_interacting_ints.size();++k) {

				if (can_interact[i]&&can_interact[j]) {
					dangle[i][j][non_interacting_ints[k]][0]=0;
					dangle[i][j][non_interacting_ints[k]][1]=0;
				}
			}
		}
	}

	//now handle the stack tables:
	//First handle dangle:
	for (int i=0;i<alphabet.size();++i) {
		for (int j=0;j<alphabet.size();++j) {
			for (int k=0; k<alphabet.size();++k) {
				for (int l=0; l<alphabet.size();++l) {

					if (!can_interact[k]||!can_interact[l]) {
						if (can_interact[i]&&can_interact[j]) {
							
							tstkh[i][j][k][l]=0;
							tstki[i][j][k][l]=0;
							tstki23[i][j][k][l]=0;
							tstki1n[i][j][k][l]=0;
							tstack[i][j][k][l]=0;
							tstkm[i][j][k][l]=0;
						}
					}
				}
			}
		}
	}


	
	//handle I, a nucleotide that is the intermolecular linker; these are listed in linker:
	//First, build a list of ints to which this applies:
	std::vector<int> linker_ints;
	linker_ints.resize(linker.size());
	for (int i=0; i<linker.size();++i) {
		linker_ints[i]=basetonum(linker[i]);
	}
	/*  This block of code is depracated because the fucntionality was moved to the class for use in multiple places
	This is now duplicate code://also, make a bool array for linkers:

	//is_linker is true for nucleotides types i that are linker nucs
	
	std::vector<bool> is_linker;
	is_linker.resize(alphabet.size());
	for (int i=0;i<is_linker.size();++i) is_linker[i]=false;
	for (int i=0;i<linker.size();++i) {
		is_linker[basetonum(linker[i])]=true;
	}
	
	*/

	//now handle individual tables:

	//First handle dangle:
	for (int i=0;i<alphabet.size();++i) {
		for (int j=0;j<alphabet.size();++j) {
			for (int k=0; k<linker_ints.size();++k) {

				if (can_interact[i]&&can_interact[j]) {
					dangle[i][j][linker_ints[k]][0]=0;
					dangle[i][j][linker_ints[k]][1]=0;
				}
			}
		}
	}

#ifdef COUNTING
    int tmp = auend.get;
#endif //COUNTING

	//now handle the stack tables:
	//Loop over the size of the alphabet for the 4 table positions
	for (int i=0;i<alphabet.size();++i) {
		for (int j=0;j<alphabet.size();++j) {
			for (int k=0; k<alphabet.size();++k) {
				for (int l=0; l<alphabet.size();++l) {

					
					if (can_interact[i]&&can_interact[j]) {
						
						
						if (isLinker(k)||isLinker(l)) {
							
							tstkh[i][j][k][l]=0;
							tstki[i][j][k][l]=0;
							tstki23[i][j][k][l]=0;
							tstki1n[i][j][k][l]=0;
							
							
							if (isLinker(k)&&isLinker(l)) {
								tstack[i][j][k][l]=0;
								tstkm[i][j][k][l]=0;

							}
							else if (isLinker(k)) {
								//These cases need to be handled as dangling ends:
								tstack[i][j][k][l]=NO_COUNT(dangle[i][j][l][2]);
								tstkm[i][j][k][l]=NO_COUNT(dangle[i][j][l][2])+penalty2(i,j,this);
							}
							else {
								//isLinker(l)

								//These cases need to be handled as dangling ends:
								tstack[i][j][k][l]=NO_COUNT(dangle[i][j][k][1]);
								tstkm[i][j][k][l] = NO_COUNT(dangle[i][j][k][1])+penalty2(i,j,this);

							}
						}//end k or l is a linker
					}//end if i and j can both interact
				}
			}
		}
	}

#ifdef COUNTING
    auend.get = tmp;
#endif //COUNTING
}

// Reads in the *.dh enthalpy tables and calls dG_T to calculate the free energy tables at the specified temperature.
// \param temperature The new temperature in degrees Kelvin
// \return Returns 0 on success. On failure, returns an error code compatible with RNA::GetErrorMessage
int datatable::ScaleToTemperature(const double temperature) {
    if (!loadedTables) return 30; // error -- data table must have already been loaded.

    //temperature is altered, so the enthalpy tables need to be read:
    datatable *localenthalpy = new datatable();
    //open the enthlpy parameters and check for errors

    #if DEBUG_OPENDAT
      cout << "Loading ENTHALPY Tables. opendat(" << data_directory << "," << alphabetName << ")"  << endl;
    #endif    

    if (localenthalpy->opendat(data_directory.c_str(), this->alphabetName.c_str(), true)==0) {
        delete localenthalpy;
        return 5;//an error has occured
    }
    //using the enthalpy parameters and the folding free energy changes at 37 degrees C,
    //	calculate folding free energy changes at temp, the current temperature, and overwrite the 
    //	folding free energy changes at 37 degrees C
    dG_T((float)temperature,*this,*localenthalpy,*this);

    this->temperature = temperature;

    //remove the enthalpy parameters, there is no need to keep them
    delete localenthalpy;
    return 0; //success
}

//Return true if this integer is a base that is the intermolecular linker
bool datatable::isLinker(const int baseNumber) {

	return LinkerInts[baseNumber];

}

//Return true if this nucleotide is a base that is the intermolecular linker
bool datatable::isLinker(const char baseChar) {

	if (std::find(linker.begin(), linker.end(), baseChar) == linker.end()) {
		return true;
	}
	else return false;

}

void de_allocate (int **v,int i) {//this function deallocates the memory used
												//in an array
int a;

	for (a=0;a<i;a++) {
   	delete[] v[a];
   }
   delete[] v;
}

void de_allocate (short int **v,int i) {//this function deallocates the memory used
												//in an array
int a;

	for (a=0;a<i;a++) {
   	delete[] v[a];
   }
   delete[] v;
}

//this function deallocates the memory used in an array											
void de_allocate (bool **v,int i) 
{   
	for (int a = 0; a < i; a++)  delete[] v[a];
    delete[] v;
}

//read is used to read data from a save file
#ifdef COUNTING
void read(ifstream *out,COUNT_CLASS(short int) *i) {

	out->read((char *) i,sizeof(*i));
}

void read(ifstream *out,COUNT_CLASS(int) *i) {

	out->read((char *) i,sizeof(*i));
}

void read(ifstream *out,COUNT_CLASS(float) *i) {

	out->read((char *) i,sizeof(*i));
}

#endif //COUNTING

void read(ifstream *out,short *i) {

	out->read((char *) i,sizeof(*i));
}

void read(ifstream *out,bool *i) {

	out->read((char *) i,sizeof(*i));
}

void read(ifstream *out,int *i) {

	out->read((char *) i,sizeof(*i));
}

void read(ifstream *out,char *i) {
	int length;

	out->read((char *) (&length), sizeof(length));
	out->read(i,length);
}

void read(ifstream *out,string *i) {
	int length,it;
	char c;


	//First read the length of the string
	out->read((char *) (&length), sizeof(length));

	//Now read the string, one character at a time:
	for (it=0;it<length;++it) {

		out->read(&c,1);
		(*i)+=c;

	}
}

void readsinglechar(ifstream *out,char *i) {

	out->read(i,sizeof(char));
}

void read(ifstream *out,float *i) {

	out->read((char *) i,sizeof(*i));
}

void read(ifstream *out,double *i) {

	out->read((char *) i,sizeof(*i));
}

void read(ifstream *out,log_double *i) {
	out->read((char *) i,sizeof(*i));
}

void read(ifstream *out, datatable* data){
	//read the alphabet info
	read(out, &(data->alphabet));
	read(out, &(data->pairing));
	read(out, &(data->not_pairing));
	read(out, &(data->non_interacting));
	read(out, &(data->linker));


	//build the LinkerInts array:
	data->LinkerInts.resize(data->alphabet.size());
	for (int li=0;li<data->LinkerInts.size();++li) {
		data->LinkerInts[li]=false;

	}
	
	for (int li=0;li<data->linker.size();++li) {

		data->LinkerInts[data->basetonum(data->linker[li])]=true;

	}

	//now open the data files:
	data->allocate_data_tables();

	read(out, &data->poppen);
	read(out, &data->maxpen);
	read(out, &data->eparam);
	read(out, &data->inter);
	read(out, &data->bulge);
	read(out, &data->hairpin);
	read(out, &data->dangle);
	read(out, &data->stack);
	read(out, &data->tstkh);
	read(out, &data->tstki);
	read(out, &data->coax);
	read(out, &data->tstackcoax);
	read(out, &data->coaxstack);
	read(out, &data->tstack);
	read(out, &data->tstkm);
	read(out, &data->tstki23);
	read(out, &data->tstki1n);
	read(out, &data->iloop11);

	for (int i=0;i<data->alphabet.size();i++) {
		for (int j=0;j<data->alphabet.size();j++) {
			for (int k=0;k<data->alphabet.size();k++) {
				for (int l=0;l<data->alphabet.size();l++) {
					for (int m=0;m<data->alphabet.size();m++) {
						for (int n=0;n<data->alphabet.size();n++) {
							for (int o=0;o<data->alphabet.size();o++) {
								if (data->pairing[i][j]&&data->pairing[n][o]) read(out,&(data->iloop21[i][j][k][l][m][n][o]));
								else data->iloop21[i][j][k][l][m][n][o]=INFINITE_ENERGY;
								for (int p=0;p<data->alphabet.size();p++) {
									if (data->pairing[i][k]&&data->pairing[j][l])
										read(out,&(data->iloop22[i][j][k][l][m][n][o][p]));
									else data->iloop22[i][j][k][l][m][n][o][p]=INFINITE_ENERGY;
								}
							}
						}
					}
				}
			}
		}
	}

	read(out, &data->numoftloops);
	read(out, &data->tloop);
	read(out, &data->numoftriloops);
	read(out, &data->triloop);
	read(out, &data->numofhexaloops);
	read(out, &data->hexaloop);

	read(out,&(data->auend));
	read(out,&(data->AUappliestoGU));
	read(out,&(data->gubonus));
	read(out,&(data->cint));
	read(out,&(data->cslope));
	read(out,&(data->c3));
	read(out,&(data->efn2a));
	read(out,&(data->efn2b));
	read(out,&(data->efn2c));
	read(out,&(data->init));
	read(out,&(data->mlasym));
	read(out,&(data->strain));
	read(out,&(data->prelog));
	read(out,&(data->singlecbulge));
}

//write is used to write data to a save file
#ifdef COUNTING
void write(ofstream *out,COUNT_CLASS(short int) *i) {

	out->write((char *) i,sizeof(*i));
}

void write(ofstream *out,COUNT_CLASS(int) *i) {

	out->write((char *) i,sizeof(*i));
}

void write(ofstream *out,COUNT_CLASS(float) *i) {

	out->write((char *) i,sizeof(*i));
}

#endif //COUNTING

void write(ofstream *out,short *i) {

	out->write((char *) i,sizeof(*i));
}

void write(ofstream *out,bool *i) {

	out->write((char *) i,sizeof(*i));
}

void write(ofstream *out,int *i) {

	out->write((char *) i,sizeof(*i));
}

void write(ofstream *out,char *i) {
	int length;

	length = (int) strlen(i)+1;

	out->write((char *) (&length), sizeof(length));
	out->write(i,strlen(i)+1);
}

void write(ofstream *out,string *i) {
	int length;

	length = (int) i->size();
	out->write((char *) (&length), sizeof(length));
	out->write(i->c_str(),length);
}


void write(ofstream *out,float *i) {

	out->write((char *) i,sizeof(*i));
}

void write(ofstream *out,double *i) {

	out->write((char *) i,sizeof(*i));
}

void write(ofstream *out,log_double *i) {

	out->write((char *) i,sizeof(*i));
}

void writesinglechar(ofstream *out, char *i) {
	out->write(i,sizeof(char));

}

void write(ofstream* out, datatable* data){

	write(out, &(data->alphabet));
	write(out, &(data->pairing));
	write(out, &(data->not_pairing));
	write(out, &(data->non_interacting));
	write(out, &(data->linker));

	write(out, &data->poppen);
	write(out, &data->maxpen);
	write(out, &data->eparam);
	write(out, &data->inter);
	write(out, &data->bulge);
	write(out, &data->hairpin);
	write(out, &data->dangle);
	write(out, &data->stack);
	write(out, &data->tstkh);
	write(out, &data->tstki);
	write(out, &data->coax);
	write(out, &data->tstackcoax);
	write(out, &data->coaxstack);
	write(out, &data->tstack);
	write(out, &data->tstkm);
	write(out, &data->tstki23);
	write(out, &data->tstki1n);
	write(out, &data->iloop11);

	for (int i=0;i<data->alphabet.size();i++) {
		for (int j=0;j<data->alphabet.size();j++) {
			for (int k=0;k<data->alphabet.size();k++) {
				for (int l=0;l<data->alphabet.size();l++) {
					for (int m=0;m<data->alphabet.size();m++) {
						for (int n=0;n<data->alphabet.size();n++) {
							for (int o=0;o<data->alphabet.size();o++) {
								if (data->pairing[i][j]&&data->pairing[n][o]) write(out,&(data->iloop21[i][j][k][l][m][n][o]));
								for (int p=0;p<data->alphabet.size();p++) {
									if (data->pairing[i][k]&&data->pairing[j][l])
										write(out,&(data->iloop22[i][j][k][l][m][n][o][p]));
								}
							}
						}
					}
				}
			}
		}
	}
	write(out, &data->numoftloops);
	write(out, &data->tloop);
	write(out, &data->numoftriloops);
	write(out, &data->triloop);
	write(out, &data->numofhexaloops);
	write(out, &data->hexaloop);

	write(out,&(data->auend));
	write(out,&(data->AUappliestoGU));
	write(out,&(data->gubonus));
	write(out,&(data->cint));
	write(out,&(data->cslope));
	write(out,&(data->c3));
	write(out,&(data->efn2a));
	write(out,&(data->efn2b));
	write(out,&(data->efn2c));
	write(out,&(data->init));
	write(out,&(data->mlasym));
	write(out,&(data->strain));
	write(out,&(data->prelog));
	write(out,&(data->singlecbulge));
}

//Calculate free energy tables at temperature T(dg),with dg 37(data) and dh 37(dhdata).
//T is in Kelvins.
void dG_T(float T, datatable &data, datatable &dhdata, datatable &dg)
{
	int a,b,c,d,e,f,g,h;


	//calculated with T
	dg.prelog=data.prelog*T/((float) 310.15);
    dg.RT=data.RT*T/((float) 310.15);
	 for(a=0;a<data.numoftloops;a++)
   	 (dg.tloop[a][1])=Tscale(T,data.tloop[a][1],dhdata.tloop[a][1]);
	 for(a=0;a<data.numoftriloops;a++)
	 (dg.triloop[a][1])=Tscale(T,data.triloop[a][1],dhdata.triloop[a][1]);
	 for(a=0;a<data.numofhexaloops;a++)
 	 (dg.hexaloop[a][1])=Tscale(T,data.hexaloop[a][1],dhdata.hexaloop[a][1]);


	for(a=1;a<5;a++)
			dg.poppen[a]=Tscale(T,data.poppen[a],dhdata.poppen[a]);

			dg.maxpen=Tscale(T,data.maxpen,dhdata.maxpen);
	for(a=1;a<7;a++)
		        dg.eparam[a]=Tscale(T,data.eparam[a],dhdata.eparam[a]);
		        dg.eparam[10]=Tscale(T,data.eparam[10],dhdata.eparam[10]);



	for(a=1;a<31;a++)

	{
		dg.inter[a]=Tscale(T,data.inter[a],dhdata.inter[a]);
		dg.bulge[a]=Tscale(T,data.bulge[a],dhdata.bulge[a]);
		dg.hairpin[a]=Tscale(T,data.hairpin[a],dhdata.hairpin[a]);


	}

	dg.auend=Tscale(T,data.auend,dhdata.auend);
	dg.gubonus=Tscale(T,data.gubonus,dhdata.gubonus);
	dg.cint=Tscale(T,data.cint,dhdata.cint);
	dg.cslope=Tscale(T,data.cslope,dhdata.cslope);
	dg.c3=Tscale(T,data.c3,dhdata.c3);
	dg.efn2a=Tscale(T,data.efn2a,dhdata.efn2a);
	dg.efn2b=Tscale(T,data.efn2b,dhdata.efn2b);
	dg.efn2c=Tscale(T,data.efn2c,dhdata.efn2c);
	dg.init=Tscale(T,data.init,dhdata.init);
	dg.mlasym=Tscale(T,data.mlasym,dhdata.mlasym);
	dg.strain=Tscale(T,data.strain,dhdata.strain);
	dg.singlecbulge=Tscale(T,data.singlecbulge,dhdata.singlecbulge);


	for(a=0;a<6;a++){
		for(b=0;b<6;b++){
		for(c=0;c<6;c++){
		for(d=0;d<6;d++){
		for(e=0;e<6;e++){
		for(f=0;f<6;f++){
			(dg.iloop11[a][b][c][d][e][f])=Tscale(T,data.iloop11[a][b][c][d][e][f],dhdata.iloop11[a][b][c][d][e][f]);
		for(g=0;g<6;g++){
			(dg.iloop21[a][b][c][d][e][f][g])=Tscale(T,data.iloop21[a][b][c][d][e][f][g],dhdata.iloop21[a][b][c][d][e][f][g]);

		for(h=0;h<6;h++){
			(dg.iloop22[a][b][c][d][e][f][g][h])=Tscale(T,data.iloop22[a][b][c][d][e][f][g][h],dhdata.iloop22[a][b][c][d][e][f][g][h]);
			}
			}
			}
			}
		(dg.tstki1n[a][b][c][d])=Tscale(T,data.tstki1n[a][b][c][d],dhdata.tstki1n[a][b][c][d]);
		(dg.tstki23[a][b][c][d])=Tscale(T,data.tstki23[a][b][c][d],dhdata.tstki23[a][b][c][d]);
		(dg.tstkm[a][b][c][d])=Tscale(T,data.tstkm[a][b][c][d],dhdata.tstkm[a][b][c][d]);
		(dg.tstkh[a][b][c][d])=Tscale(T,data.tstkh[a][b][c][d],dhdata.tstkh[a][b][c][d]);
		(dg.tstack[a][b][c][d])=Tscale(T,data.tstack[a][b][c][d],dhdata.tstack[a][b][c][d]);
		(dg.coaxstack[a][b][c][d])=Tscale(T,data.coaxstack[a][b][c][d],dhdata.coaxstack[a][b][c][d]);
		(dg.tstackcoax[a][b][c][d])=Tscale(T,data.tstackcoax[a][b][c][d],dhdata.tstackcoax[a][b][c][d]);
		(dg.coax[a][b][c][d])=Tscale(T,data.coax[a][b][c][d],dhdata.coax[a][b][c][d]);
		(dg.tstki[a][b][c][d])=Tscale(T,data.tstki[a][b][c][d],dhdata.tstki[a][b][c][d]);
		(dg.stack[a][b][c][d])=Tscale(T,data.stack[a][b][c][d],dhdata.stack[a][b][c][d]);
		if(d==1||d==2)
		(dg.dangle[a][b][c][d])=Tscale(T,data.dangle[a][b][c][d],dhdata.dangle[a][b][c][d]);
			}
			}
			}
			}


}

//This needs to be depracated, the extended alphabet needs to be handles by structure, which reads the information from disk.
/*short int tonumi(char *base)	{
short int	a;
if (!strcmp(base,"A")||!strcmp(base,"B")) (a = 1);
else if(!strcmp(base,"C")||!strcmp(base,"Z")) (a = 2);
else if(!strcmp(base,"G")||!strcmp(base,"H")) (a = 3);
else if(!strcmp(base,"U")||!strcmp(base,"V")) (a = 4);
else if(!strcmp(base,"T")||!strcmp(base,"W")) (a = 4);
else (a=0);  //this is for others, like X
return a;
}*/

//this function calculates whether a terminal pair i,j requires the end penalty
integersize penalty2(int i,int j, datatable *data) {


	if (std::find(data->alphabet[i].begin(), data->alphabet[i].end(), 'U') != data->alphabet[i].end()) {
		return data->auend;
	}
	else if (std::find(data->alphabet[j].begin(), data->alphabet[j].end(), 'U') != data->alphabet[j].end()) {
		return data->auend;
	}
	else return 0;

   //if (i==4||j==4)

   	//return data->auend;
   //else return 0;//no end penalty


}

#define DATAPATH_DEBUG false // show debug messages when setting datapath
#define DATAPATH_SUBDIR "data_tables" // default directory name (usually under RNAstructure) that holds the thermodynamic data files.
#define DATAPATH_AUTODETECT_DIRS "./data_tables", "../data_tables", "../../data_tables", ".", "..", "../.."  // look for parameter files in e.g. "./data_tables" "../data_tables" as well as ./ and ../  etc.
#define DATAPATH_AUTODETECT_FILES "rna.specification.dat", "dna.specification.dat", "autodetect.dat" // In order to auto-detect the location of the parameter files, at least one of these files must exist in DATAPATH.
#define DATAPATH_DEFAULT "." // path to use if all else fails.
#define DATAPATH_URL "RNAstructure/manual/Text/Thermodynamics.html" // path to information about DATAPATH
#define ARR_COUNT(var) (sizeof(var)/sizeof(*var))
static char CachedDataPath[maxfil]; // the cached locaton of the parameter files.
// Defines some warning messages that should be shown for various DATAPATH errors.
enum DataPathWarning { DP_WARN_NONE=0, DP_WARN_AUTODETECT_OK, DP_WARN_AUTODETECT_FAIL, DP_WARN_ENV_FAIL_DIR, DP_WARN_ENV_FAIL_FILE };

// determine if the specified directory contains the datatables. specFile can be NULL or the name of a specific file to search for.
bool testDataPath(const char* dir, const char* const alphabetName) {
    // if an alphabet is specified, verify DATAPATH by searching for the required specifications file.
    if (!is_blank(alphabetName) && fileExists(string(dir)+"/"+alphabetName+"." F_spec EXT_SPEC))  // "<alphabet>.specification.dat"
        return true;
    const char* const files[] = { DATAPATH_AUTODETECT_FILES };
    for(int f=0; f<ARR_COUNT(files); f++)
        if (fileExists(dir, files[f])) return true;
    return false;
}

// Show a warning message to the user for various DATAPATH errors.
void showDataPathWarning(const DataPathWarning warning, const char*const path) {
        // If a warning is appropriate, show the warning message.
    if (warning == DP_WARN_AUTODETECT_OK)
        cerr<<"\nUsing auto-detected DATAPATH: \""<<path<<"\" (set DATAPATH to avoid this warning)."<<endl;
    else if (warning != DP_WARN_NONE) {
        cerr<<"\nThe thermodynamic parameter files could not be located! ";
        switch(warning) {
            case DP_WARN_AUTODETECT_FAIL: cerr<<"(auto-detection failed)."; break;
            case DP_WARN_ENV_FAIL_DIR:    cerr<<"(DATAPATH is set to an invalid directory)."; break;
            case DP_WARN_ENV_FAIL_FILE:   cerr<<"(DATAPATH is set, but could not be verified)."; break;
        }
        cerr<<"\nSet the DATAPATH environment variable (e.g. [path/to/RNAstructure]/data_tables) to avoid this error.\nSee " DATAPATH_URL ".)"<<endl; 
    }
}

//!Get the directory location of thermodynamic parameter files.
// In order of priority, this will be:
//  1. Any non-NULL value set by a prior call to setDataPath.
//  2. The location pointed to by the DATAPATH environment variable.
//  3. Any of several standard locations (e.g. ./data_tables or ../data_tables etc) if any are valid.
//  4. The value of DATAPATH_DEFAULT (usually ".").
//
// This function calls setDataPath with whatever value is returned, so that so 
//   that all subsequent calls will return the cached value immediately.
//
// If alphabetName is not NULL, it is used to validate potential 
//   directories during auto-detection. (i.e. the directory must contain the 
//   file "<alphabet>.specification.dat".
const char* getDataPath(const char* const alphabetName) {
	// if it was previously set, return immediately.
    if (!is_blank(CachedDataPath)) return CachedDataPath;

    // Test whether the DATAPATH environment variable is set.
	const char* const env = getenv(DATAPATH_ENV_VAR);
	if (!is_blank(env)) {
        if (!dirExists(env)) // Determine whether or not it is a valid directory.
            showDataPathWarning(DP_WARN_ENV_FAIL_DIR, env);
        else if (!testDataPath(env, alphabetName))  // The directory exists. Verify it by searching for one of the required data files.
            showDataPathWarning(DP_WARN_ENV_FAIL_FILE, env);
		return setDataPath(env); // Always honor the path set by the user, regardless of whether it passed verification.
    }

	// Only continue if the DATAPATH environment variable was not set.
	const char* const dirs[] =  { DATAPATH_AUTODETECT_DIRS };
    for(int d=0; d<ARR_COUNT(dirs); d++)
        if (testDataPath(dirs[d], alphabetName)) {
            showDataPathWarning(DP_WARN_AUTODETECT_OK, dirs[d]); // verification was successful. Show a warning anyway, to let the user know which path was auto-detected.
            return setDataPath(dirs[d]); // Even if it failed verification, honor the path set by the user.
        }

    // DATAPATH was not set and could not be auto-detected.
    showDataPathWarning(DP_WARN_AUTODETECT_FAIL, DATAPATH_DEFAULT);
	return setDataPath(DATAPATH_DEFAULT); // set path to default to avoid repeated auto-detection attempts.
}

// Set the location where programs should look for the thermodynamic parameter files.
// Programs that call getDataPath will receive the new location.
// (This also sets the DATAPATH environment variable for backwards compatibility.)
// This is primarily necessary for the Java interface or other high-level clients
//   that might want to set the location, but are not able to modify the environment 
//   variable directly.
// the function creates a copy of the character data pointed to by the path argument (to
// prevent unintentional modification or deletion by the caller) and the function returns 
// a pointer to the copyied value.
const char* setDataPath(const char* const path) {
	if (DATAPATH_DEBUG) cout<<"setDataPath: "<<path<<endl;
    if (strcmp(CachedDataPath,path)==0) return CachedDataPath; //already equal. Just return it.
	
	// make a copy in case the caller accidentally frees/deletes *path
	strncpy(CachedDataPath, path, sizeof(CachedDataPath)-1); //min(strlen(path)+1,

	// Set the value of the evironment variable (as a backup for legacy code 
    // that still checks DATAPATH).
	#ifdef _WIN32
	_putenv_s(DATAPATH_ENV_VAR, CachedDataPath);
	#else
	setenv(DATAPATH_ENV_VAR, CachedDataPath, 0);
	#endif
	if (DATAPATH_DEBUG) cout<<"new DATAPATH: "<<CachedDataPath<<endl;
	return CachedDataPath;
}

int ReadRestraints(vector<double> &v, const char* SHAPEFile, double maximum) {
  // read data from shape file
  int pos;
  double data;
  
  // Set the default value of each position if the array is not empty.
  fill(v.begin(), v.end(), DEFAULT_RESTRAINT_VALUE); // DEFAULT_RESTRAINT_VALUE is something like -999.0
  
  ifstream infile(SHAPEFile);
  if (!infile.good()) return ERR_BAD_RESTRAINT_FILE; // file could not be opened
  while(infile >> pos >> data) {
	  // The base index (pos) should never be less than 1 because nucleobase numbering starts at 1.
	  // A base index greater than 500000 probably indicates that the file format is bad, because we probably can't accomodate this many bases.
	  if (pos < 1 || pos > 500000) return ERR_BAD_RESTRAINT_NUC_POS; 
	  if (data > -500.0) { // just ignore values less than this
		  if (pos > v.size()) {
			  // If the array is too small, resize it by adding elements with the default value.
			  v.resize(pos, DEFAULT_RESTRAINT_VALUE);
		  }
		  // Set the value at the specified position.
		  if (maximum > 0 && data > maximum) data = maximum;
		  v[pos-1] = data;
	  }
  }
  return 0;
}

// Write restraints to a file
int WriteRestraints(const vector<double> &v, const string& outfile, const bool append) {
    ofstream out(outfile.c_str(), append ? std::ios_base::app : std::ios_base::trunc);
    if (!out.good()) return ERR_BAD_RESTRAINT_FILE; 
    for(int i=0;i<v.size();i++)
        out << i << "\t" << v[i] << endl;
    out.close();
    return 0;
}

#ifdef COUNTING
/*
 * Function: get_4D_table_counts
 * 
 * This function returns that parameter usages for the elements of a 4D data table matrix,
 * such as stack or coax.
 * 
 * ifname:  The name of the file in the data tables that the data matrix was filled from.
 *          The file is used to determine the order that the data table elements are 
 *          read when filling the vector of parameter usage counts.
 * 
 * matrix:  A reference to the data table variable in the datatable object.
 */
std::vector<double> datatable::get_4D_table_counts(char* ifname, shortVector4D &matrix) {
    //use old stack reading order
    //5'a-c3'
    //3'b-d5'
    //loop over the vector v_datatable_lines which holds the datatable lines.
    std::vector<double> data_counts;

    std::vector<std::string> v_datatable_lines;

    //read in lines from datatable file, return false if file does not exist
    read_datatable_lines(ifname, v_datatable_lines);

    //loop over the vector v_datatable_lines which holds the datatable lines.
    int a, b, c, d;

    for (int i=0; i<v_datatable_lines.size(); ) {
        //read the first headline(e.g.:"AX"), and get a.
        a = basetonum(v_datatable_lines[i][0]);
        i++;
        //read the second headline(e.g.:"AY"), and get c.
        b = basetonum(v_datatable_lines[i][0]);      
        i++;

        //declare a vector storing the bases (conveted to num) at Y position, 
        //using basetonum().
        std::vector <int> v;
        //read the third headline (e.g.:"Y" title line "A   C   G   U").
        std::istringstream yline(v_datatable_lines[i]);
        //loop over the yline items(bases), push_back the number of base to 
        //vector v. the v.size determines the next loop size for reading the data.
        char base_y;
        while (yline >> base_y) {
            v.push_back(basetonum(base_y));
        }
        i++;
        for (int x=0; x<v.size(); ++x) {
            char base_x;
            std::istringstream xline(v_datatable_lines[i]);           
            //read each xline's first entry which is the base_x, and get b.
            xline >> base_x; 
            c = basetonum(base_x);
            i++;             
            //read the rest of xline (data) into matrix of v.size()*v.size()
            for (int y=0; y<v.size(); ++y) {
                //get d
                d = v[y];
                std::string temp;
                xline >> temp;

                if (matrix[a][b][c][d]._value!=INFINITE_ENERGY)
                    data_counts.push_back(matrix[a][b][c][d].get);

            }
        }
    }
    return data_counts;
}

/*
 * Function: get_gata_dangle_counts
 * 
 * This function returns that parameter usages for the elements of dangling end matrix.
 * 
 * ifname:  The name of the file in the data tables that the data matrix was filled from.
 *          The file is used to determine the order that the data table elements are 
 *          read when filling the vector of parameter usage counts.
 * 
 * matrix:  A reference to the data table variable in the datatable object.
 */
std::vector<double> datatable::get_data_dangle_counts(char* ifname, shortVector4D &matrix) {
    // Declare the vector that will contain the usage counts
    std::vector<double> data_counts;

    //declare a vector that holds the lines from the datatable file.
    std::vector<std::string> v_datatable_lines;
    //read in lines from datatable file, return false if file does not exist
    read_datatable_lines(ifname, v_datatable_lines);
    
    //from a to d, "a" changes slowest, and "d" changes fastest.
    //5'a- c3'   
    //     d    d indicates the location of c, "1" is at 3' end, "2" is at 5' end.
    //3'b- c5'
    //loop over the vector v_datatable_lines which holds the datatable lines.
    int a, b, c, d;

    // Process the data table file lines
    for (int i=0; i<v_datatable_lines.size(); ) {
        //read the first headline(e.g.:"AX") and get a.
        a = basetonum(v_datatable_lines[i][0]);
        //check length of the first headline, and get d.
        //if it is "AX" format size=2: d=1; else: d=2., old d index starts from 1
        if (v_datatable_lines[i].size() == 2) {
			//also make sure the second character is X (not a /r or other whitespace)
			if (v_datatable_lines[i][1]=='X') {
				 d = 1;
			}
			else d = 2;
        }
        else {
            d = 2;
        }
        i++;
        //read the second headline(e.g.:"A"), and get b.
        b = basetonum(v_datatable_lines[i][0]);      
        i++;
        //declare a vector storing the bases (conveted to num) at X position, 
        //using basetonum().
        std::vector <int> v;
        //read the third head_line (e.g.: x line "A   C   G   U ").
        std::istringstream xline(v_datatable_lines[i]);
        i++;
        //loop over the xline items(bases), push_back the number of base to vecor v.
        char base_x;
        while (xline >> base_x) {
            v.push_back(basetonum(base_x));
        }
        //read the data line (e.g.: data line "-0.8 -0.5 -0.8 -0.6").
        std::istringstream data_line(v_datatable_lines[i]);
        i++;
        for (int x=0; x<v.size(); ++x) {
            //get c, the number of base X in vector v. 
            c = v[x];
            std::string temp;
            data_line >> temp;
            // If the energy is not infinite, add the usage count to the count vector
            if (matrix[a][b][c][d]._value != INFINITE_ENERGY)
                data_counts.push_back(matrix[a][b][c][d].get);
        }
    }
    // Return the usage counts
    return data_counts;
}

/*
 * Function: get_data_xloop_counts
 * 
 * This function returns that parameter usages for the elements of a data table variable
 * that contains the energies of specific hairpin loops (triloops, tloops, and hexaloops).
 * 
 * ifname:  The name of the file in the data tables that the data matrix was filled from.
 *          The file is used to determine the order that the data table elements are 
 *          read when filling the vector of parameter usage counts.
 * 
 * matrix:  A reference to the data table variable in the datatable object.
 */
std::vector<double> datatable::get_data_xloop_counts(char* ifname, intVector2D &matrix) {
    // Declare a vector to store the data usage counts
    std::vector<double> data_counts;

    std::string seq;
    std::string temp;
    //declare a vector that holds the lines from the datatable file.
    std::vector<std::string> v_datatable_lines;
    //read in lines from datatable file, return false if file does not exist
    read_datatable_lines(ifname, v_datatable_lines);

    //loop over the vector v_datatable_lines which holds the datatable lines.
    for (int i=0; i<v_datatable_lines.size(); i++) {
        std::istringstream line(v_datatable_lines[i]);
        //
        line >> seq >> temp;

        // If the energy is not infinite, stoe the usage count in the vector
        if (matrix[i][1]._value!=INFINITE_ENERGY)
            data_counts.push_back(matrix[i][1].get);
    }

    // Return the parameter usage counts
    return data_counts;
}

/*
 * Function: get_data_loop_counts
 * 
 * This function returns that parameter usages for the elements of the data table variables
 * that contains the loop initiation energies.
 * 
 * ifname:  The name of the file in the data tables that the data matrix was filled from.
 *          The file is used to determine the order that the data table elements are 
 *          read when filling the vector of parameter usage counts.
 * 
 * inter:   A reference to the internal loop table variable in the datatable object.
 * bulge:   A reference to the bulge loop table variable in the datatable object.
 * hairpin: A reference to the hairpin loop table variable in the datatable object.
 */
std::vector<double> datatable::get_data_loop_counts(char* ifname, std::vector<COUNT_CLASS(short int) > &inter,
                   std::vector<COUNT_CLASS(short int) > &bulge, std::vector<COUNT_CLASS(short int) > &hairpin) {
    // Declare the vector that will store the usage counts                       
    std::vector<double> data_counts;

    std::string temp0, temp1, temp2, temp3;
    //declare a vector that holds the lines from the datatable file.
    std::vector<std::string> v_datatable_lines;
    //read in lines from datatable file, return false if file does not exist
    read_datatable_lines(ifname, v_datatable_lines);
    
    //loop over the vector v_datatable_lines which holds the datatable lines.
    for (int i=0; i<v_datatable_lines.size(); i++) {
        std::istringstream line(v_datatable_lines[i]);
        //get thre values from each line, and send to temp1, temp2, temp3
        line >> temp0 >> temp1 >> temp2 >> temp3;
        
        // Add the usage counts to the count vector
        if (inter[i+1]._value != INFINITE_ENERGY)
            data_counts.push_back(inter[i+1].get);
        
        if (bulge[i+1]._value != INFINITE_ENERGY)
            data_counts.push_back(bulge[i+1].get);
        
        if (hairpin[i+1]._value != INFINITE_ENERGY)
            data_counts.push_back(hairpin[i+1].get);
    }

    // Return the parameter usage counts
    return data_counts;
}

/*
 * Function: get_6D_table_counts
 * 
 * This function returns that parameter usages for the elements of a 6D data table matrix,
 * i.e. 1x1 internal loop energy table.
 * 
 * ifname:  The name of the file in the data tables that the data matrix was filled from.
 *          The file is used to determine the order that the data table elements are 
 *          read when filling the vector of parameter usage counts.
 * 
 * matrix:  A reference to the data table variable in the datatable object.
 */
std::vector<double> datatable::get_6D_table_counts(char* ifname, shortVector6D &matrix) {
    std::vector<double> data_counts;

    //declare a vector that holds the lines from the datatable file.
    std::vector<std::string> v_datatable_lines;
    //read in lines from datatable file, return false if file does not exist
    read_datatable_lines(ifname, v_datatable_lines);

    //use old int11 reading order.
    //5'abc3'
    //3'def5'
    //loop over the vector v_datatable_lines which holds the datatable lines.
    int a, b, c, d, e, f;

    for (int i=0; i<v_datatable_lines.size(); ) {
        //skip the headline1 as 'X'
        i++;
        //read the first headline2(e.g.:"A A") of the datatable.
        std::istringstream headline2(v_datatable_lines[i]);
        char temp1; 
        char temp2;
        headline2 >> temp1 >> temp2;
        //get a and c from headline2
        a = basetonum(temp1);
        c = basetonum(temp2);
        i++;
        //read the second headline2(e.g.:"U U") of the datatable.
        std::istringstream headline3(v_datatable_lines[i]);
        headline3 >> temp1 >> temp2;
        ///get d and f from headline3
        d = basetonum(temp1); 
        f = basetonum(temp2);      
        //skip the headline4 as 'Y'
        i++;
        //declare a vector storing the number of Y position bases, 
        //with basetonum function.
        i++;
        std::vector <int> v;
        //read the third line (e.g.:"Y" title line) from v_datatable_lines
        std::istringstream yline(v_datatable_lines[i]);
        //loop over the yline items(bases), push_back the number of base to 
        //vector v. the v.size determines the next loop size for reading the data.
        char base_y;
        while (yline >> base_y) {
            v.push_back(basetonum(base_y));
        }
        i++;
        for (int x=0; x<v.size(); ++x) {
            char base_x;
            std::istringstream xline(v_datatable_lines[i]);           
            //read each xline's first entry which is the base_x
            xline >> base_x; 
            //get b from xline fist char, base_x
            b = basetonum(base_x);
            i++;             
            //read the data into matrix of v.size()*v.size()
            for (int y=0; y<v.size(); ++y) {
                //get e from base_y position.
                e = v[y];
                std::string temp;
                xline >> temp;
                //read data into matrix
                //floor_entry_to_short(const char* inum), convert a entry char number to a short int by 
                //conversionfactor 10.
                if (matrix[a][b][c][d][e][f]._value!=INFINITE_ENERGY)
                    data_counts.push_back(matrix[a][b][c][d][e][f].get);
            }
        }
    }
    return data_counts;
}

/*
 * Function: get_7D_table_counts
 * 
 * This function returns that parameter usages for the elements of a 7D data table matrix,
 * i.e. 2x1 internal loop energy table.
 * 
 * ifname:  The name of the file in the data tables that the data matrix was filled from.
 *          The file is used to determine the order that the data table elements are 
 *          read when filling the vector of parameter usage counts.
 * 
 * matrix:  A reference to the data table variable in the datatable object.
 */
std::vector<double> datatable::get_7D_table_counts(char* ifname, shortVector7D &matrix) {
    std::vector<double> data_counts;

    //declare a vector that holds the lines from the datatable file.
    std::vector<std::string> v_datatable_lines;
    //read in lines from datatable file, return false if file does not exist
    read_datatable_lines(ifname, v_datatable_lines);

    //use old int21 reading order.
    //5'ac f3'
    //3'bdeg5'
    //   
    //loop over the vector v_datatable_lines which holds the datatable lines.
    int a, b, c, d, e, f, g;

    for (int i=0; i<v_datatable_lines.size(); ) {
        //skip the healdline1 as 'X'
        i++;
        //read the first headline2(e.g.:"A A") of the datatable.
        std::istringstream headline2(v_datatable_lines[i]);
        char temp1; 
        char temp2;
        headline2 >> temp1 >> temp2;
        //get a and f from headline2
        a = basetonum(temp1);
        f = basetonum(temp2);
        i++;
        //read the second headline3(e.g.:"U U") of the datatable.
        std::istringstream headline3(v_datatable_lines[i]);
        headline3 >> temp1 >> temp2;
        //get b and g from headline3
        b = basetonum(temp1); 
        g = basetonum(temp2);      
        //read the headline4 as 'YA'
        i++;
        // get e from headline4
        e = basetonum(v_datatable_lines[i][1]);
        //declare a vector storing the number of Y position bases, 
        //with basetonum function.
        i++;
        std::vector <int> v;
        //read the third line (e.g.:"Y" title line) from v_datatable_lines
        std::istringstream yline(v_datatable_lines[i]);
        //loop over the yline items(bases), push_back the number of base to 
        //vector v. the v.size determines the next loop size for reading the data.
        char base_y;
        while (yline >> base_y) {
            //  << base << endl;
            v.push_back(basetonum(base_y));
        }
        i++;
        for (int x=0; x<v.size(); ++x) {
            char base_x;
            std::istringstream xline(v_datatable_lines[i]);           
            //read each xline's first entry which is the base_x
            xline >> base_x; 
            //get c from xline's first char, base_x
            c = basetonum(base_x);
            i++;  
                       
            //read the data into matrix of v.size()*v.size()
            for (int y=0; y<v.size(); ++y) {
                //get d from base_y location
                d = v[y];
                std::string temp;
                xline >> temp;              
                //read data into matrix
                //floor_entry_to_short(const char* inum), convert a entry char number to a short int by 
                //conversionfactor 10.
                if (matrix[a][b][c][d][e][f][g]._value != INFINITE_ENERGY)
                    data_counts.push_back(matrix[a][b][c][d][e][f][g].get);
            }
        }
    }
    
    return data_counts;
}

/*
 * Function: get_8D_table_counts
 * 
 * This function returns that parameter usages for the elements of a 8D data table matrix,
 * i.e. 2x2 internal loop energy table.
 * 
 * ifname:  The name of the file in the data tables that the data matrix was filled from.
 *          The file is used to determine the order that the data table elements are 
 *          read when filling the vector of parameter usage counts.
 * 
 * matrix:  A reference to the data table variable in the datatable object.
 */
std::vector<double> datatable::get_8D_table_counts(char* ifname, shortVector8D &matrix) {
    std::vector<double> data_counts;

    //declare a vector that holds the lines from the datatable file.
    std::vector<std::string> v_datatable_lines;
    //read in lines from datatable file, return false if file does not exist
    read_datatable_lines(ifname, v_datatable_lines);
            
    //use olf int22 reading order
    //5'a e f b3'
    //3'c g h d5'
    //loop over the vector v_datatable_lines which holds the datatable lines.
    int a, b, c, d, e, f, g, h;

    for (int i=0; i<v_datatable_lines.size(); ) {
        //read the first line(e.g.:"A X1 Y1 A") of the datatable.
        std::istringstream headline1(v_datatable_lines[i]);
        std::vector <std::string> v_header1;
        std::string temp1;
        while (headline1 >> temp1) {
             v_header1.push_back(temp1);
        }
        //get a and b based on headline1's first base and last base.
        char base_a = v_header1[0][0]; //get fist char from a string.
        char base_b = v_header1[3][0];
        //get a and b from 
        a = basetonum(base_a);
        b = basetonum(base_b);
        i++;
        //read the second line(e.g.:"AY") of the datatable.
        std::istringstream headline2(v_datatable_lines[i]);
        std::vector <std::string> v_header2;
        std::string temp2;
        while (headline2 >> temp2) {
            v_header2.push_back(temp2);
        }
        //get c and d based on headline2's first base and last base.
        char base_c = v_header2[0][0];
        char base_d = v_header2[3][0];
        c = basetonum(base_c);
        d = basetonum(base_d);    

        i++;
        //declare a vector storing the number of Y1 position bases, 
        //with basetonum function.
        std::vector <int> v_y1;
        //read the third line (e.g.:"Y1" line) from v_datatable_lines
        std::istringstream yline1(v_datatable_lines[i]);
        //loop over the yline items(bases), push_back the number of base to 
        //vector v. the v.size determines the next loop size for reading the data.
        char base_y1;
        while (yline1 >> base_y1) {
            v_y1.push_back(basetonum(base_y1));
        }
        i++;
        //declare a vector storing the number of Y2 position bases, 
        //with basetonum function.
        std::vector <int> v_y2;
        //read the fourth line (e.g.:"Y2" line) from v_datatable_lines
        std::istringstream yline2(v_datatable_lines[i]);
        //loop over the yline items(bases), push_back the number of base to 
        //vector v. the v.size determines the next loop size for reading the data.
        char base_y2;
        while (yline2 >> base_y2) {
            v_y2.push_back(basetonum(base_y2));
        }

        i++;
        for (int x=0; x<v_y1.size(); ++x) {
            std::string base_x;
            char base_x1;
            char base_x2;
            std::istringstream xline(v_datatable_lines[i]);           
            //read each xline's first entry which is the base_x
            xline >> base_x; 
            //get e and g from xline first and scond bases
            base_x1 = base_x[0];
            base_x2 = base_x[1];
            e = basetonum(base_x1);
            g = basetonum(base_x2);
            i++;             
            //read the data into matrix of v.size()*v.size()
            for (int y=0; y<v_y1.size(); ++y) {
                //get f and h from bases of first and second yline 
                f = v_y1[y];
                h = v_y2[y];
                std::string temp;
                xline >> temp;
                //read data into matrix
                //floor_entry_to_short(const char* inum), convert a entry char number to a short int by 
                //conversionfactor 10.
                if (matrix[a][b][c][d][e][f][g][h]._value != INFINITE_ENERGY)
                    data_counts.push_back(matrix[a][b][c][d][e][f][g][h].get);  
            }
        }
    }

    return data_counts;
}

/*
 * Function: get_miscloop_counts
 * 
 * This function returns that parameter usages for the datatable variables whose values are stored in the 
 * miscloop data file. 
 * 
 * ifname:  The name of the file in the data tables that the data matrix was filled from.
 *          The file is used to determine the order that the data table elements are 
 *          read when filling the vector of parameter usage counts.
 * 
 * The other arguments are references to the datatable object variables that are defined in miscloop
 */
std::vector<double> datatable::get_miscloop_counts(char* ifname, COUNT_CLASS(float) &prelog, COUNT_CLASS(short int) &maxpen,
                   COUNT_CLASS(short int) &efn2a, COUNT_CLASS(short int) &efn2b, COUNT_CLASS(short int) &efn2c, 
                   COUNT_CLASS(short int) &strain, COUNT_CLASS(short int) &auend, COUNT_CLASS(short int) &gubonus, 
                   COUNT_CLASS(short int) &cslope, COUNT_CLASS(short int) &cint, COUNT_CLASS(short int) &c3, COUNT_CLASS(short int) &init,
                   COUNT_CLASS(short int) &singlecbulge, std::vector<COUNT_CLASS(short int) > &poppen,
                   std::vector<COUNT_CLASS(short int) > &eparam){    

    std::vector<double> data_counts;

    std::string temp;

    //pass a pointer to the vector to the function
    //declare a vector that holds the lines from the datatable file.
    std::vector<std::string> v_datatable_lines;
    //read in lines from datatable file, return false if file does not exist
    read_datatable_lines(ifname, v_datatable_lines);

    //loop over the vector v_datatable_lines which holds the datatable lines.
    for (int i=0; i<v_datatable_lines.size(); ) {

        //read the prelog.
        // temp = atof(v_datatable_lines[i].c_str());
        // prelog = temp * conversionfactor;
        temp = v_datatable_lines[i];
        data_counts.push_back(prelog.get);
        i++;
    
        //read maxpen
        temp = v_datatable_lines[i];
        // maxpen = (short) (temp * CONVERSION_FACTOR + .5);
        // maxpen = floor_to_short(temp);
        data_counts.push_back(maxpen.get);
        i++;

        //read poppen
        std::istringstream poppen_line(v_datatable_lines[i]);
        i++;
        for (int k=1; k<=4; k++){
            poppen_line >> temp;
            // poppen[k] = (short) (temp * CONVERSION_FACTOR + .5);
            data_counts.push_back(poppen[k].get);
        }

        //read eparam
        std::istringstream eparam_line(v_datatable_lines[i]);
        i++;
        eparam_line >> temp;
        // eparam[5] = (short) (floor(temp * CONVERSION_FACTOR + .5));  //constant multi-loop penalty
        data_counts.push_back(eparam[5].get);
        eparam_line >> temp;
        // eparam[6] = (short) (floor(temp * CONVERSION_FACTOR + .5));  //constant multi-loop penalty
        data_counts.push_back(eparam[6].get);
        eparam_line >> temp;
        // eparam[10] = (short) (floor(temp * CONVERSION_FACTOR + .5));
        data_counts.push_back(eparam[10].get);
        
        //read efn
        std::istringstream efn_line(v_datatable_lines[i]);
        i++;
        
        efn_line >> temp;
        data_counts.push_back(efn2a.get);
        
        efn_line >> temp;
        data_counts.push_back(efn2b.get);
        
        efn_line >> temp;
        data_counts.push_back(efn2c.get);
        
        //read mlasym
        temp = v_datatable_lines[i];
        i++;
        data_counts.push_back(mlasym.get);
        
        //read strain
        temp = v_datatable_lines[i];
        i++;
        data_counts.push_back(strain.get);
        
        //read the terminal AU penalty
        temp = v_datatable_lines[i];
        i++;
        data_counts.push_back(auend.get);

		//read whether the terminal AU penalty applies to GU, but do not pushback -- no need to count
		temp = v_datatable_lines[i];
		i++;
        
        //read the GGG hairpin bonus
        temp = v_datatable_lines[i];
        i++;
        data_counts.push_back(gubonus.get);
        
        //read the poly c hairpin penalty slope
        temp = v_datatable_lines[i];
        i++;
        data_counts.push_back(cslope.get);
        
        //read the poly c hairpin penalty intercept
        temp = v_datatable_lines[i];
        i++;
        data_counts.push_back(cint.get);
        
        //read the poly c penalty for a loop of 3
        temp = v_datatable_lines[i];
        i++;
        data_counts.push_back(c3.get);
        
        //read int
        temp = v_datatable_lines[i];
        i++;
        data_counts.push_back(init.get);
        
        //read the single C bulge bonus
        temp = v_datatable_lines[i];
        i++;
        data_counts.push_back(singlecbulge.get);
    }
    
    return data_counts;
}

// Writes a vector of parameter usage count to disk.  Might not be needed anymore
bool datatable::write_count_vector(ofstream &ofile, std::vector<double> counts){
    for (int i=0; i<counts.size(); i++){
        ofile << setprecision(outprecision) << fixed << counts[i] << endl;
    }
    return 1;
}

// Writes a file of all parameter usage counts to disk.  Might not be needed anymore
bool datatable::export_data_counts(string CountFile) {
    ofstream ofile;

    ofile.open(CountFile.c_str());

    datfile dat(getDataPath(),"rna",false);

    //01.rna.stack
    ofile <<"#01.rna.stack.dg" << endl;
    write_count_vector(ofile, get_4D_table_counts(dat.get(F_stack), stack));

    //02.rna.dangle
    ofile << endl <<  "#02.rna.dangle.dg" << endl;
    write_count_vector(ofile, get_data_dangle_counts(dat.get(F_dangle), dangle));
    
    //03.rna.tstack
    ofile << endl <<  "#03.rna.tstack.dg" << endl;
    write_count_vector(ofile, get_4D_table_counts(dat.get(F_tstack), tstack));
    
    //04.rna.tstackcoax
    ofile << endl <<  "#04.rna.tstackcoax.dg" << endl;
    write_count_vector(ofile, get_4D_table_counts(dat.get(F_tstackcoax), tstackcoax));
    
    //05.rna.tstackm
    ofile << endl <<  "#05.rna.tstackm.dg" << endl;
    write_count_vector(ofile, get_4D_table_counts(dat.get(F_tstackm), tstkm));
    
    //06.rna.coaxial
    ofile << endl <<  "#06.rna.coaxial.dg" << endl;
    write_count_vector(ofile, get_4D_table_counts(dat.get(F_coax), coax));
    
    //07.rna.miscloop
    ofile << endl <<  "#07.rna.miscloop.dg" << endl;
    write_count_vector(ofile, get_miscloop_counts(dat.get(F_miscloop), prelog, maxpen, efn2a, efn2b, efn2c, strain,
						auend, gubonus, cslope, cint, c3, init, singlecbulge, poppen,
						eparam));

    //08.rna.loop
    ofile << endl <<  "#08.rna.loop.dg" << endl;
    write_count_vector(ofile, get_data_loop_counts(dat.get(F_loop), inter, bulge, hairpin));
    
    //09.rna.triloop
    ofile << endl <<  "#09.rna.triloop.dg" << endl;
    write_count_vector(ofile, get_data_xloop_counts(dat.get(F_triloop), triloop));
       
    //10.rna.tloop
    ofile << endl <<  "#10.rna.tloop.dg" << endl;
    write_count_vector(ofile, get_data_xloop_counts(dat.get(F_tloop), tloop));

    //11.rna.hexaloop
    ofile << endl <<  "#11.rna.hexaloop.dg" << endl;
    write_count_vector(ofile, get_data_xloop_counts(dat.get(F_hexaloop), hexaloop));
    
    //12.rna.tstackh
    ofile << endl <<  "#12.rna.tstackh.dg" << endl;
    write_count_vector(ofile, get_4D_table_counts(dat.get(F_tstackh), tstkh));
    
    //13.rna.tstacki
    ofile << endl <<  "#13.rna.tstacki.dg" << endl;
    write_count_vector(ofile, get_4D_table_counts(dat.get(F_tstacki), tstki));
    
    //14.rna.tstacki1n
    ofile << endl <<  "#14.rna.tstacki1n.dg" << endl;
    write_count_vector(ofile, get_4D_table_counts(dat.get(F_tstacki1n), tstki1n));
    
    //15.rna.tstacki23
    ofile << endl <<  "#15.rna.tstacki23.dg" << endl;
    write_count_vector(ofile, get_4D_table_counts(dat.get(F_tstacki23), tstki23));

    //16.rna.coaxstack
    ofile << endl <<  "#16.rna.coaxstack.dg" << endl;
    write_count_vector(ofile, get_4D_table_counts(dat.get(F_coaxstack), coaxstack));
    
    //17.rna.int11
    ofile << endl <<  "#17.rna.int11.dg" << endl;
    write_count_vector(ofile, get_6D_table_counts(dat.get(F_int11), iloop11));
    
    //18.rna.int21
    ofile << endl <<  "#18.rna.int21.dg" << endl;
    write_count_vector(ofile, get_7D_table_counts(dat.get(F_int21), iloop21));

    //19.rna.int22
    ofile << endl <<  "#19.rna.int22.dg" << endl;
    write_count_vector(ofile, get_8D_table_counts(dat.get(F_int22), iloop22));
    
    ofile.close();

    return true;
}

// Returns a vector containing all of the parameter usage counts
std::vector<double> datatable::get_data_counts() {
    if (!loadedCov){
        // Read in the variance-covariance matrix and parameter map
        datfile dat(getDataPath(),"rna",false);
        read_cov(dat.get(F_cov), cov);
        read_param_map(dat.get(F_pmap), param_map);

        // Read in enthalpy parameter covariance if the temperature is not 37C
        if (false){
            datfile dat_h(getDataPath(),"rna",true);
            read_cov(dat_h.get(F_cov), cov_h);
        }

        // Flip the flag
        loadedCov = true;
    }

    // Declare a vector to store all of the parameter usage counts
    std::vector<double> data_counts;

    // Declare a vector that will temporarily store the usage counts for a subset of the parameters before being added to data_counts
    std::vector<double> data_sub_counts;

    // Define the datfile directory.  Currently, only rna parameters at 37 degrees are supported
    datfile dat(getDataPath(),"rna",false);

    // Get parameter usage counts for parameters stored in rna.stack
    data_counts = get_4D_table_counts(dat.get(F_stack), stack);

    // Get parameter usage counts for parameters stored in rna.dangle
    data_sub_counts = get_data_dangle_counts(dat.get(F_dangle), dangle);
    // Append the data_sub_counts to data_counts
    data_counts.insert( data_counts.end(), data_sub_counts.begin(), data_sub_counts.end() );

    // Get parameter usage counts for parameters stored in rna.tstack
    data_sub_counts = get_4D_table_counts(dat.get(F_tstack), tstack);
    // Append the data_sub_counts to data_counts
    data_counts.insert( data_counts.end(), data_sub_counts.begin(), data_sub_counts.end() );
    
    // Get parameter usage counts for parameters stored in rna.tstackcoax
    data_sub_counts = get_4D_table_counts(dat.get(F_tstackcoax), tstackcoax);
    // Append the data_sub_counts to data_counts
    data_counts.insert( data_counts.end(), data_sub_counts.begin(), data_sub_counts.end() );
    
    // Get parameter usage counts for parameters stored in rna.tstackm
    data_sub_counts = get_4D_table_counts(dat.get(F_tstackm), tstkm);
    // Append the data_sub_counts to data_counts
    data_counts.insert( data_counts.end(), data_sub_counts.begin(), data_sub_counts.end() );
    
    // Get parameter usage counts for parameters stored in rna.coaxial
    data_sub_counts = get_4D_table_counts(dat.get(F_coax), coax);
    // Append the data_sub_counts to data_counts
    data_counts.insert( data_counts.end(), data_sub_counts.begin(), data_sub_counts.end() );
    
    // Get parameter usage counts for parameters stored in rna.miscloop
    data_sub_counts = get_miscloop_counts(dat.get(F_miscloop), prelog, maxpen, efn2a, efn2b, efn2c, strain,
						auend, gubonus, cslope, cint, c3, init, singlecbulge, poppen,
						eparam);
    // Append the data_sub_counts to data_counts
    data_counts.insert( data_counts.end(), data_sub_counts.begin(), data_sub_counts.end() );
    
    // Get parameter usage counts for parameters stored in rna.loop
    data_sub_counts = get_data_loop_counts(dat.get(F_loop), inter, bulge, hairpin);
    // Append the data_sub_counts to data_counts
    data_counts.insert( data_counts.end(), data_sub_counts.begin(), data_sub_counts.end() );
    
    // Get parameter usage counts for parameters stored in rna.triloop
    data_sub_counts = get_data_xloop_counts(dat.get(F_triloop), triloop);
    // Append the data_sub_counts to data_counts   
    data_counts.insert( data_counts.end(), data_sub_counts.begin(), data_sub_counts.end() );
    
    // Get parameter usage counts for parameters stored in rna.tloop
    data_sub_counts = get_data_xloop_counts(dat.get(F_tloop), tloop);
    // Append the data_sub_counts to data_counts
    data_counts.insert( data_counts.end(), data_sub_counts.begin(), data_sub_counts.end() );
    
    // Get parameter usage counts for parameters stored in rna.hexaloop
    data_sub_counts = get_data_xloop_counts(dat.get(F_hexaloop), hexaloop);
    // Append the data_sub_counts to data_counts
    data_counts.insert( data_counts.end(), data_sub_counts.begin(), data_sub_counts.end() );
    
    // Get parameter usage counts for parameters stored in rna.tstackh
    data_sub_counts = get_4D_table_counts(dat.get(F_tstackh), tstkh);
    // Append the data_sub_counts to data_counts
    data_counts.insert( data_counts.end(), data_sub_counts.begin(), data_sub_counts.end() );
    
    // Get parameter usage counts for parameters stored in rna.tstacki
    data_sub_counts = get_4D_table_counts(dat.get(F_tstacki), tstki);
    // Append the data_sub_counts to data_counts
    data_counts.insert( data_counts.end(), data_sub_counts.begin(), data_sub_counts.end() );
    
    // Get parameter usage counts for parameters stored in rna.tstacki1n
    data_sub_counts = get_4D_table_counts(dat.get(F_tstacki1n), tstki1n);
    // Append the data_sub_counts to data_counts
    data_counts.insert( data_counts.end(), data_sub_counts.begin(), data_sub_counts.end() );
    
    // Get parameter usage counts for parameters stored in rna.tstacki23
    data_sub_counts = get_4D_table_counts(dat.get(F_tstacki23), tstki23);
    // Append the data_sub_counts to data_counts
    data_counts.insert( data_counts.end(), data_sub_counts.begin(), data_sub_counts.end() );
    
    // Get parameter usage counts for parameters stored in rna.coaxstack
    data_sub_counts = get_4D_table_counts(dat.get(F_coaxstack), coaxstack);
    // Append the data_sub_counts to data_counts
    data_counts.insert( data_counts.end(), data_sub_counts.begin(), data_sub_counts.end() );
    
    // Get parameter usage counts for parameters stored in rna.int11
    data_sub_counts = get_6D_table_counts(dat.get(F_int11), iloop11);
    // Append the data_sub_counts to data_counts
    data_counts.insert( data_counts.end(), data_sub_counts.begin(), data_sub_counts.end() );
    
    // Get parameter usage counts for parameters stored in rna.int21
    data_sub_counts = get_7D_table_counts(dat.get(F_int21), iloop21);
    // Append the data_sub_counts to data_counts
    data_counts.insert( data_counts.end(), data_sub_counts.begin(), data_sub_counts.end() );
    
    // Get parameter usage counts for parameters stored in rna.int22
    data_sub_counts = get_8D_table_counts(dat.get(F_int22), iloop22);
    // Append the data_sub_counts to data_counts
    data_counts.insert( data_counts.end(), data_sub_counts.begin(), data_sub_counts.end() );
    
    std::vector<double> tmp1;

    // Return the final vector of parameter usage counts
    for( int i =0; i< param_map.size(); i++){
        tmp1.push_back(param_map[i]*data_counts);
    }
    return tmp1;
}

// Clears the parameter usage data for a 4D datatable variable
void datatable::clear_4D_usage_data(shortVector4D &matrix){
    for(int i=0; i<alphabet.size(); i++){
        for(int j=0; j<alphabet.size(); j++){
            for(int k=0; k<alphabet.size(); k++){
                for(int l=0; l<alphabet.size(); l++){
                    matrix[i][j][k][l].get = 0;
                }
            }
        }
    }
}

// Clears the usage data for loop initiations
void datatable::clear_loop_usage_data(std::vector<COUNT_CLASS(short int) > &matrix){
    for( int i=0; i<matrix.size(); i++)
        matrix[i].get = 0;
}

// Clears the usage data for specific hairpin loops
void datatable::clear_xloop_usage_data(intVector2D &matrix){
    for( int i=0; i<matrix.size(); i++)
        matrix[i][1].get = 0;
}

// Clears the parameter usage data for all parameters
void datatable::clear_parameter_usage_data(){

    // Clear 4D data tables
    clear_4D_usage_data(stack);
    clear_4D_usage_data(tstack);
    clear_4D_usage_data(tstackcoax);
    clear_4D_usage_data(tstkm);
    clear_4D_usage_data(coax);
    clear_4D_usage_data(tstkh);
    clear_4D_usage_data(tstki);
    clear_4D_usage_data(tstki1n);
    clear_4D_usage_data(tstki23);
    clear_4D_usage_data(coaxstack);

    // Reset dangling end usage counts
    for(int i=0; i<alphabet.size(); i++){
        for(int j=0; j<alphabet.size(); j++){
            for(int k=0; k<alphabet.size(); k++){
                dangle[i][j][k][1].get = 0;
                dangle[i][j][k][2].get = 0;
            }
        }
    }  
    
    //07.rna.miscloop
        prelog.get = 0;

    maxpen.get = 0;
    
    for (int k=1; k<=4; k++){
        poppen[k].get = 0;
    }

    eparam[5].get = 0;
    eparam[6].get = 0;
    eparam[10].get = 0;
    
    efn2a.get = 0;
    efn2b.get = 0;
    efn2c.get = 0;
    
    mlasym.get = 0;
    strain.get = 0;
    auend.get = 0;
    gubonus.get = 0;
    cslope.get = 0;
    cint.get = 0;
    c3.get = 0;
    init.get = 0;
    singlecbulge.get = 0;

    //08.rna.loop
    clear_loop_usage_data(inter);
    clear_loop_usage_data(bulge);
    clear_loop_usage_data(hairpin);
    
    //09.rna.triloop
    clear_xloop_usage_data(triloop);
       
    //10.rna.tloop
    clear_xloop_usage_data(tloop);

    //11.rna.hexaloop
    clear_xloop_usage_data(hexaloop);

    //17.rna.int11
    for(int i=0; i<alphabet.size(); i++){
        for(int j=0; j<alphabet.size(); j++){
            clear_4D_usage_data(iloop11[i][j]);
            // for(int k=0; k<alphabet.size(); k++){
            //     for(int l=0; l<alphabet.size(); l++){
            //         for(int m=0; m<alphabet.size(); m++){
            //             for(int n=0; n<alphabet.size(); n++){
            //                 iloop11[i][j][k][l][m][n].get = 0;
            //             }
            //         }
            //     }
            // }
        }
    }
    
    //18.rna.int21
    for(int i=0; i<alphabet.size(); i++){
        for(int j=0; j<alphabet.size(); j++){
            for(int k=0; k<alphabet.size(); k++){
                clear_4D_usage_data(iloop21[i][j][k]);
            //     for(int l=0; l<alphabet.size(); l++){
            //         for(int m=0; m<alphabet.size(); m++){
            //             for(int n=0; n<alphabet.size(); n++){
            //                 iloop11[i][j][k][l][m][n].get = 0;
            //             }
            //         }
            //     }
            }
        }
    }    

    //19.rna.int22
    for(int i=0; i<alphabet.size(); i++){
        for(int j=0; j<alphabet.size(); j++){
            for(int k=0; k<alphabet.size(); k++){
               for(int l=0; l<alphabet.size(); l++){
                   clear_4D_usage_data(iloop22[i][j][k][l]);
            //         for(int m=0; m<alphabet.size(); m++){
            //             for(int n=0; n<alphabet.size(); n++){
            //                 iloop11[i][j][k][l][m][n].get = 0;
            //             }
            //         }
               }
            }
        }
    }
}

/*
 * Reads in the variance-covariance matrix from the DATAPATH directory
 */
bool datatable::read_cov(char* ifname, std::vector<std::vector<double> > &cov){
    //declare a vector that holds the lines from the datatable file.
    std::vector<std::string> v_datatable_lines;
    //read in lines from datatable file, return false if file does not exist
    if(!read_datatable_lines(ifname, v_datatable_lines)){
        return false;
    }

    //loop over the vector v_datatable_lines which holds the datatable lines.
    for (int i=0; i<v_datatable_lines.size(); i++) {
        std::istringstream line(v_datatable_lines[i]);
        
        std::vector <double> v;
        std::string cov_xy;

        while (line >> cov_xy) {
            v.push_back(atof(cov_xy.c_str()));
        }

        cov.push_back(v);
    }
    return true;
}

/*
 * Reads in the parameter map matrix from the DATAPATH directory.  This matrix is used to 
 * convert the data table parameter usage counts to usage counts of the free parameters that
 * compose the nearest neighbor energy model.
 */
bool datatable::read_param_map(char* ifname, std::vector<std::vector<double> > &param_map){
    //declare a vector that holds the lines from the datatable file.
    std::vector<std::string> v_datatable_lines;
    //read in lines from datatable file, return false if file does not exist
    if(!read_datatable_lines(ifname, v_datatable_lines)){
        return false;
    }

    //loop over the vector v_datatable_lines which holds the datatable lines.
    for (int i=0; i<v_datatable_lines.size(); i++) {
        std::istringstream line(v_datatable_lines[i]);
        std::string map_xy;
        std::vector <double> v;

        while (line >> map_xy) {
            v.push_back(atof(map_xy.c_str()));
        }

        param_map.push_back(v);
    }
    return true;
}

/*
 * Calculates an experimental uncertainty based on the parameter usage counts.
 */
double datatable::calculate_uncertainty(){
    // Check to see if the variance-covariance matrix and parameter map have been imported from disk
    if (!loadedCov){
        // Read in the variance-covariance matrix and parameter map
        datfile dat(getDataPath(),"rna",false);
        read_cov(dat.get(F_cov), cov);
        read_param_map(dat.get(F_pmap), param_map);

        // Read in enthalpy parameter covariance if the temperature is not 37C
        if (temperature != TEMP_37C){
            datfile dat_h(getDataPath(),"rna",true);
            read_cov(dat_h.get(F_cov), cov_h);
        }

        // Flip the flag
        loadedCov = true;
    }

    // Declare the data counts vector and two vectors that will store intermediate calculations
    std::vector<double> data_counts, tmp1, tmp2;

    // Get all of the parameter usage counts
    tmp1 = get_data_counts();

    // Convert the parameters counts to free parameter counts
    // for( int i =0; i< param_map.size(); i++){
    //     tmp1.push_back(param_map[i]*data_counts);
    // }

    // Calculate product of the variance-covariance matrix and the free parameter usage counts vector
    for( int i =0; i< cov.size(); i++){
        tmp2.push_back(cov[i]*tmp1);
    }

    // Calculate and return the uncertainty
    if (temperature == TEMP_37C){
        return sqrt(tmp1*tmp2);
    }
    else{
        // Variance in delta G at 37C
        double dG_err = tmp1*tmp2;
        
        std::vector<double> tmp3;

        for( int i =0; i< cov.size(); i++){
            tmp3.push_back(cov_h[i]*tmp1);
        }

        // Variance in delta H
        double dH_err = tmp1*tmp3;
        // Variance in delta S
        double dS_err = pow(((sqrt(dH_err)*DH_DS_CORR-sqrt(dH_err*pow(DH_DS_CORR,2)-dH_err+dG_err))/TEMP_37C),2);

        // Compute the variance in delta G at the new temperature (propagating error through dG=dH-T*dS)
        return sqrt(dH_err + pow(temperature,2)*dS_err-2*temperature*sqrt(dH_err*dS_err)*DH_DS_CORR);
    }
}
#endif //COUNTING

#undef NO_COUNT
