#include "NCM_parameters.h"
#include "io.h"
#include "zero.h"
#include <string>
#include <cstdlib>
#include <iostream>
#include <climits>
#include <cmath>
#include <cassert>
#include <type_traits>
using std::cout;
using std::string;
using std::endl;
using IO::read_datafile;
const int PRECISION = 4;

const double XIA_INTERCEPT = 6;//4.764;
const double XIA_SLOPE = -6;
const double SLOPE_SUM = -2.271-1.635;

const double RNASTRAND_SLOPE = -0.22;

bool is_hairpin_NCM(const NCM_type& n){
  return n[0] == '1';
}

NCM_id::NCM_id(NCM_type n,std::string sequence)
	: ncm(n),five(sequence),three(sequence),id(ncm+sequence)
{
//    assert(valid_NCM());
}


NCM_id::NCM_id(const int i,const int j,NCM_type n,const std::string& RNA)
	: 	ncm(n),
		five(calcFive(i,j,n,RNA)),
	    three(calcThree(i,j,n,RNA)),
		id(ncm+five+three)
{
//    assert(valid_NCM());
}

bool NCM_id::valid_NCM(){
    return(five+three).length()==(size_t)(NCM::fivep_length(ncm)+NCM::threep_length(ncm));
}


//I need a new constructor to represent a coaxial stack
std::string NCM_id::calcFive(const int i,const int j,NCM_type n,const std::string& RNA){
    if (j>i){
        return RNA.substr(i,NCM::fivep_length(n));
    }
    else {
        return RNA.substr(j-NCM::fivep_length(n)+1,
                NCM::fivep_length(n));
    }
}

std::string NCM_id::calcThree(const int i,const int j,NCM_type n,const std::string& RNA){
    if (j>i){
        return RNA.substr(j-NCM::threep_length(n)+1,
                NCM::threep_length(n));
    }
    else {
        return RNA.substr(i,NCM::threep_length(n));
    }
}

std::string NCM_id::getID() const {
	return id;
}

bool NCM_id::operator==(const NCM_id& rhs) const {
    return getID() == rhs.getID();
}

std::string NCM_id::getType() const {
	return ncm;
}

std::string NCM_id::getOuterNucs() const {
	char nucs[3];
	nucs[0] = five.front();
	nucs[1] = three.back();
	nucs[2] = '\0';
	return std::string(nucs);
}

std::string NCM_id::getInnerNucs() const {
    assert(ncm[0] != '1');//a hairpin has no inner nucleotides
	char nucs[3];
	nucs[0] = five.back();
	nucs[1] = three.front();
	nucs[2] = '\0';
	return std::string(nucs);
}

connect_id::connect_id(NCM_type i,NCM_type o, nucs n)
	:id(i+o+n),inner(i),outer(o),_nucs(n){
	}

std::string connect_id::getID() const {
	return id;
}

bool connect_id::operator==(const connect_id& rhs) const {
    return getID() == rhs.getID();
}

junction_id::junction_id(NCM_type i,NCM_type o)
	:id(i+o),inner(i),outer(o){
	}

std::string junction_id::getID() const {
	return id;
}

bool junction_id::operator==(const junction_id& rhs) const {
    return getID() == rhs.getID();
}

hinge_id::hinge_id(NCM_type i, NCM_type o, pairtype p)
	:id(i+o+p),inner(i),outer(o),pair(p){}

std::string hinge_id::getID() const {
	return id;
}

bool hinge_id::operator==(const hinge_id& rhs) const {
    return getID() == rhs.getID();
}

pair_id::pair_id(pairtype p, nucs n)
	:id(p+n),pair(p),seq(n)
	{}

std::string pair_id::getID() const {
	return id;
}

template<typename T>
parameters<T>::parameters(const options o)
: 	op(o),
	seqs(read_seqs()),
	junctions(read_junctions()),
	hinges(read_hinges()),
	pairs(read_pairs()),
	aalberts(read_aalberts()),
	connect(read_connect()),
	normalization(read_normalization()),
    mfe(std::is_integral<T>::value)
{
	/*
	for(auto entry : seqs){
		cout<< entry.first.getID()<<" "<<entry.second<<endl;
	}
	for(auto entry : junctions){
		cout<< entry.first.getID()<<" "<<entry.second<<endl;
	}
	for(auto entry : hinges){
		cout<< entry.first.getID()<<" "<<entry.second<<endl;
	}
	for(auto entry : pairs){
		cout<< entry.first.getID()<<" "<<entry.second<<endl;
	}
	*/
}

bool pair_id::operator==(const pair_id& rhs) const {
    return getID() == rhs.getID();
}

//given two NCMs, n which is being added and m which is being added onto
//give the nucleotides in the pair
std::string pn(NCM_id n,NCM_id m,bool exterior){
    if(!exterior){
        return m.getOuterNucs();
    }
    else return n.getOuterNucs();
}

std::string pn(NCM_id n, bool exterior){
    if(!exterior){
        return n.getOuterNucs();
    }
    else return n.getInnerNucs();
}

template<typename T>
T parameters<T>::energy(NCM_id theta,NCM_id phi,bool exterior, bool debug) const
{//theta is the one that's getting added on
//at some point this needs to be switched around to be consistent
	T NCM = seqs.at(theta);
	/*
	T junc = junctions.at(junction_id(phi.getType(), theta.getType()));
	T hinge = hinges.at(hinge_id(theta.getType(), phi.getType(),p));
	*/
  //string closing_nucs = exterior? theta.getInnerNucs() : theta.getOuterNucs();
  string closing_nucs = theta.getInnerNucs(); 
  //T pair = pairs.at(pair_id(string(""),closing_nucs));
	//T conn = connect.at(connect_id(theta.getType(),phi.getType(),theta.getInnerNucs()));//pn(theta,phi,exterior)));
	T conn = connect.at(connect_id(theta.getType(),phi.getType(),pn(theta,phi,exterior)));
    T e;
    if(std::is_integral<T>::value){//mfe calculation
		if(debug){
			cout<< "added NCM|seqs "<< theta.getID()<< " "<<NCM<<endl;
		}
        e = NCM + conn;
    }
    else {
        e = NCM * conn;//partition function
    }
    return e;
}

template<typename T>
T parameters<T>::junction_energy(NCM_id theta, NCM_id phi, bool exterior) const
{
	/*
	T junc = junctions.at(junction_id(phi.getType(), theta.getType()));
	T hinge = hinges.at(hinge_id(theta.getType(), phi.getType(),p));
	*/
	T conn = connect.at(connect_id(theta.getType(),phi.getType(),pn(theta,phi,exterior)));
	//T conn = connect.at(connect_id(theta.getType(),phi.getType(),theta.getInnerNucs()));
    T e;
    if(std::is_integral<T>::value){//mfe calculation
        e = conn;
    }
    else {
        e = conn;//partition function
    }
    return e;
}

/*
template<typename T>
T parameters<T>::energy(NCM_id theta,NCM_type phi,pairtype p,bool exterior) const
{//theta is the one that's getting added on
//at some point this needs to be switched around to be consistent
	T NCM = seqs.at(theta);
	T junc = junctions.at(junction_id(phi, theta.getType()));
	T hinge = hinges.at(hinge_id(theta.getType(), phi,p));
	T pair = pairs.at(pair_id(p,pn(theta,exterior)));
    T e;
    if(std::is_integral<T>::value){
		const T NORMALIZATION = (T) (4.764 * (pow(10,PRECISION)));
        e = NCM + junc + hinge + pair - NORMALIZATION;
    }
    else {
        e = NCM*junc*hinge*pair;
    }
    return e;
}
*/
template<typename T>
T parameters<T>::energy(NCM_id theta, bool debug) const
{//theta is a single-stranded NCM (a hairpin loop)
    if(std::is_integral<T>::value){
		if(debug){
			cout << "added ncm|seqs "<<seqs.at(theta)<<endl;
		}
    //T pair = pairs.at(pair_id("",theta.getOuterNucs())) +
    //  (is_hairpin_NCM(theta.getType())? T(0) : pairs.at(pair_id("",theta.getInnerNucs())));
		//return seqs.at(theta)+pair;//+aalberts_entropy(theta.getType());//-NORMALIZATION;
		return seqs.at(theta);
	}
	else {
    //T pair = pairs.at(pair_id("",theta.getOuterNucs())) *
    //  (is_hairpin_NCM(theta.getType())? T(1.0) : pairs.at(pair_id("",theta.getInnerNucs())));
		//return seqs.at(theta)*pair;
		return seqs.at(theta);
	}
}

template<typename T>
T parameters<T>::pair_energy(NCM_id theta, bool debug) const
{
  return pairs.at(pair_id("",theta.getOuterNucs()));
}


template<typename T>
T parameters<T>::pair_bonus() const{
	return this->normalization;
}

template<typename T>
T parameters<T>::mb_nuc_penalty() const{
	return 1.0;
}

template<typename T>
T parameters<T>::mb_stem_penalty() const{
	return -0.0;
}

template<typename T>
T parameters<T>::mb_closure_penalty() const{
	return -0.0;
}

template<typename T>
bool parameters<T>::mfe_calc() const {
    return mfe;
}

template<typename T>
T parameters<T>::aalberts_entropy(NCM_type n) const {
    return aalberts.at(n);
}



int NCM::fivep_length(NCM_type n)
{//char-'0'gives the int value of a numeral char
    if(n==std::string("ext")) return 0;
	return n[1]-'0';
}

int NCM::threep_length(NCM_type n)
{
	if(n[0]=='2') return n[2]-'0';
	return n[1]-'0';
}


//an integral type means this is an ENERGY
//otherwise, it's an EQUILIBIRIUM CONSTANT
template<typename T>
T read(std::string s, int precision=PRECISION){
    const double val = std::stod(s);
    if(std::is_integral<T>::value){
		if(s == string("0.0") || s==string("0.000000"))
			//TODO use the standard infinite energy value
			return IDENTITY;
        return T( -1 * floor( log(val) * (pow((double) 10, (double)precision)) + 0.5 ) );
    }
    else {
        return T(val);
    }
}

template<typename T>
std::unordered_map<NCM_id,T> parameters<T>::read_seqs(){
//input file looks like:
	//AAUU\t2.004\n
	//AAUC\t1.712\n
//procedure:
//for each allowed NCM
	//read the file
	//for line in file
		//generate xxx_ID object (the key)
		//put (key,value) into the map
//return the map
	std::unordered_map<NCM_id,T> temp_seqs;
	const std::vector<NCM_type> allowed_NCMs = {"13","14","15","16",
	"212", "213","221", "231","222","223","224","225","232","233","234","235","242",
	"243","244","252","253"};
	for(auto it=allowed_NCMs.begin();it!=allowed_NCMs.end();++it){
		std::vector <std::vector<std::string> > data =
			read_datafile(*it+"_seqs.txt",op);
		for(auto it2=data.begin();it2!=data.end();++it2){
			NCM_id key = NCM_id(*it,(*it2)[0]);
			T value = read<T>((*it2)[1].c_str());
			temp_seqs[key] = value;
		}
	}
	return temp_seqs;
}

template<typename T>
std::unordered_map<junction_id,T> parameters<T>::read_junctions(){
	std::unordered_map<junction_id,T> junctions;
	std::vector <std::vector<std::string> > data =
			read_datafile("junctions.txt",op);
	for(auto it = data.begin();it!=data.end();++it){
		junction_id key = junction_id((*it)[0],(*it)[1]);
		T value = read<T>((*it)[2].c_str());
		junctions[key] = value;
	}
	return junctions;
}

template<typename T>
std::unordered_map<connect_id,T> parameters<T>::read_connect(){
	std::unordered_map<connect_id,T> connect;
	std::vector <std::vector<std::string> > data =
			read_datafile("connect.txt",op);
	for(auto it = data.begin();it!=data.end();++it){
		connect_id key = connect_id((*it)[0],(*it)[1], (*it)[2]);
		T value = read<T>((*it)[3].c_str());
		connect[key] = value;
	}
	return connect;
}

template<typename T>
std::unordered_map<hinge_id,T> parameters<T>::read_hinges(){
	std::unordered_map<hinge_id,T> hinges;
	std::vector <std::vector<std::string> > data =
			read_datafile("hinges.txt",op);
	for(auto it = data.begin();it!=data.end();++it){
		hinge_id key = hinge_id((*it)[0],(*it)[1],(*it)[2]);
		T value = read<T>((*it)[3].c_str());
		hinges[key] = value;
	}
	return hinges;
}

template<typename T>
std::unordered_map<pair_id,T> parameters<T>::read_pairs(){
	std::unordered_map<pair_id,T> pairs;
	std::vector <std::vector<std::string> > data =
			read_datafile("pair.txt",op);
	for(auto it = data.begin();it!=data.end();++it){
		pair_id key = pair_id(string(""),(*it)[0]);
		T value = read<T>((*it)[1].c_str());
		pairs[key] = value;
	}
	return pairs;
}

template<typename T>
std::unordered_map<NCM_type,T> parameters<T>::read_aalberts(){
	std::unordered_map<NCM_type,T> aalberts;
	std::vector <std::vector<std::string> > data =
			read_datafile("aalberts.txt",op);
	for(auto it = data.begin();it!=data.end();++it){
		NCM_type key = (*it)[0];
		T value = read<T>((*it)[1].c_str());
		aalberts[key] = value;
	}
	return aalberts;
}

template<typename T>
T parameters<T>::read_normalization(){
	T value;
	std::vector <std::vector<std::string> > data =
			read_datafile("normalization.txt",op);
	value = read<T>(data[0][0]);
	return value;
}

const std::vector<NCM_type> single_NCMs = {"13","14","15","16"};

const std::vector<NCM_type> double_NCMs = {/*"212","213","221","231",*/
	"222","223","224","225","232","233","234","235","242","243",
	"244","252","253"
//    ,"ext"};
    };

const std::vector <NCM_type> single_and_double = {/*"212","213","221","231",*/
	"222","223","224","225","232","233","234","235","242","243",
	"244","252","253","13","14","15","16"
    };

const std::vector<pairtype> valid_pairtypes = {"cWW","cWH","cHW","cWS","cSW","cHH","cHS","cSS","cSH",
				"tWW","tWH","tHW","tWS","tSW","tHH","tHS","tSS","tSH"};


std::map<NCM_type,int> array_index_map(){
	std::map<NCM_type,int> m;
	std::vector<NCM_type> v = NCM::all_NCMs();
	for(unsigned int i=0;i<v.size();i++){
		m[v[i]] = i;
	}
	return m;
}
const std::map<NCM_type,int> indmap = array_index_map();

int NCM::to_int(const NCM_type& n) {
	return indmap.at(n);
}

std::vector<NCM_type> NCM::hairpins(){return single_NCMs;}
std::vector<NCM_type> NCM::stack_bulge(){return double_NCMs;}
std::vector<NCM_type> NCM::all_NCMs(){return single_and_double;}
std::vector<pairtype> NCM::all_pairtypes(){return valid_pairtypes;}


template class parameters<real_t>;
template class parameters<int>;
