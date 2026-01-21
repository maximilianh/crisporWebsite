#ifndef __NCM_PARAMETERS__
#define __NCM_PARAMETERS__

#include <map>//should this be an unordered map?
#include <unordered_map>
//if so need to figure out hash function
#include <string>
#include "constants.h"
#include "options.h"

class NCM_id{
public:
	NCM_id(NCM_type,std::string);
	NCM_id(const int i,const int j,NCM_type n,const std::string& RNA);
	std::string getID() const;
	std::string getType() const;
	std::string getOuterNucs() const;
	std::string getInnerNucs() const;
    bool operator==(const NCM_id& rhs) const;
private:
	const NCM_type ncm;
    const std::string five;
    const std::string three;
	const std::string id;
	std::string calcFive(const int i,const int j,NCM_type n,const std::string& RNA);
	std::string calcThree(const int i,const int j,NCM_type n,const std::string& RNA);
    bool valid_NCM();
};
/*
template <typename T>
inline bool operator==(const T& lhs, const T& rhs){
	return lhs.getID()==rhs.getID();
}
*/


//we have to overload std::less for the NCM_id class
//so we can use it as a key in a std::unordered_map


class junction_id{
public:
	junction_id(NCM_type,NCM_type);
	std::string getID() const;
    bool operator==(const junction_id& rhs) const;
private:
	const std::string id;
	const NCM_type inner;
	const NCM_type outer;
};

class connect_id{
public:
	connect_id(NCM_type,NCM_type,nucs);
	std::string getID() const;
    bool operator==(const connect_id& rhs) const;
private:
	const std::string id;
	const NCM_type inner;
	const NCM_type outer;
    const nucs _nucs;
};

class hinge_id{
public:
	hinge_id(NCM_type,NCM_type,pairtype);
	std::string getID() const;
    bool operator==(const hinge_id& rhs) const;
private:
	const std::string id;
	const NCM_type inner;
	const NCM_type outer;
	const pairtype pair;
};

class pair_id{
public:
	pair_id(pairtype,nucs);
	std::string getID() const;
    bool operator==(const pair_id& rhs) const;
private:
	const std::string id;
	const pairtype pair;
	const nucs seq;
};


namespace std{
    template<> struct less<NCM_id>{
       bool operator() (const NCM_id& lhs, const NCM_id& rhs) const{
           return lhs.getID() < rhs.getID();
       }
    };
    template<> struct less<junction_id>{
       bool operator() (const junction_id& lhs, const junction_id& rhs) const{
           return lhs.getID() < rhs.getID();
       }
    };
    template<> struct less<connect_id>{
       bool operator() (const connect_id& lhs, const connect_id& rhs) const{
           return lhs.getID() < rhs.getID();
       }
    };
    template<> struct less<hinge_id>{
       bool operator() (const hinge_id& lhs, const hinge_id& rhs) const{
           return lhs.getID() < rhs.getID();
       }
    };
    template<> struct less<pair_id>{
       bool operator() (const pair_id& lhs, const pair_id& rhs) const{
           return lhs.getID() < rhs.getID();
       }
    };
}


namespace std {
  template <>
  struct hash<NCM_id>
  {
    std::size_t operator()(const NCM_id& k) const {
      return std::hash<std::string>()(k.getID());
    }
  };
  template <>
  struct hash<junction_id>
  {
    std::size_t operator()(const junction_id& k) const {
      return std::hash<std::string>()(k.getID());
    }
  };
  template <>
  struct hash<connect_id>
  {
    std::size_t operator()(const connect_id& k) const {
      return std::hash<std::string>()(k.getID());
    }
  };
  template <>
  struct hash<hinge_id>
  {
    std::size_t operator()(const hinge_id& k) const {
      return std::hash<std::string>()(k.getID());
    }
  };
  template <>
  struct hash<pair_id>
  {
    std::size_t operator()(const pair_id& k) const {
      return std::hash<std::string>()(k.getID());
    }
  };

}

template<typename T>
class parameters{
public:
	parameters(const options);
	T energy(NCM_id,NCM_id,bool, bool debug=false) const;
	T energy(NCM_id,NCM_type,bool) const;
	T energy(NCM_id, bool debug=false) const;
	T junction_energy(NCM_id,NCM_id,bool) const;
	T mb_nuc_penalty() const;
	T mb_stem_penalty() const;
	T mb_closure_penalty() const;
	T aalberts_entropy(NCM_type) const;
	T pair_bonus() const;
  T pair_energy(NCM_id theta, bool debug) const;
    bool mfe_calc() const;
private:
	const options op;
	const std::unordered_map<NCM_id,T> seqs;
	std::unordered_map<NCM_id,T> read_seqs();
	const std::unordered_map<junction_id,T> junctions;
	std::unordered_map<junction_id,T> read_junctions();
	const std::unordered_map<hinge_id,T> hinges;
	std::unordered_map<hinge_id,T> read_hinges();
	const std::unordered_map<pair_id,T> pairs;
	std::unordered_map<pair_id,T> read_pairs();
	const std::unordered_map<NCM_type,T> aalberts;
	std::unordered_map<NCM_type,T> read_aalberts();
	const std::unordered_map<connect_id,T> connect;
	std::unordered_map<connect_id,T> read_connect();
	const T normalization;
	T read_normalization();
    const bool mfe;
};

namespace NCM {
	int to_int(const NCM_type& n);
	int fivep_length(NCM_type n);
	int threep_length(NCM_type n);
	std::vector<NCM_type> hairpins();
	std::vector<NCM_type> stack_bulge();
	std::vector<NCM_type> all_NCMs();
	std::vector<pairtype> all_pairtypes();
}

std::map<NCM_type,int> array_index_map();
#endif //__NCM_PARAMETERS__
