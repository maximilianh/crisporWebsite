#ifndef __SEQUENCE_H__
#define __SEQUENCE_H__
#include <string>

class sequence{
public:
	sequence(std::string t,std::string s);//call with (tag,sequence)
	int getLength() const;
	std::string toString() const;
    std::string substr(int pos, int length) const;
	std::string getTag() const;
//private:
//these cause a compiler error if they're const. but they shoud
	//not be changed.
	std::string seq;
	std::string tag;
};

#endif
