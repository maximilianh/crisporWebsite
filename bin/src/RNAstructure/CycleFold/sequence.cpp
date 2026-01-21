#include "sequence.h"
#include <algorithm>
#include <string>
#include <iostream>
using std::string;

string check_sequence(const string& s){
	string tmp; tmp.reserve(s.size());
	for (auto &c: s){
		if (::isspace(c)) continue;
		switch (c) {
			case 'A': case 'G': case 'C': case 'U':
				tmp += c; break;
			case 'a': case 'g': case 'c': case 'u':
				tmp += toupper(c); break;
			case 'T' : case 't':
				tmp += 'U'; break;
			default:
				throw string("unrecognized character '") + string(1, c) + string("' in sequence\nmaybe you need to use the -s flag to specify a .seq file?");
		}
	}
	return tmp;
}
sequence::sequence(string t,string s)
:	seq(check_sequence(s)),
	tag(t)
{}

int sequence::getLength() const
{
	return seq.size();
}

string sequence::getTag() const
{
	return tag;
}

string sequence::substr(int pos, int length) const
{
    return seq.substr(pos,length);
}

string sequence::toString() const
{
	return seq;
}

