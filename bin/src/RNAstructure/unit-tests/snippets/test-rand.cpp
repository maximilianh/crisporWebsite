#include "../src/random.cpp"
#include <iostream>
#include <vector>

using namespace std;

int main() {
	rand64 r(234);
	RandomOutput::Int<rand64> n(r);

	vector<string> s;
	s.push_back("hello world");
	s.push_back("nice day!");
	s.push_back("popcorn!!!");

	for(int i = 0; i < 100; i++) {
		cout << r(s)  << endl;
	}
	for(int i = 0; i < 100; i++) {
		cout << n()  << endl;
	}
	
	return 0;
}