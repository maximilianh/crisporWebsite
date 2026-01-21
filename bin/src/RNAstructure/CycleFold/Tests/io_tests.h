#include "../io.h"
using std::vector;
using std::string;

BOOST_AUTO_TEST_CASE(test_split){
	string text = string("this is a test");
	vector<string> test = split_on_delim(text, std::string(" "));
	vector<string> expect = {string("this"),string("is"),string("a"),string("test")};
	BOOST_CHECK_EQUAL_COLLECTIONS(test.begin(),test.end(), expect.begin(),expect.end());


	vector<string> test2 = split_on_delim(string("this       is a test"), string(" "));
	BOOST_CHECK_EQUAL_COLLECTIONS(test2.begin(),test2.end(), expect.begin(),expect.end());

	vector<string> test3 = split_on_delim(string(" this       is a test   "), string(" "));
	BOOST_CHECK_EQUAL_COLLECTIONS(test3.begin(),test3.end(), expect.begin(),expect.end());

	vector<string> test4 = split_on_delim(string(" this    \n   is a test   "), string(" \n"));
	BOOST_CHECK_EQUAL_COLLECTIONS(test4.begin(),test4.end(), expect.begin(),expect.end());
}

BOOST_AUTO_TEST_CASE(test_trim){
	string expect = string("this is a test");
	BOOST_CHECK_EQUAL(trim(string("  this is a test\n")), expect);
}
