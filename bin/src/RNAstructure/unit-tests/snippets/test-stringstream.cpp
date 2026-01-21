#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>

using namespace std;

template<class T>
T get()  {
  T setting;
  return setting;
}


void show(string input) {
  cout << "starting '" << input << "'" << endl;
  string param;
  stringstream ss(input);
  while(ss.good()) {
    getline( ss, param, ';' );
    cout << "'" << param << "' -- " << " g" << ss.good() << " e" << ss.eof() << " f" << ss.fail() << " b" << ss.bad() << endl;
  }
}

void to_str(const string& s) {
   cout << "s = '" << s << "'" << endl;
}

void getnums(string input) {
  const char* sep = "\t";
  stringstream ss(input);
  double d;
  int i=0;
  cout << endl << "\"" << input << "\"" << endl;
  while(!ss.fail()) {
    d=7.77;
    cout << setw(4) << ++i;
    ss >> d;
    cout << sep << " d=" << d;
    cout << sep << " good=" << ss.good();
    cout << sep << " fail=" << ss.fail();
    cout << sep << " bad=" << ss.bad();
    cout << sep << " eof=" << ss.eof();
    cout << endl;
  }
}

int main() {

   to_str(string("hi") + " byte!");
 
int i = get<int>();
string s = get<string>();
cout << "int " << i << endl;
cout << "str '" << s << "'" << endl;


//   to_str( ((stringstream("starting: ")) << "frogs" ).str());

   show("");
   show("   ");
   show("  ;");
   show("  ;  ");
   show("a");
   show("a;");
   show("a;   ");
   show("  a;   ");
   show("  a   ;   ");
   show("a  ;");
   show("a;   b   ;\r\nc; d");
   show("a; \r\n  b   ;\r\nc; d;");
   getnums("1");
   getnums("1 ");
   getnums("1 3    ");
   getnums("   1 3  ");
   getnums("   1 3");
   getnums("1 blarg");
   getnums("1 blarg ");

   return 0;
}
