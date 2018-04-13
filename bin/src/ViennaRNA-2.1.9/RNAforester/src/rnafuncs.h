#ifndef _RNA_FUNCS_H
#define _RNA_FUNCS_H

#include <deque>
#include <list>
#include <string>

#include "types.h"
#include "rnaforester_options.h"
#include "rna_profile_alignment.h"

#ifdef HAVE_LIBXMLPLUSPLUS
#ifdef HAVE_LIBXML2
#include <libxml++/libxml++.h>
#endif
#endif

using namespace std;
	
class RNAFuncs
{
 public:
  struct SquigglePlotOptions
  {
	  bool hideBaseNumbers;
	  Uint baseNumInterval;
	  bool greyColors;
	  bool generatePNG;
	  bool generateJPG;
          bool generateFIG;
	  double scale;
  };

  typedef struct {
    map<int,string> idmapping;
    map<int,string> comments;
    map<int,string> descriptions;
    map<int,string> names;
    map<int,string> synonyms;
    Uint xbasepos,ybasepos;
    bool xmlInput;
  } AddXmlInfos;

  static bool isRNAString(const string &str);
  static bool isViennaString(const string &str, Ulong &basePairCount, Ulong &maxDepth);
  static void drawRNAStructure(const string &seq, const string &structure, const string &filename_prefix, const string &structname, const list<pair<Uint,Uint> > &regions, const SquigglePlotOptions &options);
  static void drawRNAAlignment(const string &structure, const string &altStructure,  const string &seq1, const string &seq2, const string &strname1, const string &strname2, const string &filename_prefix, const bool atX, const SquigglePlotOptions &options);
  static void generateRNAAlignmentXML(const string &structure, const string &altStructure, const string &seq1, const string &seq2, const string &strname1, const string &strname2, ostream &s);
  static void printAli(const string &name1, const string &id2, const string &seq1, const string &seq2, const string &str1, const string &str3);
  
  static Uint treeSize(const string &viennaStr);
	
#ifdef HAVE_LIBXMLPLUSPLUS
#ifdef HAVE_LIBXML2
	static string getXSDURL();
	static void printMAliXML(deque<pair<double,RNAProfileAlignment*> > &resultList, const RNAforesterOptions &options,double &minPairProb,AddXmlInfos &xmlInfos,const string &outputFile);
	static void printPAliXML(const string &id1, const string &name2, const string &seq1, const string &seq2, const string &str1, const string &str3, double &score, const RNAforesterOptions &options,AddXmlInfos &xmlInfos,const string &outputFile);
	static void printMapping(map<int,string> &mapping);
#endif	
#endif
	static string UpperCase(const string &str);
};

#endif




