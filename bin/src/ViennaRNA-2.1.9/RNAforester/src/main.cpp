/*
  Copyright by Matthias Hoechsmann (C) 2002-2004
  =====================================                                   
  You may use, copy and distribute this file freely as long as you
  - do not change the file,
  - leave this copyright notice in the file,
  - do not make any profit with the distribution of this file
  - give credit where credit is due
  You are not allowed to copy or distribute this file otherwise
  The commercial usage and distribution of this file is prohibited
  Please report bugs and suggestions to <mhoechsm@TechFak.Uni-Bielefeld.DE>
*/

#include <deque>
#include <functional>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <list>
#include <sstream>
#include <string>
#include <cstring>
#include <climits>
#include <map>

//#include <sys/timeb.h>
#include <sys/times.h>
#include <unistd.h>
#include <string.h>

#ifndef WIN32
#include "config.h"
#endif

#include "Arguments.h"
#include "alignment.h"
#include "debug.h"

//#include "global_alignment.h"
#include "treeedit.h"

#ifdef HAVE_LIBXMLPLUSPLUS
#ifdef HAVE_LIBXML2
#include <libxml++/libxml++.h>
#include <libxml/xmlschemas.h>
#include "xsd.h"
#endif
#endif

#include "misc.h"
#include "progressive_align.h"
#include "rna_alignment.h"
#include "rnaforest.h"
#include "rnaforestsz.h"
#include "rnafuncs.h"
#include "rna_profile_alignment.h"
#include "rna_algebra.h"
#include "rnaforester_options.h"

#include "alignment.t.cpp"
//#include "global_alignment.t.cpp"
#include "treeedit.t.cpp"
//#include "ppforest.t.cpp"

using namespace std;

/* ****************************************** */
/*          Definitions and typedefs          */
/* ****************************************** */

struct ToLower : public unary_function<char,char> {
	char operator()(char a)
	{
		return tolower(a);
	};
};

template <class L>
void makeDotFileAli(const PPForest<L> &ppf, const RNAforesterOptions &options)
{
	if(options.has(RNAforesterOptions::OutputAlignmentDotFormat))
	{
		string filename;
		options.get(RNAforesterOptions::OutputAlignmentDotFormat,filename,string("ali.dot"));
		ofstream s(filename.c_str());
		ppf.printDot(s);
	}
}

template <class L>
void makeDotFileInp(const PPForest<L> &ppf, const RNAforesterOptions &options, Uint count)
{
	if(options.has(RNAforesterOptions::MakeDotForInputTrees))
	{
		ostringstream ss;
		ofstream os;
		ss << "input" << count << ".dot";
		os.open(ss.str().c_str());
		ppf.printDot(os);
		os.close();
	}
}

static const string RNAFORESTER_VERSION = "1.5";
static const string PROMPT = "Input string (upper or lower case); & to end for multiple alignments, @ to quit\n";
static const string SCALE = "....,....1....,....2....,....3....,....4....,....5....,....6....,....7....,....8\n";

void alignMultiple(deque<RNAProfileAlignment*> &alignList, Score &score,const RNAforesterOptions &options,RNAFuncs::AddXmlInfos &xmlInfos);
void alignPairwise(deque<RNAForest*> &inputListPW,Score &score,const RNAforesterOptions &options,RNAFuncs::AddXmlInfos &xmlInfos);
void cutAfterChar(string &s,char c);

void editPairwise(list<RNAForestSZ*> &inputListSZ,Score &score,RNAforesterOptions &options);
void alignPairwiseSimple(deque<RNAForest*> &inputListPW,Score &score,RNAforesterOptions &options);

#ifdef HAVE_LIBXMLPLUSPLUS
#ifdef HAVE_LIBXML2
extern "C" {
  bool validateXSD(string filename);
}
#endif
#endif

static void showversion(const char *prog)
{
	cout << prog << ", version " << RNAFORESTER_VERSION << endl;
	cout << "Copyright Matthias Hoechsmann 2001-2004," << endl << "mhoechsm@techfak.uni-bielefeld.de" << endl;
}

int main(int argc, const char **argv)
{	
	string buffer;
	string baseStr,viennaStr,nameStr;
	Ulong basePairCount,maxDepth;
	deque<RNAForest*> inputListPW;
	deque<RNAProfileAlignment*> alignList;
	bool showScale=true,multipleAlign=false;
	istream *inputStream=NULL;
	istringstream inputString;
	string tmpInput = "";
	ifstream *inputFile=NULL;
	Uint structure_count=1;
	int suboptPercent=100;
	double minPairProb=0.25;
#ifdef HAVE_LIBXMLPLUSPLUS
#ifdef HAVE_LIBXML2
	xmlpp::Document* xmlDoc;
#endif
#endif
	string xmlOrig;
	map<int,string> comments;
	map<int,string> idMapping;
	map<int,string> descriptions;
	map<int,string> names;
	map<int,string> synonyms;
	int seqID = 0;
	stringstream ss;
	RNAFuncs::AddXmlInfos xmlInfos;
	
	list<RNAForestSZ*> inputListSZ;

	try
	{
		RNAforesterOptions options(argc,argv);

		// check options
		if(options.has(RNAforesterOptions::Help))
		{
			options.help();
			exit(EXIT_SUCCESS);
		}

		if(options.has(RNAforesterOptions::SecretHelp))
		{
			options.secretHelp();
			exit(EXIT_SUCCESS);
		}

		if(options.has(RNAforesterOptions::Version))
		{
			showversion(argv[0]);
			exit(EXIT_SUCCESS);      
		}

		// read score values
		Score score(options);
		if(!options.has(RNAforesterOptions::ShowOnlyScore))
			score.print();

		// check option suboptimals
		if(options.has(RNAforesterOptions::LocalSubopts))
		{
			options.get(RNAforesterOptions::LocalSubopts,suboptPercent,100);
			if(suboptPercent<0 || suboptPercent>100)
			{
				cerr << "error: value for parameter --subopts must be in range from 0 to 100" << endl;
				exit(EXIT_FAILURE);
			}
		}

#ifdef HAVE_LIBRNA  // This features require the ViennaRNA library		
		options.get(RNAforesterOptions::PredictMinPairProb,minPairProb,0.25);
		if(options.has(RNAforesterOptions::PredictProfile))
		  cout << "Minimum required basepair probability (-pmin): " << minPairProb << endl;
#endif
		
		// show if suboptimals
		if(!options.has(RNAforesterOptions::ShowOnlyScore))
			if(options.has(RNAforesterOptions::LocalSimilarity) && options.has(RNAforesterOptions::LocalSubopts))
			  cout << "calculate suboptimals within " << suboptPercent << "% of global optimum" << endl << endl;

		// profile search
		if(options.has(RNAforesterOptions::ProfileSearch))
		  {
		    string filename;
		    
		    options.get(RNAforesterOptions::ProfileSearch,filename,string(""));
		    if(filename=="")
		      {
			cerr << "no profile filename" << endl;
			exit(EXIT_FAILURE);
		      }
		    
		    RNAProfileAlignment *rnaProfileAli=new RNAProfileAlignment(filename);
		    alignList.push_back(rnaProfileAli);
		  }

		if(options.has(RNAforesterOptions::ReadFromFile))
		{
			string filename;
			options.get(RNAforesterOptions::ReadFromFile,filename,string(""));
			inputFile=new ifstream(filename.c_str());
			if(inputFile->fail())
			{
				cerr << "cannot open file: \"" << filename << "\"" << endl;
				exit(EXIT_FAILURE);
			}
#ifdef HAVE_LIBXMLPLUSPLUS
#ifdef HAVE_LIBXML2 							
			// getline extracts all characters of the first line until '\n' to check whether inputFile is an xml file
			getline(*inputFile,buffer);
			
			// input file is an xml file 
			if(buffer.find("<?xml",0)==0){
			  buffer="";
			  
			  // validation
			  //if(!validateXSD(filename)){
			  //  exit(EXIT_FAILURE);
			  //}
			  
			  // create Dom parser
			  xmlpp::DomParser domParser;
			  domParser.parse_file(filename);
			 			  
			  xmlDoc = domParser.get_document();
			  xmlpp::Element* elemRoot = xmlDoc->get_root_node();
			  xmlpp::Node::NodeList nodeListRnaStructure = elemRoot->get_children();
			  
			  // for each rnastructure element
			  xmlpp::Node::NodeList::iterator it1;
			  xmlpp::Node::NodeList::iterator it2;
			  xmlpp::Node::NodeList::iterator it3;
			  for(it1=nodeListRnaStructure.begin();it1!=nodeListRnaStructure.end();it1++){
			    if((*it1)->get_name() == "rnastructure"){
			      seqID++;
			      tmpInput.append(">");
			      ss << seqID;
			      tmpInput.append(ss.str());
			      ss.str("");
			      tmpInput.append("\n");
			    }
			   			       
			    xmlpp::Node::NodeList nodeListSeq = (*it1)->get_children();
			    for(it2=nodeListSeq.begin();it2!=nodeListSeq.end();it2++){
			      if((*it2)->get_name()=="sequence"){
				xmlpp::Element* elemSeq = (xmlpp::Element*)(*it2);
				
				// map internal id to seqID :: JK debug
				cout << "seqID: " << seqID << endl;
				cout << "idMapping[seqID]: " << elemSeq->get_attribute("seqID")->get_value() << endl;
				//
				idMapping[seqID] = elemSeq->get_attribute("seqID")->get_value();

				xmlpp::Node::NodeList nodeListName = elemSeq->get_children();
				for(it3=nodeListName.begin();it3!=nodeListName.end();it3++){
				  if((*it3)->get_name()=="name" && ((xmlpp::Element*)(*it3))->get_child_text()){
				    // map the name to seqID
				    names[seqID] = ((xmlpp::Element*)(*it3))->get_child_text()->get_content();
				  }
				  if((*it3)->get_name()=="synonyms" && ((xmlpp::Element*)(*it3))->get_child_text()){
				    // map synonyms to seqID
				    synonyms[seqID] = ((xmlpp::Element*)(*it3))->get_child_text()->get_content();
				  }
				  if((*it3)->get_name()=="description" && ((xmlpp::Element*)(*it3))->get_child_text()){
				    // map the description to seqID
				    descriptions[seqID] = ((xmlpp::Element*)(*it3))->get_child_text()->get_content();
				  }
				  if((*it3)->get_name()=="freeSequence" || (*it3)->get_name()=="nucleicAcidSequence"){
				    xmlpp::Element* elemFreeSequence = (xmlpp::Element*)(*it3);
				    tmpInput.append(elemFreeSequence->get_child_text()->get_content());
				    tmpInput.append("\n");
				  }
				} // end for it3
			      }
			      
			      // get comment if available
			      if((*it2)->get_name()=="comment" && ((xmlpp::Element*)(*it2))->get_child_text()){
				comments[seqID] = ((xmlpp::Element*)(*it2))->get_child_text()->get_content();
			      }
			      			    			      
			      // get structure
			      if((*it2)->get_name()=="structure" && ((xmlpp::Element*)(*it2))->get_child_text()){
				xmlpp::Element* elemStr = (xmlpp::Element*)(*it2);
				tmpInput.append(elemStr->get_child_text()->get_content());
				tmpInput.append("\n");
			      }
			    } // end for it2
			    
			  } // end for it1
			  xmlOrig = xmlDoc->write_to_string();

			  // struct with additional xml infos
			  xmlInfos.idmapping = idMapping;
			  xmlInfos.comments = comments;
			  xmlInfos.descriptions = descriptions;
			  xmlInfos.names = names;
			  xmlInfos.synonyms = synonyms;
			  xmlInfos.xmlInput = true;
				  
			  inputString.str(tmpInput.c_str());
			  inputStream=&inputString;
			  
			} else {

			  xmlInfos.xmlInput = false;
			  // no xml file write the first line from buffer back to stream
			  for(int i=buffer.size();i>=0;i--){
			    inputFile->putback(buffer[i]);
			  }
#endif
#endif
			  inputStream=inputFile;

#ifdef HAVE_LIBXMLPLUSPLUS	
#ifdef HAVE_LIBXML2	  
			}
#endif
#endif
		} else {
		  inputStream=&cin;
		}
		
		if(showScale)
		  {
		    if(!options.has(RNAforesterOptions::NoScale) && !options.has(RNAforesterOptions::ReadFromFile))
		      cout << endl << PROMPT << SCALE;
		    showScale=false;
		  }

		for(;;)
		{
		  getline(*inputStream,buffer);
		  
			if(inputStream->eof())
			  {	
			    if(options.has(RNAforesterOptions::Multiple) && !options.has(RNAforesterOptions::ProfileSearch))
				buffer="&";
			    else
			      exit(EXIT_SUCCESS);
			  }

			if(buffer.empty())
				continue;

			// quit if character is @
			if(buffer[0]=='@')
			  break;

			// delete '\r' at line end from non unix files
			if(buffer[buffer.size()-1]=='\r')
				buffer.erase(buffer.size()-1);

			// check for name of structure
			if(buffer[0]=='>')
			{
			  /*if(buffer.find(" ",0) != string::npos){
			    int t = buffer.find(" ",0);
			    nameStr=&buffer[t];
			  } else {
			    nameStr=&buffer[1];
			    }*/
			  nameStr=&buffer[1];
			  continue;
			}

			// cut after blank
			cutAfterChar(buffer,' ');

			// check for aligning multiple structures
			// if input is read from file the eof has the same meaning as &
			if( buffer[0]=='&')
				multipleAlign=true;
			else
			{
				// check for base string
			  
			  if(RNAFuncs::isRNAString(buffer))
				{
				  baseStr=buffer;
				  // convert to small letters
				  transform(baseStr.begin(),baseStr.end(),baseStr.begin(),ToLower());
				  // t -> u
				  replace(baseStr.begin(),baseStr.end(),'t','u');
				  // delete '.'  and '-' from alignment files
				  remove(baseStr.begin(),baseStr.end(),'.');
				  remove(baseStr.begin(),baseStr.end(),'-');

#ifdef HAVE_LIBRNA  // This features require the ViennaRNA library				  
				  if(options.has(RNAforesterOptions::PredictProfile))
				    {
				      ostringstream ss;		 
				      string constraint;

				      // if there is no name given (> ...) use a counter
				      if(nameStr=="")
					ss << "> " << structure_count;
				      else
					ss << nameStr;	
				      
				      cout << "Predicting structure profile for sequence: " << ss.str() << endl;
				      RNAProfileAlignment *rnaProfileAli=new RNAProfileAlignment(baseStr,ss.str(),constraint,minPairProb);
				      alignList.push_back(rnaProfileAli);
				      makeDotFileInp(*rnaProfileAli,options,structure_count);
				      structure_count++;
				    }
#endif				  

				  //				  continue;
				}
			  else
			    {

				// check for vienna string	
			  if(RNAFuncs::isViennaString(buffer,basePairCount,maxDepth))
			    {

#ifdef HAVE_LIBRNA  // This features require the ViennaRNA library			      
			      // skip structure lines if structures are predicted
			      if(options.has(RNAforesterOptions::PredictProfile))
				{
				  cout << "ignoring structure: " << buffer << endl;
				  continue;
				}
#endif
			      
			     viennaStr=buffer;	       
			    }
			  else   
				{ 
					cerr << "The input sequence is neither an RNA/DNA string nor in vienna format." << endl;
					cerr << "line: " << buffer << endl;
					showScale=true;
					exit(EXIT_FAILURE);
				}

				// add structure to input list
				if(options.has(RNAforesterOptions::Multiple))
				{
					ostringstream ss;		 

					// if there is no name given (> ...) use a counter
					if(nameStr=="")
					  ss << "> " << structure_count;
					else
					  ss << nameStr;	


					RNAProfileAlignment *rnaProfileAli=new RNAProfileAlignment(baseStr,viennaStr,ss.str());
					makeDotFileInp(*rnaProfileAli,options,structure_count);					
					alignList.push_back(rnaProfileAli); 	
				}
				else
				{
				  if(options.has(RNAforesterOptions::TreeEdit))
				    {
				      RNAForestSZ *rnaForestSZ=new RNAForestSZ(baseStr,viennaStr,nameStr);
				      inputListSZ.push_back(rnaForestSZ);
				    }
				  else
				    {
				      RNAForest *rnaForest=new RNAForest(baseStr,viennaStr,nameStr);			
				      nameStr="";
				      makeDotFileInp(*rnaForest,options,structure_count);
				      inputListPW.push_back(rnaForest);
				    }
				}

				structure_count++;      
				showScale=true;
			    }
			}

			// ***** multiple alignment
			if((options.has(RNAforesterOptions::Multiple) && multipleAlign) || (options.has(RNAforesterOptions::ProfileSearch) && alignList.size()==2))
			{
				alignMultiple(alignList,score,options,xmlInfos);
				multipleAlign=false;
				structure_count=1;
				
				if(options.has(RNAforesterOptions::ProfileSearch))
				  {
				    string filename;
		    
				    options.get(RNAforesterOptions::ProfileSearch,filename,string(""));		    		    
				    RNAProfileAlignment *rnaProfileAli=new RNAProfileAlignment(filename);
				    alignList.push_back(rnaProfileAli);
				  }
				else
				  break;

				//				break;
			}

			// ***** pairwise alignment
			if(inputListPW.size()==2)
			{
			  if(options.has(RNAforesterOptions::GlobalAlignment))
			  {
			    alignPairwiseSimple(inputListPW,score,options);
			  } else 
			    alignPairwise(inputListPW,score,options,xmlInfos);
			  break;
			}
			
			if(inputListSZ.size()==2)
			{
			  editPairwise(inputListSZ,score,options);
			  break;
			}

		}

		// free dynamic allocated memory
		deque<RNAForest*>::const_iterator it;
		for(it = inputListPW.begin(); it!=inputListPW.end(); it++){
		  delete *it;
		}

		DELETE(inputFile);

		return (0);
	}
	catch(RNAforesterOptions::IncompatibleException e)
	{
		e.showError();
		return(EXIT_FAILURE);
	} 
	catch(RNAforesterOptions::RequiresException e)
	{
		e.showError();
		return(EXIT_FAILURE);
	}
#ifdef HAVE_LIBXMLPLUSPLUS
#ifdef HAVE_LIBXML2
	catch(xmlpp::validity_error ve) 
	{
	  cout << ve.what() << endl;
	//return(EXIT_FAILURE);
	} 
	catch(xmlpp::parse_error pe)
	{
	  cout << pe.what() << endl;
	  //return(EXIT_FAILURE);
	} 
	catch(xmlpp::exception e)
	{
	  cout << e.what() << endl;
	  //return(EXIT_FAILURE);
	} 
#endif
#endif
	  	
}



void cutAfterChar(string &s,char c){
	string::size_type pos=s.find(c);
	if(pos!=string::npos)
		s.erase(pos);
}

void alignMultiple(deque<RNAProfileAlignment*> &alignList, Score &score,const RNAforesterOptions &options,RNAFuncs::AddXmlInfos &xmlInfos)
{
	DoubleScoreProfileAlgebraType *alg;
	deque<pair<double,RNAProfileAlignment*> > resultList;
//	double optScore;
	Uint clusterNr=1;
	double minPairProb;
	string outputFile;

	options.get(RNAforesterOptions::ConsensusMinPairProb,minPairProb,0.5);
	
	// distance or similarity
	if(options.has(RNAforesterOptions::CalculateDistance))
		alg=new DoubleDistProfileAlgebra(score);
	else
		alg=new DoubleSimiProfileAlgebra(score);

	cout << endl;	

	progressiveAlign(alignList,resultList,alg,options);

	cout << endl;
	cout << "*** Results ***" << endl << endl; 
	cout << "Minimum basepair probability for consensus structure (-cmin): " << minPairProb << endl << endl;

	deque<pair<double,RNAProfileAlignment*> >::const_iterator it;
	for(it=resultList.begin();it!=resultList.end();it++)
	  {
	    cout << "RNA Structure Cluster Nr: " << clusterNr << endl;
	    cout << "Score: " << it->first << endl;
	    cout << "Members: " << it->second->getNumStructures() << endl << endl;
	    if(options.has(RNAforesterOptions::FastaOutput))
	      {
		it->second->printFastaAli(false);
		cout << endl;
	      }
	    else
	      {  

		it->second->printSeqAli();

#ifdef HAVE_LIBRNA  // This features require the ViennaRNA library		
		// print alignment
		it->second->printStrAli();
		cout << endl;
#endif		
	      }

	    // save profile
	    if(options.has(RNAforesterOptions::SaveProfile))
	    {
	      string filename;
	      filename=options.generateFilename(RNAforesterOptions::SaveProfile,".pro", "rna.pro",clusterNr);
      it->second->save(filename);
	    }
	    

		// print consensus structure
		cout << "Consensus sequence/structure:" << endl;
		it->second->printConsensus(minPairProb);
		cout << endl;

#ifdef HAVE_LIBG2  // This features require the g2 library		
		// generate squiggle plots
			if(options.has(RNAforesterOptions::MakeSquigglePlot))
			{
				RNAProfileAlignment::SquigglePlotOptions sqOptions;
				string filename;

				// plot showing full base information
				filename=options.generateFilename(RNAforesterOptions::MakeSquigglePlot,".ps", "rnaprofile.ps",clusterNr);
				sqOptions.greyColors=options.has(RNAforesterOptions::SquiggleGreyColors);
				sqOptions.minPairProb=minPairProb;
				sqOptions.mostLikelySequence=false;
				it->second->squigglePlot(filename,sqOptions);

				// plot showing consensus sequence
				filename=options.generateFilename(RNAforesterOptions::MakeSquigglePlot,"_cons.ps", "rnaprofile_cs.ps",clusterNr);
				sqOptions.mostLikelySequence=true;
				it->second->squigglePlot(filename,sqOptions);
			}
#endif		      

		clusterNr++;
	  }
#ifdef HAVE_LIBXMLPLUSPLUS
#ifdef HAVE_LIBXML2
	if(options.has(RNAforesterOptions::XmlOutputFile)){
	  options.get(RNAforesterOptions::XmlOutputFile,outputFile,string(""));
	} else {
	  outputFile = "result.xml";
	}

	if(options.has(RNAforesterOptions::GenerateXML))
        { 
		RNAFuncs::printMAliXML(resultList,options,minPairProb,xmlInfos,outputFile);	
	}
#endif
#endif

	/*
	if(!options.has(RNAforesterOptions::ShowOnlyScore))
	{
		cout << "global optimal score: ";
	}
	cout << optScore << endl;

	list<RNAProfileAlignment*>::iterator it=inputMapProfile.begin();
	RNAProfileAlignment *ppf=it;
	if(!options.has(RNAforesterOptions::ShowOnlyScore))
	{
		// generate dot file
		makeDotFileAli(*ppf,options);

		ppf->print();
		cout << endl;
		ppf->printConsensus();
	}
	*/

	// save profile alignment to binary file
	/*	      if(options.has(RNAforesterOptions::SaveMultipleAliFile))
	{
	string filename;

	filename=generateFilename(options,RNAforesterOptions::SaveMultipleAliFile,".mta", "unknown.dot");
	ofstream s(filename.c_str());
	f1->save(s);
	}
	*/

	/*
	// generate squiggle plots
	if(options.has(RNAforesterOptions::MakeSquigglePlot))
	{
		RNAProfileAlignment::SquigglePlotOptions sqOptions;
		string filename;

		// plot showing full base information
		filename=options.generateFilename(RNAforesterOptions::MakeSquigglePlot,".ps", "rnaprofile.ps");
		sqOptions.greyColors=options.has(RNAforesterOptions::SquiggleGreyColors);
		sqOptions.mostLikelySequence=false;
		ppf->squigglePlot(filename,sqOptions);

		// plot showing consensus sequence
		filename=options.generateFilename(RNAforesterOptions::MakeSquigglePlot,"_cons.ps", "rnaprofile_cons.ps");
		sqOptions.greyColors=options.has(RNAforesterOptions::SquiggleGreyColors);
		sqOptions.mostLikelySequence=true;
		ppf->squigglePlot(filename,sqOptions);
	}
	*/

	
	delete alg;
}


void alignPairwise(deque<RNAForest*> &inputListPW,Score &score,const RNAforesterOptions &options,RNAFuncs::AddXmlInfos &xmlInfos)
{
	IntScoreRNA_AlgebraType *alg;
	Uint xbasepos,ybasepos,xlen,ylen;
	list<pair<Uint,Uint> > xsubopts;
	list<pair<Uint,Uint> > ysubopts;
	string seq1,seq2,str1,str2;
	char s[8];
	int suboptPercent;
	Uint count=1;
	RNAFuncs::SquigglePlotOptions sqOptions;
	tms tmsStart, tmsEnd;
	string outputFile;

	// read options
	options.get(RNAforesterOptions::LocalSubopts,suboptPercent,100);

#ifdef HAVE_LIBG2  // This features require the g2 library	
	// generate squiggle plot
	if(options.has(RNAforesterOptions::MakeSquigglePlot))
	{	      
		// get sq options
		sqOptions.hideBaseNumbers=options.has(RNAforesterOptions::SquiggleHideBaseNumbers);
		options.get(RNAforesterOptions::SquiggleBaseNumberInterval,sqOptions.baseNumInterval,(Uint)20);
		sqOptions.greyColors=options.has(RNAforesterOptions::SquiggleGreyColors);
		options.get(RNAforesterOptions::SquiggleScaleFactor,sqOptions.scale,1.0);
		sqOptions.generateFIG=options.has(RNAforesterOptions::SquiggleGenerateFIG);		
#ifdef HAVE_LIBGD
		sqOptions.generatePNG=options.has(RNAforesterOptions::SquiggleGeneratePNG);
		sqOptions.generateJPG=options.has(RNAforesterOptions::SquiggleGenerateJPG);
#endif
		
	}
#endif
	
	// distance or similarity
	if(options.has(RNAforesterOptions::CalculateDistance))
		alg=new IntDistRNA_Algebra(score);
	else
	{
		if(options.has(RNAforesterOptions::RIBOSUMScore))
			alg=new RIBOSUM8560(score);
		else
			alg=new IntSimiRNA_Algebra(score);
	}

	RNAForest *f1=inputListPW.front();
	RNAForest *f2=inputListPW.back();

	if(options.has(RNAforesterOptions::SpaceTimeInfo))
	  {
	    IntScore_AlgebraType *algGlobClassic;

	    algGlobClassic=new IntDist_Algebra(score);

	    if(!options.has(RNAforesterOptions::ShowOnlyScore))
	      {
		cout << "F1_NUMNODES" << ";";
		cout << "F2_NUMNODES" << ";";
		cout << "F1_DEGREE" << ";";
		cout << "F2_DEGREE" << ";";
		cout << "F1_LEAVES" << ";"; 
		cout << "F2_LEAVES" << ";";
		cout << "F1_DEPTH" << ";"; 
		cout << "F2_DEPTH" << ";";
		cout << "F1_NUMCSFS" << ";";
		cout << "F2_NUMCSFS" << ";";
		cout << "TABLE_SIZE_4D" << ";";
		cout << "TABLE_SIZE_2D" << ";";
		cout << "GLOBALI_TIME" << ";";
		cout << "GLOBALI_TIME_SPEEDUP" << ";";
		cout << "LOCALALI_TIME" << ";";
		cout << "LOCALALI_TIME_SPEEDUP" << ";";
		
		cout << endl;
	      }

	    cout << f1->size() << "\t";
	    cout << f2->size() << "\t"; 
	    cout << f1->maxDegree() << "\t";
	    cout << f2->maxDegree() << "\t";
	    cout << f1->numLeaves() << "\t";
            cout << f2->numLeaves() << "\t";
	    cout << f1->maxDepth() << "\t";
	    cout << f2->maxDepth() << "\t";
	    cout << f1->getNumCSFs() << "\t";
	    cout << f2->getNumCSFs() << "\t";
	    cout << f1->size()*f2->size()*f1->maxDegree()*f1->maxDegree() << "\t";
	    cout << f1->getNumCSFs()*f2->getNumCSFs() << "\t";

 	    // global alignment
	    {
	      times(&tmsStart);
	      Alignment<int,RNA_Alphabet,RNA_AlphaPair> ali(f1,f2,*algGlobClassic,false,true);
	      times(&tmsEnd);
	      cout <<((double) (tmsEnd.tms_utime - tmsStart.tms_utime))/sysconf(_SC_CLK_TCK) << "\t"; 
	      //	    cerr << "#" << ali.getGlobalOptimum() << "#";
	    }

 	    // global alignment speedup
	    {
	      times(&tmsStart);
	      Alignment<int,RNA_Alphabet,RNA_AlphaPair> ali(f1,f2,*algGlobClassic,false,false);
	      times(&tmsEnd);
	      cout <<((double) (tmsEnd.tms_utime - tmsStart.tms_utime))/sysconf(_SC_CLK_TCK) << "\t"; 
	      //	    cerr << "#" << ali.getGlobalOptimum() << "#";
	    }

 	    // local alignment
	    {
	      times(&tmsStart);
	      Alignment<int,RNA_Alphabet,RNA_AlphaPair> ali(f1,f2,*algGlobClassic,true,true);
	      times(&tmsEnd);
	      cout <<((double) (tmsEnd.tms_utime - tmsStart.tms_utime))/sysconf(_SC_CLK_TCK) << "\t"; 
	      //	    cerr << "#" << ali.getLocalOptimum() << "#";
	    }

 	    // local alignment speedup
	    {
	      times(&tmsStart);
	      Alignment<int,RNA_Alphabet,RNA_AlphaPair> ali(f1,f2,*algGlobClassic,true,false);
	      times(&tmsEnd);
	      cout <<((double) (tmsEnd.tms_utime - tmsStart.tms_utime))/sysconf(_SC_CLK_TCK) << "\t"; 
	      //	    cerr << "#" << ali.getLocalOptimum() << "#";
	    }
	    

	    cout << endl;

	    exit(EXIT_SUCCESS);
	  }
	
	Alignment<int,RNA_Alphabet,RNA_AlphaPair> ali(f1,f2,*alg);
	RNA_Alignment ppfali;
	ppfali.setStructureNames(f1->getName(),f2->getName());
	double optScore;


	if(!options.has(RNAforesterOptions::ShowOnlyScore))
	{
	  if(options.has(RNAforesterOptions::SmallInLarge))
	    {
	      cout << "small-in-large optimal score: ";
	    }
	  else
	    {
	      if(options.has(RNAforesterOptions::LocalSimilarity))
		cout << "local optimal score: ";
	      else
		cout << "global optimal score: ";
	    }
	}

	if(options.has(RNAforesterOptions::SmallInLarge))
	  {
	    cout << ali.getSILOptimum() << endl;
	    optScore = ali.getSILOptimum();
	  }
	else
	  {
	    if(options.has(RNAforesterOptions::LocalSimilarity))
	      {
		cout << ali.getLocalOptimum() << endl;
		optScore = ali.getLocalOptimum();
	      }  
	    else
	      {
		cout << ali.getGlobalOptimum() << endl;
		optScore = ali.getGlobalOptimum();
		if(options.has(RNAforesterOptions::RelativeScore))
		  cout << ali.getGlobalOptimumRelative() << endl;
	      }
	  }


	if(!options.has(RNAforesterOptions::ShowOnlyScore))
	{
	  if(options.has(RNAforesterOptions::SmallInLarge))
	    {
	      ali.getOptSILAlignment(ppfali,ybasepos);
	      cout << "starting at position: " << ybasepos << endl << endl; 
	    }
	  else
	    {	      
	      if(options.has(RNAforesterOptions::LocalSimilarity))
		{
		  ali.resetOptLocalAlignment(suboptPercent);
		  ali.getOptLocalAlignment(ppfali,xbasepos,ybasepos);		      		    
		  cout << "starting at positions: " << xbasepos << "," << ybasepos << endl << endl; 
		  xmlInfos.xbasepos = xbasepos;
		  xmlInfos.ybasepos = ybasepos;
		}
	      else		  
		ali.getOptGlobalAlignment(ppfali);		  
	    }

		// generate dot file
		makeDotFileAli(ppfali,options);

		// print alignment
		ppfali.getSequenceAlignments(seq1,seq2);				
		ppfali.getStructureAlignment(str1,true);
		ppfali.getStructureAlignment(str2,false);

		if(options.has(RNAforesterOptions::FastaOutput)){
		  cout << ppfali.getStructureNameX() << endl;
		  cout << seq1 << endl;
		  cout << str1 << endl;
		  cout << ppfali.getStructureNameY() << endl;
		  cout << seq2 << endl;
		  cout << str2 << endl;
		  cout << endl;
		}
		else {
		  RNAFuncs::printAli(ppfali.getStructureNameX(),ppfali.getStructureNameY(),seq1,seq2,str1,str2);
		}

		xlen=seq1.size();
		ylen=seq2.size();
		if(options.has(RNAforesterOptions::LocalSimilarity))
		{
			xsubopts.push_back(pair<Uint,Uint>(xbasepos,xlen));
			ysubopts.push_back(pair<Uint,Uint>(ybasepos,ylen));
		}

#ifdef HAVE_LIBG2  // This features require the g2 library		
		if(options.has(RNAforesterOptions::MakeSquigglePlot))
		  {
		    sprintf(s,"%d",count);
		    ppfali.squigglePlot(s,sqOptions);      
		  }
#endif	       

		// suboptimal alignments
		if(options.has(RNAforesterOptions::LocalSubopts))
		{
			while(ali.nextLocalSuboptimum())
			{
				count++;
				cout << "local optimal score: ";
				cout << ali.getLocalOptimum() << endl;
				ali.getOptLocalAlignment(ppfali,xbasepos,ybasepos);
				cout << "starting at positions: " << xbasepos << "," << ybasepos << endl << endl; 

				// print alignment
				ppfali.getSequenceAlignments(seq1,seq2);
				ppfali.getStructureAlignment(str1,true);
				ppfali.getStructureAlignment(str2,false);
				RNAFuncs::printAli(ppfali.getStructureNameX(),ppfali.getStructureNameY(),seq1,seq2,str1,str2);  
				xlen=seq1.size();
				ylen=seq2.size();
				xsubopts.push_back(pair<Uint,Uint>(xbasepos,xlen));
				ysubopts.push_back(pair<Uint,Uint>(ybasepos,ylen));	

#ifdef HAVE_LIBG2  // This features require the g2 library		
		                if(options.has(RNAforesterOptions::MakeSquigglePlot))		  
				{
				  sprintf(s,"%d",count);				
				  ppfali.squigglePlot(s,sqOptions);      	    	  
				}
#endif			       
			}
		}

	}

#ifdef HAVE_LIBRNA  // This features require the ViennaRNA library	
#ifdef HAVE_LIBXMLPLUSPLUS
#ifdef HAVE_LIBXML2
	// generate xml
	if(options.has(RNAforesterOptions::XmlOutputFile)){
	  options.get(RNAforesterOptions::XmlOutputFile,outputFile,string(""));
	} else {
	  outputFile = "result.xml";
	}
	if(options.has(RNAforesterOptions::GenerateXML))
	{
	  RNAFuncs::printPAliXML(ppfali.getStructureNameX(),ppfali.getStructureNameY(),seq1,seq2,str1,str2,optScore,options,xmlInfos,outputFile);
	}
#endif	
#endif
#endif

#ifdef HAVE_LIBG2  // This features require the g2 library
	// generate squiggle plot
	if(options.has(RNAforesterOptions::MakeSquigglePlot))
	{
		f1->plot2d("x_str", xsubopts, sqOptions);
		f2->plot2d("y_str", ysubopts, sqOptions);  
	}
#endif	

	// clear input list
	deque<RNAForest*>::const_iterator it;
	for(it = inputListPW.begin(); it!=inputListPW.end(); it++)
		delete *it;
	inputListPW.clear();
	delete alg;
}


void editPairwise(list<RNAForestSZ*> &inputListSZ,Score &score,RNAforesterOptions &options)
{
  //  timeb t1,t2;
  IntScoreSZAlgebraType *alg;

//  if(options.has(RNAforesterOptions::CalculateDistance))
    alg=new IntDistSZAlgebra(score);
//  else
//    alg=new IntSimiSZAlgebra(score);
  
  RNAForestSZ *f1=inputListSZ.front();
  RNAForestSZ *f2=inputListSZ.back();

  //  ftime(&t1);
  Mapping<int,RNA_Alphabet> mapping(f1,f2,*alg);
  //  ftime(&t2);

//	f1->printParameters("F1");
//	f2->printParameters("F2");

  cout << "Global optimum: " << mapping.getGlobalOptimum() << endl;
  //  cout << "Calculation Time ms: " << (t2.time*1000+t2.millitm) - (t1.time*1000+t1.millitm) << endl;
}



void alignPairwiseSimple(deque<RNAForest*> &inputListPW,Score &score,RNAforesterOptions &options)
{
  //  timeb t1,t2;
  IntScore_AlgebraType *alg;

  //  if(options.has(RNAforesterOptions::CalculateDistance))
    alg=new IntDist_Algebra(score);
//  else
//    alg=new ScoreAlgebraSimple(score);

  RNAForest *f1=inputListPW.front();
  RNAForest *f2=inputListPW.back();

  //  ftime(&t1);
  Alignment<int,RNA_Alphabet,RNA_AlphaPair> ali(f1,f2,*alg,false);
  //  ftime(&t2);

//	f1->printParameters("F1");
  //f2->printParameters("F2");

  cout << "Global optimum: " << ali.getGlobalOptimum() << endl;
  //  cout << "Calculation Time ms: " << (t2.time*1000+t2.millitm) - (t1.time*1000+t1.millitm) << endl;
}

#ifdef HAVE_LIBXML2
#ifdef HAVE_LIBXMLPLUSPLUS
extern "C" {
  bool validateXSD(string filename){
    xmlSchemaParserCtxtPtr ctxt;
    xmlSchemaValidCtxtPtr valCtxt;
    xmlSchemaPtr schema;
    int val;

    ctxt = xmlSchemaNewMemParserCtxt(XSD_STRING,sizeof(XSD_STRING));
        
    schema = xmlSchemaParse(ctxt);
    
    valCtxt = xmlSchemaNewValidCtxt(schema);

    val = xmlSchemaValidateFile(valCtxt,filename.c_str(),0);

    if(val==0){
      return true;
    } else {
      return false;
    }
  } 
}
#endif
#endif

