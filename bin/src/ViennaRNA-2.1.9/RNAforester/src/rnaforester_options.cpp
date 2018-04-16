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

#include <iostream>
#include <sstream>

#ifndef WIN32
#include "config.h"
#endif

#ifndef HAVE_LIBRNA
#undef HAVE_LIBG2
#endif

#include "rnaforester_options.h"

RNAforesterOptions::RNAforesterOptions(int argc, const char **argv)
: m_argv(argv)
{
	nrArgs = argc;	
	m_options=new OptionsInfo[NumberOfOptions];

	setOption(Help,                      "--help","","                    ","shows this help info",false);
	setOption(Version,                   "--version","","                 ","shows version information",false);
	setOption(OutputAlignmentDotFormat,  "-dot","=file","                 ","show alignment forest in dot format",true);
	setOption(CalculateDistance,         "-d","","                        ","calculate distance instead of similarity",false);
	setOption(RelativeScore,             "-r","","                        ","calculate relative score",false);
	setOption(LocalSimilarity,           "-l","","                        ","local similarity",false);
	setOption(LocalSubopts,		     "-so","=int","                   ","local suboptimal alignments within int%",false);
	setOption(SmallInLarge,              "-s","","                        ","small-in-large similarity",false);
	setOption(Multiple,                  "-m","","                        ","multiple alignment mode",false);
	setOption(ClusterThreshold,          "-mt","=double","                ","clustering threshold",false);
	setOption(ClusterJoinCutoff,         "-mc","=double","                ","clustering cutoff",false);
#ifdef HAVE_LIBRNA  // This features require the ViennaRNA library
	setOption(PredictProfile,            "-p","","                        ","predict structures from sequences",false);
	setOption(PredictMinPairProb,	     "-pmin","=double","              ","minimum basepair frequency for prediction",false);
#endif
	setOption(SaveProfile,	             "-sp","=file","                  ","save profile",false);
	setOption(ProfileSearch,	     "-ps","=file","                  ","profile search",false);

	setOption(TreeEdit,                  "-e","","                        ","use tree edit model to calculate distance and similarity",true);
	setOption(GlobalAlignment,           "-g","","                        ","calculate global alignment in its original version",true);
	setOption(BpRepScore,                "-pm","=int","                   ","basepair(bond) match score",false);
	setOption(BpDelScore,                "-pd","=int","                   ","basepair bond indel score",false);
	setOption(BMatchScore,               "-bm","=int","                   ","base match score",false);
	setOption(BRepScore,                 "-br","=int","                   ","base mismatch score",false);
	setOption(BDelScore,                 "-bd","=int","                   ","base indel score",false);
	setOption(RIBOSUMScore,              "--RIBOSUM","","                 ","RIBOSUM85-60 scoring matrix (base-pair substitutions)",false);
	setOption(ConsensusMinPairProb,	     "-cmin","=double","              ","minimum basepair frequency for consensus structure",false);
#ifdef HAVE_LIBG2  // This features require the g2 library
	setOption(MakeSquigglePlot,          "-2d","","                       ","generate alignment 2D plots in postscript format",false);
	setOption(SquiggleHideBaseNumbers,   "--2d_hidebasenum","","          ","hide base numbers in 2D plot",false);
	setOption(SquiggleBaseNumberInterval,"--2d_basenuminterval","=n","    ","show every n-th base number",false);
	setOption(SquiggleGreyColors,        "--2d_grey","","                 ","use only grey colors in 2D plots",false);
	setOption(SquiggleScaleFactor,       "--2d_scale","=double","         ","scale factor for the 2d plots",false);
	setOption(SquiggleGenerateFIG,       "--2d_fig","","                  ","generate additional fig file of 2d plot",false);
#ifdef HAVE_LIBGD  // This features require the gd library
	setOption(SquiggleGeneratePNG,       "--2d_png","","                  ","generate additional png file of 2d plot",false);
	setOption(SquiggleGenerateJPG,       "--2d_jpg","","                  ","generate additional jpg file of 2d plot",false);
#endif
#endif
	setOption(ReadFromFile,              "-f","=file","                   ","read input from file",false);
	//  setOption(SaveMultipleAliFile,       "-sm","=file","                  ","save multiple alignment as binary file",true);
	setOption(NoScale,                   "--noscale","","                 ","suppress output of scale",false);
	setOption(MakeDotForInputTrees,      "-idot","","                     ","make dot files for the input trees",true);
#ifdef HAVE_LIBRNA  // This features require the ViennaRNA library and the libxml++ library  
#ifdef HAVE_LIBXMLPLUSPLUS
#ifdef HAVE_LIBXML2 	
	setOption(GenerateXML,               "--xml","","                     ","generate xml file in RNAStructAlignmentML format",false);
	setOption(XmlOutputFile,             "--xml_output","","                     ","name of xml output file",false);
#endif
#endif
#endif       
	setOption(SecretHelp,                "--shelp","","                   ","shows this help info",true);
	setOption(ShowOnlyScore,             "--score","","                   ","compute only scores, no alignment",false);
	setOption(FastaOutput,               "--fasta","","                   ","generate fasta output of alignments",false);

	setOption(SpaceTimeInfo,             "--spacetime","","               ","space and time measurements",true);

	// set the arguments that can be seperated by spaces
	stringstream ss;
	for(int i=0;i<NumberOfOptions;i++)
	{
		if(!m_options[i].parameter.empty())
			ss << m_options[i].tag << "|";
	}

	Arguments::setArgumentsWithSpaces(ss.str());
	m_args = new Arguments(argc,argv);

	// check options for compatibility
	//exclude(LocalSimilarity,Multiple);
	exclude(LocalSimilarity,SmallInLarge);
	//exclude(SmallInLarge,Multiple);

#ifdef HAVE_LIBG2  // This features require the g2 library	
	exclude(SquiggleHideBaseNumbers,Multiple);
	exclude(SquiggleBaseNumberInterval,Multiple);
#endif	
	exclude(CalculateDistance,LocalSimilarity);
	exclude(CalculateDistance,RIBOSUMScore);
	exclude(CalculateDistance,RelativeScore);
	exclude(CalculateDistance,Multiple);
	exclude(CalculateDistance,SmallInLarge);
	exclude(Multiple,RIBOSUMScore);
	exclude(RIBOSUMScore,BpRepScore);
	exclude(RIBOSUMScore,BMatchScore);
	exclude(RIBOSUMScore,BDelScore);
	
	requires(LocalSubopts,LocalSimilarity);
	requires(ClusterThreshold,Multiple);
	requires(ClusterJoinCutoff,Multiple);
#ifdef HAVE_LIBRNA  // This features require the ViennaRNA library	
	requires(PredictProfile,Multiple);
	requires(PredictMinPairProb,PredictProfile);
#endif	


	//	exclude(TreeEdit,Multiple);
	//	exclude(TreeEdit,LocalSimilarity);
	//	exclude(TreeEdit,RIBOSUMScore);
	//	exclude(GlobalAlignment,Multiple);
	//	exclude(GlobalAlignment,LocalSimilarity);
	//	exclude(GlobalAlignment,RIBOSUMScore);
	//	exclude(GlobalAlignment,TreeEdit);
}

RNAforesterOptions::~RNAforesterOptions()
{
	delete[] m_options;
	delete m_args;
}

inline void RNAforesterOptions::setOption(RNAforesterOption i,string tag, string parameter, string filler, string description, bool hidden)
{
	m_options[i].tag=tag;
	m_options[i].parameter=parameter;
	m_options[i].filler=filler;
	m_options[i].description=description;
	m_options[i].hidden=hidden;
}

bool RNAforesterOptions::has(RNAforesterOption option) const
{
	return m_args->has(m_options[option].tag);
}

const char** RNAforesterOptions::getArgs() const
{
	return m_argv;
}

const unsigned int RNAforesterOptions::getNrOfOptions() const
{
	return nrArgs;
}

void RNAforesterOptions::help()
{
	cout << "Usage: " << m_argv[0] << " [options]" << endl;

	for(int i=0;i<NumberOfOptions;i++)
	{
		if(!m_options[i].hidden)
			cout << m_options[i].tag << m_options[i].parameter << m_options[i].filler << m_options[i].description << endl;
	}
}

void RNAforesterOptions::secretHelp()
{
	cout << "Usage: " << m_argv[0] << " [options]" << endl;
	cout << "help including hidden parameters" << endl;

	for(int i=0;i<NumberOfOptions;i++)
	{
		cout << m_options[i].tag << m_options[i].parameter << m_options[i].filler << m_options[i].description;
		if(m_options[i].hidden)
			cout << " *";
		cout << endl;
	}
}

string RNAforesterOptions::generateFilename(RNAforesterOption option, const string &suffix, const string &defName, Uint count) const
{
	string s;

	get(option,s,string(""));
	if(s=="")
	{
		if(has(ReadFromFile))
		{
			ostringstream ss;

			get(ReadFromFile,s,defName);
			ss << s;
			if(count)
				ss << "_" << count;

			ss << suffix << '\0';

			s=ss.str();
		}
		else
		{
			s=defName;
		}
	}

	return s;
}

void RNAforesterOptions::exclude(RNAforesterOption opt1, RNAforesterOption opt2)
{
	if(has(opt1) && has(opt2))
		throw IncompatibleException(m_options[opt1].tag,m_options[opt2].tag);
}

void RNAforesterOptions::requires(RNAforesterOption opt1, RNAforesterOption opt2)
{
	if(has(opt1) && !has(opt2))
		throw RequiresException(m_options[opt1].tag,m_options[opt2].tag);
}

// ************************************************

void  RNAforesterOptions::IncompatibleException::showError()
{
	cerr << "The options " << m_tag1 << " and " << m_tag2 << " exclude each other." << endl;
}

void  RNAforesterOptions::RequiresException::showError()
{
  cerr << "The options " << m_tag1 << " requires option " << m_tag2 << "." <<  endl;
}
