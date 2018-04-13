#ifndef _RNAFORESTER_OPTIONS_H
#define _RNAFORESTER_OPTIONS_H

//#include <string>
#ifndef WIN32
#include "config.h"
#endif

#ifndef HAVE_LIBRNA
#undef HAVE_LIBG2
#endif

#include "Arguments.h"
#include "types.h"

class RNAforesterOptions
{
 private:
  struct OptionsInfo
  {
    string tag;
    string parameter;
    string filler;
    string description;
    bool hidden;
  };

 public:
  enum RNAforesterOption
  {
    Help=0,
    Version,
    OutputAlignmentDotFormat,
    CalculateDistance,
    RelativeScore,
    LocalSimilarity,
    LocalSubopts,
    SmallInLarge,
    Multiple,
    ClusterThreshold,
    ClusterJoinCutoff,
#ifdef HAVE_LIBRNA  // This features require the ViennaRNA library
    PredictProfile,
    PredictMinPairProb,
#endif
    SaveProfile,
    ProfileSearch,
    TreeEdit,
    GlobalAlignment,
    BpRepScore,
    BpDelScore,
    BMatchScore,
    BRepScore,
    BDelScore,
    RIBOSUMScore,
    ConsensusMinPairProb,
#ifdef HAVE_LIBG2  // This features require the g2 library
    MakeSquigglePlot,
    SquiggleHideBaseNumbers,
    SquiggleBaseNumberInterval,
    SquiggleGreyColors,
    SquiggleScaleFactor,
    SquiggleGenerateFIG,
#ifdef HAVE_LIBGD  // This features require the gd library
    SquiggleGeneratePNG,
    SquiggleGenerateJPG,
#endif
#endif
    ShowOnlyScore,
    FastaOutput,
    ReadFromFile,
    NoScale,
    MakeDotForInputTrees,
#ifdef HAVE_LIBRNA  // This features require the ViennaRNA library
#ifdef HAVE_LIBXMLPLUSPLUS
    GenerateXML,
    XmlOutputFile,
#endif
#endif
    SpaceTimeInfo,
    SecretHelp,
    NumberOfOptions  // this must always be the last entry in the enum
  };

  class IncompatibleException
    {
    private:
      string m_tag1;
      string m_tag2;

    public:
      IncompatibleException(string tag1,string tag2)
	: m_tag1(tag1), m_tag2(tag2) {};

      void showError();
    };

    class RequiresException
    {
    private:
      string m_tag1;
      string m_tag2;

    public:
      RequiresException(string tag1,string tag2)
	: m_tag1(tag1), m_tag2(tag2) {};

      void showError();
    };

  RNAforesterOptions(int argc, const char **argv);
  ~RNAforesterOptions();

  bool has(RNAforesterOption option) const;
  const char** getArgs() const;
  const unsigned int getNrOfOptions() const;
  
  template<class T>
  void get(RNAforesterOption option, T &var, T def) const
    {
      m_args->get(m_options[option].tag,var,def);
    }

  void help();
  void secretHelp();
  string generateFilename(RNAforesterOption option, const string &suffix, const string &defName, Uint count=0) const;

private:
  OptionsInfo *m_options;
  Arguments *m_args;
  const char **m_argv;
  unsigned int nrArgs;
  
  inline void setOption(RNAforesterOption,string tag, string parameter, string filler, string description, bool hidden);
  void exclude(RNAforesterOption opt1, RNAforesterOption opt2);
  void requires(RNAforesterOption opt1, RNAforesterOption opt2);

};

#endif
