#include "rna_algebra.h"

Score::Score(RNAforesterOptions &options)
{
  m_isLocal=options.has(RNAforesterOptions::LocalSimilarity);
  m_isRIBOSUM=options.has(RNAforesterOptions::RIBOSUMScore);

  // distance or similarity ?
  if(options.has(RNAforesterOptions::CalculateDistance))
    {
      m_isDistance=true;
      m_bp_rep_score = 0;
      m_bp_del_score = 3;
      m_b_match_score = 0;
      m_b_rep_score =1;
      m_b_del_score =2;      
    }
  else
    {
      m_isDistance=false;

	  if(options.has(RNAforesterOptions::RIBOSUMScore))
	  {
		  m_bp_del_score=-100;
		  m_b_del_score=-200;
	  }
	  else
	  {
		m_bp_rep_score = 10;
		m_bp_del_score =-5;
		m_b_match_score = 1;
		m_b_rep_score = 0;
                m_b_del_score =-10;
	  }
    }
  
  // read scores
  options.get(RNAforesterOptions::BpRepScore,m_bp_rep_score,m_bp_rep_score);
  options.get(RNAforesterOptions::BpDelScore,m_bp_del_score,m_bp_del_score);
  options.get(RNAforesterOptions::BMatchScore,m_b_match_score,m_b_match_score);
  options.get(RNAforesterOptions::BRepScore,m_b_rep_score,m_b_rep_score);
  options.get(RNAforesterOptions::BDelScore,m_b_del_score,m_b_del_score);
}

Score::Score(const Score &s)
{
  // copy scores
  m_isDistance = s.m_isDistance;
  m_isLocal = s.m_isLocal;
  m_bp_rep_score = s.m_bp_rep_score;
  m_bp_del_score = s.m_bp_del_score;
  m_b_match_score = s.m_b_match_score;
  m_b_rep_score = s.m_b_rep_score;
  m_b_del_score = s.m_b_del_score;   
}

void Score::print()
{
  // show score parameters
  cout << "*** Scoring parameters ***" << endl << endl; 

      cout << "Scoring type: ";
      if(m_isDistance)
	cout << "distance" << endl;
      else
	{
	  if(m_isLocal)
	    cout << "local ";
	  
	  cout << "similarity" << endl;	 
	}
           
      cout << "Scoring parameters:" << endl;

	  if(m_isRIBOSUM)
	  {
		  cout << "RIBOSUM85-60 Scoring matrix" << endl;
		  cout << "pd:  " << m_bp_del_score << endl;
		  cout << "bd:  " << m_b_del_score << endl << endl;
	  }
	  else
	  {
		  cout << "pm:   " << m_bp_rep_score << endl;
		  cout << "pd:   " << m_bp_del_score << endl;
		  cout << "bm:   " << m_b_match_score << endl;
		  cout << "br:   " << m_b_rep_score << endl;
		  cout << "bd:   " << m_b_del_score << endl << endl;
	  }
}
