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

#include <algorithm>
#include <fstream>

#include "alignment.h"
#include "progressive_align.h"

#include "alignment.t.cpp"

using namespace std;

// !!! this operator is defined as > !!!
bool operator < (pair<double,RNAProfileAlignment*> &l, pair<double,RNAProfileAlignment*> &r)
{
  if(l.second->getNumStructures() < r.second->getNumStructures())
    return false;
  else
    return true;
}

void progressiveAlign(deque<RNAProfileAlignment*> &inputList, deque<pair<double,RNAProfileAlignment*> > &resultList, const DoubleScoreProfileAlgebraType *alg,const RNAforesterOptions &options)
{
    RNAProfileAliMapType inputMapProfile;
	deque<RNAProfileAliKeyPairType> alignList;
	RNAProfileAliMapType::iterator it,it2;
	long x,y,bestx,besty;  // keys of stuctures in inputMapProfile
	RNAProfileAlignment *ppfx=NULL,*ppfy=NULL,*ppf=NULL;
	int level=1;
	double bestScore,threshold,cutoff;
	string clusterfilename;
	ofstream s;
	Matrix<double> *score_mtrx;	
	bool local;

	local=options.has(RNAforesterOptions::LocalSimilarity);

	cout << "*** Calculation ***" << endl << endl;

	// create inputMapProfile to access a profile by an index value
	deque<RNAProfileAlignment*>::const_iterator inpIt;
	long i=1;
	for(inpIt=inputList.begin();inpIt!=inputList.end();inpIt++)
	  {
	    inputMapProfile[i]=*inpIt;
	    i++;
	  }
	inputList.clear();

	// create matrix for all against all comparison
	bestScore=alg->worst_score();
	score_mtrx=new Matrix<double>(inputMapProfile.size(),inputMapProfile.size());

	// set threshold for the clustering algorithm
	if(options.has(RNAforesterOptions::CalculateDistance))
		options.get(RNAforesterOptions::ClusterThreshold,threshold,20.0);
	else
		options.get(RNAforesterOptions::ClusterThreshold,threshold,0.7);
	cout << "clustering threshold is: " << threshold << endl;

	// set cutoff value for clustering
	if(options.has(RNAforesterOptions::CalculateDistance))
	        options.get(RNAforesterOptions::ClusterJoinCutoff,cutoff,100.0);
	else
	        options.get(RNAforesterOptions::ClusterJoinCutoff,cutoff,0.0);
	cout << "join clusters cutoff is: " << cutoff << endl << endl;

	// generate dot file
	clusterfilename=options.generateFilename(RNAforesterOptions::Help,"_cluster.dot", "cluster.dot");  // use Help as dummy
	s.open(clusterfilename.c_str());
	s << "digraph forest" << endl << "{" << endl;

	// generate nodes for the input forests
	for(it=inputMapProfile.begin();it!=inputMapProfile.end();it++) 
	{
		ppf=it->second;

		s << "\"" << ppf->getName() << "\"" << "[label=\"" << ppf->getName() << "\"]" << endl;
		s << "\"" << ppf->getName() << "\"" << "[label=\"" << ppf->getName() << "\"]" << endl;
	}


	// compute all pairwise distances
	// !! NOTE !! iterating through the map is ordered by key number
	// as i only calculate a triangle matrix this is a prerequisite
	cout << "Computing all pairwise similarities" << endl;
	for(it=inputMapProfile.begin();it!=inputMapProfile.end();it++)
	{	  
		x=it->first;
		ppfx=it->second;
		for(it2=inputMapProfile.begin();it2->first<it->first;it2++)
		{
			y=it2->first;
			ppfy=it2->second;

			Alignment<double,RNA_Alphabet_Profile,RNA_Alphabet_Profile> ali(ppfx,ppfy,*alg,local);
			if(local)
			  score_mtrx->setAt(x-1,y-1,ali.getLocalOptimum());
			else
			  score_mtrx->setAt(x-1,y-1,ali.getGlobalOptimumRelative());
			
			cout << x << "," << y << ": " << score_mtrx->getAt(x-1,y-1) << endl;

			WATCH(DBG_MULTIPLE,"progressiveAlign",x);
			WATCH(DBG_MULTIPLE,"progressiveAlign",y);
			WATCH(DBG_MULTIPLE,"progressiveAlign",ali.getGlobalOptimumRelative());
		}
	}
	cout << endl;

	while(inputMapProfile.size()>1)
	{
		// find the best score of all pairwise alignments
		bestScore=alg->worst_score();
		for(it=inputMapProfile.begin();it!=inputMapProfile.end();it++)
		{
			x=it->first;
			for(it2=inputMapProfile.begin();it2->first < it->first;it2++)
			{
				double old_bestScore=bestScore;

				y=it2->first;	      
				WATCH(DBG_MULTIPLE,"progressiveAlign",x);
				WATCH(DBG_MULTIPLE,"progressiveAlign",y);
				WATCH(DBG_MULTIPLE,"progressiveAlign",score_mtrx->getAt(x-1,y-1));

				WATCH(DBG_MULTIPLE,"progressiveAlign",bestScore);
				bestScore=alg->choice(bestScore,score_mtrx->getAt(x-1,y-1));

				if(bestScore!=old_bestScore)
				{
					bestx=it->first;
					besty=it2->first;
					old_bestScore=bestScore;
				}
			}
		}

		WATCH(DBG_MULTIPLE,"progressiveAlign",bestScore);
		cout << "joining alignments:" << endl;

		// if threshold is set generate a list of best pairs within threshold 
		if(alg->choice(bestScore,threshold)!= threshold)
		{
			Graph graph;
			int *mate;
			int maximize;

			graph=makePairsGraph(inputMapProfile,alg,score_mtrx,threshold);
			if(options.has(RNAforesterOptions::CalculateDistance))
				maximize=0;
			else
				maximize=1;

			mate = Weighted_Match(graph,1,maximize);
			for(x=1;x<=score_mtrx->xDim();x++)  // !! begins at 1 !!
			{
				if(x<mate[x])
				{		
					// if it is a best pair put it in the align list
					ppfx=inputMapProfile[x];
					inputMapProfile.erase(x);
					alignList.push_back(make_pair(x,ppfx));
					ppfy=inputMapProfile[mate[x]];
					inputMapProfile.erase(mate[x]);
					alignList.push_back(make_pair(mate[x],ppfy));
				}
			}
			free(mate);
			free(graph);
		}
		else
		{
			// if there us no pair below the threshold
			// combine those two profile forests, that produced the best score
			ppfx=inputMapProfile[bestx];
			ppfy=inputMapProfile[besty];
			inputMapProfile.erase(bestx);
			inputMapProfile.erase(besty);
			alignList.push_back(make_pair(bestx,ppfx));
			alignList.push_back(make_pair(besty,ppfy));
		}

		// align all forests in the align list. 
		while(alignList.size()>1)
		{
			string aliName;

			x=alignList.front().first;
			ppfx=alignList.front().second;
			alignList.erase(alignList.begin());
			y=alignList.front().first;
			ppfy=alignList.front().second;
			alignList.erase(alignList.begin());

			// compute alignment again 
			Alignment<double,RNA_Alphabet_Profile,RNA_Alphabet_Profile> bestali(ppfx,ppfy,*alg,local);
			if(local)
			  bestScore=bestali.getLocalOptimum();
			else
			  bestScore=bestali.getGlobalOptimumRelative();

			// test, if score is worse than the cutoff value
			if(alg->choice(bestScore,cutoff)== cutoff)
			{
				// copy involved forests to the result list
				cout << x << "," << y << ": alignment is below cutoff." << endl; 
				if(ppfx->getNumStructures()>1)
				  resultList.push_back(make_pair(ppfx->maxScore(*alg),ppfx));
				if(ppfy->getNumStructures()>1)
				  resultList.push_back(make_pair(ppfy->maxScore(*alg),ppfy));							  
			}
			else
			{
			        // calculate optimal alignment and add it to inputMapProfile

			  ppf=new RNAProfileAlignment(ppfx->getNumStructures(),ppfy->getNumStructures());

			  if(local)
			    {
			      Uint xbasepos,ybasepos;
			      bestali.getOptLocalAlignment(*ppf,xbasepos,ybasepos);
			    }
			  else
			    {
			      bestali.getOptGlobalAlignment(*ppf);
			    }
			
				ppf->addStrNames(ppfx->getStrNames());
				ppf->addStrNames(ppfy->getStrNames());
				aliName=ppfx->getName() + "." +ppfy->getName();
				ppf->setName(aliName);

				// for debug purposes
				//string dotfilename;
				//dotfilename=ppf->getName() + string(".dot");
				//				ofstream s("test.dot");
				//				ppf->printDot(s);

				// generate nodes in cluster file (dot format)
				s << "\"" << ppf->getName() << "\"" << "[shape=\"diamond\" label=\"" << bestScore << "\"]" << endl;
				s << "\"" << ppf->getName() << "\"" << "-> {\"" <<  ppfx->getName() << "\" \"" << ppfy->getName() << "\"}";
				
				cout << x << "," << y << ": " << bestScore << " -> " << min(x,y) << endl;
				
				delete ppfx;
				delete ppfy;      

				// calculate distance to all forests in the list
				cout << "Calculate similarities to other clusters" << endl;
				ppfx=ppf;
				// x remains x !!
				for(it=inputMapProfile.begin();it!=inputMapProfile.end();it++)
				  {
				    y=it->first;
				    ppfy=it->second;
				    Alignment<double,RNA_Alphabet_Profile,RNA_Alphabet_Profile> ali(ppfx,ppfy,*alg,local);
				    if(local)
				      {
					score_mtrx->setAt(min(x-1,y-1),max(x-1,y-1),ali.getLocalOptimum());  // min - max = fill the upper triangle
					cout << min(x,y) << "," << max(x,y) << ": " << ali.getLocalOptimum() <<  endl;					
				      }
				    else
				      {
					score_mtrx->setAt(min(x-1,y-1),max(x-1,y-1),ali.getGlobalOptimumRelative());  // min - max = fill the upper triangle
					cout << min(x,y) << "," << max(x,y) << ": " << ali.getGlobalOptimumRelative() <<  endl;
				      }
				  } 
				cout << endl;
				
				// ... and append it to the list
				inputMapProfile.insert(make_pair(x,ppf));
			}
		}

		level++;
	}

	assert(inputMapProfile.size()<2);

	// copy last profile to resultList
	if(inputMapProfile.size()==1)
	{
	    ppf=inputMapProfile.begin()->second;
		resultList.push_back(make_pair(ppf->maxScore(*alg),ppf));
		inputMapProfile.clear();
	}

	// end of dot file
	s << "}" << endl;

	// sort result list
	//	std::sort(resultList.begin(),resultList.end());
	
	
	delete score_mtrx;
}

Graph makePairsGraph(const RNAProfileAliMapType &inputMapProfile, const DoubleScoreProfileAlgebraType *alg, const Matrix<double> *score_mtrx, double threshold)
{
	Graph graph;
	RNAProfileAliMapType::const_iterator it,it2;
	RNAProfileAlignment *ppfx=NULL,*ppfy=NULL;

	graph = NewGraph(score_mtrx->xDim());

	for(int i=0;i<score_mtrx->xDim();i++)
	{
		Xcoord(graph,i+1) = 0;
		Ycoord(graph,i+1) = 0;
	}

	for(it=inputMapProfile.begin();it!=inputMapProfile.end();it++)
	{
		ppfx=it->second;
		for(it2=inputMapProfile.begin();it2->first<it->first;it2++)
		{
			double score;
			ppfy=it2->second;

			score=score_mtrx->getAt(it->first-1,it2->first-1);
			if(alg->choice(score,threshold) != threshold)  // is it better than the threshold ?
			{
				AddEdge (graph,it->first,it2->first,(int)(score*100.0));      
			}
		}
	}


	WriteGraph (graph,"test.out");
	return graph;
}
