/*****************************************************************************
*   Small Class to make Shape-Templated for MT2Analysis                      *
*****************************************************************************/

#ifndef MT2Shapes_HH
#define MT2Shapes_HH

#include "MT2tree.hh"
#include "helper/Utilities.hh"
#include "helper/Monitor.hh"
#include "THStack.h"
#include "TTree.h"
#include <map>
#include <sstream>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <math.h>
#include <vector>
#include <stdlib.h>
#include <stdio.h>

using namespace std;

//________________________________________________________________________________
class MT2Shapes  {

public:
	MT2Shapes();
	MT2Shapes(TString);
	MT2Shapes(TString, TString);
	MT2Shapes(TString, TString, std::ostringstream*);
	virtual ~MT2Shapes();


	void init(TString filename = "samples.dat");
	void loadSamples(const char* filename = "samples.dat");
	int  GetNSamples(){return fSamples.size();};
	int  GetNShapes(){return nShapes;};
	void SetPrintSummary(bool printsummary){fPrintSummary=printsummary;};
	void SetDraw(bool draw){fDraw=draw;};
	void SetWrite(bool write){fWrite=write;};
	void SetPileUpWeights(bool pileup){fDoPileUpWeights=pileup;};
	void SetbSFWeights(bool bSF){fbSFReWeight=bSF;};
	void SetLogStream(std::ostringstream* stream){fLogStream=stream;};
	void Print();

	struct sample{
		TString name;
		TString sname;
		TString type;
		TString shapename;
		TFile *file;
		TTree *tree;
		float xsection;
		float nevents;
		float kfact;
		float PU_avg_weight;
		float lumi;
		int color;
	};
	std::vector<sample>  fSamples;
	
	TH1D* fh_shapes[100];

	void setVerbose(int v){ fVerbose = v;};
	void setOutputDir(TString dir){ fOutputDir = Util::MakeOutputDir(dir); };
	void setOutputFile(TString filename){ fOutputFile = Util::MakeOutputFile(fOutputDir + filename); };
	void GetShapes(TString var, TString cuts, int njets, int nbjets, int nleps, TString selection_name, TString HLT,
			         TString xtitle, const int nbins, const double *bins);
	void GetShapes(TString var, TString cuts, int njets, int nbjets, int nleps, TString selection_name, TString HLT,
			         TString xtitle, const int nbins, const double min, const double max);
	void GetShapes(TString var, TString cuts, TString selection_name, TString HLT,
			         TString xtitle, const int nbins, const double *bins);
	void GetShapes(TString var, TString cuts, TString selection_name, TString HLT,
			         TString xtitle, const int nbins, const double min, const double max);
	int test;

private:

	TString fOutputDir;
	TFile *fOutputFile;
	int fVerbose;
	TString fPath;
	int nShapes;
	bool fPrintSummary;
	bool fDraw;
	bool fWrite;
	bool fCout;
	bool fDoPileUpWeights;
	bool fbSFReWeight;
	std::ostringstream* fLogStream;

	MT2tree* fMT2tree;
	TTree*   fTree;

	void DrawHisto(TH1* h_orig, TString canvname,  Option_t *drawopt);
	void FixOverAndUnderflowBins(TH1D* h, bool overflow=true, bool underflow=true);

};

#endif
