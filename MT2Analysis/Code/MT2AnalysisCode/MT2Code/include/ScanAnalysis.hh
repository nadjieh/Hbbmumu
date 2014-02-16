/*****************************************************************************
*   Small Class to make plots for MassAnalysis                               *
*****************************************************************************/

#ifndef ScanAnalysis_HH
#define ScanAnalysis_HH

#include "MT2tree.hh"
//#include "MassPlotter.hh"
#include "helper/Utilities.hh"
//#include "helper/Monitor.hh"
#include "THStack.h"
#include "TTree.h"
#include <map>
#include "TH2F.h"


using namespace std;

//________________________________________________________________________________
class ScanAnalysis {
  
public:
  ScanAnalysis();
  ScanAnalysis(TString);
  virtual ~ScanAnalysis();
  
  void init(TString filename = "samples.dat", TString xsecFilename="");
  void loadSamples(const char* filename = "samples.dat");
  void setVerbose(int v){ fVerbose = v;};
  void setOutputDir(TString dir){ fOutputDir = Util::MakeOutputDir(dir); };
  void setOutputFile(TString filename){ fOutputFile = fOutputDir + filename; };

  struct sample{
    TString name;
    TString sname;
    TString shapename;
    TString type;
    TFile *file;
    TTree *tree;
    float xsection;
    float nevents;
    float kfact;
    float lumi;
    int color;
  };
  std::vector<sample>  fSamples;
  
  void Analysis(TString outfile, TString filter="");
  
private:
  
  
  TString fOutputDir;
  TString fOutputFile;
  int fVerbose;
  TString fPath;
  
  MT2tree* fMT2tree;
  TTree*   fTree;
  
  map<int, TString> ProcessIDMap;

  bool isMT2Event();
  bool isMT2MuEvent();
  bool isMT2EleEvent();

  bool isMT2bEvent();
  bool isMT2bMuEvent();
  bool isMT2bEleEvent();
  
  map <  pair<float, float>, map<TString, float>  >  getXSec(TString filename);
  void createHistos(TString Label, map<TString, TH2F*> *h, int nHTBins, float HTBins[10], int nMT2Bins, float MT2Bins[10][10]);
  TString GetBinString(const int nHTBins, float HTBins[10], const int nMT2Bins, float MT2Bins[10][10]);

  map <  pair<float, float>, map<TString, float>  >  NLOXsecMap;
};

#endif
