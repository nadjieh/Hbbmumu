/*****************************************************************************
*   Small Class to make plots for MassAnalysis                               *
*****************************************************************************/

#ifndef ZInvEstFromW_HH
#define ZInvEstFromW_HH

#include "MT2tree.hh"
//#include "MassPlotter.hh"
#include "helper/Utilities.hh"
//#include "helper/Monitor.hh"
#include "THStack.h"
#include "TTree.h"
#include <map>


using namespace std;

//________________________________________________________________________________
class ZInvEstFromW {
  
public:
  ZInvEstFromW();
  ZInvEstFromW(TString);
  //ZInvEstFromW(TString, TString);
  virtual ~ZInvEstFromW();
  
  void init(TString filename = "samples.dat");
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

  //  void createCounterHistos(map<TString, TH1F*> *hm, map<string, string> cutLabels, TString sName);
  
  float*  Analysis(TString MT2_REGIME, bool IS_MC, TString outfile, int njets, int nele=0, int nmu=0,  float MT2CUT_LOW=-1, float MT2CUT_HIGH=999999999, TString trigger="", TString filter="", bool BENRICH=false, bool CORRECT_EFF=false, float* INPUT=NULL);
  void createHistos(map<TString, TH1F*> *histos, vector<TString> samplesName );
  void printBLeptEff(TString sample, map<TString, TH1F*> histos, map<TString, float> varFloat);

  

private:
  
  
  TString fOutputDir;
  TString fOutputFile;
  int fVerbose;
  TString fPath;
  
  MT2tree* fMT2tree;
  TTree*   fTree;
  
  
};

#endif
