#include "ScanAnalysis.hh"

#include "TLorentzVector.h"
#include "helper/Utilities.hh"

#include "TString.h"
#include "TObject.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TDirectoryFile.h"
#include "TH2.h"
#include "TH2D.h"
#include "TH3.h"
#include "TTree.h"
#include "TEventList.h"
#include "TCut.h"
#include "TTreeFormula.h"
#include "TRandom.h"
#include "TROOT.h"

#include <sstream>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <math.h>
#include <time.h> // access to date/time

#include "helper/Hemisphere.hh"


//utils
#ifndef Utilities_HH
#include "Utilities.hh"
#endif

using namespace std;

//____________________________________________________________________________
ScanAnalysis::ScanAnalysis(){
  // Default constructor, no samples are set
}

//____________________________________________________________________________

ScanAnalysis::ScanAnalysis(TString outputdir/*, TString outputfile*/){
  // Explicit constructor with output directory and output file
  setOutputDir(outputdir);
  //setOutputFile(outputfile);
}

//____________________________________________________________________________
ScanAnalysis::~ScanAnalysis(){
  //fOutputFile->Close();
  //delete fOutputFile;
}

//____________________________________________________________________________

void ScanAnalysis::init(TString filename, TString xsecFilename){
  if(fVerbose > 0) cout << "------------------------------------" << endl;
  if(fVerbose > 0) cout << "Initializing MassPlotter ... " << endl;
  //Util::SetStyle();
  loadSamples(filename);

  //defining process map
  ProcessIDMap[1] = "ng";
  ProcessIDMap[2] = "ns";
  ProcessIDMap[3] = "nn";
  ProcessIDMap[4] = "ll";
  ProcessIDMap[5] = "sb";
  ProcessIDMap[6] = "ss";
  ProcessIDMap[7] = "tb";
  ProcessIDMap[8] = "bb";
  ProcessIDMap[9] = "gg";
  ProcessIDMap[10] = "sg";

  //getting xsecs
  if(xsecFilename!=""){
    NLOXsecMap =  getXSec(xsecFilename);
  }


 }


void ScanAnalysis::Analysis( TString outputfile,  TString filter){
  
  setOutputFile(outputfile);
  TFile *fout = TFile::Open(fOutputFile,"RECREATE");
  
  //FIXME number of events per point, to be fixed
  int TotalNumEventsPoint = 10000.;
  
  // bin definition
  const int nHTBins = 2;
  float HTBins[10];
  HTBins[0] =  750; HTBins[1] = 950; HTBins[2] = 9999.;
  const int nMT2Bins = 5, nMT2bBins = 4;
  float MT2Bins[10][10], MT2bBins[10][10];
  
  MT2Bins[0][0] = 150;  MT2Bins[0][1] =  200;  MT2Bins[0][2] = 275;  MT2Bins[0][3] = 375;  MT2Bins[0][4] =  500;  MT2Bins[0][5] =  9999;
  MT2Bins[1][0] = 150;  MT2Bins[1][1] =  200;  MT2Bins[1][2] = 275;  MT2Bins[1][3] = 375;  MT2Bins[1][4] =  500;  MT2Bins[1][5] =  9999;

  MT2bBins[0][0] = 125;  MT2bBins[0][1] =  150;  MT2bBins[0][2] = 200;  MT2bBins[0][3] = 300;  MT2bBins[0][4] =  9999;
  MT2bBins[1][0] = 125;  MT2bBins[1][1] =  150;  MT2bBins[1][2] = 180;  MT2bBins[1][3] = 260;  MT2bBins[1][4] =  9999;
  
  //create histos
  map<TString, TH2F*> histos;
  createHistos("MT2", &histos, nHTBins, HTBins, nMT2Bins, MT2Bins);
  createHistos("MT2Ele", &histos, nHTBins, HTBins, nMT2Bins, MT2Bins);
  createHistos("MT2Mu", &histos, nHTBins, HTBins, nMT2Bins, MT2Bins);
  createHistos("MT2b", &histos, nHTBins, HTBins, nMT2bBins, MT2bBins);
  createHistos("MT2bEle", &histos, nHTBins, HTBins, nMT2bBins, MT2bBins);
  createHistos("MT2bMu", &histos, nHTBins, HTBins, nMT2bBins, MT2bBins);

  for(map<TString,TH2F*>::iterator h= histos.begin(); h!=histos.end();h++)    h->second->Sumw2();

  ////// Samples loop
  for(size_t i = 0; i < fSamples.size(); ++i){
    if( fSamples[i].tree->GetEntries()==0 ) continue;

    TString tsname = fSamples[i].sname;
    TString type = fSamples[i].type;
    
    fMT2tree = new MT2tree();
    fSamples[i].tree->SetBranchAddress("MT2tree", &fMT2tree);
    
    //Filtering TTree
    TString myCuts = filter;
    fSamples[i].tree->Draw(">>selList", myCuts);
    TEventList *myEvtList = (TEventList*)gDirectory->Get("selList");
    fSamples[i].tree->SetEventList(myEvtList);
    int counter=0;
    Long64_t nbytes = 0, nb = 0;
    if(fVerbose>2) cout <<  "\t\t Filtering done, size=" <<myEvtList->GetN()  << endl;


    ///////////analysis
    if(myEvtList->GetN() ==0) continue;
    while(myEvtList->GetEntry(counter++) !=-1){
      int jentry = myEvtList->GetEntry(counter-1);
      nb =  fSamples[i].tree->GetEntry(jentry);   nbytes += nb;
      
      float HT = fMT2tree->misc.HT;
      float MT2 = fMT2tree->misc.MT2;
      
      float M0 = fMT2tree->Susy.M0;
      float M12 = fMT2tree->Susy.M12;
      pair< float, float > susyPoint; susyPoint.first = M0; susyPoint.second = M12;
      
      //subprocess label
      TString procLab = ProcessIDMap[fMT2tree->GenProcessID];
      
      //xsecs are in pb I suppose, so weight is normalized to 1/pb
      float weight =  NLOXsecMap[susyPoint][procLab]/(float)TotalNumEventsPoint;
      
      //event selection
      bool isMT2 = isMT2Event();
      bool isMT2Ele = isMT2EleEvent();
      bool isMT2Mu = isMT2MuEvent();

      bool isMT2b = isMT2bEvent();
      bool isMT2bEle = isMT2bEleEvent();
      bool isMT2bMu = isMT2bMuEvent();


      //Signal regions
      if(isMT2){
	TString binString = GetBinString(nHTBins, HTBins, nMT2Bins,MT2Bins);
	//cout << "MT2-"+binString << " " << HT << " " << MT2 <<endl;
	histos["MT2-"+binString]->Fill(M0, M12, weight);
      }
      if(isMT2b){
	TString binString = GetBinString(nHTBins, HTBins, nMT2bBins,MT2bBins);
	histos["MT2b-"+binString]->Fill(M0, M12, weight);
      }

      //Control regions
      if(isMT2Ele){
	TString binString = GetBinString(nHTBins, HTBins, nMT2Bins,MT2Bins);
	histos["MT2Ele-"+binString]->Fill(M0, M12, weight);
      } 
      if(isMT2Mu){
	TString binString = GetBinString(nHTBins, HTBins, nMT2Bins,MT2Bins);
	histos["MT2Mu-"+binString]->Fill(M0, M12, weight);
      } 
      if(isMT2bEle){
	TString binString = GetBinString(nHTBins, HTBins, nMT2bBins,MT2bBins);
	histos["MT2bEle-"+binString]->Fill(M0, M12, weight);
      }
      if(isMT2bMu){
	TString binString = GetBinString(nHTBins, HTBins, nMT2bBins,MT2bBins);
	histos["MT2bMu-"+binString]->Fill(M0, M12, weight);
      }
      
    }//end event loop
  }//end samples loop
  
  fout->Write();
  fout->Close();
  
  cout << "Output file in " << outputfile << endl;
  
}


bool ScanAnalysis::isMT2Event(){
  bool passed = true;
  
  if(!(fMT2tree->NJetsIDLoose40>2 )) passed=false;
  if(!(fMT2tree->misc.MET>30 )) passed=false;
  if(!(fMT2tree->misc.Jet0Pass==1 )) passed=false;
  if(!(fMT2tree->misc.Jet1Pass==1 )) passed=false;
  if(!(fMT2tree->misc.PassJetID==1 )) passed=false;
  if(!(fMT2tree->misc.Vectorsumpt <70 )) passed=false;
  if(!(fMT2tree->misc.MinMetJetDPhi > 0.3 )) passed=false;
  if(!(fMT2tree->misc.SecondJPt > 100 )) passed=false;
  if(!(fMT2tree->misc.HT > 750 )) passed=false;
  if(!(fMT2tree->misc.MT2 > 150 )) passed=false;
  if(!(fMT2tree->NEles==0 || fMT2tree->ele[0].lv.Pt()<10)) passed=false;
  if(!(fMT2tree->NMuons==0 || fMT2tree->muo[0].lv.Pt()<10)) passed=false;
  
  return passed;
}

bool ScanAnalysis::isMT2EleEvent(){
  bool passed = true;
  
  if(!(fMT2tree->NJetsIDLoose40>2 )) passed=false;
  if(!(fMT2tree->misc.MET>30 )) passed=false;
  if(!(fMT2tree->misc.Jet0Pass==1 )) passed=false;
  if(!(fMT2tree->misc.Jet1Pass==1 )) passed=false;
  if(!(fMT2tree->misc.PassJetID==1 )) passed=false;
  if(!(fMT2tree->misc.Vectorsumpt <70 )) passed=false;
  if(!(fMT2tree->misc.MinMetJetDPhi > 0.3 )) passed=false;
  if(!(fMT2tree->misc.SecondJPt > 100 )) passed=false;
  if(!(fMT2tree->misc.HT > 750 )) passed=false;
  if(!(fMT2tree->misc.MT2 > 150 )) passed=false;
  if(!(fMT2tree->NEles==1 && fMT2tree->ele[0].lv.Pt()>10)) passed=false;
  if(!(fMT2tree->ele[0].MT < 100)) passed=false;

  if(!(fMT2tree->NMuons==0 || fMT2tree->muo[0].lv.Pt()<10)) passed=false;
  
  return passed;
}

bool ScanAnalysis::isMT2MuEvent(){
  bool passed = true;
  
  if(!(fMT2tree->NJetsIDLoose40>2 )) passed=false;
  if(!(fMT2tree->misc.MET>30 )) passed=false;
  if(!(fMT2tree->misc.Jet0Pass==1 )) passed=false;
  if(!(fMT2tree->misc.Jet1Pass==1 )) passed=false;
  if(!(fMT2tree->misc.PassJetID==1 )) passed=false;
  if(!(fMT2tree->misc.Vectorsumpt <70 )) passed=false;
  if(!(fMT2tree->misc.MinMetJetDPhi > 0.3 )) passed=false;
  if(!(fMT2tree->misc.SecondJPt > 100 )) passed=false;
  if(!(fMT2tree->misc.HT > 750 )) passed=false;
  if(!(fMT2tree->misc.MT2 > 150 )) passed=false;
  if(!(fMT2tree->NEles==0 || fMT2tree->ele[0].lv.Pt()<10)) passed=false;
  if(!(fMT2tree->NMuons==1 && fMT2tree->muo[0].lv.Pt()>10)) passed=false;
  if(!(fMT2tree->muo[0].MT < 100)) passed=false;
 
  return passed;
}

bool ScanAnalysis::isMT2bEvent(){
  bool passed = true;
  
  if(!(fMT2tree->NJetsIDLoose40>3 )) passed=false;
  if(!(fMT2tree->misc.MET>30 )) passed=false;
  if(!(fMT2tree->misc.Jet0Pass==1 )) passed=false;
  if(!(fMT2tree->misc.Jet1Pass==1 )) passed=false;
  if(!(fMT2tree->misc.PassJetID==1 )) passed=false;
  if(!(fMT2tree->misc.Vectorsumpt <70 )) passed=false;
  if(!(fMT2tree->misc.MinMetJetDPhi > 0.3 ||  
       fMT2tree->misc.MinMetJetDPhiIndex>3 )) passed=false;
  if(!(fMT2tree->misc.LeadingJPt > 150 )) passed=false;
  if(!(fMT2tree->misc.SecondJPt > 100 )) passed=false;
  if(!(fMT2tree->NBJets>0 )) passed=false;
  if(!(fMT2tree->misc.HT > 750 )) passed=false;
  if(!(fMT2tree->misc.MT2 > 125 )) passed=false;
  if(!(fMT2tree->NEles==0 || fMT2tree->ele[0].lv.Pt()<10)) passed=false;
  if(!(fMT2tree->NMuons==0 || fMT2tree->muo[0].lv.Pt()<10)) passed=false;
  
  return passed;
}


bool ScanAnalysis::isMT2bEleEvent(){
  bool passed = true;
  
  if(!(fMT2tree->NJetsIDLoose40>3 )) passed=false;
  if(!(fMT2tree->misc.MET>30 )) passed=false;
  if(!(fMT2tree->misc.Jet0Pass==1 )) passed=false;
  if(!(fMT2tree->misc.Jet1Pass==1 )) passed=false;
  if(!(fMT2tree->misc.PassJetID==1 )) passed=false;
  if(!(fMT2tree->misc.Vectorsumpt <70 )) passed=false;
  if(!(fMT2tree->misc.MinMetJetDPhi > 0.3 ||  
       fMT2tree->misc.MinMetJetDPhiIndex>3 )) passed=false;
  if(!(fMT2tree->misc.LeadingJPt > 150 )) passed=false;
  if(!(fMT2tree->misc.SecondJPt > 100 )) passed=false;
  if(!(fMT2tree->NBJets>0 )) passed=false;
  if(!(fMT2tree->misc.HT > 750 )) passed=false;
  if(!(fMT2tree->misc.MT2 > 125 )) passed=false;
  if(!(fMT2tree->NEles==1 && fMT2tree->ele[0].lv.Pt()>10)) passed=false;
  if(!(fMT2tree->NMuons==0 || fMT2tree->muo[0].lv.Pt()<10)) passed=false;
  if(!(fMT2tree->ele[0].MT < 100)) passed=false;

  return passed;
}

bool ScanAnalysis::isMT2bMuEvent(){
  bool passed = true;
  
  if(!(fMT2tree->NJetsIDLoose40>3 )) passed=false;
  if(!(fMT2tree->misc.MET>30 )) passed=false;
  if(!(fMT2tree->misc.Jet0Pass==1 )) passed=false;
  if(!(fMT2tree->misc.Jet1Pass==1 )) passed=false;
  if(!(fMT2tree->misc.PassJetID==1 )) passed=false;
  if(!(fMT2tree->misc.Vectorsumpt <70 )) passed=false;
  if(!(fMT2tree->misc.MinMetJetDPhi > 0.3 ||  
       fMT2tree->misc.MinMetJetDPhiIndex>3 )) passed=false;
  if(!(fMT2tree->misc.LeadingJPt > 150 )) passed=false;
  if(!(fMT2tree->misc.SecondJPt > 100 )) passed=false;
  if(!(fMT2tree->NBJets>0 )) passed=false;
  if(!(fMT2tree->misc.HT > 750 )) passed=false;
  if(!(fMT2tree->misc.MT2 > 125 )) passed=false;
  if(!(fMT2tree->NEles==0 || fMT2tree->ele[0].lv.Pt()<10)) passed=false;
  if(!(fMT2tree->NMuons==1 && fMT2tree->muo[0].lv.Pt()>10)) passed=false;
  if(!(fMT2tree->muo[0].MT < 100)) passed=false;

  return passed;
}

//get string for MT2Bin
TString ScanAnalysis::GetBinString(int nHTBins, float HTBins[10], int nMT2Bins, float MT2Bins[10][10]){
  TString binString="";
  float HT = fMT2tree->misc.HT;
  float MT2 = fMT2tree->misc.MT2;
  for(int HT_i=0; HT_i<nHTBins; HT_i++){
    float HT_lo = HTBins[HT_i];
    float HT_hi = HTBins[HT_i+1];
    bool passedHT = false;
    if( HT>HT_lo && HT<=HT_hi) passedHT=true;
    if(!passedHT) continue;
    
    for(int MT2_i=0; MT2_i<nMT2Bins; MT2_i++){
      float MT2_lo = MT2Bins[HT_i][MT2_i];
      float MT2_hi = MT2Bins[HT_i][MT2_i+1];
      if(MT2>MT2_lo && MT2<=MT2_hi){
	ostringstream labStream;
	labStream << "HT_" << HT_lo << "to"<<HT_hi << "-MT2_"<<MT2_lo << "to"<< MT2_hi;
	binString = labStream.str();
      }
    }
    
  }
  return binString;
}



//histograms
void ScanAnalysis::createHistos(TString Label, map<TString, TH2F*> *histos, int nHTBins, float HTBins[10], int nMT2Bins, float MT2Bins[10][10]){
  TString binString="";
  for(int HT_i=0; HT_i<nHTBins; HT_i++){
    float HT_lo = HTBins[HT_i];
    float HT_hi = HTBins[HT_i+1];
    for(int MT2_i=0; MT2_i<nMT2Bins; MT2_i++){
      float MT2_lo = MT2Bins[HT_i][MT2_i];
      float MT2_hi = MT2Bins[HT_i][MT2_i+1];
      
      ostringstream labStream;
      labStream << "HT_" << HT_lo << "to"<<HT_hi << "-MT2_"<<MT2_lo << "to"<< MT2_hi;
      binString = labStream.str();

      TString hName = Label+"-"+binString;
      (*histos)[hName] = new TH2F(hName,"",300,0,3000,100,0,1000); 
    }
  }
  
 
}





map <  pair<float, float>, map<TString, float>  >  ScanAnalysis::getXSec(TString filename){
  map <  pair<float, float>, map<TString, float>  > xsec;
  std::ifstream inFile(filename);
  std::string sLine;
  while(std::getline(inFile, sLine)) {
    //Do here whatever you need to do
    int init_m0 = sLine.find("m0", 0);
    if(init_m0 == string::npos) continue;
    
    //file structure
    //|(scale=1.0) m0=1000, m1/2=100, tanbeta=3, A0=0, sign(mu)=+ve| 0.224611 | 0.0101205 | 53.1719 | 3.5998e-06 | 0.00592 | 0.0386 | 0.024632 | 0.001203 | 66.1 | 2.3 |
    //| Sub-processes | ng | ns | nn | ll | sb | ss | tb | bb | gg | sg |
    
    int m0, m12, tanb;
    string str, str2;
    float ng,ns,nn,ll,sb,ss,tb,bb,gg,sg;
    sscanf(sLine.c_str(), " |%*s m0=%*d, m1/2=%*d, %*s A0=0, sign(mu)=+ve| %f | %e | %e | %e | %e | %e | %e | %e | %e | %e |",
	   &ng,&ns,&nn,&ll,&sb,&ss,&tb,&bb,&gg,&sg);
    sscanf(sLine.c_str(), " |%*s m0=%d, m1/2=%d, %*s |",
	   &m0, &m12);
    
    pair< float, float > susyPoint;
    map< TString, float > subXsec;
    susyPoint.first = m0;         susyPoint.second = m12;
    subXsec["ng"] = ng;
    subXsec["ns"] = ns;
    subXsec["nn"] = nn;
    subXsec["ll"] = ll;
    subXsec["sb"] = sb;
    subXsec["ss"] = ss;
    subXsec["tb"] = tb;
    subXsec["bb"] = bb;
    subXsec["gg"] = gg;
    subXsec["sg"] = sg;
    xsec[ susyPoint ] = subXsec;
  }
  inFile.close();
  
  //access example
  pair< float, float > susyPoint; susyPoint.first = 990; susyPoint.second = 500;
  map< TString, float > subXsec = xsec[susyPoint];
  return xsec;
}




//____________________________________________________________________________
void ScanAnalysis::loadSamples(const char* filename){
  fSamples.clear();
  char buffer[200];
  ifstream IN(filename);
  
  char /*ParName[100],*/ StringValue[1000];
  float ParValue;

  if(fVerbose > 0) cout << "------------------------------------" << endl;
  if(fVerbose > 0) cout << "Sample File  " << filename << endl;
  int counter(0);

  while( IN.getline(buffer, 200, '\n') ){
    // ok = false;
    if (buffer[0] == '#') {
      continue; // Skip lines commented with '#'
    }
    if( !strcmp(buffer, "GENERAL") ){
      IN.getline(buffer, 200, '\n');
      sscanf(buffer, "Path\t%s", StringValue);
      fPath = StringValue;
      cout << fPath << endl;

      if(fVerbose >0){
	cout << " ----  " << endl;
	cout << "  Path " << fPath << endl;
      }

    }
    if( !strcmp(buffer, "SAMPLE")){

      sample s;
      IN.getline(buffer, 200, '\n');
      sscanf(buffer, "Name\t%s", StringValue);
      s.name = TString(StringValue);

      IN.getline(buffer, 200, '\n');
      sscanf(buffer, "SName\t%s", StringValue);
      s.sname = TString(StringValue);

      IN.getline(buffer, 200, '\n');
      sscanf(buffer, "ShapeName\t%s", StringValue);
      s.shapename = TString(StringValue);

      IN.getline(buffer, 200, '\n');
      sscanf(buffer, "File\t%s", StringValue);
      TString file =fPath+StringValue;
      TFile *f = TFile::Open(file);
      s.file = f;
      s.tree = (TTree*)f->Get("MassTree");

      IN.getline(buffer, 200, '\n');
      sscanf(buffer, "Xsection\t%f", &ParValue);
      s.xsection = ParValue;

      /*
      IN.getline(buffer, 200, '\n');
      sscanf(buffer, "Nevents\t%f", &ParValue);
      s.nevents = ParValue;
      */
      IN.getline(buffer, 200, '\n');
      sscanf(buffer, "Kfact\t%f", &ParValue);
      s.kfact = ParValue;

      IN.getline(buffer, 200, '\n');
      sscanf(buffer, "Lumi\t%f", &ParValue);
      s.lumi = ParValue;

      IN.getline(buffer, 200, '\n');
      sscanf(buffer, "Type\t%s", StringValue);
      s.type = StringValue;

      IN.getline(buffer, 200, '\n');
      sscanf(buffer, "Color\t%f", &ParValue);
      s.color = ParValue;

      if(fVerbose > 0){
	cout << " ---- " << endl;
	cout << "  New sample added: " << s.name << endl;
	cout << "   Sample no.      " << counter << endl;
	cout << "   Short name:     " << s.sname << endl;
	cout << "   File:           " << (s.file)->GetName() << endl;
 	//cout << "   Events:         " << s.nevents  << endl;
	cout << "   Events in tree: " << s.tree->GetEntries() << endl;
	cout << "   Xsection:       " << s.xsection << endl;
	cout << "   Lumi:           " << s.lumi << endl;
	cout << "   kfactor:        " << s.kfact << endl;
	cout << "   type:           " << s.type << endl;
	cout << "   Color:          " << s.color << endl;
      }
      fSamples.push_back(s);
      counter++;
    }
  }
  if(fVerbose > 0) cout << "------------------------------------" << endl;
}
