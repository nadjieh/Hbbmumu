#include "TEventList.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TString.h"
#include "TObject.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TH2.h"
#include "TH3.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TCut.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TH3F.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TH3D.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TAxis.h"
#include "TSystem.h"
#include "TMap.h"
#include <sstream>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <math.h>
#include <vector>
#include <stdlib.h>
#include <stdio.h>
#include <map>
#include <string>
#include <cmath>
#include "/dataLOCAL/paktinat/Projects/MT2/MT2_V02-03-02/MT2Analysis/Code/MT2AnalysisCode/MT2Code/include/MT2tree.hh"
#include "/dataLOCAL/paktinat/Projects/MT2/MT2_V02-03-02/MT2Analysis/Code/MT2AnalysisCode/TESCO/include/helper/Utilities.hh"

using namespace std;

void load(const char* filename = "/dataLOCAL/paktinat/Projects/MT2/MT2_V02-03-02/MT2Analysis/Code/MT2AnalysisCode/MT2Code/samples/samplesMine.dat");
void GetBEfficiency_InputHistograms();

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
	float PU_avg_weight;
	float lumi;
	int color;
};

vector < sample >  fSamples;
const int fVerbose = 3;
TString fPath;

//GetBEfficiency input histograms for pt binned only (1d), pt and eta binned (2d), pt binned in eta slices (1d)
void GetBEfficiency_InputHistograms(){
	gROOT->ProcessLine(".x SetStyle_PRD.C");

	gROOT->SetStyle("Plain");
	gStyle->SetOptStat(0);

	bool runQCD   = false;//set false if want to save some time //useful if look at high MT2
	bool onlyQCD  = false;//set false if want to save some time //useful if look at high MT2
	bool calcsusy = true;//at the moment has no genflavour // skip it
	bool runData  = false;//set false if want to save some time //at the moment not used, only for xchecks

	bool fHT      = true;//to change  htskim
	bool fMET     = false; //to change metskim
	bool fChange  = true; //use this only if you want to change the skims by hand

	//not implemented here, yet
	bool ttbarxsecup   = false;
	bool ttbarxsecdown = false;
	bool wjetsxsecup   = false;
	bool wjetsxsecdown = false;

//	TString samples         = "/shome/haweber/MT2Analysis_8TeV/Code/MT2AnalysisCode/RootMacros/samples/samples_Stop_53X_0bMET30_corr.dat";//Dileptons
//	TString samples         = "/shome/haweber/MT2Analysis_8TeV/Code/MT2AnalysisCode/RootMacros/samples/samples_type1MET_CHS_53X_MET.dat";//MET
//	TString samples         = "/shome/haweber/MT2Analysis_8TeV/Code/MT2AnalysisCode/RootMacros/samples/samples_type1MET_CHS_53X_HT.dat";//HT
//	TString samples         = "/shome/haweber/MT2Analysis_8TeV/Code/MT2AnalysisCode/RootMacros/samples/samples_type1MET_CHS_53X_MET_newTTbar.dat";//MET
	TString samples         = "/dataLOCAL/paktinat/Projects/MT2/MT2_V02-03-02/MT2Analysis_8TeV/Code/MT2AnalysisCode/RootMacros/samples/samplesMine.dat";//HT
//	TString samples         = "/shome/haweber/MT2Analysis_8TeV/Code/MT2AnalysisCode/RootMacros/samples/samples_type1MET_CHS_53X_METorHT.dat";//HT
	TString outputdir 	= "/dataLOCAL/paktinat/Projects/MT2/MT2_V02-03-02/MT2Analysis/Code/Efficiencies/InputHistos/";//if not empty end with slash
//	TString outputdir 	= "/shome/haweber/MT2Analysis_8TeV/Code/Efficiencies/test/InputHistos/";//if not empty end with slash
	Util::MakeOutputDir(outputdir);
//	TString outputname 	= "BEfficiencies_histos_DileptonicStop_CSVM.root";
//	TString outputname 	= "BEfficiencies_histos_ge2jets_0l_CSVM_HT750orMET200_MinDPhiPt40_noQCD.root";
//	TString outputname 	= "BEfficiencies_histos_ge2jets_0l_CSVM_HT750orMET200_MinDPhi4_noQCD.root";
//	TString outputname 	= "BEfficiencies_histos_ge2jets_CSVM_HT750orMET200_MinDPhiPt40_noQCD.root";
//	TString outputname 	= "BEfficiencies_histos_ge2jets_CSVM_HT750orMET200_MinDPhi4.root";
//	TString outputname 	= "BEfficiencies_histos_ge2jets_CSVM_HT450andMET200_MinDPhi4_newTTbar.root";
	TString outputname 	= "BEfficiencies_histos_ge2jets_CSVM_HT750andMET30_MinDPhi4_newTTbar_noQCD.root";

//	TString outputname 	= "BEfficiencies_histos_HT_ge3jets_0l_CSVM_allMT2_MinDPhi40.root";
//	TString outputname 	= "BEfficiencies_histos_HT_ge3jets_0l_CSVM_allMT2.root";
//	TString outputname 	= "BEfficiencies_histos_HT_ge3jets_ge1l_CSVM_allMT2.root";
//	TString outputname 	= "BEfficiencies_histos_HT_ge6jets_0l_CSVM_allMT2.root";
//	TString outputname 	= "BEfficiencies_histos_HT_ge6jets_ge1l_CSVM_allMT2.root";
//	TString outputname 	= "BEfficiencies_histos_HT_eq2jets_0l_CSVM_allMT2.root";
//	TString outputname 	= "BEfficiencies_histos_HT_eq2jets_ge1l_CSVM_allMT2.root";
//	TString outputname 	= "BEfficiencies_histos_HT_ge3jets_0l_CSVM_MT2lt125.root";
//	TString outputname 	= "BEfficiencies_histos_HT_ge3jets_ge1l_CSVM_MT2lt125.root";
//	TString outputname 	= "BEfficiencies_histos_HT_ge6jets_0l_CSVM_MT2lt125.root";
//	TString outputname 	= "BEfficiencies_histos_HT_ge6jets_ge1l_CSVM_MT2lt125.root";
//	TString outputname 	= "BEfficiencies_histos_HT_eq2jets_0l_CSVM_MT2lt125.root";
//	TString outputname 	= "BEfficiencies_histos_HT_eq2jets_ge1l_CSVM_MT2lt125.root";
//	TString outputname 	= "BEfficiencies_histos_HT_ge3jets_0l_CSVM_MT2ge125.root";
//	TString outputname 	= "BEfficiencies_histos_HT_ge3jets_ge1l_CSVM_MT2ge125.root";
//	TString outputname 	= "BEfficiencies_histos_HT_ge6jets_0l_CSVM_MT2ge125.root";
//	TString outputname 	= "BEfficiencies_histos_HT_ge6jets_ge1l_CSVM_MT2ge125.root";
//	TString outputname 	= "BEfficiencies_histos_HT_eq2jets_0l_CSVM_MT2ge125.root";
//	TString outputname 	= "BEfficiencies_histos_HT_eq2jets_ge1l_CSVM_MT2ge125.root";
//	TString outputname 	= "BEfficiencies_histos_MET_ge3jets_0l_CSVM_allMT2.root";
//	TString outputname 	= "BEfficiencies_histos_MET_ge3jets_ge1l_CSVM_allMT2.root";
//	TString outputname 	= "BEfficiencies_histos_MET_ge6jets_0l_CSVM_allMT2.root";
//	TString outputname 	= "BEfficiencies_histos_MET_ge6jets_ge1l_CSVM_allMT2.root";
//	TString outputname 	= "BEfficiencies_histos_MET_eq2jets_0l_CSVM_allMT2.root";
//	TString outputname 	= "BEfficiencies_histos_MET_eq2jets_ge1l_CSVM_allMT2.root";
//	TString outputname 	= "BEfficiencies_histos_MET_ge3jets_0l_CSVM_MT2lt125.root";
//	TString outputname 	= "BEfficiencies_histos_MET_ge3jets_ge1l_CSVM_MT2lt125.root";
//	TString outputname 	= "BEfficiencies_histos_MET_ge6jets_0l_CSVM_MT2lt125.root";
//	TString outputname 	= "BEfficiencies_histos_MET_ge6jets_ge1l_CSVM_MT2lt125.root";
//	TString outputname 	= "BEfficiencies_histos_MET_eq2jets_0l_CSVM_MT2lt125.root";
//	TString outputname 	= "BEfficiencies_histos_MET_eq2jets_ge1l_CSVM_MT2lt125.root";
//	TString outputname 	= "BEfficiencies_histos_MET_ge3jets_0l_CSVM_MT2ge125.root";
//	TString outputname 	= "BEfficiencies_histos_MET_ge3jets_ge1l_CSVM_MT2ge125.root";
//	TString outputname 	= "BEfficiencies_histos_MET_ge6jets_0l_CSVM_MT2ge125.root";
//	TString outputname 	= "BEfficiencies_histos_MET_ge6jets_ge1l_CSVM_MT2ge125.root";
//	TString outputname 	= "BEfficiencies_histos_MET_eq2jets_0l_CSVM_MT2ge125.root";
//	TString outputname 	= "BEfficiencies_histos_MET_eq2jets_ge1l_CSVM_MT2ge125.root";

	int btaggername = 4; //0: TCHE, 1: TCHP, 2: SSVHE, 3: SSVHP, 4: CSV, 5: JP
	//CSV: L: 0.244, M: 0.679, T: 0.898
	//JP:  L: 0.275, M: 0.545, T: 0.790
	//TCHPT: 3.41, SSVHEM: 1.74, SSVHPT: 2.00
	double bdiscr_value = 0.679;
	int WP = 1;// WP: 0: loose, 1: medium, 2: tight
	double bpt = 20.;
	double beta = 2.4;
	const int Nptbins   = 17; double ptbins[Nptbins+1]     = {20, 30, 40, 50, 60, 70, 80, 100, 120, 160, 210, 260, 320, 400, 500, 600, 800, 1000};
//	const int Nptbins   = 16; double ptbins[Nptbins+1]     = {20, 30, 40, 50, 60, 70, 80, 100, 120, 160, 210, 260, 320, 400, 500, 670, 750};
	const int NetabinsL =  4; double etabinsL[NetabinsL+1] = {0.0, 0.5, 1.0, 1.5, 2.4};
	const int NetabinsM =  3; double etabinsM[NetabinsM+1] = {0.0, 0.8, 1.6, 2.4};
	const int NetabinsT =  1; double etabinsT[NetabinsT+1] = {0.0, 2.4};

	map<string, TH2D*> histos2d;
	map<string, TH1D*> histospt;

	const int sampletypesize = 10;
	string sample_type[sampletypesize] = {"QCD", "WJets", "ZJets", "TTbar", "SingleTop", "Other", "mc", "mcNOqcd", "susy", "data"};//same as type in samples.dat

	const int prefixsize = 10;
	string prefix[10] = {"B", "C", "L", "NonB", "All", "BTaggedB", "BTaggedC", "BTaggedL", "BTaggedNonB", "BTaggedAll"};


	for(int is = 0; is<sampletypesize; ++is){
	   for(int js = 0; js<10; ++js){
		string hs = string("_") + sample_type[is];
		string mapname;
		mapname = prefix[js] + "Jet" + hs;
		if(WP==0 && histos2d.count(mapname) == 0 ) histos2d[mapname] = new TH2D((mapname+string("_2d")).c_str(), "eff", Nptbins, ptbins, NetabinsL, etabinsL);
		if(WP==1 && histos2d.count(mapname) == 0 ) histos2d[mapname] = new TH2D((mapname+string("_2d")).c_str(), "eff", Nptbins, ptbins, NetabinsM, etabinsM);
		if(WP==2 && histos2d.count(mapname) == 0 ) histos2d[mapname] = new TH2D((mapname+string("_2d")).c_str(), "eff", Nptbins, ptbins, NetabinsT, etabinsT);
		if(         histospt.count(mapname) == 0 ) histospt[mapname] = new TH1D( mapname.c_str(),                "eff", Nptbins, ptbins);
	   }
	}
	
	for(map<string,TH2D*>::iterator h=histos2d.begin(); h!=histos2d.end();++h){
		//cout << "Load " << h->first << endl;
		h->second->Sumw2();}
	for(map<string,TH1D*>::iterator h=histospt.begin(); h!=histospt.end();++h){
		//cout << "Load " << h->first << endl;
		h->second->Sumw2();}

	cout << "Histograms loaded" << endl;
	
	std::ostringstream cutStream;
	cutStream       << " " 	  
     << "(misc.ProcessID!=6||(misc.Event!=1689009&&misc.Event!=2275452&&misc.Event!=1946785&&misc.Event!=1936763&&misc.Event!=1890738&&misc.Event!=1757319))" << "&&" //summer QCD 53X
     << "(misc.ProcessID!=10||(Susy.MassChi==350&&Susy.MassLSP==50))" << "&&" // T2bb M_stop = 350, M_lsp = 50
 //    << "(misc.ProcessID!=10||(Susy.MassChi==600&&Susy.MassLSP==50))" << "&&" // T2bb M_stop = 600, M_lsp = 50

//	  << "misc.MT2 >=0"                                      << "&&" 
//	  << "misc.MT2 < 125"                                    << "&&" 
//	  << "misc.MT2 >=125"                                    << "&&" 
	  << "misc.MET>=30"                                      << "&&"
	  << "misc.HT>=450"                                      << "&&"
	  << "misc.Jet0Pass ==1"                                 << "&&"
	  << "misc.Jet1Pass ==1"                                 << "&&"
//	  << "misc.SecondJPt  >100"                              << "&&"
	  << "misc.PassJetID ==1"                                << "&&"
	  << "misc.Vectorsumpt < 70"                             << "&&"
//	  << "(NJetsIDLoose40 +NEles +NMuons)>=3"                << "&&"
	  << "misc.MinMetJetDPhi4>0.3"                           << "&&"
//	  << "misc.MinMetJetDPhiPt40>0.3"                        << "&&"
//        << "NBJetsCSVM     >=1"                                << "&&"
//        JetSkim
	  << "NJetsIDLoose40 >=2"                                << "&&"
//	  << "NJetsIDLoose40 >=3"                                << "&&"
//	  << "NJetsIDLoose40 >=6"                                << "&&"
//	  << "NJetsIDLoose40 ==2"                                << "&&"
//        HTskim
//	  << "misc.HT >=750"                                     << "&&" 
//        METskim
//	  << "misc.MET >=200"                                    << "&&"
//	  LeptonSkim
//	  << "(NEles+NMuons) ==0"                                << "&&"
//	  << "(NEles+NMuons) >=1"                                << "&&"
// Stop Skim -------------------------
//	  << "NJetsIDLoose40 >=2"                                << "&&"
//        << "(NEles+NMuons) >=2"                                << "&&"
//        << "NBJetsCSVM     >=1"                                << "&&"
//	  << "NJetsIDLoose40 >=3"                                << "&&" 
// Noise -- first 6 are official filters
	  << "( misc.ProcessID==10 || misc.HBHENoiseFlag == 0)"  << "&&"
	  << "misc.CSCTightHaloIDFlag == 0"                      << "&&"
	  << "(misc.ProcessID==10||misc.hcalLaserEventFlag==0)"  << "&&"
//	  << "(misc.isFastSim || misc.HBHENoiseFlag == 0)"       << "&&"
//	  << "(misc.isFastSim || misc.CSCTightHaloIDFlag == 0)"  << "&&"
//	  << "(misc.isFastSim || misc.hcalLaserEventFlag == 0)"  << "&&"
	  << "misc.trackingFailureFlag == 0"                     << "&&"
	  << "misc.eeBadScFlag == 0"                             << "&&"
	  << "misc.EcalDeadCellTriggerPrimitiveFlag == 0";

	TString cuts = cutStream.str().c_str();
/*
	std::ostringstream triggerStreamEE;
	triggerStreamEE << "( "
	  << "(trigger.HLT_DiElectrons==1)" << " )";
	TString triggerEE = triggerStreamEE.str().c_str();
	std::ostringstream triggerStreamEMu;
	triggerStreamEMu << "( "
	  << "(trigger.HLT_EMu==1)" << " )";
	TString triggerEMu = triggerStreamEMu.str().c_str();
	std::ostringstream triggerStreamMuMu;
	triggerStreamMuMu << "( "
	  << "(trigger.HLT_DiMuons==1)" << " )";
	TString triggerMuMu = triggerStreamMuMu.str().c_str();
*//*
	std::ostringstream triggerStream;
	triggerStream << "( "
			<< "trigger.HLT_PFHT650_v5 == 1 || trigger.HLT_PFHT650_v6 == 1 || trigger.HLT_PFHT650_v7 == 1 || "
			<< "trigger.HLT_PFHT650_v8 == 1 || trigger.HLT_PFHT650_v9 == 1 || "
			<< "trigger.HLT_PFNoPUHT650_v1 == 1 || trigger.HLT_PFNoPUHT650_v3 == 1)";
	
	TString trigger = triggerStream.str().c_str();
*/
	std::ostringstream triggerStream;
	triggerStream << "( "
			<< "trigger.HLT_PFMET150_v2 == 1 || trigger.HLT_PFMET150_v3 == 1 || trigger.HLT_PFMET150_v4 == 1 || "
			<< "trigger.HLT_PFMET150_v5 == 1 || trigger.HLT_PFMET150_v6 == 1 || trigger.HLT_PFMET150_v7 == 1 )";
	TString trigger = triggerStream.str().c_str();
  

	load(samples.Data());

	for(size_t i = 0; i < fSamples.size(); ++i){

   	    if(runData==false && fSamples[i].type=="data") continue;
	    if(calcsusy==false && fSamples[i].type=="susy") continue;
	    if(runQCD==false   && fSamples[i].sname=="QCD") continue;
	    if(onlyQCD==true   && fSamples[i].sname!="QCD" && runQCD) continue;

	    string sampletype = (string)fSamples[i].type;
	    bool METskim = false;
	    bool HTskim  = false;
	    TString path = fSamples[i].file->GetName();
	    if(path.Contains("METskim")) METskim = true;
	    if(path.Contains("HTskim"))  HTskim  = true;
	    if(fChange){
		if(fHT)  { HTskim = true;  METskim = false; }
		if(fMET) { HTskim = false; METskim = true;  }
	    	if(!(fMET || fHT)) { HTskim = false; METskim = false; }
	    }
//		cout << true << " " << path << " " << METskim << " " << HTskim << endl;
	    if(!(METskim || HTskim)) continue;
	    if(sampletype==(string)"mc"){
		if(fSamples[i].sname=="QCD") sampletype = (string)fSamples[i].sname;
		else if(fSamples[i].sname=="Wtolnu") sampletype = (string)"WJets";
		else if(fSamples[i].sname=="DY")     sampletype = (string)"ZJets";
		else if(fSamples[i].name=="TTbar")   sampletype = (string)"TTbar";
		else if(fSamples[i].name=="TTbar_Had")   sampletype = (string)"TTbar";
		else if(fSamples[i].name=="TTbar_SemiLep")   sampletype = (string)"TTbar";
		else if(fSamples[i].name=="TTbar_FullyLep")   sampletype = (string)"TTbar";
		else if(fSamples[i].sname=="TTbarV") sampletype = (string)"TTbar";
		else if(fSamples[i].sname=="Top")    sampletype = (string)"SingleTop";//no ttbar
		//else if(fSamples[i].sname=="VV" || fSamples[i].sname=="VVV") sampletype = (string)"VV/VVV";
		else sampletype = (string)"Other";
	    }
	//    if(sampletype==(string)"susy") sampletype=(string)"Stop";

	    Double_t sample_weight = fSamples[i].xsection * fSamples[i].kfact * fSamples[i].lumi / (fSamples[i].nevents);
	    if(fVerbose>2) cout << "B-Efficiency: looping over " << fSamples[i].name << endl;
	    if(fVerbose>2) cout << "              sample has weight " << sample_weight << " and " << fSamples[i].tree->GetEntries() << " entries" << endl; 
	    if(fVerbose>2 && fSamples[i].tree->GetEntries()==0) cout << "skip sample, has no entries" << endl;
            if(fSamples[i].tree->GetEntries()==0) continue;
	    MT2tree* fMT2tree = new MT2tree();
	    fSamples[i].tree->SetBranchAddress("MT2tree", &fMT2tree);
	    Long64_t nentries =  fSamples[i].tree->GetEntries();
	    Long64_t nbytes = 0, nb = 0;
	    int nev =0;

	    TString myCuts = cuts;
	    if(HTskim)  myCuts = myCuts + " && misc.HT>750";
	    if(METskim) myCuts = myCuts + " && misc.HT<=750 && misc.MET>200";//METskim excluding HTskim

	/*
	    if( fSamples[i].type=="data" && fSamples[i].sname=="EE-Data") { myCuts += " && " + triggerEE; }//cuts to be aplied only on data
	    else if( fSamples[i].type=="data" && fSamples[i].sname=="EMu-Data") { myCuts += " && " + triggerEMu; }//cuts to be aplied only on data
	    else if( fSamples[i].type=="data" && fSamples[i].sname=="MuMu-Data") { myCuts += " && " + triggerMuMu; } //cuts to be aplied only on data
	    else if(fSamples[i].type=="data") {cout << "data not usuable" << " type " << fSamples[i].type << " Sname " << fSamples[i].sname << endl; continue; }//not 
	*/
	    if( fSamples[i].type=="data") { myCuts += " && " + trigger; }

   	    cout << "Cuts for Flow: " << myCuts << endl;
   	    fSamples[i].tree->Draw(">>selList", myCuts);

	    TEventList *myEvtList = (TEventList*)gDirectory->Get("selList");

	    fSamples[i].tree->SetEventList(myEvtList);

	    int counter=0;
	    cout << "Filtering done, size=" <<myEvtList->GetN()  << endl;
	    if(myEvtList->GetSize()==0) continue;
        while(myEvtList->GetEntry(counter++) !=-1){	
      		int jentry = myEvtList->GetEntry(counter-1);
            
            nb =  fSamples[i].tree->GetEntry(jentry);   nbytes += nb;
            fSamples[i].tree->SetBranchAddress("MT2tree", &fMT2tree);
            
            if ( fVerbose>2 && counter % 50000 == 0 )  cout << "+++ Proccessing event " << counter << endl;


	    //if((fSamples[i].sname!="VVV"&&fSamples[i].sname!="TTbarV" && fSamples[i].type!="susy"&&sampletype!="Stop")  && fMT2tree->misc.HBHENoiseFlag!=0     ) continue;//test
	    //if((fSamples[i].sname!="VVV"&&fSamples[i].sname!="TTbarV" && fSamples[i].type!="susy"&&sampletype!="Stop")  && fMT2tree->misc.hcalLaserEventFlag!=0) continue;//test
            Double_t weight = sample_weight;
            if (!fMT2tree->misc.isData ) weight = weight * fMT2tree->pileUp.Weight; // pile-up reweighting for MC 
 /*           Bool_t recoedee   = false;// exact 2 ele, 0 muo
            Bool_t recoedemu  = false;// exact 1 ele, 1 muo
            Bool_t recoedmumu = false;// exact 0 ele, 2 muo
		//change this - might have some inefficiencies - three leptons with one not passing IDMedium==1
		//20-20 selection --> third lepton is 10 GeV 
            if(fMT2tree->NEles>=2 &&fMT2tree->ele[0].lv.Pt()>20&&fMT2tree->ele[1].lv.Pt()>20&&fMT2tree->ele[2].lv.Pt()<20&&(fMT2tree->NMuons==0||fMT2tree->muo[0].lv.Pt()<10) && fMT2tree->ele[0].IDMedium==1 && fMT2tree->ele[1].IDMedium==1) recoedee   = true;
            if(fMT2tree->NMuons>=2&&fMT2tree->muo[0].lv.Pt()>20&&fMT2tree->muo[1].lv.Pt()>20&&fMT2tree->muo[2].lv.Pt()<10&&(fMT2tree->NEles ==0||fMT2tree->ele[0].lv.Pt()<10)) recoedmumu = true;
            if(fMT2tree->NMuons>=1&&fMT2tree->muo[0].lv.Pt()>20&&fMT2tree->muo[1].lv.Pt()<10&&fMT2tree->NEles>=1 &&fMT2tree->ele[0].lv.Pt()>20&&fMT2tree->ele[1].lv.Pt()<10 && fMT2tree->ele[0].IDMedium==1) recoedemu  = true;//XXX change emu to 20/20 selection, left ee,mumu at 20/10 for now
	    //id cuts?

            Bool_t recoedosee   = false;// opposite sign
            Bool_t recoedosemu  = false;// opposite sign
            Bool_t recoedosmumu = false;// opposite sign
	    if(recoedee   && (fMT2tree->ele[0].Charge)*(fMT2tree->ele[1].Charge)==(-1)) recoedosee   = true;
	    if(recoedemu  && (fMT2tree->ele[0].Charge)*(fMT2tree->muo[0].Charge)==(-1)) recoedosemu  = true;
	    if(recoedmumu && (fMT2tree->muo[0].Charge)*(fMT2tree->muo[1].Charge)==(-1)) recoedosmumu = true;
            Bool_t recoedeenZ   = false;// off Z
            Bool_t recoedemunZ  = false;// off Z -> is not requirement, so this variable should be useless
            Bool_t recoedmumunZ = false;// off Z
	    if(recoedee   && ((fMT2tree->ele[0].lv + fMT2tree->ele[1].lv).M()<81 || (fMT2tree->ele[0].lv + fMT2tree->ele[1].lv).M()>101) ) recoedeenZ   = true;
	    if(recoedemu  && ((fMT2tree->ele[0].lv + fMT2tree->muo[0].lv).M()<81 || (fMT2tree->ele[0].lv + fMT2tree->muo[0].lv).M()>101) ) recoedemunZ  = true;
	    if(recoedmumu && ((fMT2tree->muo[0].lv + fMT2tree->muo[1].lv).M()<81 || (fMT2tree->muo[0].lv + fMT2tree->muo[1].lv).M()>101) ) recoedmumunZ = true;

	    if(!(recoedee)   && !(recoedemu)   && !(recoedmumu)  ) continue;//require dilepton
	    if(!(recoedosee) && !(recoedosemu) && !(recoedosmumu) ) continue;//require os-dilepton, maybe comment this for background est.

	    if(recoedee   && !(fMT2tree->misc.isData) ) weight = 0.960 * weight * 0.953366;//from SSb Florida
	    if(recoedemu  && !(fMT2tree->misc.isData) ) weight = 0.934 * weight;//from SSb Florida
	    if(recoedmumu && !(fMT2tree->misc.isData) ) weight = 0.875 * weight * 1.031213;//from SSb Florida
*/
	   	if(wjetsxsecup  ==true && fSamples[i].sname=="Wtolnu")	weight = weight*0.7;
	   	if(wjetsxsecdown==true && fSamples[i].sname=="Wtolnu")	weight = weight*1.3;
	   	if(ttbarxsecup  ==true && fSamples[i].name== "TTbar" )	weight = weight*0.85;
	   	if(ttbarxsecdown==true && fSamples[i].name== "TTbar" )	weight = weight*1.15;

		for(int k = 0; k<fMT2tree->NJets; ++k){
			if(fMT2tree->jet[k].isPFIDLoose==false) continue;
			if(fMT2tree->jet[k].Flavour<=-7777 && fMT2tree->misc.isData==false) continue;//MC sample veto if no btagging information
			if(abs(fMT2tree->jet[k].Flavour>100)) cout << "jet " << k << "has flavour " << fMT2tree->jet[k].Flavour << endl;
			float btagdiscr=-1;
			if(btaggername==0)      btagdiscr = fMT2tree->jet[k].bTagProbTCHE;
			else if(btaggername==1) btagdiscr = fMT2tree->jet[k].bTagProbTCHP;
			else if(btaggername==2) btagdiscr = fMT2tree->jet[k].bTagProbSSVHE;
			else if(btaggername==3) btagdiscr = fMT2tree->jet[k].bTagProbSSVHP;
			else if(btaggername==4) btagdiscr = fMT2tree->jet[k].bTagProbCSV;
			else if(btaggername==5) btagdiscr = fMT2tree->jet[k].bTagProbJProb;
			else                    btagdiscr = fMT2tree->jet[k].bTagProbCSV;
			float jetpt = fMT2tree->jet[k].lv.Pt();
			float jeteta= fabs(fMT2tree->jet[k].lv.Eta() );
			//float jetMT2= fMT2tree->misc.MT2;
			int jetabsflavour = abs(fMT2tree->jet[k].Flavour);
			if(jetpt<bpt)   continue;
			if(jeteta>beta) continue;
			histospt[string("AllJet_") + sampletype]->Fill(jetpt, weight );
			histos2d[string("AllJet_") + sampletype]->Fill(jetpt, jeteta, weight );
			if(jetabsflavour==5){
				histospt[string("BJet_") + sampletype]->Fill(jetpt, weight );
				histos2d[string("BJet_") + sampletype]->Fill(jetpt, jeteta, weight );
			}
			else if(jetabsflavour==4){
				histospt[string("CJet_") + sampletype]->Fill(jetpt, weight );
				histos2d[string("CJet_") + sampletype]->Fill(jetpt, jeteta, weight );
			}
			else{
				histospt[string("LJet_") + sampletype]->Fill(jetpt, weight );
				histos2d[string("LJet_") + sampletype]->Fill(jetpt, jeteta, weight );
			}
			if(btagdiscr>=bdiscr_value){
				histospt[string("BTaggedAllJet_") + sampletype]->Fill(jetpt, weight );
				histos2d[string("BTaggedAllJet_") + sampletype]->Fill(jetpt, jeteta, weight );
				if(jetabsflavour==5){
					histospt[string("BTaggedBJet_") + sampletype]->Fill(jetpt, weight );
					histos2d[string("BTaggedBJet_") + sampletype]->Fill(jetpt, jeteta, weight );
				}
				else if(jetabsflavour==4){
					histospt[string("BTaggedCJet_") + sampletype]->Fill(jetpt, weight );
					histos2d[string("BTaggedCJet_") + sampletype]->Fill(jetpt, jeteta, weight );
				}
				else{
					histospt[string("BTaggedLJet_") + sampletype]->Fill(jetpt, weight );
					histos2d[string("BTaggedLJet_") + sampletype]->Fill(jetpt, jeteta, weight );
				}
			}
		}//for(int k = 0; k<fMT2tree->NJets; ++k)
	   }
	   delete fMT2tree;
	}//for(size_t i = 0; i < fSamples.size(); ++i)

	cout << "add overflow bins" << endl;
	for(map<string,TH2D*>::iterator h=histos2d.begin(); h!=histos2d.end();++h){
		//there should not be any underflow, just look at overflow
		int nbinsx = h->second->GetNbinsX();
		int nbinsy = h->second->GetNbinsY();
		for(int ix = 1; ix<= nbinsx;++ix){
			h->second->SetBinContent(ix, nbinsy, 
						h->second->GetBinContent(ix,nbinsy) +
						h->second->GetBinContent(ix, nbinsy+1));
			h->second->SetBinError(ix, nbinsy, 
						sqrt(h->second->GetBinError(ix,nbinsy)*h->second->GetBinError(ix,nbinsy) +
						h->second->GetBinError(ix, nbinsy+1)*h->second->GetBinError(ix, nbinsy+1)) );
		}
		for(int iy = 1; iy<= nbinsy;++iy){
			h->second->SetBinContent(nbinsx, iy, 
						h->second->GetBinContent(nbinsx,iy) +
						h->second->GetBinContent(nbinsx+1, iy));
			h->second->SetBinError(nbinsx, iy, 
						sqrt(h->second->GetBinError(nbinsx,iy)*h->second->GetBinError(nbinsx,iy) +
						h->second->GetBinError(nbinsx+1, iy)*h->second->GetBinError(nbinsx+1, iy)) );
		}
		h->second->SetBinContent(nbinsx, nbinsy, 
						h->second->GetBinContent(nbinsx,nbinsy) +
						h->second->GetBinContent(nbinsx+1, nbinsy+1));
		h->second->SetBinError(nbinsx, nbinsy, 
					sqrt(h->second->GetBinError(nbinsx,nbinsy)*h->second->GetBinError(nbinsx,nbinsy) +
					h->second->GetBinError(nbinsx+1, nbinsy+1)*h->second->GetBinError(nbinsx+1, nbinsy+1)) );
	}
	for(map<string,TH1D*>::iterator h=histospt.begin(); h!=histospt.end();++h){
		//there should not be any underflow, just look at overflow
		int nbinsx = h->second->GetNbinsX();
		h->second->SetBinContent(nbinsx, h->second->GetBinContent(nbinsx) + h->second->GetBinContent(nbinsx+1));
		h->second->SetBinError(  nbinsx, h->second->GetBinError(nbinsx  )*h->second->GetBinError(nbinsx  ) +
                                                 h->second->GetBinError(nbinsx+1)*h->second->GetBinError(nbinsx+1) );
	}


	cout << "merge mc histograms into one" << endl;
	//fill mc histogram from the samples histograms
	for(int ists = 0; ists<sampletypesize; ++ists){
		if((sample_type[ists])==("data")) continue;
		if((sample_type[ists])==("susy")) continue;
		if((sample_type[ists])==("mc") ) continue;//have all mc subsamples only
		string helptype = string("_") + sample_type[ists];
		string helpmc = "_mc";
		string helpmcnoqcd = "_mcNOqcd";
		histospt[string("BJet")   + helpmc]->Add(histospt[string("BJet")   + helptype], 1.);
		histospt[string("CJet")   + helpmc]->Add(histospt[string("CJet")   + helptype], 1.);
		histospt[string("LJet")   + helpmc]->Add(histospt[string("LJet")   + helptype], 1.);
		histospt[string("AllJet") + helpmc]->Add(histospt[string("AllJet") + helptype], 1.);
		histos2d[string("BJet")   + helpmc]->Add(histos2d[string("BJet")   + helptype], 1.);
		histos2d[string("CJet")   + helpmc]->Add(histos2d[string("CJet")   + helptype], 1.);
		histos2d[string("LJet")   + helpmc]->Add(histos2d[string("LJet")   + helptype], 1.);
		histos2d[string("AllJet") + helpmc]->Add(histos2d[string("AllJet") + helptype], 1.);
		histospt[string("BTaggedBJet")   + helpmc]->Add(histospt[string("BTaggedBJet")   + helptype], 1.);
		histospt[string("BTaggedCJet")   + helpmc]->Add(histospt[string("BTaggedCJet")   + helptype], 1.);
		histospt[string("BTaggedLJet")   + helpmc]->Add(histospt[string("BTaggedLJet")   + helptype], 1.);
		histospt[string("BTaggedAllJet") + helpmc]->Add(histospt[string("BTaggedAllJet") + helptype], 1.);
		histos2d[string("BTaggedBJet")   + helpmc]->Add(histos2d[string("BTaggedBJet")   + helptype], 1.);
		histos2d[string("BTaggedCJet")   + helpmc]->Add(histos2d[string("BTaggedCJet")   + helptype], 1.);
		histos2d[string("BTaggedLJet")   + helpmc]->Add(histos2d[string("BTaggedLJet")   + helptype], 1.);
		histos2d[string("BTaggedAllJet") + helpmc]->Add(histos2d[string("BTaggedAllJet") + helptype], 1.);
		if(sample_type[ists]==string("QCD")) continue;
		histospt[string("BJet")   + helpmcnoqcd]->Add(histospt[string("BJet")   + helptype], 1.);
		histospt[string("CJet")   + helpmcnoqcd]->Add(histospt[string("CJet")   + helptype], 1.);
		histospt[string("LJet")   + helpmcnoqcd]->Add(histospt[string("LJet")   + helptype], 1.);
		histospt[string("AllJet") + helpmcnoqcd]->Add(histospt[string("AllJet") + helptype], 1.);
		histos2d[string("BJet")   + helpmcnoqcd]->Add(histos2d[string("BJet")   + helptype], 1.);
		histos2d[string("CJet")   + helpmcnoqcd]->Add(histos2d[string("CJet")   + helptype], 1.);
		histos2d[string("LJet")   + helpmcnoqcd]->Add(histos2d[string("LJet")   + helptype], 1.);
		histos2d[string("AllJet") + helpmcnoqcd]->Add(histos2d[string("AllJet") + helptype], 1.);
		histospt[string("BTaggedBJet")   + helpmcnoqcd]->Add(histospt[string("BTaggedBJet")   + helptype], 1.);
		histospt[string("BTaggedCJet")   + helpmcnoqcd]->Add(histospt[string("BTaggedCJet")   + helptype], 1.);
		histospt[string("BTaggedLJet")   + helpmcnoqcd]->Add(histospt[string("BTaggedLJet")   + helptype], 1.);
		histospt[string("BTaggedAllJet") + helpmcnoqcd]->Add(histospt[string("BTaggedAllJet") + helptype], 1.);
		histos2d[string("BTaggedBJet")   + helpmcnoqcd]->Add(histos2d[string("BTaggedBJet")   + helptype], 1.);
		histos2d[string("BTaggedCJet")   + helpmcnoqcd]->Add(histos2d[string("BTaggedCJet")   + helptype], 1.);
		histos2d[string("BTaggedLJet")   + helpmcnoqcd]->Add(histos2d[string("BTaggedLJet")   + helptype], 1.);
		histos2d[string("BTaggedAllJet") + helpmcnoqcd]->Add(histos2d[string("BTaggedAllJet") + helptype], 1.);

	}
	cout << "merge C and L histograms to NonB" << endl;
	for(int ists = 0; ists<sampletypesize; ++ists){
		string helptype = string("_") + sample_type[ists];
		histospt[string("NonBJet")   + helptype]->Add(histospt[string("CJet")   + helptype], 1.);
		histospt[string("NonBJet")   + helptype]->Add(histospt[string("LJet")   + helptype], 1.);
		histos2d[string("NonBJet")   + helptype]->Add(histos2d[string("CJet")   + helptype], 1.);
		histos2d[string("NonBJet")   + helptype]->Add(histos2d[string("LJet")   + helptype], 1.);
		histospt[string("BTaggedNonBJet")   + helptype]->Add(histospt[string("BTaggedCJet")   + helptype], 1.);
		histospt[string("BTaggedNonBJet")   + helptype]->Add(histospt[string("BTaggedLJet")   + helptype], 1.);
		histos2d[string("BTaggedNonBJet")   + helptype]->Add(histos2d[string("BTaggedCJet")   + helptype], 1.);
		histos2d[string("BTaggedNonBJet")   + helptype]->Add(histos2d[string("BTaggedLJet")   + helptype], 1.);
	}
	cout << "save histogram" << endl;
	TFile *fsavefile = new TFile(outputdir+outputname,"RECREATE");
	fsavefile->cd();
	for(map<string,TH2D*>::iterator h=histos2d.begin(); h!=histos2d.end();++h){
		h->second->Write();}
	for(map<string,TH1D*>::iterator h=histospt.begin(); h!=histospt.end();++h){
		h->second->Write();}
	fsavefile->Close();

	cout << "Saved histograms in " << outputdir << outputname << endl;

}//void GetBEfficiency_InputHistograms()




void load(const char* filename){
	
	fSamples.clear();
	char buffer[200];
	ifstream IN(filename);

	char ParName[100], StringValue[1000];
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

			IN.getline(buffer, 200, '\n');			//new
			sscanf(buffer, "ShapeName\t%s", StringValue);	//new
			s.shapename = TString(StringValue);		//new

			IN.getline(buffer, 200, '\n');
			sscanf(buffer, "File\t%s", StringValue);
			TString file =fPath+StringValue;
			TFile *f = TFile::Open(file);
			s.file = f;
			s.tree = (TTree*)f->Get("MassTree");
			
			IN.getline(buffer, 200, '\n');
			sscanf(buffer, "Xsection\t%f", &ParValue);
			s.xsection = ParValue;
			
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

			TH1F *h_PUWeights = (TH1F*) s.file->Get("h_PUWeights");
			TH1F *h_Events    = (TH1F*) s.file->Get("h_Events");
			if(h_PUWeights==0 || h_Events==0){
				cout << "ERROR: sample " << (s.file)->GetName() << " does not have PU and NEvents histos! " << endl;
				exit(1);
			}
			s.type!="data" ? s.PU_avg_weight = h_PUWeights->GetMean()    : s.PU_avg_weight =1;
			s.type!="data" ? s.nevents       = h_Events   ->GetEntries() : s.nevents       =1;
			delete h_PUWeights;
			delete h_Events;
			if(fVerbose > 0){
				cout << " ---- " << endl;
				cout << "  New sample added: " << s.name << endl;
				cout << "   Sample no.      " << counter << endl;
				cout << "   Short name:     " << s.sname << endl;
				cout << "   File:           " << (s.file)->GetName() << endl;
				cout << "   Events:         " << s.nevents  << endl;
				cout << "   Events in tree: " << s.tree->GetEntries() << endl; 
				cout << "   Xsection:       " << s.xsection << endl;
				cout << "   Lumi:           " << s.lumi << endl;
				cout << "   kfactor:        " << s.kfact << endl;
				cout << "   avg PU weight:  " << s.PU_avg_weight << endl;
				cout << "   type:           " << s.type << endl;
				cout << "   Color:          " << s.color << endl;
			}
			fSamples.push_back(s);
			counter++;
		}
	}
	if(fVerbose > 0) cout << "------------------------------------" << endl;
}
