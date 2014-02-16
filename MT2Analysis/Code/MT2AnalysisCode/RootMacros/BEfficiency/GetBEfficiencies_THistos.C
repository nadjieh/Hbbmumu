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
#include "TEfficiency.h"
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
#include "/shome/haweber/MT2Analysis_8TeV/Code/MT2AnalysisCode/MT2Code/include/MT2tree.hh"


void GetBEfficiencies_THistos();

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

void GetBEfficiencies_THistos(){
	gROOT->ProcessLine(".x SetStyle_PRD.C");

	gROOT->SetStyle("Plain");
	gStyle->SetOptStat(0);

	TString inputdir 	= "/shome/haweber/MT2Analysis_8TeV/Code/Efficiencies/InputHistos/";//if not empty end with slash
	TString outputdir 	= "/shome/haweber/MT2Analysis_8TeV/Code/Efficiencies/";//if not empty end with slash
	TString inputname 	= "BEfficiencies_histos_ge2jets_CSVM_HT750andMET30_or_HT450andMET200_MinDPhi4_newTTbar.root";
	TString outputname 	= "BEfficiencies_ge2jets_CSVM_HT750andMET30_or_HT450andMET200_MinDPhi4_newTTbar.root";


//	TString outputdir 	= "/shome/haweber/MT2Analysis_8TeV/Code/Efficiencies/test/";//if not empty end with slash
//	TString inputdir 	= "/shome/haweber/MT2Analysis_8TeV/Code/Efficiencies/test/InputHistos/";//if not empty end with slash
//	TString outputname 	= "BEfficiencies_HT_ge3jets_0l_CSVM_allMT2_MinDPhi40.root";
//	TString inputname 	= "BEfficiencies_histos_HT_ge3jets_0l_CSVM_allMT2_MinDPhi40.root";
//	TString outputname 	= "BEfficiencies_HT_ge3jets_0l_CSVM_allMT2.root";
//	TString inputname 	= "BEfficiencies_histos_HT_ge3jets_0l_CSVM_allMT2.root";
//	TString outputname 	= "BEfficiencies_HT_ge3jets_ge1l_CSVM_allMT2.root";
//	TString inputname 	= "BEfficiencies_histos_HT_ge3jets_ge1l_CSVM_allMT2.root";
//	TString outputname 	= "BEfficiencies_HT_ge6jets_0l_CSVM_allMT2.root";
//	TString inputname 	= "BEfficiencies_histos_HT_ge6jets_0l_CSVM_allMT2.root";
//	TString outputname 	= "BEfficiencies_HT_ge6jets_ge1l_CSVM_allMT2.root";
//	TString inputname 	= "BEfficiencies_histos_HT_ge6jets_ge1l_CSVM_allMT2.root";
//	TString outputname 	= "BEfficiencies_HT_eq2jets_0l_CSVM_allMT2.root";
//	TString inputname 	= "BEfficiencies_histos_HT_eq2jets_0l_CSVM_allMT2.root";
//	TString outputname 	= "BEfficiencies_HT_eq2jets_ge1l_CSVM_allMT2.root";
//	TString inputname 	= "BEfficiencies_histos_HT_eq2jets_ge1l_CSVM_allMT2.root";
//	TString outputname 	= "BEfficiencies_HT_ge3jets_0l_CSVM_MT2lt125.root";
//	TString inputname 	= "BEfficiencies_histos_HT_ge3jets_0l_CSVM_MT2lt125.root";
//	TString outputname 	= "BEfficiencies_HT_ge3jets_ge1l_CSVM_MT2lt125.root";
//	TString inputname 	= "BEfficiencies_histos_HT_ge3jets_ge1l_CSVM_MT2lt125.root";
//	TString outputname 	= "BEfficiencies_HT_ge6jets_0l_CSVM_MT2lt125.root";
//	TString inputname 	= "BEfficiencies_histos_HT_ge6jets_0l_CSVM_MT2lt125.root";
//	TString outputname 	= "BEfficiencies_HT_ge6jets_ge1l_CSVM_MT2lt125.root";
//	TString inputname 	= "BEfficiencies_histos_HT_ge6jets_ge1l_CSVM_MT2lt125.root";
//	TString outputname 	= "BEfficiencies_HT_eq2jets_0l_CSVM_MT2lt125.root";
//	TString inputname 	= "BEfficiencies_histos_HT_eq2jets_0l_CSVM_MT2lt125.root";
//	TString outputname 	= "BEfficiencies_HT_eq2jets_ge1l_CSVM_MT2lt125.root";
//	TString inputname 	= "BEfficiencies_histos_HT_eq2jets_ge1l_CSVM_MT2lt125.root";
//	TString outputname 	= "BEfficiencies_HT_ge3jets_0l_CSVM_MT2ge125.root";
//	TString inputname 	= "BEfficiencies_histos_HT_ge3jets_0l_CSVM_MT2ge125.root";
//	TString outputname 	= "BEfficiencies_HT_ge3jets_ge1l_CSVM_MT2ge125.root";
//	TString inputname 	= "BEfficiencies_histos_HT_ge3jets_ge1l_CSVM_MT2ge125.root";
//	TString outputname 	= "BEfficiencies_HT_ge6jets_0l_CSVM_MT2ge125.root";
//	TString inputname 	= "BEfficiencies_histos_HT_ge6jets_0l_CSVM_MT2ge125.root";
//	TString outputname 	= "BEfficiencies_HT_ge6jets_ge1l_CSVM_MT2ge125.root";
//	TString inputname 	= "BEfficiencies_histos_HT_ge6jets_ge1l_CSVM_MT2ge125.root";
//	TString outputname 	= "BEfficiencies_HT_eq2jets_0l_CSVM_MT2ge125.root";
//	TString inputname 	= "BEfficiencies_histos_HT_eq2jets_0l_CSVM_MT2ge125.root";
//	TString outputname 	= "BEfficiencies_HT_eq2jets_ge1l_CSVM_MT2ge125.root";
//	TString inputname 	= "BEfficiencies_histos_HT_eq2jets_ge1l_CSVM_MT2ge125.root";
//	TString outputname 	= "BEfficiencies_MET_ge3jets_0l_CSVM_allMT2.root";
//	TString inputname 	= "BEfficiencies_histos_MET_ge3jets_0l_CSVM_allMT2.root";
//	TString outputname 	= "BEfficiencies_MET_ge3jets_ge1l_CSVM_allMT2.root";
//	TString inputname 	= "BEfficiencies_histos_MET_ge3jets_ge1l_CSVM_allMT2.root";
//	TString outputname 	= "BEfficiencies_MET_ge6jets_0l_CSVM_allMT2.root";
//	TString inputname 	= "BEfficiencies_histos_MET_ge6jets_0l_CSVM_allMT2.root";
//	TString outputname 	= "BEfficiencies_MET_ge6jets_ge1l_CSVM_allMT2.root";
//	TString inputname 	= "BEfficiencies_histos_MET_ge6jets_ge1l_CSVM_allMT2.root";
//	TString outputname 	= "BEfficiencies_MET_eq2jets_0l_CSVM_allMT2.root";
//	TString inputname 	= "BEfficiencies_histos_MET_eq2jets_0l_CSVM_allMT2.root";
//	TString outputname 	= "BEfficiencies_MET_eq2jets_ge1l_CSVM_allMT2.root";
//	TString inputname 	= "BEfficiencies_histos_MET_eq2jets_ge1l_CSVM_allMT2.root";
//	TString outputname 	= "BEfficiencies_MET_ge3jets_0l_CSVM_MT2lt125.root";
//	TString inputname 	= "BEfficiencies_histos_MET_ge3jets_0l_CSVM_MT2lt125.root";
//	TString outputname 	= "BEfficiencies_MET_ge3jets_ge1l_CSVM_MT2lt125.root";
//	TString inputname 	= "BEfficiencies_histos_MET_ge3jets_ge1l_CSVM_MT2lt125.root";
//	TString outputname 	= "BEfficiencies_MET_ge6jets_0l_CSVM_MT2lt125.root";
//	TString inputname 	= "BEfficiencies_histos_MET_ge6jets_0l_CSVM_MT2lt125.root";
//	TString outputname 	= "BEfficiencies_MET_ge6jets_ge1l_CSVM_MT2lt125.root";
//	TString inputname 	= "BEfficiencies_histos_MET_ge6jets_ge1l_CSVM_MT2lt125.root";
//	TString outputname 	= "BEfficiencies_MET_eq2jets_0l_CSVM_MT2lt125.root";
//	TString inputname 	= "BEfficiencies_histos_MET_eq2jets_0l_CSVM_MT2lt125.root";
//	TString outputname 	= "BEfficiencies_MET_eq2jets_ge1l_CSVM_MT2lt125.root";
//	TString inputname 	= "BEfficiencies_histos_MET_eq2jets_ge1l_CSVM_MT2lt125.root";
//	TString outputname 	= "BEfficiencies_MET_ge3jets_0l_CSVM_MT2ge125.root";
//	TString inputname 	= "BEfficiencies_histos_MET_ge3jets_0l_CSVM_MT2ge125.root";
//	TString outputname 	= "BEfficiencies_MET_ge3jets_ge1l_CSVM_MT2ge125.root";
//	TString inputname 	= "BEfficiencies_histos_MET_ge3jets_ge1l_CSVM_MT2ge125.root";
//	TString outputname 	= "BEfficiencies_MET_ge6jets_0l_CSVM_MT2ge125.root";
//	TString inputname 	= "BEfficiencies_histos_MET_ge6jets_0l_CSVM_MT2ge125.root";
//	TString outputname 	= "BEfficiencies_MET_ge6jets_ge1l_CSVM_MT2ge125.root";
//	TString inputname 	= "BEfficiencies_histos_MET_ge6jets_ge1l_CSVM_MT2ge125.root";
//	TString outputname 	= "BEfficiencies_MET_eq2jets_0l_CSVM_MT2ge125.root";
//	TString inputname 	= "BEfficiencies_histos_MET_eq2jets_0l_CSVM_MT2ge125.root";
//	TString outputname 	= "BEfficiencies_MET_eq2jets_ge1l_CSVM_MT2ge125.root";
//	TString inputname 	= "BEfficiencies_histos_MET_eq2jets_ge1l_CSVM_MT2ge125.root";

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

        TFile *infile = TFile::Open(inputdir + inputname);
	map<string, TH2D*> histos2d;
	map<string, TH1D*> histospt;

	map<string, TH2F*> effhistos2d;
	map<string, TH1F*> effhistospt;
	map<string, TEfficiency*> effteffpt;
	map<string, TEfficiency*> effteff2d;


	const int sampletypesize = 10;
	string sample_type[sampletypesize] = {"QCD", "WJets", "ZJets", "TTbar", "SingleTop", "Other", "mc", "mcNOqcd", "susy", "data"};//same as type in samples.dat

	string prefix[5] = {"B", "C", "L", "NonB", "All"};

	for(int is = 0; is<sampletypesize; ++is){
	   for(int js = 0; js<5; ++js){
		string hs = string("_") + sample_type[is];
		string mapname;
		mapname = prefix[js] + "Jet" + hs;
		if(histospt.count(mapname) == 0 ) histospt[mapname] = (TH1D*)infile->Get(mapname.c_str());
		if(histos2d.count(mapname) == 0 ) histos2d[mapname] = (TH2D*)infile->Get((mapname+string("_2d")).c_str());
		mapname = "BTagged" + prefix[js] + "Jet" + hs;
		if(histospt.count(mapname) == 0 ) histospt[mapname] = (TH1D*)infile->Get(mapname.c_str());
		if(histos2d.count(mapname) == 0 ) histos2d[mapname] = (TH2D*)infile->Get((mapname+string("_2d")).c_str());
		if(js>=5) continue;
		mapname = prefix[js] + "Efficiency" + hs;
		if(WP==0 && effhistos2d.count(mapname) == 0 ) effhistos2d[mapname] = new TH2F((mapname+string("_2d")).c_str(), "eff", Nptbins, ptbins, NetabinsL, etabinsL);
		if(WP==1 && effhistos2d.count(mapname) == 0 ) effhistos2d[mapname] = new TH2F((mapname+string("_2d")).c_str(), "eff", Nptbins, ptbins, NetabinsM, etabinsM);
		if(WP==2 && effhistos2d.count(mapname) == 0 ) effhistos2d[mapname] = new TH2F((mapname+string("_2d")).c_str(), "eff", Nptbins, ptbins, NetabinsT, etabinsT);
		if(         effhistospt.count(mapname) == 0 ) effhistospt[mapname] = new TH1F( mapname.c_str(),                "eff", Nptbins, ptbins);
		mapname = "TEff_"+mapname;
		string passname = "BTagged" + prefix[js] + "Jet" + hs;
		string totname  = prefix[js] + "Jet" + hs;
		if(effteffpt.count(mapname) == 0 ) effteffpt[mapname] = new TEfficiency((*(histospt[passname])), (*(histospt[totname])) );
		effteffpt[mapname]->SetNameTitle((mapname).c_str(), "eff");
		if(effteff2d.count(mapname) == 0 ) effteff2d[mapname] = new TEfficiency((*(histos2d[passname])), (*(histos2d[totname])) );
		effteff2d[mapname]->SetNameTitle((mapname+string("_2d")).c_str(), "eff");
	   }
	}
	for(int is = 0; is<sampletypesize; ++is){
	   for(int js = 0; js<5; ++js){
		string hs = string("_") + sample_type[is];
		string mapname  = "TEff_" + prefix[js] + "Efficiency" + hs;
		string mapnameh = prefix[js] + "Efficiency" + hs;
		string passname = "BTagged" + prefix[js] + "Jet" + hs;
		string totname  = prefix[js] + "Jet" + hs;
		for(int i = 1; i<=histospt[ totname]->GetNbinsX(); ++i){
			int bin = i;
			double eff    = effteffpt[mapname]->GetEfficiency(bin);
			double errup  = effteffpt[mapname]->GetEfficiencyErrorUp(bin);
			double errlow = effteffpt[mapname]->GetEfficiencyErrorLow(bin);
			double err = errup;
			if(errlow>errup) err = errlow;
			if(eff==1||eff==0||err==1||err==0){
				double err1, err2;
				double eff1 = histospt[passname]->IntegralAndError(1,histospt[ totname]->GetNbinsX(), err1);
				double eff2 = histospt[ totname]->IntegralAndError(1,histospt[ totname]->GetNbinsX(), err2);
				if(eff2!=0){
				if(eff==1||eff==0) eff = eff1/eff2;
				if(err==1||err==0) err = sqrt(pow(err1/eff2,2)+pow((err2*eff1)/(eff2*eff2),2));
				}
			}
			//if(eff==1||eff==0) eff = histospt[passname]->Integral() / histospt[totname]->Integral();
			effhistospt[mapnameh]->SetBinContent(bin, eff);
			effhistospt[mapnameh]->SetBinError(  bin, err);
		}
		for(int i = 1; i<=histos2d[ totname]->GetNbinsX(); ++i){
		   for(int j = 1; j<=histos2d[ totname]->GetNbinsY(); ++j){
			int bin = histos2d[ totname]->GetBin(i,j);
			double eff    = effteff2d[mapname]->GetEfficiency(bin);
			double errup  = effteff2d[mapname]->GetEfficiencyErrorUp(bin);
			double errlow = effteff2d[mapname]->GetEfficiencyErrorLow(bin);
			double err = errup;
			if(errlow>errup) err = errlow;
			if(eff==1||eff==0||err==1||err==0){
				double err1, err2;
				double eff1 = histos2d[passname]->IntegralAndError(1,histos2d[ totname]->GetNbinsX(), 1,histos2d[ totname]->GetNbinsY(), err1);
				double eff2 = histos2d[ totname]->IntegralAndError(1,histos2d[ totname]->GetNbinsX(), 1,histos2d[ totname]->GetNbinsY(), err2);
				if(eff2!=0){
				if(eff==1||eff==0) eff = eff1/eff2;
				if(err==1||err==0) err = sqrt(pow(err1/eff2,2)+pow((err2*eff1)/(eff2*eff2),2));
				}
			}
			//if(eff==1||eff==0) eff = histos2d[passname+"_2d"]->Integral() / histospt[totname+"_2d"]->Integral();
			effhistos2d[mapnameh]->SetBinContent(bin, eff);
			effhistos2d[mapnameh]->SetBinError(  bin, err);
		   }
		}
		//effhistos2d[mapnameh]->Divide(histos2d[passname],histos2d[ totname],1,1,"B");
	   }
	}

	TFile *fsavefile = new TFile(outputdir+outputname,"RECREATE");
	fsavefile->cd();
	for(map<string,TH2F*>::iterator h=effhistos2d.begin(); h!=effhistos2d.end();++h){
		h->second->Write();}
	for(map<string,TH1F*>::iterator h=effhistospt.begin(); h!=effhistospt.end();++h){
		h->second->Write();}
	for(map<string,TEfficiency*>::iterator h=effteff2d.begin(); h!=effteff2d.end();++h){
		h->second->Write();}
	for(map<string,TEfficiency*>::iterator h=effteffpt.begin(); h!=effteffpt.end();++h){
		h->second->Write();}
	fsavefile->Close();

	cout << "Saved histograms in " << outputdir << outputname << endl;

}//void GetBEfficiencies_THistos()



