#include "TFile.h"
#include "TLegend.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1D.h"
#include <sstream>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <math.h>
#include <vector>
#include <stdlib.h>
#include <stdio.h>
#include <map>
#include "../TESCO/include/helper/Utilities.hh"

using namespace std;

class Sample {
public:
	Sample(string f, float n, float x);
	~Sample();
	float nevents;
	float xsection;
	TFile* file;
	TTree* t;
	TH1D*  h1;
};

Sample::Sample(string f, float n, float x){
	nevents  = n;
	xsection = x;

	file     = TFile::Open(f.c_str());
	t        = (TTree*)file->Get("MassTree");
}

Sample::~Sample(){
	if(t!=0) delete t;
	if(file!=0) delete file;
}

TH1D* GetHisto(TTree* chain, TString var, TString basecut, TString name, double weight, bool normalized, int nbins, float min, float max ){
	TH1D* h = new TH1D(name,"", nbins, min, max );
	h->Sumw2();
	basecut =TString::Format("(%.10f)*(%s)", weight, basecut.Data());
	TString varname = var + h->GetName();
	cout << "+++ drawing " << var << " with cut " << basecut << endl; 
	int n = chain->Draw(varname,basecut,"goff");
	cout << name << " " << h->Integral() << endl;
	if(normalized) h->Scale(1./h->Integral());
	return h;
}

void run_PhotonsVsZll() {

	// User input -----------------------------------------------------
	const float lumi=4600;
	const int maxnsamples= 10;
	Sample *photons[10];
	Sample *Zll[10];
	int fVerbose=0;
	Bool_t fGenLevelMCGamma=false;
	Bool_t fGenLevelMCZll  =false;
	TString saveDir        ="GammaJetsPrediction/20111212_ZllVsPhotonComparison_test/";
	TString saveName       ="20111212_ZllVsPhotonComparison";
	Util::MakeOutputDir(saveDir);
	// User input -----------------------------------------------------

	TString RECOorGEN="";
	if(      fGenLevelMCGamma && fGenLevelMCZll ){RECOorGEN="_GammaGen_ZllGen";}
	else if(!fGenLevelMCGamma && fGenLevelMCZll ){RECOorGEN="_GammaReco_ZllGen";}
	else if( fGenLevelMCGamma && !fGenLevelMCZll){RECOorGEN="_GammaGen_ZllReco";}
	else if(!fGenLevelMCGamma && !fGenLevelMCZll){RECOorGEN="_GammaReco_ZllReco";}
	saveName+=RECOorGEN;

	// MC histos -----------------------------------------------
	photons[0] =new Sample("/shome/pnef/MT2Analysis/MT2trees/MT2_V01-03-06/20111212_ZllandGJets_nocuts/GenVBosonPtgt100_Skim/GJets_TuneZ2_40_HT_100_7TeV-madgraph-Summer11-PU_S4.root" , 12730972, 25690); 
	photons[1] =new Sample("/shome/pnef/MT2Analysis/MT2trees/MT2_V01-03-06/20111212_ZllandGJets_nocuts/GenVBosonPtgt100_Skim/GJets_TuneZ2_100_HT_200_7TeV-madgraph-Summer11-PU_S4.root",  1536287,  5213); 
	photons[2] =new Sample("/shome/pnef/MT2Analysis/MT2trees/MT2_V01-03-06/20111212_ZllandGJets_nocuts/GenVBosonPtgt100_Skim/GJets_TuneZ2_200_HT_inf_7TeV-madgraph-Summer11-PU_S4.root",  9377170,   798.3);
	Zll[0]   =new Sample("/shome/pnef/MT2Analysis/MT2trees/MT2_V01-03-06/20111212_ZllandGJets_nocuts/GenVBosonPtgt100_Skim/DYJetsToLL_PtZ-100_TuneZ2_7TeV-madgraph-tauola-Summer11-PU_S4.root",1133490, 25.1); // comparing LO xsection  
	
	std::ostringstream cutStream1, cutStream2;
	if(fGenLevelMCZll){
	cutStream1 << " " 
		   << "misc.MET >=0"                                      << "&&"
		   << "GenDiLeptPt(20,2.4,71,111,true)>149"               << "&&"
                   << "fabs(GenDiLeptRapidity(20,2.4,71,111,true))<1.4"    << "&&"
	           << "misc.CrazyHCAL==0";
	} else{
	cutStream1 << " " 
		   << "misc.MET >=0"                                      << "&&"
		   << "RecoOSDiLeptPt(20,2.4,71,111)>149"                 << "&&"
                   << "fabs(RecoOSDiLeptRapidity(20,2.4,71,111))<1.4"     << "&&"
	           << "misc.CrazyHCAL==0";
	}
	
	if(fGenLevelMCGamma){
	cutStream2 << " " 
		   << "misc.MET >=0"                     << "&&"
		   << "GenPhoton[0].Pt()>149"            << "&&"
                   << "fabs(GenPhoton[0].Rapidity())<1.4" << "&&"
	           << "misc.CrazyHCAL==0";
	}else{
	cutStream2 << " " 
		   << "misc.MET >=0"                     << "&&"
		   << "photon[0].lv.Pt()>149"            << "&&"
                   << "fabs(photon[0].lv.Rapidity())<1.4"   << "&&"
		   << "photon[0].isEGMtightID==1"                         << "&&"
	           << "misc.CrazyHCAL==0";
	}

  	TString cuts1   = cutStream1.str().c_str();
  	TString cuts2   = cutStream2.str().c_str();

	TH1D* h_photon0;
	TH1D* h_photon1;
	TH1D* h_photon2;
	TH1D* h_Zll0;

	if(fGenLevelMCGamma){
	h_photon0 =  GetHisto(photons[0]->t, "GenPhoton[0].Pt()>>",                  cuts2 , "h_photon0",lumi*photons[0]->xsection/photons[0]->nevents, false, 100, 100, 800);
	h_photon1 =  GetHisto(photons[1]->t, "GenPhoton[0].Pt()>>",                  cuts2 , "h_photon1",lumi*photons[1]->xsection/photons[1]->nevents, false, 100, 100, 800);
	h_photon2 =  GetHisto(photons[2]->t, "GenPhoton[0].Pt()>>",                  cuts2 , "h_photon2",lumi*photons[2]->xsection/photons[2]->nevents, false, 100, 100, 800);
	}else{
	h_photon0 =  GetHisto(photons[0]->t, "photon[0].lv.Pt()>>",                  cuts2 , "h_photon0",lumi*photons[0]->xsection/photons[0]->nevents, false, 100, 100, 800);
	h_photon1 =  GetHisto(photons[1]->t, "photon[0].lv.Pt()>>",                  cuts2 , "h_photon1",lumi*photons[1]->xsection/photons[1]->nevents, false, 100, 100, 800);
	h_photon2 =  GetHisto(photons[2]->t, "photon[0].lv.Pt()>>",                  cuts2 , "h_photon2",lumi*photons[2]->xsection/photons[2]->nevents, false, 100, 100, 800);
	}
	if(fGenLevelMCZll){
	h_Zll0    =  GetHisto(Zll[0]->t,     "GenDiLeptPt(20,2.4,71,111,true)>>",    cuts1 , "h_Zll0",   lumi*Zll[0]->xsection/Zll[0]->nevents,         false, 100, 100, 800);
	}else{
	h_Zll0    =  GetHisto(Zll[0]->t,     "RecoOSDiLeptPt(20,2.4,71,111)>>"  ,    cuts1 , "h_Zll0",   lumi*Zll[0]->xsection/Zll[0]->nevents,         false, 100, 100, 800);
	}
	
	// data histos --------------------------------------------
	
	TFile* f  = new TFile("/shome/pnef/MT2Analysis/Code/Misc/PhotonsVsZll/Kostas_bosonPt_4.6fb-1.root","OPEN");
	TH1F*  h_DataOrigZll   = (TH1F*) f->Get("zpt");
	TH1F*  h_DataOrigGamma = (TH1F*) f->Get("gpt");
	
	TH1F*  h_DataZll   = h_DataOrigZll  ->Clone("h_DataZll");
	TH1F*  h_DataGamma = h_DataOrigGamma->Clone("h_DataGamma");
	// Start plotting from 149 GeV:
	for(int i=0; i<h_DataZll->GetNbinsX(); ++i){ // start plotting from 150
		if(h_DataZll  ->GetBinLowEdge(i)<149) {h_DataZll  ->SetBinContent(i,0); h_DataZll  ->SetBinError(i,0);}
		if(h_DataGamma->GetBinLowEdge(i)<149) {h_DataGamma->SetBinContent(i,0); h_DataGamma->SetBinError(i,0);}
	}
	
	h_DataZll  ->Rebin(2); // get 100 bins from 100, to 800
	h_DataGamma->Rebin(2);

//	h_photon0->SetLineColor(kCyan);
//	h_photon1->SetLineColor(kCyan+1);
//	h_photon2->SetLineColor(kCyan+2);
//	h_Zll0->SetLineColor(kOrange);
//	
//	TLegend* leg = new TLegend(0.2, 0.6, 0.5, 0.9);	
//	leg->AddEntry(h_photon0,h_photon0->GetName(),"l" );
//	leg->AddEntry(h_photon1,h_photon1->GetName(),"l" );
//	leg->AddEntry(h_photon2,h_photon2->GetName(),"l" );
//	leg->AddEntry(h_Zll0 ,h_Zll0 ->GetName(),"l" );
//	leg -> SetFillColor(0);
//	leg -> SetBorderSize(0);
//
//	TCanvas *c1 = new TCanvas("c1", "c1", 500, 500);
//	c1->cd();
//	h_photon0->Draw("hist");
//	h_photon1->Draw("histsame");
//	h_photon2->Draw("histsame");
//	h_Zll0-> Draw("histsame");
//	leg->Draw();

	// get summed histograms and make ratio
	TCanvas *c2 = new TCanvas("Zll_vs_Photons", "", 500, 500);
	TH1D* h_photon = (TH1D*) h_photon0->Clone("h_photon");
	h_photon->Add(h_photon1);
	h_photon->Add(h_photon2);
	TH1D* h_Zll = (TH1D*) h_Zll0->Clone("h_Zll");
	
	// normalize Photons MC to data
	h_photon->Scale(h_DataGamma->Integral()/h_photon->Integral());
	h_Zll   ->Scale(h_DataZll  ->Integral()/h_Zll->Integral());
	
	h_Zll->SetMarkerStyle(22);
	h_Zll->SetLineColor(kOrange);
	h_Zll->SetMarkerColor(kOrange);
	h_Zll->SetMarkerSize(1.2);
	h_DataZll->SetMarkerStyle(26);
	h_DataZll->SetMarkerColor(kBlack);
	h_DataZll->SetLineColor(kBlack);
	h_photon->SetLineColor(kMagenta+2);
	h_photon->SetMarkerColor(kMagenta+2);
	h_photon->SetMarkerStyle(20);
	h_photon->SetMarkerSize(1.2);
	h_DataGamma->SetMarkerStyle(4);
	h_DataGamma->SetMarkerColor(kBlack);
	h_DataGamma->SetLineColor(kBlack);
	gPad->SetLogy(1);

	
	TLegend* leg2 = new TLegend(0.2, 0.6, 0.5, 0.9);	
	leg2->AddEntry(h_photon,"Gamma+Jets MC","p" );
	leg2->AddEntry(h_DataGamma,"Photon Data","p" );
	leg2->AddEntry(h_Zll ,"Zll MC","p" );
	leg2->AddEntry(h_DataZll,"Zll Data","p" );
	leg2 -> SetFillColor(0);
	leg2 -> SetBorderSize(0);
//	
	h_photon->SetXTitle("boson p_{T} (GeV)");
	h_photon->SetYTitle("Events / 7 GeV");
	h_photon->Draw("EX");
	h_DataGamma->Draw("EXsame");
	h_Zll ->Draw("EXsame");
	h_DataZll ->Draw("EXsame");
	leg2->Draw();

	TCanvas *c3 = new TCanvas("Zll_over_photon_ratio", "", 500, 500);
	TH1D* h_ratio     = (TH1D*) h_Zll     ->Clone("h_ratio");	
	TH1D* h_ratioData = (TH1D*) h_DataZll ->Clone("h_ratioData");	
	h_ratio->SetYTitle("#frac{Z(ll)}{#gamma}");
	h_ratio->SetXTitle("boson p_{T} (GeV)");
	h_ratio    ->Divide(h_photon);
	h_ratioData->Divide(h_DataGamma);
	h_ratio    ->SetMarkerStyle(20);
	h_ratio    ->SetMarkerSize(1.2);
	h_ratio    ->SetLineColor(kMagenta+2);
	h_ratio    ->SetMarkerColor(kMagenta+2);
	h_ratioData->SetMarkerStyle(26);
	h_ratioData->SetMarkerColor(kBlack);
	h_ratio    ->Draw("EX0");
	h_ratioData->Draw("EX0same");
	
	TLegend* leg3 = new TLegend(0.2, 0.6, 0.5, 0.9);	
	leg3->AddEntry(h_ratioData,"Zll/photon Data","p" );
	leg3->AddEntry(h_ratio    ,"Zll/photon MC","p" );
	leg3 -> SetFillColor(0);
	leg3 -> SetBorderSize(0);
	leg3 ->Draw();
	
	// save all
	TString fileName = saveDir+"/"+saveName+".root";
	TFile *file = new TFile(fileName.Data(), "RECREATE");
	c2->Write();
	c3->Write();
	h_ratio->Write();
	h_ratioData->Write();
	h_photon->Write();
	h_Zll->Write();
	h_DataGamma->Write();
	h_DataZll->Write();
	file->Close();
	TString plotname      =saveName;
	TString plotnameRatio =saveName+"Ratio";
	Util::Print(c2, plotname     , saveDir, 0);
	Util::Print(c3, plotnameRatio, saveDir, 0);
}

