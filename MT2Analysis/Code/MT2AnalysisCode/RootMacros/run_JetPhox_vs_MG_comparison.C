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

void run_JetPhox_vs_MG_comparison() {

	int fVerbose =1;
	const float lumi=1;
	const int maxnsamples= 10;
	Sample *photons[10];
	float fMinPt=50;
	float fUncertainty=0.20; // with in percent of error band for ratio plit
	
	photons[0]      =new Sample("/shome/pnef/MT2Analysis/MT2trees/MT2_V01-03-02/20111122_MC_nocuts_removedPhoton/GenPhotonPtgt50/GJets_TuneZ2_40_HT_100_7TeV-madgraph-Summer11-PU_S4.root" , 12730972, 25690);
	photons[1]      =new Sample("/shome/pnef/MT2Analysis/MT2trees/MT2_V01-03-02/20111122_MC_nocuts_removedPhoton/GenPhotonPtgt50/GJets_TuneZ2_100_HT_200_7TeV-madgraph-Summer11-PU_S4.root",  1536287,  5213);
	photons[2]      =new Sample("/shome/pnef/MT2Analysis/MT2trees/MT2_V01-03-02/20111122_MC_nocuts_removedPhoton/GenPhotonPtgt50/GJets_TuneZ2_200_HT_inf_7TeV-madgraph-Summer11-PU_S4.root",  9377170,   798.3*1.1); // upscaled by 10% to get a smooth intercalibration in partonic HT

	TFile* fJetPhox1     = new TFile("/shome/pnef/MT2Analysis/Code/Misc/JetPhoxPhotons/jetphox_dironef_MT2cut.root"      ,"OPEN");
	TFile* fJetPhox2     = new TFile("/shome/pnef/MT2Analysis/Code/Misc/JetPhoxPhotons/jetphox_dironef_MT2cut_pt200.root","OPEN");
	TFile* fJetPhox3     = new TFile("/shome/pnef/MT2Analysis/Code/Misc/JetPhoxPhotons/jetphox_dironef_MT2cut_pt400.root","OPEN");
	TH1F*  hJetPhoxOrig1  = (TH1F*) fJetPhox1 ->Get("h21");
	TH1F*  hJetPhoxOrig2  = (TH1F*) fJetPhox2 ->Get("h21");
	TH1F*  hJetPhoxOrig3  = (TH1F*) fJetPhox3 ->Get("h21");

	// GetJetPhox Pt Spectrum starting from 50 GeV
	TH1F* h_JetPhox1 = new TH1F("h_JetPhox",       "", 100, 0, 1000); h_JetPhox1->Sumw2();
	TH1F* h_JetPhox2 = new TH1F("h_JetPhox_Pt200", "", 100, 0, 1000); h_JetPhox2->Sumw2();
	TH1F* h_JetPhox3= hJetPhoxOrig3->Clone("h_JetPhox_Pt400");  h_JetPhox3->Sumw2();
	for(int i=1; i<h_JetPhox1->GetNbinsX(); ++i){
		if(hJetPhoxOrig1->GetBinLowEdge(i)>=fMinPt && hJetPhoxOrig1->GetBinLowEdge(i)<200){
			h_JetPhox1->SetBinContent(i,hJetPhoxOrig1->GetBinContent(i));
			h_JetPhox1->SetBinError  (i,hJetPhoxOrig1->GetBinError(i));
		}
		if(hJetPhoxOrig2->GetBinLowEdge(i)>=200 && hJetPhoxOrig2->GetBinLowEdge(i)<400){
			h_JetPhox2->SetBinContent(i,hJetPhoxOrig2->GetBinContent(i));
			h_JetPhox2->SetBinError  (i,hJetPhoxOrig2->GetBinError(i));
		}
	}

	std::ostringstream cutStream1, cutStream2;
	cutStream1 << " " 
                   << "GenPhoton[0].Pt()>" << fMinPt       << "&&"
	           << "misc.CrazyHCAL==0";
	
  	TString cuts1   = cutStream1.str().c_str();

	TH1D* h_photon0 =  GetHisto(photons[0]->t, "GenPhoton[0].Pt()>>", cuts1 , "h_photon0",lumi*photons[0]->xsection/photons[0]->nevents, false, 100, 0, 1000);
	TH1D* h_photon1 =  GetHisto(photons[1]->t, "GenPhoton[0].Pt()>>", cuts1 , "h_photon1",lumi*photons[1]->xsection/photons[1]->nevents, false, 100, 0, 1000);
	TH1D* h_photon2 =  GetHisto(photons[2]->t, "GenPhoton[0].Pt()>>", cuts1 , "h_photon2",lumi*photons[2]->xsection/photons[2]->nevents, false, 100, 0, 1000);
	
	h_photon0->SetLineColor(kViolet+2);
	h_photon1->SetLineColor(kViolet+2);
	h_photon2->SetLineColor(kViolet+2);
	h_JetPhox1 ->SetMarkerStyle(20);
	h_JetPhox2 ->SetMarkerStyle(20);
	h_JetPhox3 ->SetMarkerStyle(20);
	
	TLegend* leg = new TLegend(0.2, 0.6, 0.5, 0.9);	
	leg->AddEntry(h_photon0,h_photon0->GetName(),"l" );
	leg->AddEntry(h_photon1,h_photon1->GetName(),"l" );
	leg->AddEntry(h_photon2,h_photon2->GetName(),"l" );
	leg->AddEntry(h_JetPhox1,"JetPhox NLO"      , "l" );
	leg->AddEntry(h_JetPhox2,"JetPhox NLO pt >200"  , "l" );
	leg->AddEntry(h_JetPhox3,"JetPhox NLO pt >400"  , "l" );
	leg -> SetFillColor(0);
	leg -> SetBorderSize(0);

	TCanvas *c1 = new TCanvas("c1", "c1", 500, 500);
	c1->cd();
	h_photon0->Draw("hist");
	h_photon1->Draw("histsame");
	h_photon2->Draw("histsame");
	h_JetPhox1->Draw("EX0same");
	h_JetPhox2->Draw("EX0same");
	h_JetPhox3->Draw("EX0same");
	leg->Draw();

	// get summed histograms and make ratio
	TCanvas *c2 = new TCanvas("c2", "c2", 500, 500);
	TH1D* h_photon = (TH1D*) h_photon0->Clone("h_photon");
	h_photon->Add(h_photon1);
	h_photon->Add(h_photon2);
	
	TH1F* h_JetPhox = (TH1F*) h_JetPhox1->Clone("h_JetPhox");
	h_JetPhox->Add(h_JetPhox2);
	h_JetPhox->Add(h_JetPhox3);

	TLegend* leg2 = new TLegend(0.2, 0.6, 0.5, 0.9);	
	leg2->AddEntry(h_photon,"MG Gamma+Jets","p" );
	leg2->AddEntry(h_JetPhox ,"JetPhox","p" );
	leg2 -> SetFillColor(0);
	leg2 -> SetBorderSize(0);


	cout << "******************************************" << endl;
	cout << "h_JetPhox->Integral() " << h_JetPhox->Integral() << endl;
	cout << "h_photon-> Integral() " << h_photon->Integral() << endl;
	cout << "MG to JetPhox ratio   " << h_photon->Integral()/h_JetPhox->Integral() << endl;
	cout << "******************************************" << endl;


	h_photon ->Scale(h_JetPhox->Integral()/h_photon->Integral()); // scale MG to JetPhox

	h_photon->SetXTitle("#gamma-Pt (GeV)");
	h_photon->SetYTitle("#frac{d#sigma}{dp_{T}}");
	h_photon->SetMarkerStyle(20);
	h_photon->SetMarkerStyle(20);
	h_photon->SetMarkerColor(kMagenta+2);
	h_photon->SetLineColor(kMagenta+2);
	h_photon->Draw("EX0");
	h_JetPhox->SetMarkerStyle(4);
	h_JetPhox->Draw("EX0same");
	gPad->SetLogy(1);
	leg2->Draw();


	//Ratio plot
	TCanvas *c3 = new TCanvas("c3", "c3", 500, 500);
	TH1D* h_ratio = (TH1D*) h_JetPhox ->Clone("h_ratio");	
	h_ratio->Divide(h_photon);
	h_ratio->Draw("EX0");
	

	// Variable bin ratio plot
	TCanvas *c4 = new TCanvas("c4", "c4", 500, 500);
	const int nbins = 9;
	double  bins[nbins+1]   = {0, 50, 60, 80, 120, 200, 280, 400, 600, 800};

	TH1D* h_photon_variable  = new TH1D("h_photon_variable" , "", nbins, bins);
	TH1D* h_jetphox_variable = new TH1D("h_jetphox_variable", "", nbins, bins);
	for(int i=1; i<=h_photon->GetNbinsX(); ++i){
		int binnumber=1;
		for(int b=1; b<=h_photon_variable->GetNbinsX(); ++b){
			if(h_photon->GetBinLowEdge(i)>=h_photon_variable->GetBinLowEdge(b)+h_photon_variable->GetBinWidth(b)){
				binnumber=b+1;
			}else if (h_photon->GetBinLowEdge(i)>=bins[nbins]){
				binnumber=0;
			} else{
				break;
			}		
		}
		if(fVerbose>2)cout << "variable bin number " << binnumber << " variable low edge " << h_photon_variable->GetBinLowEdge(binnumber) << endl;
		if(fVerbose>2)cout << "uniform  bin number " << i         << " uniform  low edge " << h_photon->GetBinLowEdge(i) << endl;
		double bincont   = h_photon_variable->GetBinContent(binnumber);
		double binerror  = h_photon_variable->GetBinError(binnumber);
		bincont  = bincont + h_photon->GetBinContent(i);
		binerror = sqrt(binerror*binerror +  h_photon->GetBinError(i)*h_photon->GetBinError(i));
		h_photon_variable->SetBinContent(binnumber,bincont);
		h_photon_variable->SetBinError  (binnumber,binerror);
		if(fVerbose>2) cout << "filling " << bincont  << " to bin " << binnumber << endl; 
		bincont   = h_jetphox_variable->GetBinContent(binnumber);
		binerror  = h_jetphox_variable->GetBinError(binnumber);
		bincont  = bincont + h_JetPhox->GetBinContent(i);
		binerror = sqrt(binerror*binerror +  h_JetPhox->GetBinError(i)*h_JetPhox->GetBinError(i));
		h_jetphox_variable->SetBinContent(binnumber,bincont);
		h_jetphox_variable->SetBinError  (binnumber,binerror);

	}

	h_jetphox_variable->Divide(h_photon_variable);	


	TH1D* hErrorbar = new TH1D("hErrorbar", "", 1,0, bins[nbins]);
	hErrorbar->SetYTitle("(JetPhox NLO)/(MadGraph)");
	hErrorbar->SetXTitle("#gamma-Pt (GeV)");
	hErrorbar->SetBinContent(1,1);
	hErrorbar->SetBinError(1,fUncertainty);
	hErrorbar->SetFillColor(5); 
	hErrorbar->SetFillStyle(3001);
	hErrorbar->SetMinimum(0);
	hErrorbar->SetMaximum(2);
	hErrorbar->Draw("e2");
	h_jetphox_variable->SetMarkerStyle(4);
	h_jetphox_variable->Draw("EXsame");



	

}

