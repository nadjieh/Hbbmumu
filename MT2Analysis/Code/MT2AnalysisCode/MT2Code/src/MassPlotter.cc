/*****************************************************************************
*   Small Class to make plots for MassAnalysis                               *
*****************************************************************************/

#include "MassPlotter.hh"

#include "TLorentzVector.h"
#include "TGraphAsymmErrors.h"
#include "helper/Utilities.hh"
#include "helper/Monitor.hh"

#include "TLatex.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TLine.h"
#include "TArrow.h"
#include "TBox.h"
#include "TMath.h"

#include "TString.h"
#include "TObject.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TH2.h"
#include "TH3.h"
#include "TF1.h"
#include "THStack.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TEventList.h"
#include "TCut.h"
#include "TTreeFormula.h"
#include "TStyle.h"
#include "TRandom.h"
#include "TROOT.h"
#include "TVirtualPad.h"
#include "TClonesArray.h"

#include "TGraph2D.h"
#include "TMatrixD.h"
#include "TMatrixT.h"
#include "TVectorD.h"
#include "TPaveStats.h"
#include <sstream>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <math.h>
#include <vector>
#include <stdlib.h>
#include <stdio.h>
#include <map>
#include <time.h> // access to date/time

using namespace std;

 Objectproperties::Objectproperties(TString type_):type(type_), nBins(10),lowVal(-1000.),highVal(0.){
		if(type == "pt"){	
			nBins = 100;
			lowVal = -2.;
			highVal = 1000.;
		}
                else if(type == "phi"){
                        nBins = 70;
                        lowVal = -3.5;
                        highVal = 3.5;
                }
                else if(type == "mult"){
                        nBins = 20;
                        lowVal = -0.5;
                        highVal = 19.5;
                }
		else std::cout<<"Bad Type!!!"<<std::endl;


	}
 Objectproperties::Objectproperties():type(""), nBins(10),lowVal(-1000.),highVal(0.){}
//____________________________________________________________________________
MassPlotter::MassPlotter(){
	fSave=true;
	fisPhoton=false;
	fPUReweight=true;
	fbSFReWeight=true;
// Default constructor, no samples are set
  cout.precision(2);
  cout.setf(ios::fixed,ios::floatfield);

}

//____________________________________________________________________________
MassPlotter::MassPlotter(TString outputdir){
// Explicit constructor with output directory
	setOutputDir(outputdir);
	fSave=true;
	fisPhoton=false;
	fPUReweight=true;
	fbSFReWeight=true;
  cout.precision(2);
  cout.setf(ios::fixed,ios::floatfield);

}

//____________________________________________________________________________
MassPlotter::MassPlotter(TString outputdir, TString outputfile){
// Explicit constructor with output directory and output file
	setOutputDir(outputdir);
	setOutputFile(outputfile);
	fSave=true;
	fisPhoton=false;
	fPUReweight=true;
	fbSFReWeight=true;
  cout.precision(2);
  cout.setf(ios::fixed,ios::floatfield);

}

//____________________________________________________________________________
MassPlotter::~MassPlotter(){
	fOutputFile->Close();
	delete fOutputFile;
}

//____________________________________________________________________________
void MassPlotter::init(TString filename){
	if(fVerbose > 0) cout << "------------------------------------" << endl;
	if(fVerbose > 0) cout << "Initializing MassPlotter ... " << endl;
	//Util::SetStyle();
	loadSamples(filename);

	// map for Z->nunu plots
	// to be used for Z->nunu           to be used for Z->ll
	RemoveLeptMap["misc.MET"]       = "Znunu.METplusLeptsPtReco";
}

//___________________________________________________________________________
void MassPlotter::makeSmallCopy(int nevents, int sample, TString cuts, TString trigger){
  // 	if(fSamples.size()<sample) return;


  cout<<" trigger "<<trigger<<endl;
  cout<<" cuts "<<cuts<<endl;

  for(int ii = 0; ii < fSamples.size(); ii++){
    if(fSamples.size() > sample && ii != sample)
      continue;
    //  if(fSamples[ii].type != "data" || fSamples[ii].name != "Run2012D1F")
    // continue;

    cout << "making small copy of tree " << fSamples[ii].name << endl;
 //      TFile *newfile = new TFile(fOutputDir+"/"+fSamples[ii].name+"_small.root","recreate");
    TString fileName = fSamples[ii].file->GetName();

    fileName =  fileName.ReplaceAll("_parkedsmallb.root", "_parkedverysmall.root");

    TFile *newfile = new TFile(fOutputDir+"/"+fileName,"recreate");
    TTree *newtree = fSamples[ii].tree->CloneTree(0);
    TH1F *h_PUWeights = (TH1F*) fSamples[ii].file->Get("h_PUWeights");
    TH1F *h_Events    = (TH1F*) fSamples[ii].file->Get("h_Events");
    
    if(fSamples[ii].file->Get("h_SMSEvents")){
      h_SMSEvents    = (TH2F*) fSamples[ii].file->Get("h_SMSEvents");
      h_mSugraEvents = (TH2F*) fSamples[ii].file->Get("h_mSugraEvents");
    }

    TString myCuts = cuts;
    if( fSamples[ii].type=="data") myCuts += " && " + trigger;
 

    fSamples[ii].tree->Draw(">>selList", myCuts);
    TEventList *myEvtList = (TEventList*)gDirectory->Get("selList");
    fSamples[ii].tree->SetEventList(myEvtList);


    Int_t nentries = myEvtList->GetN();

    for (Int_t i=0;i<nentries; i++) {
      if(i > nevents) continue;

      if ( i % 100000 == 0 ){ 	
	fprintf(stdout, "\rProcessed events: %6d of %6d ", i + 1, nentries);
	fflush(stdout);
      }
 
      fSamples[ii].tree->GetEntry(myEvtList->GetEntry(i));
      newtree->Fill();
   
    }
    cout << endl;
	
    h_PUWeights->Write();
    h_Events->Write();
    if(fSamples[ii].file->Get("h_SMSEvents")){
      h_SMSEvents->Write();
      h_mSugraEvents->Write();
    }

    newfile->Write();
	
    delete newfile;
  }
}

//___________________________________________________________________________
void MassPlotter::makePlots(){
	double dPhisplit[4]={0., 1.0, 1.5, 3.142};
	MakeMT2PredictionAndPlots(false, dPhisplit, 2.5);  // not cleaned, fudgefactor 2
}


//__________________________________________________________________________
void MassPlotter::MakeMT2PredictionAndPlots(bool cleaned , double dPhisplit[], double fudgefactor){
	TString cleanflag;
	if(cleaned) cleanflag="cleaned";
	else        cleanflag="notcleaned";
	std::ostringstream o0, o1, o2, o3;
	o0 << dPhisplit[0];
	o1 << dPhisplit[1];
	o2 << dPhisplit[2];
	o3 << dPhisplit[3];
	TString dPhirange1 = (TString) o0.str() + "to" +(TString) o1.str();
	TString dPhirange2 = (TString) o2.str() + "to" +(TString) o3.str();
	TString dPhirange  = (TString) o0.str()+(TString) o1.str() + (TString) o2.str() +(TString) o3.str();
	
	std::vector<sample>  QCDSamples;
	std::vector<sample>  Samples_NOTData;
	std::vector<sample>  QCDandDataSamples;
	for(int i=0; i< fSamples.size(); ++i){
		if(fSamples[i].sname=="QCD"){
			QCDSamples.push_back(fSamples[i]);
			QCDandDataSamples.push_back(fSamples[i]);
		}
		if(fSamples[i].type!="data"){
			Samples_NOTData.push_back(fSamples[i]);
		}
		if(fSamples[i].type=="data"){
			QCDandDataSamples.push_back(fSamples[i]);
		}
	}

	double ZerotoPi[4]    ={0., 3.142, 3.142, 3.142};
	double controlsplit[4]={dPhisplit[2], dPhisplit[3], dPhisplit[3], dPhisplit[3]};

	std::ostringstream cutStream;
	cutStream  
		  << "misc.MET > 30"                                << "&&"
		  << "misc.HT > 300"                                << "&&"
		  << "misc.Jet0Pass == 1"                           << "&&"
		  << "misc.Jet1Pass == 1"                           << "&&"
		  << "misc.PassJetID == 1"                          << "&&"
		  << "misc.Vectorsumpt<70"                          << "&&"
		  << "misc.MinMetJetDPhi>0.3"                       << "&&" 
		  << "misc.EcalDeadCellBEFlag==1"                   << "&&"
		  << "misc.HBHENoiseFlag == 1"                   ;

	TString cuts = cutStream.str().c_str();
	
	//                 variable                            cuts njets  nlepts HLT title     bins               flip_order  log  composite    ratio  stacked overlay 
	//	MakePlot(fSamples,"misc.MT2" ,                         cuts, -3,   0,    "",  "MT2" , gNMT2bins, gMT2bins , false,  true ,  true,      true,  true,  false);
//	MakePlot(fSamples,"hemi[0].lv1.M()" ,                  cuts, -3,   0,    "",  "h1.M", 30,       0, 1000  , false,  true ,  true,      true,  true,  false);
//	MakePlot(fSamples,"hemi[0].lv2.M()" ,                  cuts, -3,   0,    "",  "h1.M", 30,       0, 1000  , false,  true ,  true,      true,  true,  false);

}

//________________________________________________________________________

void MassPlotter::makeplot(TString var, TString maincuts, TString basecuts, int njets, int nbjets, int nleps, TString HLT,  TString xtitle,
			   const int nbins, const double min, const double max,
			   bool flip_order, bool logflag, bool composited, bool ratio, bool stacked, 
			   bool overlaySUSY, float overlayScale, bool add_underflow, bool saveHistos, TString saveMacro){

  MakePlot(fSamples, var, maincuts, basecuts, njets, nbjets, nleps, HLT, xtitle, nbins, min, max, flip_order, logflag, composited, ratio, stacked, overlaySUSY, overlayScale, add_underflow, saveHistos, saveMacro);
  
}

void MassPlotter::makePlot(TString var, TString cuts, int njets, int nbjets, int nleps, TString HLT,  TString xtitle,
			   const int nbins, const double min, const double max,
			   bool flip_order, bool logflag, bool composited, bool ratio, bool stacked, 
			   bool overlaySUSY, float overlayScale, bool add_underflow, bool saveHistos, TString saveMacro){

  MakePlot(fSamples, var, cuts, njets, nbjets, nleps, HLT, xtitle, nbins, min, max, flip_order, logflag, composited, ratio, stacked, overlaySUSY, overlayScale, add_underflow, saveHistos, saveMacro);
  
}

//________________________________________________________________________
void MassPlotter::Makeplot(TString var, TString maincuts, TString basecuts, int njets, int nbjets, int nleps, TString HLT, TString xtitle, 
			   const int nbins,  const double *bins,
			   bool flip_order, bool logflag, bool composited, bool ratio, bool stacked, 
			   bool overlaySUSY, float overlayScale, bool add_underflow, bool saveHistos, TString saveMacro){

	  MakePlot(fSamples, var, maincuts, basecuts, njets, nbjets, nleps, HLT, xtitle, nbins, bins, flip_order, logflag, composited, ratio, stacked, overlaySUSY, overlayScale, add_underflow, saveHistos, saveMacro);

}

void MassPlotter::MakePlot(TString var, TString cuts, int njets, int nbjets, int nleps, TString HLT, TString xtitle, 
			   const int nbins,  const double *bins,
			   bool flip_order, bool logflag, bool composited, bool ratio, bool stacked, 
			   bool overlaySUSY, float overlayScale, bool add_underflow, bool saveHistos, TString saveMacro){

	  MakePlot(fSamples, var, cuts, njets, nbjets, nleps, HLT, xtitle, nbins, bins, flip_order, logflag, composited, ratio, stacked, overlaySUSY, overlayScale, add_underflow, saveHistos, saveMacro);

}

// ________________________________________________________________________
void MassPlotter::PrintWEfficiency(int sample_index ,TString process,  std::string lept, Long64_t nevents, bool includeTaus){
	sample Sample = fSamples[sample_index];

	std::cout << setfill('=') << std::setw(70) << "" << std::endl;
	cout << "printing W efficiency: \n"
	     << Sample.name << endl;	
	std::cout << setfill('-') << std::setw(70) << "" << std::endl;

	enum counters_t { count_begin, all=count_begin, presel, NJetsIDLoose , PassJetID, MinMetJetDPhi, VectorSumPt, MT2,  count_end };
	Monitor counters[count_end];
	TString lablx[count_end] = {"all events", "presel", "NJetsIDLoose" , "PassJetID", "MinMetJetDPhi", "VectorSumPt", "MT2"};

	string to_measure;
	if     (lept == "ele" && process=="W" )       {to_measure = "W->enu  recoed";} 
	else if(lept == "muo" && process=="W" )       {to_measure = "W->munu recoed";} 
	else if(lept == "ele" && process=="Top")      {to_measure = "Top W->enu  recoed";}
	else if(lept == "muo" && process=="Top")      {to_measure = "Top W->munu  recoed";}

	fMT2tree = new MT2tree();
	Sample.tree->SetBranchAddress("MT2tree", &fMT2tree);
	Long64_t nentries =  Sample.tree->GetEntries();
	Long64_t nbytes = 0, nb = 0;
	for (Long64_t jentry=0; jentry<min(nentries, nevents);jentry++) {
		nb = Sample.tree->GetEntry(jentry);   nbytes += nb;
		Sample.tree->SetBranchAddress("MT2tree", &fMT2tree);
      
		if ( fVerbose>2 && jentry % 50000 == 0 )  cout << "+++ Proccessing event " << jentry << endl;

		bool leptfound(false);
		bool eventgood(false);
		bool acceptance(false);
		if(process=="W"){
			if(lept=="ele" && fMT2tree->GenLeptFromW(11, 0 , 1000,includeTaus)==true)                                              eventgood =true;
			if(lept=="muo" && fMT2tree->GenLeptFromW(13, 0 , 1000,includeTaus)==true)                                              eventgood =true;
			if(lept=="ele" && fMT2tree->GenLeptFromW(11,10,  2.4 ,includeTaus)==true)                                              acceptance=true;
			if(lept=="muo" && fMT2tree->GenLeptFromW(13,10,  2.4 ,includeTaus)==true)                                              acceptance=true;
			if(lept=="ele" && fMT2tree->NEles==1 && fMT2tree->NMuons==0 && fMT2tree->GenLeptFromW(11, 0, 1000,includeTaus)==true && fMT2tree->ele[0].lv.Pt()>10)  leptfound =true;
			if(lept=="muo" && fMT2tree->NMuons==1&& fMT2tree->NEles== 0 && fMT2tree->GenLeptFromW(13, 0, 1000,includeTaus)==true && fMT2tree->muo[0].lv.Pt()>10)  leptfound =true;
		}else if(process=="Top"){ 
			if(lept=="ele" && fMT2tree->TopDecayModeResult(11)==true)                                                 eventgood =true;
			if(lept=="muo" && fMT2tree->TopDecayModeResult(13)==true)                                                 eventgood =true;
			if(lept=="ele" && fMT2tree->TopDecayModeResult(11)==true && fMT2tree->SLTopAccept(10, 2.4)==true )        acceptance =true;
			if(lept=="muo" && fMT2tree->TopDecayModeResult(13)==true && fMT2tree->SLTopAccept(10, 2.4)==true )        acceptance =true;
			if(lept=="ele" && fMT2tree->NEles==1 && fMT2tree->NMuons==0 && eventgood &&fMT2tree->ele[0].lv.Pt()>10)             leptfound =true;
			if(lept=="muo" && fMT2tree->NMuons==1&& fMT2tree->NEles== 0 && eventgood &&fMT2tree->muo[0].lv.Pt()>10)             leptfound =true;
		}

		Double_t weight = fMT2tree->pileUp.Weight;

		// all events
		if(eventgood)   counters[all].fill("all events", weight);
		if(acceptance)  counters[all].fill("acceptance", weight);
		if(leptfound)   counters[all].fill(to_measure, weight);
		
		// presel
		if(fMT2tree->misc.MET                     < 30    )  continue;
		if(fMT2tree->misc.HT                      < 300   )  continue;
		if(fMT2tree->misc.Jet0Pass                ==0     )  continue;
		if(fMT2tree->misc.Jet1Pass                ==0     )  continue;
//		if(fMT2tree->misc.LeadingJPt              < 150   )  continue;
		if(fMT2tree->misc.SecondJPt               < 100   )  continue;
//		if(fMT2tree->NBJets                       < 1     )  continue;
//		if(fMT2tree->misc.HBHENoiseFlag           ==0     )  continue;
		if(fMT2tree->misc.CrazyHCAL               ==1     )  continue;
		if(eventgood)   counters[presel].fill("presel", weight);
		if(acceptance)  counters[presel].fill("acceptance", weight);
		if(leptfound)   counters[presel].fill(to_measure, weight);
		
		
		// NJetsIDLoose
		if(fMT2tree->NJetsIDLoose                 <3     ) continue;  
		if(eventgood)   counters[NJetsIDLoose].fill("NJetsIDLoose", weight);
		if(acceptance)  counters[NJetsIDLoose].fill("acceptance", weight);
		if(leptfound)   counters[NJetsIDLoose].fill(to_measure, weight);

		// PassJetID
		if(fMT2tree->misc.PassJetID               ==0     ) continue;
		if(eventgood)   counters[PassJetID].fill("PassJetID", weight);
		if(acceptance)  counters[PassJetID].fill("acceptance", weight);
		if(leptfound)   counters[PassJetID].fill(to_measure, weight);
		
		// MinMetJetDPhi
		if(fMT2tree->misc.MinMetJetDPhi <0.3 )       continue;
		if(eventgood)   counters[MinMetJetDPhi].fill("MinMetJetDPhi", weight);
		if(acceptance)  counters[MinMetJetDPhi].fill("acceptance", weight);
		if(leptfound)   counters[MinMetJetDPhi].fill(to_measure, weight);
		
		// VectorSumPt
		if(fMT2tree->misc.Vectorsumpt             >70     )  continue;
		if(eventgood)   counters[VectorSumPt].fill("VectorSumPt", weight);
		if(acceptance)  counters[VectorSumPt].fill("acceptance", weight);
		if(leptfound)   counters[VectorSumPt].fill(to_measure, weight);
		
		// MT2
//		if(fMT2tree->misc.MT2             <400     )  continue;
//		if(fMT2tree->misc.MT2             <150     )  continue;
		if(fMT2tree->misc.MT2 <200 || fMT2tree->misc.MT2>400 )  continue;
//		if(fMT2tree->misc.MT2 <150 || fMT2tree->misc.MT2>300 )  continue;
//		if(fMT2tree->misc.MT2 <100 || fMT2tree->misc.MT2>150  )  continue;
		if(eventgood)   counters[MT2].fill("MT2", weight);
		if(acceptance)  counters[MT2].fill("acceptance", weight);
		if(leptfound)   counters[MT2].fill(to_measure, weight);
		
	}

	// print stats     
	std::cout << setfill('=') << std::setw(70) << "" << std::endl;
	std::cout << "Statistics" << std::endl;
	std::cout << setfill('-') << std::setw(70) << "" << setfill(' ') << std::endl;
	for ( counters_t iCount=count_begin; iCount<count_end; iCount = counters_t(iCount+1) ) {
		counters[iCount].print();
		std::cout << setfill('-') << std::setw(70) << "" << setfill(' ') << std::endl;  
	}
	std::cout << setfill('=') << std::setw(70) << "" << std::endl;

	// fill histo
	TH1D* h_all = new TH1D("h_all", "", count_end, 0., (double) count_end );
	TH1D* h_acc = new TH1D("h_acc", "", count_end, 0., (double) count_end );
	TH1D* h_rec = new TH1D("h_rec", "", count_end, 0., (double) count_end );
	TH1D* h_1   = new TH1D("h_1"  , "", count_end, 0., (double) count_end );
	TH1D* h_2   = new TH1D("h_2"  , "", count_end, 0., (double) count_end );
	TH1D* h_3   = new TH1D("h_3"  , "", count_end, 0., (double) count_end );
	for(int i=0; i<count_end; ++i){
		h_1    ->GetXaxis()->SetBinLabel(i+1, lablx[i]);
		h_2    ->GetXaxis()->SetBinLabel(i+1, lablx[i]);
		h_3    ->GetXaxis()->SetBinLabel(i+1, lablx[i]);
		h_all->SetBinContent(i+1,counters[i].counts((string) lablx[i]));
		h_acc->SetBinContent(i+1,counters[i].counts("acceptance"));
		h_rec->SetBinContent(i+1,counters[i].counts(to_measure));
	}
	h_1    ->Sumw2();
	h_2    ->Sumw2();
	h_3    ->Sumw2();
	h_1    ->Divide(h_acc,h_all);	
	h_2    ->Divide(h_rec,h_acc);	
	h_3    ->Divide(h_rec,h_all);	
	TCanvas *col  = new TCanvas("Wevents", "", 0, 0, 900, 700);
	col -> cd();
	h_1    ->SetLineColor(kBlue);
	h_1    ->SetMarkerStyle(20);
	h_1    ->SetMinimum(0);
	h_1    ->SetMaximum(1);
	h_1    ->SetDrawOption("E");
	h_1    ->Draw();
	string name  ="Wevents_"+lept+"_acceptance_"+(string) Sample.name; 
	h_2    ->SetLineColor(kBlue);
	h_2    ->SetMinimum(0);
	h_2    ->SetMaximum(1);
	h_2    ->SetMarkerStyle(20);
	h_2    ->SetDrawOption("E");
	h_2    ->Draw();
	name  ="Wevents_"+lept+"_reco_"+(string) Sample.name; 
	h_3    ->SetLineColor(kBlue);
	h_3    ->SetMinimum(0);
	h_3    ->SetMaximum(1);
	h_3    ->SetMarkerStyle(20);
	h_3    ->SetDrawOption("E");
	h_3    ->Draw();
	name  ="Wevents_"+lept+"_prob_"+(string) Sample.name; 

	cout << "lept " << lept << " process " << process << endl;	
	//__________________________________
	if        (lept=="ele" && process=="W") {
		fWpred.Wenu_acc        =h_1->GetBinContent(count_end);
		fWpred.Wenu_acc_err    =h_1->GetBinError(count_end); 
		fWpred.Wenu_rec        =h_2->GetBinContent(count_end);
		fWpred.Wenu_rec_err    =h_2->GetBinError(count_end);  
		fWpred.Wenu_prob       =h_3->GetBinContent(count_end);
		fWpred.Wenu_prob_err   =h_3->GetBinError(count_end);  
	}else if  (lept=="ele" && process=="Top"){
		fWpred.TopWenu_acc     =h_1->GetBinContent(count_end);
		fWpred.TopWenu_acc_err =h_1->GetBinError(count_end); 
		fWpred.TopWenu_rec     =h_2->GetBinContent(count_end);
		fWpred.TopWenu_rec_err =h_2->GetBinError(count_end);  
		fWpred.TopWenu_prob    =h_3->GetBinContent(count_end);
		fWpred.TopWenu_prob_err=h_3->GetBinError(count_end);  
	}else if  (lept=="muo" && process=="W") {
		fWpred.Wmunu_acc       =h_1->GetBinContent(count_end);
		fWpred.Wmunu_acc_err   =h_1->GetBinError(count_end); 
		fWpred.Wmunu_rec       =h_2->GetBinContent(count_end);
		fWpred.Wmunu_rec_err   =h_2->GetBinError(count_end); 
		fWpred.Wmunu_prob      =h_3->GetBinContent(count_end);
		fWpred.Wmunu_prob_err  =h_3->GetBinError(count_end);  
	}else if  (lept=="muo" && process=="Top") {
		fWpred.TopWmunu_acc    =h_1->GetBinContent(count_end);
		fWpred.TopWmunu_acc_err=h_1->GetBinError(count_end); 
		fWpred.TopWmunu_rec    =h_2->GetBinContent(count_end);
		fWpred.TopWmunu_rec_err=h_2->GetBinError(count_end); 
		fWpred.TopWmunu_prob    =h_3->GetBinContent(count_end);
		fWpred.TopWmunu_prob_err=h_3->GetBinError(count_end);  
	}

	//_________________________________
	delete h_1;
	delete h_2;
	delete h_3;
	delete h_all;
	delete h_acc;
	delete h_rec;
	delete col;
}

//________________________________________________________________________

void MassPlotter::PrintCutFlow(int njets, int nleps, TString trigger, TString cuts){
  
  Monitor counters[fSamples.size()];

  const int   nProc      = 17;  
  TString  cnames[nProc] = {"QCD", "W+jets", "Z+jets", "Top","Other", "Total Bkg.", "data", "LM1", "LM2", "LM3", "LM4", "LM5", "LM8", "LM9", "LM11", "LM12", "LM13"};
  Monitor  ccount[nProc], ccount_100[nProc];

  for(size_t i = 0; i < fSamples.size(); ++i){

    Double_t sample_weight =0;
    if(fPUReweight) sample_weight= fSamples[i].xsection * fSamples[i].kfact * fSamples[i].lumi / (fSamples[i].nevents * fSamples[i].PU_avg_weight );
    else            sample_weight= fSamples[i].xsection * fSamples[i].kfact * fSamples[i].lumi / (fSamples[i].nevents );
    if(fVerbose>2) cout << "PrintCutFlow: looping over " << fSamples[i].name << endl;
    if(fVerbose>2) cout << "              sample has weight " << sample_weight << " and " << fSamples[i].tree->GetEntries() << " entries" << endl; 
    if(fVerbose>2) cout << "              Average PU weight: " << fSamples[i].PU_avg_weight << endl;
    if(fVerbose>2) cout << "              PU reweightung:    " << fPUReweight << endl;
    if(fVerbose>2) cout << "              Original Entries: " << fSamples[i].nevents << endl;
    
    fMT2tree = new MT2tree();
    fSamples[i].tree->SetBranchAddress("MT2tree", &fMT2tree);
    Long64_t nentries =  fSamples[i].tree->GetEntries();
    Long64_t nbytes = 0, nb = 0;
    int nev =0;

    //leo tweak - filtering out the TTree
    TString myCuts = cuts;

    if( fSamples[i].type=="data") myCuts += " && " + trigger; //cuts to be aplied only on data
    cout << "Cuts for Flow: " << myCuts << endl;
    fSamples[i].tree->Draw(">>selList", myCuts);


    TEventList *myEvtList = (TEventList*)gDirectory->Get("selList");
    fSamples[i].tree->SetEventList(myEvtList);
    int counter=0;
    cout << "Filtering done, size=" <<myEvtList->GetN()  << endl;
    
    if(myEvtList->GetSize()==0) continue;
    
    while(myEvtList->GetEntry(counter++) !=-1){
      
      int jentry = myEvtList->GetEntry(counter-1);
      
      //for (Long64_t jentry=0; jentry<nentries;jentry++) {
      nb =  fSamples[i].tree->GetEntry(jentry);   nbytes += nb;
      fSamples[i].tree->SetBranchAddress("MT2tree", &fMT2tree);
      
      if ( fVerbose>2 && jentry % 50000 == 0 )  cout << "+++ Proccessing event " << jentry << endl;
      
      Double_t weight = sample_weight;
      if (!fMT2tree->misc.isData ) weight = weight * fMT2tree->pileUp.Weight; // pile-up reweighting for MC 

       bool isMT2gt100 = fMT2tree->misc.MT2 > 100.;

      if( fMT2tree->NJetsIDLoose < 1 || !fMT2tree->misc.Jet0Pass                                      )  continue;
      if( fMT2tree->NJetsIDLoose < 2 || !fMT2tree->misc.Jet1Pass    || fMT2tree->misc.SecondJPt < 100 )  continue;
      if( !(!fMT2tree->misc.isData   || (fMT2tree->misc.Run <162803 || fMT2tree->misc.Run >162909  )) )  continue;
      if(trigger=="HT"){
	      if( fMT2tree->misc.isData ==1 
		  && fMT2tree->trigger.HLT_PFHT650_v5 ==0 
		  && fMT2tree->trigger.HLT_PFHT650_v6 ==0 
		  && fMT2tree->trigger.HLT_PFHT650_v7 ==0 ) continue;
	      if( fMT2tree->misc.MET     < 30)  continue;
      }else if(trigger=="MHT_HT"){
	      if( fMT2tree->misc.isData ==1 
		  && fMT2tree->trigger.HLT_PFHT350_PFMET100_v3 ==0 
		  && fMT2tree->trigger.HLT_PFHT350_PFMET100_v4 ==0 
		  && fMT2tree->trigger.HLT_PFHT350_PFMET100_v5 ==0 ) continue;
	      if( fMT2tree->misc.MET       <   30) continue;
	      if( fMT2tree->misc.HT        <  450) continue;
      }

      if     (njets>1) {
	if(fMT2tree->NJetsIDLoose != njets) continue;
	TString text = TString::Format("All events (jets == %d)",njets);
	counters[i].fill(text.Data(),weight);
	FillMonitor(ccount, fSamples[i].sname, fSamples[i].type, text, weight);
	if(isMT2gt100)     FillMonitor(ccount_100, fSamples[i].sname, fSamples[i].type, text, weight);
      }
      else if(njets <-1) {
	if(fMT2tree->NJetsIDLoose < abs(njets)) continue;
	TString text = TString::Format("All events(jets >= %d)",abs(njets));
	counters[i].fill(text.Data(),weight);
	FillMonitor(ccount, fSamples[i].sname, fSamples[i].type, text, weight);
	if(isMT2gt100)     FillMonitor(ccount_100, fSamples[i].sname, fSamples[i].type, text, weight);
      }

      if(fMT2tree->PassMinMetJetDPhi03() ==0)      continue;      
      counters[i].fill("Minimum DPhi(MET,jet) > 0.3",weight);
      FillMonitor(ccount, fSamples[i].sname, fSamples[i].type, "Minimum DPhi(MET,jet) > 0.3", weight);
      if(isMT2gt100)     FillMonitor(ccount_100, fSamples[i].sname, fSamples[i].type, "Minimum DPhi(MET,jet) > 0.3", weight);

      if( fMT2tree->misc.HBHENoiseFlag != 1 )  continue;
      if( fMT2tree->misc.CrazyHCAL     != 0 )  continue;
      counters[i].fill("HBHE noise veto",weight);
      FillMonitor(ccount, fSamples[i].sname, fSamples[i].type, "HBHE noise veto", weight);
      if(isMT2gt100)     FillMonitor(ccount_100, fSamples[i].sname, fSamples[i].type, "HBHE noise veto", weight);

//      if( fMT2tree->misc.EcalDeadCellBEFlag != 1 )  continue;
//      counters[i].fill("Boundary energy veto (dead ecal)",weight);
//      FillMonitor(ccount, fSamples[i].sname, fSamples[i].type, "Boundary energy veto (dead ecal)", weight);
//      if(isMT2gt100)     FillMonitor(ccount_100, fSamples[i].sname, fSamples[i].type, "Boundary energy veto (dead ecal)", weight);

      if( fMT2tree->misc.Vectorsumpt > 70. )  continue;
      counters[i].fill("VectorSumPt < 70",weight);
      FillMonitor(ccount, fSamples[i].sname, fSamples[i].type, "VectorSumPt < 70", weight);
      if(isMT2gt100)     FillMonitor(ccount_100, fSamples[i].sname, fSamples[i].type, "VectorSumPt < 70", weight);

      if( !fMT2tree->misc.PassJetID )  continue;
      counters[i].fill("jets > 50GeV failing PFID event veto",weight);
      FillMonitor(ccount, fSamples[i].sname, fSamples[i].type, "jets > 50GeV failing PFID event veto", weight);
      if(isMT2gt100)     FillMonitor(ccount_100, fSamples[i].sname, fSamples[i].type, "jets > 50GeV failing PFID event veto", weight);

      // Only considering (so far): no lepton requirement (<0);  1 lepton (==1), lepton veto (otherwise)
      string nLeps = (std::string) TString::Format("%d",abs(nleps));
      if(nleps !=-10){ 
	      if(nleps>0){ // exactly nleps lepton
		string cut_name = "NLeps == "+nLeps;
		if( (fMT2tree->NEles + fMT2tree->NMuons) != abs(nleps))  continue;
		counters[i].fill(cut_name,weight);
		FillMonitor(ccount, fSamples[i].sname, fSamples[i].type, cut_name, weight);
		if(isMT2gt100)     FillMonitor(ccount_100, fSamples[i].sname, fSamples[i].type, cut_name, weight);
	      }
	      else if( nleps==0 ) { // lepton veto
		string cut_name = "Lepton veto";
		if( (fMT2tree->NEles + fMT2tree->NMuons) != 0 )  continue;
		counters[i].fill(cut_name,weight);
		FillMonitor(ccount, fSamples[i].sname, fSamples[i].type, cut_name, weight);
		if(isMT2gt100)     FillMonitor(ccount_100, fSamples[i].sname, fSamples[i].type, cut_name, weight);
	      }
	      else if(nleps<0){ // at least nleps lepton
		string cut_name = "NLeps >= "+nLeps;
		if( (fMT2tree->NEles + fMT2tree->NMuons) < abs(nleps))  continue;
		counters[i].fill(cut_name ,weight);
		FillMonitor(ccount, fSamples[i].sname, fSamples[i].type, cut_name, weight);
		if(isMT2gt100)     FillMonitor(ccount_100, fSamples[i].sname, fSamples[i].type, cut_name, weight);
	      }
      }
      
//      if( !(fMT2tree->GetNBtags(2,1.74) >= 1) ) continue;
//      TString text = "b-tags >=1";
//      counters[i].fill(text.Data(),weight);
//      FillMonitor(ccount, fSamples[i].sname, fSamples[i].type, text, weight);
//      if(isMT2gt100)     FillMonitor(ccount_100, fSamples[i].sname, fSamples[i].type, text, weight);
      
      if( fMT2tree->misc.MT2 < 80. )  continue;
      counters[i].fill("MT2 > 80 GeV" ,weight);
      FillMonitor(ccount, fSamples[i].sname, fSamples[i].type, "MT2 > 80 GeV", weight);
      if( fMT2tree->misc.MT2 < 100. )  continue;
      counters[i].fill("MT2 > 100 GeV",weight);
      FillMonitor(ccount, fSamples[i].sname, fSamples[i].type, "MT2 > 100 GeV", weight);
      if( fMT2tree->misc.MT2 < 120. )  continue;
      counters[i].fill("MT2 > 120 GeV",weight);
      FillMonitor(ccount, fSamples[i].sname, fSamples[i].type, "MT2 > 120 GeV", weight);
      if( fMT2tree->misc.MT2 < 135. )  continue;
      counters[i].fill("MT2 > 135 GeV",weight);
      FillMonitor(ccount, fSamples[i].sname, fSamples[i].type, "MT2 > 135 GeV", weight);
      if( fMT2tree->misc.MT2 < 150. )  continue;
      counters[i].fill("MT2 > 150 GeV",weight);
      FillMonitor(ccount, fSamples[i].sname, fSamples[i].type, "MT2 > 150 GeV", weight);
      if( fMT2tree->misc.MT2 < 165. )  continue;
      counters[i].fill("MT2 > 165 GeV",weight);
      FillMonitor(ccount, fSamples[i].sname, fSamples[i].type, "MT2 > 165 GeV", weight);
      if( fMT2tree->misc.MT2 < 180. )  continue;
      counters[i].fill("MT2 > 180 GeV",weight);
      FillMonitor(ccount, fSamples[i].sname, fSamples[i].type, "MT2 > 180 GeV", weight);
      if( fMT2tree->misc.MT2 < 200. )  continue;
      counters[i].fill("MT2 > 200 GeV",weight);
      FillMonitor(ccount, fSamples[i].sname, fSamples[i].type, "MT2 > 200 GeV", weight);
      if( fMT2tree->misc.MT2 < 225. )  continue;
      counters[i].fill("MT2 > 225 GeV",weight);
      FillMonitor(ccount, fSamples[i].sname, fSamples[i].type, "MT2 > 225 GeV", weight);
      if( fMT2tree->misc.MT2 < 250. )  continue;
      counters[i].fill("MT2 > 250 GeV",weight);
      FillMonitor(ccount, fSamples[i].sname, fSamples[i].type, "MT2 > 250 GeV", weight);
      if( fMT2tree->misc.MT2 < 275. )  continue;
      counters[i].fill("MT2 > 275 GeV",weight);
      FillMonitor(ccount, fSamples[i].sname, fSamples[i].type, "MT2 > 275 GeV", weight);
      if( fMT2tree->misc.MT2 < 300. )  continue;
      counters[i].fill("MT2 > 300 GeV",weight);
      FillMonitor(ccount, fSamples[i].sname, fSamples[i].type, "MT2 > 300 GeV", weight);
      if( fMT2tree->misc.MT2 < 325. )  continue;
      counters[i].fill("MT2 > 325 GeV",weight);
      FillMonitor(ccount, fSamples[i].sname, fSamples[i].type, "MT2 > 325 GeV", weight);
      if( fMT2tree->misc.MT2 < 350. )  continue;
      counters[i].fill("MT2 > 350 GeV",weight);
      FillMonitor(ccount, fSamples[i].sname, fSamples[i].type, "MT2 > 350 GeV", weight);
      if( fMT2tree->misc.MT2 < 375. )  continue;
      counters[i].fill("MT2 > 375 GeV",weight);
      FillMonitor(ccount, fSamples[i].sname, fSamples[i].type, "MT2 > 375 GeV", weight);
      if( fMT2tree->misc.MT2 < 400. )  continue;
      counters[i].fill("MT2 > 400 GeV",weight);
      FillMonitor(ccount, fSamples[i].sname, fSamples[i].type, "MT2 > 400 GeV", weight);
      if( fMT2tree->misc.MT2 < 425. )  continue;
      counters[i].fill("MT2 > 425 GeV",weight);
      FillMonitor(ccount, fSamples[i].sname, fSamples[i].type, "MT2 > 425 GeV", weight);
      if( fMT2tree->misc.MT2 < 450. )  continue;
      counters[i].fill("MT2 > 450 GeV",weight);
      FillMonitor(ccount, fSamples[i].sname, fSamples[i].type, "MT2 > 450 GeV", weight);
      if( fMT2tree->misc.MT2 < 475. )  continue;
      counters[i].fill("MT2 > 475 GeV",weight);
      FillMonitor(ccount, fSamples[i].sname, fSamples[i].type, "MT2 > 475 GeV", weight);
      if( fMT2tree->misc.MT2 < 500. )  continue;
      counters[i].fill("MT2 > 500 GeV",weight);
      FillMonitor(ccount, fSamples[i].sname, fSamples[i].type, "MT2 > 500 GeV", weight);
      
    }
    delete fMT2tree;
  }

  std::cout << setfill('=') << std::setw(70) << "" << std::endl;
  std::cout << "Statistics - by sample" << std::endl;
  std::cout << setfill('_') << std::setw(70) << "" << setfill(' ') << std::endl;
  for ( size_t i = 0; i < fSamples.size(); ++i){
    std::cout << "++++  " << fSamples[i].name << std::endl;  
    std::cout << setfill('-') << std::setw(70) << "" << setfill(' ') << std::endl;    
    counters[i].print();
    std::cout << setfill('_') << std::setw(70) << "" << setfill(' ') << std::endl;    
  }  

  std::cout << setfill('=') << std::setw(70) << "" << std::endl;
  std::cout << "Statistics - by process" << std::endl;
  std::cout << setfill('_') << std::setw(70) << "" << setfill(' ') << std::endl;
  for ( size_t i = 0; i < nProc; ++i){
    std::cout << "++++  " << cnames[i] << std::endl;  
    std::cout << setfill('-') << std::setw(70) << "" << setfill(' ') << std::endl;    
    ccount[i].print();
    std::cout << setfill('_') << std::setw(70) << "" << setfill(' ') << std::endl;    
  }  

  std::cout << setfill('=') << std::setw(70) << "" << std::endl;
  std::cout << "Statistics - by process (MT2 > 100GeV)" << std::endl;
  std::cout << setfill('_') << std::setw(70) << "" << setfill(' ') << std::endl;
  for ( size_t i = 0; i < nProc; ++i){
    std::cout << "++++  " << cnames[i] << std::endl;  
    std::cout << setfill('-') << std::setw(70) << "" << setfill(' ') << std::endl;    
    ccount_100[i].print();
    std::cout << setfill('_') << std::setw(70) << "" << setfill(' ') << std::endl;    
  }  

}

// MT2 vs HT cutflot  --------------------------------------------------------------------
void MassPlotter::PrintCutFlowMT2vsHT(TString trigger, TString cuts){
  
  // define here the HT & MT2 sampling to be used
  int   numHTbins   =  15;
  float minHT       = 400;
  float HTgridsize  =  50;
  float minMT2cut   = 100;
  float maxMT2cut   = 700;
  float MT2gridsize =  25;
  
  const int   nProc      = 17;  
  TString  cnames[nProc] = {"QCD", "W+jets", "Z+jets", "Top", "Other", "Total Bkg.", "data", "LM1", "LM2", "LM3", "LM4", "LM5", "LM8", "LM9", "LM11", "LM12", "LM13"};
  Monitor Monitors[numHTbins][nProc];
  Monitor counters[fSamples.size()];

  // Loop over samples
  for(size_t i = 0; i < fSamples.size(); ++i){

    Double_t sample_weight =0;
    if(fPUReweight) sample_weight= fSamples[i].xsection * fSamples[i].kfact * fSamples[i].lumi / (fSamples[i].nevents * fSamples[i].PU_avg_weight );
    else            sample_weight= fSamples[i].xsection * fSamples[i].kfact * fSamples[i].lumi / (fSamples[i].nevents );
    if(fVerbose>2) cout << "PrintCutFlow: looping over " << fSamples[i].name << endl;
    if(fVerbose>2) cout << "              sample has weight " << sample_weight << " and " << fSamples[i].tree->GetEntries() << " entries" << endl; 
    if(fVerbose>2) cout << "              Average PU weight: " << fSamples[i].PU_avg_weight << endl;
    if(fVerbose>2) cout << "              PU reweightung:    " << fPUReweight << endl;
    if(fVerbose>2) cout << "              Original Entries: " << fSamples[i].nevents << endl;

    
    fMT2tree = new MT2tree();
    fSamples[i].tree->SetBranchAddress("MT2tree", &fMT2tree);
    Long64_t nentries =  fSamples[i].tree->GetEntries();
    Long64_t nbytes = 0, nb = 0;
    int nev =0;

    // don't do anything for empty trees: (maybe tree already skimmed)
    if (fSamples[i].tree->GetEntries()==0) continue;

    //leo tweak - filtering out the TTree
    TString myCuts = cuts;

    if( fSamples[i].type=="data" && trigger!="") myCuts += " && " + trigger; //cuts to be aplied only on data
    cout << "              Cuts for Flow: " << myCuts << endl;
    fSamples[i].tree->Draw(">>selList", myCuts);


    TEventList *myEvtList = (TEventList*)gDirectory->Get("selList");
    fSamples[i].tree->SetEventList(myEvtList);
    int counter=0;
    cout << "              Filtering done, size=" <<myEvtList->GetN()  << endl;
    
    if(myEvtList->GetSize()==0) continue;
    
    // loop over entried in MT2trees
    while(myEvtList->GetEntry(counter++) !=-1){
      
      int jentry = myEvtList->GetEntry(counter-1);
      
      nb =  fSamples[i].tree->GetEntry(jentry);   nbytes += nb;
      fSamples[i].tree->SetBranchAddress("MT2tree", &fMT2tree);
      
      if ( fVerbose>2 && jentry % 50000 == 0 )  cout << "+++ Proccessing event " << jentry << endl;
      
      Double_t weight = sample_weight;
      if (!fMT2tree->misc.isData ) weight = weight * fMT2tree->pileUp.Weight; // pile-up reweighting for MC 

	counters[i].fill("All events",weight); // fill all events that pass selection, without MT2 cut
	float MT2cut =minMT2cut; // initial cut on MT2
	while(MT2cut <=maxMT2cut ){
		if(fMT2tree->misc.MT2 >MT2cut){
			TString cutname = TString::Format("MT2 > %.0f GeV",MT2cut);
			counters[i].fill(cutname.Data(),weight);
			for(int HTbin=0; HTbin< numHTbins; ++HTbin){
				float HTcut = minHT + HTbin*HTgridsize;
				if( fMT2tree->misc.HT > HTcut)   FillMonitor(Monitors[HTbin], fSamples[i].sname, fSamples[i].type,cutname.Data() , weight);
			
			}
		}
		MT2cut+=MT2gridsize;
	}
    }
    delete fMT2tree;
  }

  std::ostringstream logStream;
  // by sample, no HT cuts
  logStream << setfill('=') << std::setw(70) << "" << std::endl;
  logStream << "Statistics - by sample, no HT cut" << std::endl;
  logStream << setfill('_') << std::setw(70) << "" << setfill(' ') << std::endl;
  for ( size_t i = 0; i < fSamples.size(); ++i){
    logStream << "++++  " << fSamples[i].name << std::endl;  
    logStream << setfill('-') << std::setw(70) << "" << setfill(' ') << std::endl;    
    logStream << counters[i];
    logStream << setfill('_') << std::setw(70) << "" << setfill(' ') << std::endl;    
  }  
  ofstream f_log ("Cutflow_bySample_defaultcuts.log", ios::app);
  f_log << logStream.str();
  logStream.str(""); // clear the ostringstream

  for ( int HTbin =0; HTbin < numHTbins; ++HTbin){
	  logStream << setfill('=') << std::setw(70) << "" << std::endl;
	  TString HTcut= TString::Format("%.0f", minHT+HTbin*HTgridsize);
	  logStream << setfill('_') << std::setw(70) << "" << setfill(' ') << std::endl;
	  logStream << "the following cuts were applied: " << std::endl;
	  logStream << cuts                                << std::endl;
 	  logStream << "with pfHT > " << HTcut	<< " GeV " << std::endl;  
	  logStream << setfill('_') << std::setw(70) << "" << setfill(' ') << std::endl;
	  logStream << "Statistics - by process (HT > " << HTcut << " GeV)" << std::endl;
	  logStream << setfill('_') << std::setw(70) << "" << setfill(' ') << std::endl;
	  for ( size_t i = 0; i < nProc; ++i){
	    logStream << "++++  " << cnames[i] << std::endl;  
	    logStream << setfill('-') << std::setw(70) << "" << setfill(' ') << std::endl;    
	    logStream << Monitors[HTbin][i];
	    logStream << setfill('_') << std::setw(70) << "" << setfill(' ') << std::endl;    
	  }  

	TString cutflowlog = "Cutflow_HT" +HTcut + "_byProcess.log";
	ofstream f_log2 (cutflowlog.Data(), ios::app);
	f_log2 << logStream.str();
	logStream.str(""); // clear the ostringstream
  }

}

//________________________________________________________________________

void MassPlotter::FillMonitor(Monitor *count, TString sname, TString type, TString cut, double weight){
  if     (sname == "QCD"    ) 	count[ 0].fill(cut.Data(),weight);
  else if(sname == "Wtolnu" ) 	count[ 1].fill(cut.Data(),weight);
  else if(sname == "DY"     ) 	count[ 2].fill(cut.Data(),weight);
  else if(sname == "Top"    )	count[ 3].fill(cut.Data(),weight);
  else if(sname == "VV" || sname == "VGamma" || sname == "Photons") 	
                                count[ 4].fill(cut.Data(),weight);
  if     (type  == "mc"     )   count[ 5].fill(cut.Data(),weight);
  else if(type  == "data"   )	count[ 6].fill(cut.Data(),weight);
  else if(sname == "LM1"    )	count[ 7].fill(cut.Data(),weight);
  else if(sname == "LM2"    )	count[ 8].fill(cut.Data(),weight);
  else if(sname == "LM3"    )	count[ 9].fill(cut.Data(),weight);
  else if(sname == "LM4"    )	count[10].fill(cut.Data(),weight);
  else if(sname == "LM5"    )	count[11].fill(cut.Data(),weight);
  else if(sname == "LM8"    )	count[12].fill(cut.Data(),weight);
  else if(sname == "LM9"    )	count[13].fill(cut.Data(),weight);
  else if(sname == "LM11"    )	count[14].fill(cut.Data(),weight);
  else if(sname == "LM12"    )	count[15].fill(cut.Data(),weight);
  else if(sname == "LM13"    )	count[16].fill(cut.Data(),weight);
}

//________________________________________________________________________

void MassPlotter::plotSig(TString var, TString cuts, TString xtitle, int nbins, double min, double max, bool flip_order, int type ){

        TString varname = Util::removeFunnyChar(var.Data());

	TH1D*    h_mc_sum    = new TH1D   (varname+"mc_sum", "", nbins, min, max );
	TH1D*    h_susy      = new TH1D   (varname+"susy"  , "", nbins, min, max );	
	h_mc_sum    ->Sumw2();
	h_susy      ->Sumw2();
	// vector of all histos
	vector<TH1D*> h_samples;

	for(size_t i = 0; i < fSamples.size(); ++i){
	        if(!(fSamples[i].type=="susy") && !(fSamples[i].type=="mc")) continue;

		h_samples.push_back(new TH1D(varname+"_"+fSamples[i].name, "", nbins, min, max));
		h_samples[i] -> Sumw2();

		Double_t weight = fSamples[i].xsection * fSamples[i].kfact * fSamples[i].lumi / (fSamples[i].nevents);

		if(fVerbose>2) cout << "MakePlot: looping over " << fSamples[i].sname << endl;
		if(fVerbose>2) cout << "           sample has weight " << weight << " and " << fSamples[i].tree->GetEntries() << " entries" << endl; 
	
		TString variable  = TString::Format("%s>>%s",var.Data(),h_samples[i]->GetName());
		TString selection = TString::Format("(%f) * (%s)",weight,cuts.Data());
		  
		if(fVerbose>2) cout << "+++++ Drawing " << variable  << endl
				    << "\twith cuts: "  << selection << endl;

		int nev = fSamples[i].tree->Draw(variable.Data(),selection.Data(),"goff");
		
		if(fVerbose>2) cout << "\tevents found : "  <<  nev << endl
				    << "\t->Integral() : "  <<  h_samples[i]->Integral() << endl;
		
		if( fSamples[i].type=="mc"        )   h_mc_sum->Add(h_samples[i]);
		else if( fSamples[i].type=="susy" )   h_susy  ->Add(h_samples[i]);
	}
	Float_t  x[nbins], y[nbins];
	for (int i = 1; i <=nbins+1; i++){
	  x[i-1] = h_mc_sum->GetBinLowEdge(i);
	  float s = h_susy  ->Integral(i,nbins+1);
	  float b = h_mc_sum->Integral(i,nbins+1);
	  switch (type){
	  case 0:
	    y[i-1] = s/sqrt(b);
	    break;
	  case 1:
	    y[i-1] = s/sqrt(s+b);
	    break;
	  case 2:
	    y[i-1] = s/b;
	    break;
	  }
	}
	TGraph *sig = new TGraph(nbins+1,x,y);
	sig->SetTitle("");
	sig->GetXaxis()->SetTitle(xtitle);
	sig->SetMarkerStyle(20);
	switch (type){
	case 0:
	  sig->GetYaxis()->SetTitle("S/#sqrt{B}");
	  break;
	case 1:
	  sig->GetYaxis()->SetTitle("S/#sqrt{S+B}");
	  break;
	case 2:
	  sig->GetYaxis()->SetTitle("S/B");
	  break;
	}
	sig->Draw("ACP");
}
//__________________________________________________________________________
void MassPlotter::CompSamples(TString var, TString cuts, TString optcut, bool RemoveLepts, 
		              TString xtitle, const int nbins, const double *bins, bool add_underflow, bool logflag, double scale_factor, bool normalize){

  	CompSamples(fSamples, var, cuts, optcut, RemoveLepts, xtitle, nbins, bins, add_underflow, logflag, scale_factor, normalize);
}
//___________________________________________________________________________________________________________________
void MassPlotter::compSamples(TString var, TString cuts, TString optcut, bool RemoveLepts,
			   TString xtitle, const int nbins, const double min, const double max, bool add_underflow, bool logflag, double scale_factor, bool normalize){

  	CompSamples(fSamples, var, cuts, optcut, RemoveLepts, xtitle, nbins, min, max, add_underflow, logflag, scale_factor, normalize);

}
// _______________________________________________________________________
void MassPlotter::CompSamples(std::vector<sample> Samples, TString var, TString cuts, TString optcut, bool RemoveLepts,
			   TString xtitle, const int nbins, const double min, const double max, bool add_underflow, bool logflag, double scale_factor, bool normalize){

  	double bins[nbins];
  	bins[0] = min;
  	for(int i=1; i<=nbins; i++)
    	bins[i] = min+i*(max-min)/nbins;
  	CompSamples(Samples, var, cuts, optcut, RemoveLepts, xtitle, nbins, bins, add_underflow, logflag, scale_factor, normalize);

}
// __________________________________________________________________________
void MassPlotter::CompSamples(std::vector<sample> Samples, TString var, TString cuts, TString optcut, bool RemoveLepts,
			   TString xtitle, const int nbins, const double *bins, bool add_underflow, bool logflag, double scale_factor, bool normalize){

        TString varname = Util::removeFunnyChar(var.Data());
	TH1D*    histos[2];
	cout << Samples.size() << endl;
	histos[0] = new TH1D(varname+Samples[0].sname  , "", nbins, bins);
	histos[1] = new TH1D(varname+Samples[1].sname  , "", nbins, bins);

	histos[0] -> Sumw2(); 
	histos[1] -> Sumw2();
	histos[0] -> SetLineColor(kBlue); 
	histos[1] -> SetLineColor(kRed); 
	histos[0] -> SetLineWidth(4);
	histos[1] -> SetLineWidth(4);

	// legend
	TLegend* Legend1 = new TLegend(.71,.54,.91,.92);

	for(int i=0; i<2; ++i){
		Double_t weight = scale_factor * Samples[i].xsection * Samples[i].kfact * Samples[i].lumi / (Samples[i].nevents);
		if(fVerbose>2) cout << "GetHisto: looping over " << Samples[i].sname << endl;
		if(fVerbose>2) cout << "           sample has weight " << weight << " and " << Samples[i].tree->GetEntries() << " entries" << endl; 
		// exchange variable to plot with the corresponding one with removed leptons;	
		TString variable; 
		TString selection; 
		if(RemoveLepts && (Samples[i].sname == "DYToLL" || Samples[i].sname=="QCD_2")){
			MapType::iterator iter = RemoveLeptMap.begin();
			iter = RemoveLeptMap.find(var);
			if (iter != RemoveLeptMap.end() ){
				variable = iter -> second;
			} else {cout << "found RemoveLepts==true, but no corresponding map" << endl; variable = var;}
		}else {
			variable = var;
		}
		if(optcut!="_" && Samples[i].sname == "DYToLL"){
			TString newcuts = cuts + "&&" + optcut;  
			if(Samples[i].type!="data" && fPUReweight==true ) selection      = TString::Format("(%.15f*pileUp.Weight) * (%s)",weight,newcuts.Data());
			else                                              selection      = TString::Format("(%.15f) * (%s)"              ,weight,newcuts.Data()); 
		} else {
			if(Samples[i].type!="data" && fPUReweight==true ) selection      = TString::Format("(%.15f*pileUp.Weight) * (%s)",weight,cuts.Data());
			else                                              selection      = TString::Format("(%.15f) * (%s)"              ,weight,cuts.Data()); 
		}

		variable          = TString::Format("%s>>%s",variable.Data(),histos[i]->GetName());
		
			  
		if(fVerbose>2) cout << "+++++ Drawing " << variable  << endl
				    << "\twith cuts: "  << selection << endl;

		int nev = Samples[i].tree->Draw(variable.Data(),selection.Data(),"goff");

		// Add underflow & overflow bins
		// This failed for older ROOT version when the first(last) bin is empty
		// and there are underflow (overflow) events --- must check whether this 
		// is still the case
		if(add_underflow) {
			histos[i]->SetBinContent(1,
				    histos[i]->GetBinContent(0) + histos[i]->GetBinContent(1));
			histos[i]->SetBinError(1,
				  sqrt(histos[i]->GetBinError(0)*histos[i]->GetBinError(0)+
				       histos[i]->GetBinError(1)*histos[i]->GetBinError(1) ));
		}
		histos[i]->SetBinContent(histos[i]->GetNbinsX(),
					    histos[i]->GetBinContent(histos[i]->GetNbinsX()  )+ 
					    histos[i]->GetBinContent(histos[i]->GetNbinsX()+1) );
		histos[i]->SetBinError(histos[i]->GetNbinsX(),
					  sqrt(histos[i]->GetBinError(histos[i]->GetNbinsX()  )*
					       histos[i]->GetBinError(histos[i]->GetNbinsX()  )+
					       histos[i]->GetBinError(histos[i]->GetNbinsX()+1)*
					       histos[i]->GetBinError(histos[i]->GetNbinsX()+1)  ));
		
		if(fVerbose>2) cout << "\tevents found : "  <<  nev << endl
				    << "\t->Integral() : "  <<  histos[i]->Integral() << endl;
	     	Legend1 ->AddEntry(histos[i], Samples[i].sname, "l"); 
	}
	
	plotRatio(histos[0], histos[1], logflag, normalize, varname+"_"+Samples[0].sname+"_"+Samples[1].sname, Legend1 , xtitle, "");
	
	
}

//________________________________________________________________________

void MassPlotter::MakePlot(std::vector<sample> Samples, TString var, TString maincuts, TString basecuts, int njets, int nbjets, int nleps, TString HLT,
			   TString xtitle, const int nbins, const double min, const double max,
			   bool flip_order, bool logflag, bool composited, bool ratio, 
			   bool stacked, bool overlaySUSY, float overlayScale, bool add_underflow, bool saveHistos, TString saveMacro){

  double bins[nbins];
  bins[0] = min;
  for(int i=1; i<=nbins; i++)
    bins[i] = min+i*(max-min)/nbins;
  MakePlot(Samples, var, maincuts, basecuts, njets, nbjets, nleps, HLT, xtitle, nbins, bins, flip_order, logflag, composited, ratio, stacked, overlaySUSY, overlayScale, add_underflow, saveHistos, saveMacro);

}

void MassPlotter::MakePlot(std::vector<sample> Samples, TString var, TString cuts, int njets, int nbjets, int nleps, TString HLT,
			   TString xtitle, const int nbins, const double min, const double max,
			   bool flip_order, bool logflag, bool composited, bool ratio, 
			   bool stacked, bool overlaySUSY, float overlayScale, bool add_underflow, bool saveHistos, TString saveMacro){

  double bins[nbins];
  bins[0] = min;
  for(int i=1; i<=nbins; i++)
    bins[i] = min+i*(max-min)/nbins;
  MakePlot(Samples, var, cuts, njets, nbjets, nleps, HLT, xtitle, nbins, bins, flip_order, logflag, composited, ratio, stacked, overlaySUSY, overlayScale, add_underflow, saveHistos, 	saveMacro);

}

//________________________________________________________________________

void MassPlotter::MakePlot(std::vector<sample> Samples, TString var, TString cuts, int njets, int nbjets, int nleps, TString HLT,
			   TString xtitle, const int nbins, const double *bins, 
			   bool flip_order, bool logflag, bool composited, bool ratio, 
			   bool stacked, bool overlaySUSY, float overlayScale, bool add_underflow, bool saveHistos, TString saveMacro){

  TString basecuts = "";
  TString maincuts = cuts;
  MakePlot(Samples, var, maincuts, basecuts, njets, nbjets, nleps, HLT, xtitle, nbins, bins, flip_order, logflag, composited, ratio, stacked, overlaySUSY, overlayScale, add_underflow, saveHistos, saveMacro);

}

void MassPlotter::MakePlot(std::vector<sample> Samples, TString var, TString maincuts, TString basecuts, int njets, int nbjets, int nleps, TString HLT,
			   TString xtitle, const int nbins, const double *bins, 
			   bool flip_order, bool logflag, bool composited, bool ratio, 
			   bool stacked, bool overlaySUSY, float overlayScale, bool add_underflow, bool saveHistos, TString saveMacro){
        TString varname = Util::removeFunnyChar(var.Data());

	TString nMinCut = var;
	nMinCut += TString::Format(">=%g", bins[0]);

	TString nJets, nJetsVar = "NJetsIDLoose40";
	if (njets>=10) {
	  nJets =  "(" + nJetsVar + TString::Format(">=%d",njets/10);
	  nJets += "&&"+ nJetsVar + TString::Format("<=%d",njets%10)+")";
	}
	else{
	  nJets = nJetsVar + (njets < 0 ? ">=" : "==");
	  nJets = nJets + TString::Format("%d",abs(njets));
	}

	if(njets == -10)
	  nJets = nJetsVar + ">= 0";


	//saeid
// 	TString nBJets = "NBJets40CSVM";    // nbjets = -10  --> >=0 b-tags
	TString nBJets = "NBJetsCSVT";    // nbjets = -10  --> >=0 b-tags
	//saeid
	nBJets += nbjets < 0 ? ">=" : "==";
	nBJets += nbjets==-10 ? "0" : TString::Format("%d",abs(nbjets));

	TString  nLeps;
	if     (nleps < 0 )  nLeps = " && (NEles + NMuons) >=";
	else if(nleps >=0  ) nLeps = " && (NEles + NMuons) ==";
	nLeps += TString::Format("%d",abs(nleps));
	if     (nleps ==-10) nLeps = " "; 
	if     (nleps ==-11) nLeps = " && NEles ==1 && NMuons ==0"; 
	if     (nleps ==-13) nLeps = " && NEles ==0 && NMuons ==1"; 

	TString tmpname = varname+"_CUT_"+nJets+"_"+nBJets+"_"+nLeps;
	tmpname.ReplaceAll(">=" ,".ge");
	tmpname.ReplaceAll("<=" ,".le");
	tmpname.ReplaceAll(">" ,".gt");
	tmpname.ReplaceAll("<" ,".lt");
	tmpname.ReplaceAll("==",".eq");
	tmpname.ReplaceAll("!=",".ne");
	tmpname.ReplaceAll("&&","_");
	tmpname.ReplaceAll("||","_");
	tmpname.ReplaceAll("misc.","");
	TString varname2 = Util::removeFunnyChar(tmpname.Data());
	
	THStack* h_stack     = new THStack(varname2, "");
  	TH1D*    h_data      = new TH1D   (varname2+"data"  , "", nbins, bins );
	TH1D*    h_mc_sum    = new TH1D   (varname2+"mc_sum", "", nbins, bins );
	TH1D*    h_susy      = new TH1D   (varname2+"susy"  , "", nbins, bins );	

	// h_data
	h_data -> Sumw2();
	h_data -> SetMarkerStyle(20);
	h_data -> SetMarkerColor(kBlack);
	h_data -> SetLineColor(kBlack);
	h_data -> SetStats(false);

	h_mc_sum -> SetFillStyle(3004);
	h_mc_sum -> SetFillColor(kBlack);
	h_mc_sum -> SetStats(0);
	
	h_mc_sum    ->Sumw2();
	h_susy      ->Sumw2();


	// vector of all histos
	vector<TH1D*> h_samples;

	// arrays of composited stuff
	TH1D    *h_composited[7];
	//TString  cnames[5] = {"QCD", "W/Z/#gamma production", "Top production", "susy", "data"};
	//int      ccolor[5] = {401, 418, 602, 0, 632};
	TString  cnames[7] = {"QCD", "W+jets", "Z+jets", "Top", "Other", "susy", "data"};
	int      ccolor[7] = { 401,   417,       419,      855,    603,     0,     632};
	vector<TH1D*> h_signals;
	int n_sig=0;
	TString  snames[5] = {"signal1","signal2","signal3","signal4","signal5"};
	int      scolor[5] = {     kRed+3,   kOrange-3,   kViolet-6,   kBlack,   kRed-4};
	for (int i=0; i<7; i++){
	  h_composited[i] = new TH1D(varname2+"_"+cnames[i], "", nbins, bins);
	  h_composited[i] -> Sumw2();
	  h_composited[i] -> SetFillColor  (stacked ? ccolor[i] : 0);
	  h_composited[i] -> SetLineColor  (ccolor[i]);
	  if (!stacked) {
	    h_composited[i] -> SetLineWidth(4);
	  }
	  h_composited[i] -> SetMarkerColor(ccolor[i]);
	  h_composited[i] -> SetStats(false);
	}
	if(fisPhoton){
		h_composited[4] -> SetFillColor(kViolet-3); 
		h_composited[4] -> SetLineColor(kViolet-3);
		cnames[4]="Photons";
	}

	// legend
	TLegend* Legend1 = new TLegend(.71,.68,.91,.92);
	Legend1->SetName(varname2+"_legend");
	//TLegend* Legend1 = new TLegend(.3,.5,.6,.88);

	for(size_t i = 0; i < Samples.size(); ++i){
		h_samples.push_back(new TH1D(varname2+"_"+Samples[i].name, "", nbins, bins));
		h_samples[i] -> Sumw2();
		h_samples[i] -> SetFillColor(stacked ? Samples[i].color : 0);
		h_samples[i] -> SetLineColor(Samples[i].color);
		if (!stacked) {
		  h_samples[i] -> SetLineWidth(4);
		}
		h_samples[i] -> SetMarkerColor(Samples[i].color);
		h_samples[i] -> SetStats(false);
		if(Samples[i].type == "susy" ){
			h_samples[i] -> SetLineColor(kBlack);
			h_samples[i] -> SetLineStyle(kDotted);
			h_composited[5] -> SetLineColor(kBlack);
			h_composited[5] -> SetLineStyle(kDotted);
		}
	
		Double_t weight=0;
		if(fPUReweight) weight = Samples[i].xsection * Samples[i].kfact * Samples[i].lumi / (Samples[i].nevents*Samples[i].PU_avg_weight);
		else            weight = Samples[i].xsection * Samples[i].kfact * Samples[i].lumi / (Samples[i].nevents);

		if(fVerbose>2) cout << "===========================================================" <<  endl;
		if(fVerbose>2) cout << "MakePlot: looping over " << Samples[i].sname << "-----------------------------------" <<  endl;
		if(fVerbose>2) cout << "  +++++++ xsection:    "    << Samples[i].xsection << " k-fact " << Samples[i].kfact << endl;
		if(fVerbose>2) cout << "  +++++++ tot events:  "  << Samples[i].nevents  << " avg pu weight " << Samples[i].PU_avg_weight << endl;
		if(fVerbose>2) cout << "  +++++++ PU reweight: "  << fPUReweight << endl;
		if(fVerbose>2) cout << "  +++++++ weight:      " << weight << " and " << Samples[i].tree->GetEntries() << " entries" << endl; 
	
		TString variable  = TString::Format("%s>>%s",var.Data(),h_samples[i]->GetName());
		TString theCuts = nJets + "&&" + nBJets + nLeps;
		if(!add_underflow) theCuts = theCuts + "&&" + nMinCut;
		theCuts = theCuts + "&&" + maincuts;
		if(basecuts!="") theCuts = theCuts + "&&" + basecuts;
		if(Samples[i].type=="data" && HLT!="") theCuts += " &&("+HLT+")"; // triggers for data

		TString btagweight = "1.00"; //stored btag weights up to >=3, default is weight==1 to avoid non-existing weights
		if(nbjets>=0 && nbjets<=3) btagweight = TString::Format("SFWeight.BTagCSV40eq%d",abs(nbjets));
		else if(nbjets>=-3)        btagweight = TString::Format("SFWeight.BTagCSV40ge%d",abs(nbjets));

		TString selection;
		if(     Samples[i].type!="data" && fPUReweight && fbSFReWeight) selection = TString::Format("(%.15f*pileUp.Weight*%s) * (%s)",weight, btagweight.Data(), theCuts.Data());
		else if(Samples[i].type!="data" && fPUReweight                ) selection = TString::Format("(%.15f*pileUp.Weight) * (%s)",   weight,                    theCuts.Data());
		else if(Samples[i].type!="data" &&                fbSFReWeight) selection = TString::Format("(%.15f*%s) * (%s)",              weight, btagweight.Data(), theCuts.Data());
		else                                                            selection = TString::Format("(%.15f) * (%s)",                 weight,                    theCuts.Data()); 

		if(     Samples[i].type=="susy")  selection = TString::Format("(%.15f) * (%s)",weight, theCuts.Data());

	//	TString selection;
	//	if(Samples[i].type!="data" && fPUReweight) selection      = TString::Format("(%.15f*pileUp.Weight) * (%s)",weight,theCuts.Data());
	//	else                                       selection      = TString::Format("(%.15f) * (%s)"              ,weight,theCuts.Data()); 
		  
		if(fVerbose>2) cout << "  +++++++ Drawing      " << variable  << endl
				    << "  +++++++ with cuts:   " << setw(40)  << selection << endl;

		int nev = Samples[i].tree->Draw(variable.Data(),selection.Data(),"goff");

		// Add underflow & overflow bins
		// This failed for older ROOT version when the first(last) bin is empty
		// and there are underflow (overflow) events --- must check whether this 
		// is still the case
		if(add_underflow) {
		  h_samples[i]->SetBinContent(1,
					      h_samples[i]->GetBinContent(0) + h_samples[i]->GetBinContent(1));
		  h_samples[i]->SetBinError(1,
					    sqrt(h_samples[i]->GetBinError(0)*h_samples[i]->GetBinError(0)+
						 h_samples[i]->GetBinError(1)*h_samples[i]->GetBinError(1) ));
		}
		h_samples[i]->SetBinContent(h_samples[i]->GetNbinsX(),
					    h_samples[i]->GetBinContent(h_samples[i]->GetNbinsX()  )+ 
					    h_samples[i]->GetBinContent(h_samples[i]->GetNbinsX()+1) );
		h_samples[i]->SetBinError(h_samples[i]->GetNbinsX(),
					  sqrt(h_samples[i]->GetBinError(h_samples[i]->GetNbinsX()  )*
					       h_samples[i]->GetBinError(h_samples[i]->GetNbinsX()  )+
					       h_samples[i]->GetBinError(h_samples[i]->GetNbinsX()+1)*
					       h_samples[i]->GetBinError(h_samples[i]->GetNbinsX()+1)  ));
		
		if(fVerbose>2) cout << "  +++++++ MC   events : "  <<  nev << endl;

		/// event count with errors
		TH1F * clone = (TH1F*)h_samples[i]->Clone();
		clone->Rebin(clone->GetNbinsX());
		if(fVerbose>2) cout << "  +++++++ Events:       " << clone->GetBinContent(1) << " +- " << clone->GetBinError(1) << endl;
		
		if (Samples[i].sname.Contains("QCD")) {
		  h_composited[0]->Add(h_samples[i]);
		}
		else if(Samples[i].sname.Contains("Top") ){
		  h_composited[3]->Add(h_samples[i]);
		}
		else if(Samples[i].sname=="Wtolnu") {
		  h_composited[1]->Add(h_samples[i]);
		}
		else if(Samples[i].sname=="DY") {
		  h_composited[2]->Add(h_samples[i]);
		}
		else if(Samples[i].sname=="VV" || Samples[i].sname=="VGamma" || Samples[i].sname=="Photons") {
		  h_composited[4]->Add(h_samples[i]);
		}
		else if(Samples[i].type.Contains("susy")){
		  h_composited[5]->Add(h_samples[i]);
		  cnames[5] = Samples[i].sname;
		  snames[n_sig] = Samples[i].sname;
		  char buf[10];
		  sprintf(buf,"no%d", n_sig);
		  //itoa(n_sig, buf, 10);
		  TString numstring = TString(buf);
		  h_signals.push_back(new TH1D(varname2+"_"+snames[n_sig]+"_"+numstring, "", nbins, bins));
		  h_signals[n_sig] -> Sumw2();
		  h_signals[n_sig] -> SetFillColor  (stacked ? scolor[n_sig] : 0);
		  h_signals[n_sig] -> SetLineColor  (scolor[n_sig]);
		  h_signals[n_sig] -> SetLineWidth  (3);
		  if (!stacked) {
		    h_signals[n_sig] -> SetLineWidth(4);
		  }
		  h_signals[n_sig] -> SetMarkerColor(scolor[n_sig]);
		  h_signals[n_sig] -> SetStats(false);

		  h_signals[n_sig]->Add(h_samples[i]);
		  n_sig++;
		}
		else if(Samples[i].type.Contains("data")){
		  h_composited[6]->Add(h_samples[i]);
		}

	}

	if (composited){
	  if(flip_order){
	    for (int i=5; i>=0; --i){
	      if (!stacked)  {
	        h_composited[i]->Scale(1./h_composited[i]->Integral());
	        h_composited[i] -> SetMinimum(0.00005);
	        //if (fVerbose>2) cout << " h_composited[" << i<< "]->Integral = " << h_composited[i]->Integral()  << endl;
	      }
	      if( i!=5 || !overlaySUSY) h_stack  -> Add(h_composited[i]);
	      if( i==5)                 h_susy   -> Add(h_composited[i]);
	      else        	        h_mc_sum -> Add(h_composited[i]);
	      if( i==5 && overlaySUSY)
	        //Legend1 ->AddEntry(h_composited[i], cnames[i]  , "f");  
		for (int iii=0; iii<n_sig; iii++) Legend1 ->AddEntry(h_signals[iii], snames[iii]  , "f");  
	      else
	        Legend1 ->AddEntry(h_composited[i], cnames[i], stacked ? "f" : "l");
	    }
	  }
	  else{
	    for (int i=0; i<6; ++i){
	      if (!stacked)  {
	        h_composited[i]->Scale(1./h_composited[i]->Integral());
	        h_composited[i] -> SetMinimum(0.00005);
	        //if (fVerbose>2) cout << " h_composited[" << i<< "]->Integral = " << h_composited[i]->Integral()  << endl;
	      }
	      if( i!=5 || !overlaySUSY) h_stack  -> Add(h_composited[i]);
	      if( i==5)                 h_susy   -> Add(h_composited[i]);
	      else        	      h_mc_sum -> Add(h_composited[i]);
	      if( i==5 && overlaySUSY){
	        //if(h_composited[i]->GetEntries()>0) Legend1 ->AddEntry(h_composited[i], cnames[i] + (overlayScale ? overlayScale==1? "":
		//  					      TString::Format(" x %.0f",overlayScale) :
		//  					      " scaled to data")   , "f");  
		for (int iii=0; iii<n_sig; iii++) Legend1 ->AddEntry(h_signals[iii], snames[iii]+  (overlayScale ? overlayScale==1? "":
		  					      TString::Format(" x %.0f",overlayScale) :
		  					      " scaled to data")   , "f");  
	      }else{
	        if(h_composited[i]->GetEntries()>0) Legend1 ->AddEntry(h_composited[i], cnames[i], stacked ? "f" : "l");
	      }
	    }
	  }
	  h_data->Add(h_composited[6]);
	  if(h_data->Integral()>0 && stacked){
	    Legend1     ->AddEntry(h_data, "data", "p");
	  }
	  if(fVerbose > 2 && composited) {
		           cout << "------------------------------------"                << endl
	                        << "QCD Integral:      " << h_composited[0]->Integral()  << endl
			        << "W+jets Integral:   " << h_composited[1]->Integral()  << endl
			        << "Z+jets Integral:   " << h_composited[2]->Integral()  << endl
			        << "Top Integral:      " << h_composited[3]->Integral()  << endl
			        << "Other Integral:    " << h_composited[4]->Integral()  << endl
				<< "TOTAL BG:          " << h_composited[0]->Integral()+h_composited[1]->Integral()+h_composited[2]->Integral()+h_composited[3]->Integral()+h_composited[4]->Integral() <<endl
			        << "SUSY:              " << h_composited[5]->Integral()  << endl
			        << "Data:              " << h_composited[6]->Integral()  << endl;

			   //Same, but with errors
			   string OverallSamples[8] = {"QCD","W+jets","Z+Jets","Top","Other","SUSY","Data","TOTAL BG"};
			   for(int os=0; os<8; os++){
			     string tSample = OverallSamples[os];
			     if(tSample!="TOTAL BG"){
			       TH1F* clone = (TH1F*) h_composited[os]->Clone();
			       clone->Rebin(clone->GetNbinsX());
			       if(fVerbose>2) cout << tSample << ": " << clone->GetBinContent(1) << " +- " << clone->GetBinError(1) << endl;
			     }
			     else{
			       TH1F* clone = (TH1F*) h_composited[0]->Clone();
			       clone->Add(h_composited[1]);
			       clone->Add(h_composited[2]);
                               clone->Add(h_composited[3]);
			       clone->Add(h_composited[4]);
			       clone->Rebin(clone->GetNbinsX());
                               if(fVerbose>2) cout << tSample << ": " << clone->GetBinContent(1) << " +- " << clone->GetBinError(1) << endl;
			     }
			   }
			  // save nevents for W background prediction
			  if(nleps==-11){
			  	fWpred.QCD_bg_e   = h_composited[0]->Integral();	
				fWpred.W_bg_e     = h_composited[1]->Integral();
				fWpred.Z_bg_e     = h_composited[2]->Integral();
				fWpred.Top_bg_e   = h_composited[3]->Integral();
				fWpred.Other_bg_e = h_composited[4]->Integral();
			  }else if(nleps==-13){
			  	fWpred.QCD_bg_mu   = h_composited[0]->Integral();	
				fWpred.W_bg_mu     = h_composited[1]->Integral();
				fWpred.Z_bg_mu     = h_composited[2]->Integral();
				fWpred.Top_bg_mu   = h_composited[3]->Integral();
				fWpred.Other_bg_mu = h_composited[4]->Integral();
			  } 
	  }
	}
	else{
	  for(int i=0; i<Samples.size(); ++i){
	    if (!stacked) {
	      h_samples[i]->Scale(1./h_samples[i]->Integral());
	      h_samples[i] -> SetMinimum(0.0001);
	    }
	    if((Samples[i].sname=="QCD" )){
	      h_samples[i] -> SetMinimum(0.01);
	      h_stack->Add(h_samples[i]);
	      h_mc_sum->Add(h_samples[i]);
	      Legend1 ->AddEntry(h_samples[i], Samples[i].sname, "f");
	    }
	    if(Samples[i].sname!="QCD" && Samples[i].type!="data"){
	      if(Samples[i].type=="susy") h_susy->Add(h_samples[i]);
	      if(Samples[i].type!="susy" || !overlaySUSY) h_stack->Add(h_samples[i]);
	      if(Samples[i].type!="susy") h_mc_sum->Add(h_samples[i]);
	      if(Samples[i].type=="susy" && overlaySUSY){
	         if(h_samples[i]->Integral()>0) Legend1 ->AddEntry(h_samples[i], Samples[i].sname + (overlayScale ? 
								       TString::Format(" x %.0f",overlayScale) :
								       " scaled to data")   , "f") ;
	      }else{
		 if(h_samples[i]->Integral()>0) Legend1 ->AddEntry(h_samples[i], Samples[i].sname, stacked ? "f" : "l");
	      }	
	    }
	    if(Samples[i].type == "data"){
	      h_data -> Add(h_samples[i]);
	    }
	    
	  }
	  if(h_data->Integral()>0){
	    Legend1     ->AddEntry(h_data, "data", "l");
	  }
	}

	TString ytitle = "Events";

	TString oname = var+ "_CUT_";
	oname += njets < 0 ? TString::Format("ge%dJets",abs(njets)) : TString::Format("%dJets",abs(njets));
	oname += nbjets == -10 ? "" : nbjets < 0 ? TString::Format("_ge%dBtag",abs(nbjets)) : TString::Format("_%dBtag",abs(nbjets));
	oname += nleps == 0 ? "_0Lep_" : nleps == -10 ? "_anyLep_" : nleps == -11 ? "_1Ele_" : nleps == -13 ? "_1Muo_" :
			     nleps < 0 ? TString::Format("_ge%dLeps_",abs(nleps)) : TString::Format("_%dLeps_",abs(nleps));
	oname += maincuts;
	oname.ReplaceAll(">=" ,".ge");
	oname.ReplaceAll("<=" ,".le");
	oname.ReplaceAll(">" ,".gt");
	oname.ReplaceAll("<" ,".lt");
	oname.ReplaceAll("==",".eq");
	oname.ReplaceAll("!=",".ne");
	oname.ReplaceAll("&&","_");
	oname.ReplaceAll("||","_");
	oname.ReplaceAll("misc.","");
	oname.ReplaceAll("LeadingJPt","J1Pt");
	oname.ReplaceAll("SecondJPt","J2Pt");
	oname.ReplaceAll("LeptConfig","LepCfg");
	oname.ReplaceAll("Vectorsumpt","VSPT");
	oname.ReplaceAll("hcalLaserEventFlag","hcalFlg");
	oname.ReplaceAll("trackingFailureFlag","trkFailFlg");
	oname.ReplaceAll("eeBadScFlag","eeBadFlg");
	oname.ReplaceAll("EcalDeadCellTriggerPrimitiveFlag","EcalTPFlg");
	oname.ReplaceAll("EcalDeadCellBEFlag","BEFlg");
	oname.ReplaceAll("HBHENoiseFlag","HBHEFlg");
	oname.ReplaceAll("CSCTightHaloIDFlag","CSCFlg");
	oname.ReplaceAll("NJetsIDLoose40","NJIDLoose40");
	oname.ReplaceAll("NJetsIDLoose","NJIDLoose");
	oname.ReplaceAll("isPFIDLoose","isJLoose");
	oname.ReplaceAll("IsGoodPFJet","IsGoodPFJ");
	oname.ReplaceAll("MinMetJetDPhi","MinDPhi");
	oname.ReplaceAll("trigger.HLT_HT260_MHT60_v2","HLT_HT260_MHT60_v2");
	oname.ReplaceAll("trigger.HLT_HT250_MHT60_v3","HLT_HT250_MHT60_v3");
	oname.ReplaceAll("trigger.HLT_HT250_MHT60_v2","HLT_HT250_MHT60_v2");
	oname.ReplaceAll("trigger.HLT_HT440_v2","HLT_HT440_v2");
	oname.ReplaceAll("trigger.HLT_HT450_v2","HLT_HT450_v2");
	oname.ReplaceAll("trigger.HLT_HT500_v3","HLT_HT500_v3");
	oname.ReplaceAll(",","-");
        TString outname = Util::removeFunnyChar(oname.Data());
	outname.ReplaceAll("_ProcessID.ne10_Susy.MassChi.eq350_Susy.MassLSP.eq50","");
	outname.ReplaceAll("_ProcessID.eq10","");
	outname.ReplaceAll("NMuons.eq0_muo0.lv.Pt.lt10_NEles.eq0_ele0.lv.Pt.lt10","noLepPt10");
	outname.ReplaceAll("NMuons.gt0_muo0.lv.Pt.gt10_NEles.gt0_ele0.lv.Pt.gt10","gt1lepPt10");
	outname.ReplaceAll("Jet0Pass.eq1","J0Pass");
	outname.ReplaceAll("Jet1Pass.eq1","J1Pass");
	outname.ReplaceAll("PassJetID.eq1","PassJID");
	outname.ReplaceAll("Flg.eq0","Flg");
	outname.ReplaceAll("_ProcessID.ne6_Event.ne1689009_Event.ne2275452_Event.ne1946785_Event.ne1936763_Event.ne1890738_Event.ne1757319_","_");
	outname.ReplaceAll("ProcessID.ne6_Event.ne4160010_Event.ne5022935_Event.ne1323244_Event.ne1305531_Event.ne4053304_Event.ne5630056_Event.ne3539458_ProcessID.ne10_Susy.MassChi.eq450_Susy.MassLSP.eq150","");
	outname.ReplaceAll("_ProcessID.ne6_Event.ne4160010_Event.ne5022935_Event.ne1323244_Event.ne1305531_Event.ne4053304_Event.ne5630056_Event.ne3539458","");
	outname.ReplaceAll("_HBHEFlg_CSCFlg_hcalFlg_trkFailFlg_eeBadFlg_EcalTPFlg_CrazyHCAL.eq0","");

	if(saveHistos){
		TString fileName = fOutputDir;
		if(!fileName.EndsWith("/")) fileName += "/";
		Util::MakeOutputDir(fileName);
		fileName = fileName + outname + ".root";
		TFile *savefile = new TFile(fileName.Data(), "RECREATE");
		savefile ->cd();
		for(unsigned int nh = 0; nh<7; ++nh)                h_composited[nh]->Write();
		for(unsigned int nh = 0; nh<h_signals.size(); ++nh) h_signals[nh]->Write();
		h_stack->Write();
		h_mc_sum->Write();
		h_data->Write();
		h_susy->Write();
		Legend1->Write();
		savefile->Close();
		std::cout << "Saved histograms in " << savefile->GetName() << std::endl;
	}

	if(!stacked) {
	  outname = outname + (flip_order ? "_flipped" : "") + (logflag ? "_log" : "") + "_shape";
	  printHisto(h_stack, h_data, h_mc_sum, Legend1 , outname, "hist", logflag, xtitle, ytitle, njets, nleps, false);

	}
	else if (!overlaySUSY){
	  outname = outname + (flip_order ? "_flipped" : "") + (logflag ? "_log" : "") + (composited ? "_comp" : "");
	  printHisto(h_stack, h_data, h_mc_sum, Legend1 , outname, "hist", logflag, xtitle, ytitle, njets, nleps);
	  if (ratio) 
	    plotRatioStack(h_stack,  h_mc_sum, h_data,  logflag, false, outname, Legend1, xtitle, ytitle, njets, nleps);
	}
	else {
	  outname =  outname + (flip_order ? "_flipped" : "") + (logflag ? "_log" : "") + (composited ? "_comp" : "") + "_overlay";	  
	  //printHisto(h_stack, h_data, h_mc_sum, h_susy, Legend1 , outname, "hist", logflag, xtitle, ytitle, njets, nbjets, nleps, overlayScale);
	  printHisto(h_stack, h_data, h_mc_sum, h_signals, n_sig, Legend1 , outname, "hist", logflag, xtitle, ytitle, njets, nbjets, nleps, overlayScale);
	  if (ratio) 
	    plotRatioStack(h_stack,  h_mc_sum, h_data, h_susy,  logflag, false, outname, Legend1, xtitle, ytitle, njets, nbjets, nleps, overlayScale, saveMacro);

	}
	

// 	for(int i=0; i<h_samples.size(); ++i){
// 		delete h_samples[i];
// 	}
// 	delete h_mc_sum;
// 	delete Legend1;
// 	delete h_data;
// 	delete h_susy;
	
}




//__________________________________________________________________________
void MassPlotter::abcd_MT2(TString var, TString basecut, TString upper_cut, TString lower_cut, const int nbins,const double min, const double max, double fit_min, double fit_max){

  double bins[nbins];
  bins[0] = min;
  for(int i=1; i<=nbins; i++)
    bins[i] = min+i*(max-min)/nbins;
  ABCD_MT2(var, basecut, upper_cut, lower_cut, nbins, bins, fit_min, fit_max);
  
}

//__________________________________________________________________________
void MassPlotter::ABCD_MT2(TString var, TString basecut, TString upper_cut, TString lower_cut, const int nbins, const double *bins, double fit_min, double fit_max){
  TH2D *h_ABCD_MT2_qcd           = new TH2D("ABCD_MT2__qcd"           , "",  100, 0, 600, 180, 0, TMath::Pi());  h_ABCD_MT2_qcd          ->Sumw2();
  TH2D *h_ABCD_MT2_susy          = new TH2D("ABCD_MT2__susy"          , "",  100, 0, 600, 180, 0, TMath::Pi());  h_ABCD_MT2_susy         ->Sumw2();
  TH1D *h_ABCD_lower_y_band_susy = new TH1D("ABCD_lower_y_band__susy" , "",  nbins, bins);		        h_ABCD_lower_y_band_susy->Sumw2();
  TH1D *h_ABCD_upper_y_band_susy = new TH1D("ABCD_upper_y_band__susy" , "",  nbins, bins);		        h_ABCD_upper_y_band_susy->Sumw2();
  TH1D *h_ABCD_lower_y_band_data = new TH1D("ABCD_lower_y_band__data" , "",  nbins, bins);		        h_ABCD_lower_y_band_data->Sumw2();
  TH1D *h_ABCD_upper_y_band_data = new TH1D("ABCD_upper_y_band__data" , "",  nbins, bins);		        h_ABCD_upper_y_band_data->Sumw2();
  TH1D *h_ABCD_lower_y_band_qcd  = new TH1D("ABCD_lower_y_band__qcd"  , "",  nbins, bins);		        h_ABCD_lower_y_band_qcd ->Sumw2();
  TH1D *h_ABCD_upper_y_band_qcd  = new TH1D("ABCD_upper_y_band__qcd"  , "",  nbins, bins);		        h_ABCD_upper_y_band_qcd ->Sumw2();
  TH1D *h_ABCD_lower_y_band_mc   = new TH1D("ABCD_lower_y_band__mc"   , "",  nbins, bins);		        h_ABCD_lower_y_band_mc  ->Sumw2();
  TH1D *h_ABCD_upper_y_band_mc   = new TH1D("ABCD_upper_y_band__mc"   , "",  nbins, bins);		        h_ABCD_upper_y_band_mc  ->Sumw2();
  TH1D *h_ABCD_lower_y_band_sub  = new TH1D("ABCD_lower_y_band__sub"  , "",  nbins, bins);		        h_ABCD_lower_y_band_sub ->Sumw2();
  TH1D *h_ABCD_upper_y_band_sub  = new TH1D("ABCD_upper_y_band__sub"  , "",  nbins, bins);		        h_ABCD_upper_y_band_sub ->Sumw2();
  TH1D *ratio_data               = new TH1D("ratio__data"             , "",  nbins, bins);		        ratio_data              ->Sumw2();
  TH1D *ratio_qcd                = new TH1D("ratio__qcd"              , "",  nbins, bins);		        ratio_qcd               ->Sumw2();
  TH1D *ratio_sub                = new TH1D("ratio__sub"              , "",  nbins, bins);		        ratio_sub               ->Sumw2();
  
  for(size_t i = 0; i < fSamples.size(); ++i){
    
    Long64_t nentries =  fSamples[i].tree->GetEntries();
    Long64_t nbytes = 0, nb = 0;
    
    Double_t weight = fSamples[i].xsection * fSamples[i].kfact * fSamples[i].lumi / (fSamples[i].nevents);
    if(fVerbose>2) cout << "ABCD: looping over " << fSamples[i].name << endl;
    if(fVerbose>2) cout << "      sample has weight " << weight << endl; 
    
//     TEntryList *elist = new TEntryList("elist","elist");
//     int nev = fSamples[i].tree->Draw(">>elist",basecut.Data(),"entrylist");
//     fSamples[i].tree->SetEntryList(elist);
    
    TString  weights   = (fSamples[i].type!="data" ?  TString::Format("(%.15f*pileUp.Weight)",weight) : TString::Format("(%.15f)",weight));
    TString selection = TString::Format("(%s) * (%s)"      ,weights.Data(),basecut.Data());
    TString sel_up    = TString::Format("(%s) * (%s && %s)",weights.Data(),basecut.Data(),upper_cut.Data());
    TString sel_lo    = TString::Format("(%s) * (%s && %s)",weights.Data(),basecut.Data(),lower_cut.Data());

    // Remove jets intentionaly!!!
    // RemoveAndRecalcMT2() has to be called before any selection to make the internal changes effective!
    TString removeJet = "1";
    //TString removeJet = "SmearAndRecalcMT2()>-5.";
    //TString removeJet = "RemoveAndRecalcMT2(0,0.000001)>-5.";
    //TString removeJet = "RemoveAndRecalcMT2(1,0.005)>-5.";
    TString selection_QCD = TString::Format("(%s) * (%s && %s)"      ,weights.Data(),removeJet.Data(),basecut.Data());
    TString sel_up_QCD    = TString::Format("(%s) * (%s && %s && %s)",weights.Data(),removeJet.Data(),basecut.Data(),upper_cut.Data());
    TString sel_lo_QCD    = TString::Format("(%s) * (%s && %s && %s)",weights.Data(),removeJet.Data(),basecut.Data(),lower_cut.Data());

    TString variable;
    TString var2 = "misc.MT2";
    if (fSamples[i].type == "data"){
      variable  = TString::Format("%s>>+%s",var2.Data(),h_ABCD_lower_y_band_data->GetName());
      int nevL = fSamples[i].tree->Draw(variable,sel_lo,"goff");
      variable  = TString::Format("%s>>+%s",var2.Data(),h_ABCD_upper_y_band_data->GetName());
      int nevU = fSamples[i].tree->Draw(variable,sel_up,"goff");
      if(fVerbose>2) cout << "\tData lower, events found : "  <<  nevL << endl
			  << "\t->Integral() : "  <<  h_ABCD_lower_y_band_data->Integral() + h_ABCD_lower_y_band_data->GetBinContent(nbins+1) << endl;
      if(fVerbose>2) cout << "\tData upper, events found : "  <<  nevU << endl
			  << "\t->Integral() : "  <<  h_ABCD_upper_y_band_data->Integral() + h_ABCD_upper_y_band_data->GetBinContent(nbins+1) << endl;
    }
    else if (fSamples[i].type == "mc" && fSamples[i].name == "QCD_Pt_120to170"){   // WARNING: checking statistical fluctuation in this pt bin!!!
      variable  = TString::Format("%s:misc.MT2>>+%s",var.Data(),h_ABCD_MT2_qcd->GetName());
      TString sel_u    = /*sel_up_QCD; /*/TString::Format("(%s) * (%s && misc.MT2<80&& %s)",weights.Data(),basecut.Data(),upper_cut.Data());
      TString sel_l    = /*sel_lo_QCD; /*/TString::Format("(%s) * (%s && misc.MT2<80&& %s)",weights.Data(),basecut.Data(),lower_cut.Data());
      int nev2d= fSamples[i].tree->Draw(variable,selection_QCD,"goff");
      variable  = TString::Format("%s>>+%s",var2.Data(),h_ABCD_lower_y_band_qcd->GetName());
      int nevL = fSamples[i].tree->Draw(variable,sel_l,"goff");
      variable  = TString::Format("%s>>+%s",var2.Data(),h_ABCD_upper_y_band_qcd->GetName());
      int nevU = fSamples[i].tree->Draw(variable,sel_u,"goff");
      if(fVerbose>2) cout << "\t" << fSamples[i].name << " lower, events found : "  <<  nevL << endl
			  << "\t->Integral() : "  <<  h_ABCD_lower_y_band_qcd->Integral() + h_ABCD_lower_y_band_qcd->GetBinContent(nbins+1) << endl;
      if(fVerbose>2) cout << "\t" << fSamples[i].name << " upper, events found : "  <<  nevU << endl
			  << "\t->Integral() : "  <<  h_ABCD_upper_y_band_qcd->Integral() + h_ABCD_upper_y_band_qcd->GetBinContent(nbins+1) << endl;
    }
    else if (fSamples[i].type == "mc" && fSamples[i].name == "QCD_Pt_170to300"){   // WARNING: checking statistical fluctuation in this pt bin!!!
      variable  = TString::Format("%s:misc.MT2>>+%s",var.Data(),h_ABCD_MT2_qcd->GetName());
      TString sel_u    = /*sel_up_QCD; /*/TString::Format("(%s) * (%s && misc.MT2<140&& %s)",weights.Data(),basecut.Data(),upper_cut.Data());
      TString sel_l    = /*sel_lo_QCD; /*/TString::Format("(%s) * (%s && misc.MT2<170&& %s)",weights.Data(),basecut.Data(),lower_cut.Data());
      int nev2d= fSamples[i].tree->Draw(variable,selection_QCD,"goff");
      variable  = TString::Format("%s>>+%s",var2.Data(),h_ABCD_lower_y_band_qcd->GetName());
      int nevL = fSamples[i].tree->Draw(variable,sel_l,"goff");
      variable  = TString::Format("%s>>+%s",var2.Data(),h_ABCD_upper_y_band_qcd->GetName());
      int nevU = fSamples[i].tree->Draw(variable,sel_u,"goff");
      if(fVerbose>2) cout << "\t" << fSamples[i].name << " lower, events found : "  <<  nevL << endl
			  << "\t->Integral() : "  <<  h_ABCD_lower_y_band_qcd->Integral() + h_ABCD_lower_y_band_qcd->GetBinContent(nbins+1) << endl;
      if(fVerbose>2) cout << "\t" << fSamples[i].name << " upper, events found : "  <<  nevU << endl
			  << "\t->Integral() : "  <<  h_ABCD_upper_y_band_qcd->Integral() + h_ABCD_upper_y_band_qcd->GetBinContent(nbins+1) << endl;
    }
    else if (fSamples[i].type == "mc" && fSamples[i].name == "QCD_Pt_300to470"){   // WARNING: checking statistical fluctuation in this pt bin!!!
      variable  = TString::Format("%s:misc.MT2>>+%s",var.Data(),h_ABCD_MT2_qcd->GetName());
      TString sel_u    = /*sel_up_QCD; /*/TString::Format("(%s) * (%s && misc.MT2<200&& %s)",weights.Data(),basecut.Data(),upper_cut.Data());
      TString sel_l    = /*sel_lo_QCD; /*/TString::Format("(%s) * (%s && %s)",weights.Data(),basecut.Data(),lower_cut.Data());
      int nev2d= fSamples[i].tree->Draw(variable,selection_QCD,"goff");
      variable  = TString::Format("%s>>+%s",var2.Data(),h_ABCD_lower_y_band_qcd->GetName());
      int nevL = fSamples[i].tree->Draw(variable,sel_l,"goff");
      variable  = TString::Format("%s>>+%s",var2.Data(),h_ABCD_upper_y_band_qcd->GetName());
      int nevU = fSamples[i].tree->Draw(variable,sel_u,"goff");
      if(fVerbose>2) cout << "\t" << fSamples[i].name << " lower, events found : "  <<  nevL << endl
			  << "\t->Integral() : "  <<  h_ABCD_lower_y_band_qcd->Integral() + h_ABCD_lower_y_band_qcd->GetBinContent(nbins+1) << endl;
      if(fVerbose>2) cout << "\t" << fSamples[i].name << " upper, events found : "  <<  nevU << endl
			  << "\t->Integral() : "  <<  h_ABCD_upper_y_band_qcd->Integral() + h_ABCD_upper_y_band_qcd->GetBinContent(nbins+1) << endl;
    }
    else if (fSamples[i].type == "mc" && fSamples[i].sname == "QCD"){
      variable  = TString::Format("%s:misc.MT2>>+%s",var.Data(),h_ABCD_MT2_qcd->GetName());
      int nev2d= fSamples[i].tree->Draw(variable,selection_QCD,"goff");
      variable  = TString::Format("%s>>+%s",var2.Data(),h_ABCD_lower_y_band_qcd->GetName());
      int nevL = fSamples[i].tree->Draw(variable,sel_lo_QCD,"goff");
      variable  = TString::Format("%s>>+%s",var2.Data(),h_ABCD_upper_y_band_qcd->GetName());
      int nevU = fSamples[i].tree->Draw(variable,sel_up_QCD,"goff");
      if(fVerbose>2) cout << "\t" << fSamples[i].name << " lower, events found : "  <<  nevL << endl
			  << "\t->Integral() : "  <<  h_ABCD_lower_y_band_qcd->Integral() + h_ABCD_lower_y_band_qcd->GetBinContent(nbins+1) << endl;
      if(fVerbose>2) cout << "\t" << fSamples[i].name << " upper, events found : "  <<  nevU << endl
			  << "\t->Integral() : "  <<  h_ABCD_upper_y_band_qcd->Integral() + h_ABCD_upper_y_band_qcd->GetBinContent(nbins+1) << endl;
    }
    else if (fSamples[i].type == "mc"){
      variable  = TString::Format("%s>>+%s",var2.Data(),h_ABCD_lower_y_band_mc->GetName());
      int nevL = fSamples[i].tree->Draw(variable,sel_lo,"goff");
      variable  = TString::Format("%s>>+%s",var2.Data(),h_ABCD_upper_y_band_mc->GetName());
      int nevU = fSamples[i].tree->Draw(variable,sel_up,"goff");
      if(fVerbose>2) cout << "\t" << fSamples[i].name << " lower, events found : "  <<  nevL << endl
			  << "\t->Integral() : "  <<  h_ABCD_lower_y_band_mc->Integral() + h_ABCD_lower_y_band_mc->GetBinContent(nbins+1) << endl;
      if(fVerbose>2) cout << "\t" << fSamples[i].name << " upper, events found : "  <<  nevU << endl
			  << "\t->Integral() : "  <<  h_ABCD_upper_y_band_mc->Integral() + h_ABCD_upper_y_band_mc->GetBinContent(nbins+1) << endl;
    }
    else if (fSamples[i].type == "susy"){
      variable  = TString::Format("%s:misc.MT2>>+%s",var.Data(),h_ABCD_MT2_susy->GetName());
      int nev2d= fSamples[i].tree->Draw(variable,selection,"goff");
      variable  = TString::Format("%s>>+%s",var2.Data(),h_ABCD_lower_y_band_susy->GetName());
      int nevL = fSamples[i].tree->Draw(variable,sel_lo,"goff");
      variable  = TString::Format("%s>>+%s",var2.Data(),h_ABCD_upper_y_band_susy->GetName());
      int nevU = fSamples[i].tree->Draw(variable,sel_up,"goff");
      if(fVerbose>2) cout << "\t" << fSamples[i].name << " lower, events found : "  <<  nevL << endl
			  << "\t->Integral() : "  <<  h_ABCD_lower_y_band_susy->Integral() + h_ABCD_lower_y_band_susy->GetBinContent(nbins+1) << endl;
      if(fVerbose>2) cout << "\t" << fSamples[i].name << " upper, events found : "  <<  nevU << endl
			  << "\t->Integral() : "  <<  h_ABCD_upper_y_band_susy->Integral() + h_ABCD_upper_y_band_susy->GetBinContent(nbins+1) << endl;
    }
  }
  
  
  h_ABCD_lower_y_band_sub->Add(h_ABCD_lower_y_band_data,h_ABCD_lower_y_band_mc,1,-1);
  h_ABCD_upper_y_band_sub->Add(h_ABCD_upper_y_band_data,h_ABCD_upper_y_band_mc,1,-1);

  for (int i= 0; i<=nbins; i++){
    if (h_ABCD_lower_y_band_sub->GetBinContent(i)<0)  h_ABCD_lower_y_band_sub->SetBinContent(i, 0.);
    if (h_ABCD_upper_y_band_sub->GetBinContent(i)<0)  h_ABCD_upper_y_band_sub->SetBinContent(i, 0.);
  }

  ratio_qcd ->Divide(h_ABCD_upper_y_band_qcd ,h_ABCD_lower_y_band_qcd );
  ratio_data->Divide(h_ABCD_upper_y_band_data,h_ABCD_lower_y_band_data);
  ratio_sub ->Divide(h_ABCD_upper_y_band_sub ,h_ABCD_lower_y_band_sub );
  
  h_ABCD_lower_y_band_data->SetLineColor(1);  h_ABCD_upper_y_band_data->SetLineColor(2);
  h_ABCD_lower_y_band_qcd ->SetLineColor(1);  h_ABCD_upper_y_band_qcd ->SetLineColor(2);
  h_ABCD_lower_y_band_sub ->SetLineColor(1);  h_ABCD_upper_y_band_sub ->SetLineColor(2);
  h_ABCD_lower_y_band_sub ->SetLineStyle(3);  h_ABCD_upper_y_band_sub ->SetLineStyle(3);   ratio_sub->SetLineStyle(3);
  ratio_qcd ->SetMarkerStyle(24);   ratio_qcd ->SetMarkerSize(0.8);
  ratio_data->SetMarkerStyle(24);   ratio_data->SetMarkerSize(0.8);
  ratio_sub ->SetMarkerStyle(24);   ratio_sub ->SetMarkerSize(0.8);
  ratio_qcd ->SetXTitle("M_{T2}");   ratio_qcd ->SetYTitle("ratio");
  ratio_sub ->SetXTitle("M_{T2}");   ratio_sub ->SetYTitle("ratio");

  //TF1 *f_qcd = new TF1("f_qcd","pol0(0)+expo(1)",bins[0], bins[nbins]);   f_qcd->SetLineColor(8);
  //TF1 *f_sub = new TF1("f_sub","pol0(0)+expo(1)",bins[0], bins[nbins]);   f_sub->SetLineColor(8);
  //TF1 *f_qcd = new TF1("f_qcd","exp([0]+300.*[1])+expo(0)",bins[0], bins[nbins]);   f_qcd->SetLineColor(8);
  //TF1 *f_sub = new TF1("f_sub","exp([0]+300.*[1])+expo(0)",bins[0], bins[nbins]);   f_sub->SetLineColor(8);
  TF1 *f_qcd1 = new TF1("f_qcd1","expo(0)",bins[0], bins[nbins]);   f_qcd1->SetLineColor(8);
  TF1 *f_sub1 = new TF1("f_sub1","expo(0)",bins[0], bins[nbins]);   f_sub1->SetLineColor(8);
  TF1 *f_qcd2 = new TF1("f_qcd2","exp([0]+200.*[1])+expo(0)",bins[0], bins[nbins]);   f_qcd2->SetLineColor(9);
  TF1 *f_sub2 = new TF1("f_sub2","exp([0]+200.*[1])+expo(0)",bins[0], bins[nbins]);   f_sub2->SetLineColor(8);
  TF1 *f_qcd3 = new TF1("f_qcd3","expo(0)+[2]"              ,bins[0], bins[nbins]);   f_qcd3->SetLineColor(50);

  gStyle->SetPalette(1);
  gStyle->SetOptFit (0);
  gStyle->SetOptStat(0);

  TCanvas *can = new TCanvas("can", "2d qcd", 0, 0, 900, 700);
  can->SetLogz(1);
  h_ABCD_MT2_qcd->SetMinimum(0.0001);
  h_ABCD_MT2_qcd->Draw("colz");
  TCanvas *ccan = new TCanvas("ccan", "2d susy", 0, 0, 900, 700);
  ccan->SetLogz(1);
  h_ABCD_MT2_susy->SetMinimum(0.0001);
  h_ABCD_MT2_susy->Draw("colz");
  TCanvas *can2 = new TCanvas("can2", "qcd", 0, 0, 900, 700);
  can2->SetLogy(1);
  h_ABCD_lower_y_band_qcd ->Draw();  
  h_ABCD_upper_y_band_qcd ->Draw("same");
  TCanvas *can3 = new TCanvas("can3", "data w/ and w/o substraction", 0, 0, 900, 700);
  can3->SetLogy(1);
  h_ABCD_lower_y_band_data->Draw();  
  h_ABCD_upper_y_band_data->Draw("same");
  h_ABCD_lower_y_band_sub ->Draw("same");  
  h_ABCD_upper_y_band_sub ->Draw("same");
  TCanvas *can4 = new TCanvas("can4", "ratio qcd", 0, 0, 900, 700);
  can4->SetLogy(1);
  ratio_qcd->Draw("E1");
  //f_qcd->FixParameter(0,0.006);
  ratio_qcd->SetTitle("QCD");
  ratio_qcd->Fit("f_qcd1","0","",fit_min,fit_max);
  f_qcd1->Draw("same");
  f_qcd2->SetParameter(0,f_qcd1->GetParameter(0));
  f_qcd2->SetParameter(1,f_qcd1->GetParameter(1));
  f_qcd2->Draw("same");
  f_qcd3->SetParameter(0,f_qcd1->GetParameter(0));
  f_qcd3->SetParameter(1,f_qcd1->GetParameter(1));
  f_qcd3->SetParameter(2,0.002);
  ratio_qcd->Fit("f_qcd3","0","",fit_min,bins[nbins]);
  f_qcd3->Draw("same");
  TLine *l_qcd_min = new TLine(fit_min,f_qcd1->GetYaxis()->GetXmin(),fit_min,f_qcd1->GetMaximum()); l_qcd_min->SetLineStyle(2); l_qcd_min->Draw("same");
  TLine *l_qcd_max = new TLine(fit_max,f_qcd1->GetYaxis()->GetXmin(),fit_max,f_qcd1->GetMaximum()); l_qcd_max->SetLineStyle(2); l_qcd_max->Draw("same");

  TCanvas *can5 = new TCanvas("can5", "ratio data", 0, 0, 900, 700);
  can5->SetLogy(1);
  ratio_sub ->Draw("PE1");
  ratio_data->Draw("sameE1");
  //float p0 = f_qcd->GetParameter(0);
  //f_sub->FixParameter(0,0.006);
  //ratio_data->Fit("f_sub1","0","",fit_min,fit_max);  // fit w/o subtracting background
  ratio_sub->Fit("f_sub1","0","",fit_min,fit_max);   // background subtracted fit
  f_sub1->Draw("same");
  f_sub2->SetParameter(0,f_sub1->GetParameter(0));
  f_sub2->SetParameter(1,f_sub1->GetParameter(1));
  f_sub2->Draw("same");
  TLine *l_sub_min = new TLine(fit_min,f_sub1->GetYaxis()->GetXmin(),fit_min,f_sub1->GetMaximum()); l_sub_min->SetLineStyle(2); l_sub_min->Draw("same");
  TLine *l_sub_max = new TLine(fit_max,f_sub1->GetYaxis()->GetXmin(),fit_max,f_sub1->GetMaximum()); l_sub_max->SetLineStyle(2); l_sub_max->Draw("same");
  TLegend *leg = new TLegend(.55,.5,.85,.75);
  leg->SetFillColor(0);
  leg->AddEntry(ratio_data,"data","le");
  leg->AddEntry(ratio_sub ,"data (EWK subtracted)","le");
  leg->Draw("same");

  PrintABCDPredictions(var, basecut, upper_cut, lower_cut, f_qcd1,f_sub1,f_qcd3);
  
}

//_________________________________________________________________________________
void MassPlotter::PrintABCDPredictions(TString var, TString basecut, TString upper_cut, TString lower_cut, TF1* func_qcd, TF1* func_sub, TF1* func_qcd_model){  
  int nbins=180;
  float min=100., max=1000.;
  TH1D *h_pred_lower_y_band_data = new TH1D("pred_lower_y_band__data" , "", nbins,min,max);		        h_pred_lower_y_band_data->Sumw2();
  TH1D *h_pred_upper_y_band_data = new TH1D("pred_upper_y_band__data" , "", nbins,min,max);		        h_pred_upper_y_band_data->Sumw2();
  TH1D *h_pred_lower_y_band_qcd  = new TH1D("pred_lower_y_band__qcd"  , "", nbins,min,max);		        h_pred_lower_y_band_qcd ->Sumw2();
  TH1D *h_pred_upper_y_band_qcd  = new TH1D("pred_upper_y_band__qcd"  , "", nbins,min,max);		        h_pred_upper_y_band_qcd ->Sumw2();
  TH1D *h_pred_lower_y_band_mc   = new TH1D("pred_lower_y_band__mc"   , "", nbins,min,max);		        h_pred_lower_y_band_mc  ->Sumw2();
  TH1D *h_pred_upper_y_band_mc   = new TH1D("pred_upper_y_band__mc"   , "", nbins,min,max);		        h_pred_upper_y_band_mc  ->Sumw2();
  TH1D *h_pred_lower_y_band_sub  = new TH1D("pred_lower_y_band__sub"  , "", nbins,min,max);		        h_pred_lower_y_band_sub ->Sumw2();
  TH1D *h_pred_upper_y_band_sub  = new TH1D("pred_upper_y_band__sub"  , "", nbins,min,max);		        h_pred_upper_y_band_sub ->Sumw2();
  TH1D *h_pred_lower_y_band_susy = new TH1D("pred_lower_y_band__susy" , "", nbins,min,max);		        h_pred_lower_y_band_susy->Sumw2();
  TH1D *h_pred_upper_y_band_susy = new TH1D("pred_upper_y_band__susy" , "", nbins,min,max);		        h_pred_upper_y_band_susy->Sumw2();
  //TF1 *f_qcd = new TF1("f_qcd","pol0(0)+expo(1)",100,1000);
  //TF1 *f_sub = new TF1("f_sub","pol0(0)+expo(1)",100,1000);
//   TF1 *f_qcd = new TF1("f_qcd","exp([0]+300.*[1])+expo(0)",100,1000);
//   TF1 *f_sub = new TF1("f_sub","exp([0]+300.*[1])+expo(0)",100,1000);
  TF1 *f_qcd = new TF1("f_qcd","expo(0)",min,max);
  TF1 *f_sub = new TF1("f_sub","expo(0)",min,max);
  f_qcd->SetParameters(func_qcd->GetParameter(0),func_qcd->GetParameter(1));//,func_qcd->GetParameter(2));
  f_sub->SetParameters(func_sub->GetParameter(0),func_sub->GetParameter(1));//,func_sub->GetParameter(2));
  TF1 *f_qcd2 = new TF1("f_qcd2","exp([0]+200.*[1])+expo(0)",min,max);
  TF1 *f_sub2 = new TF1("f_sub2","exp([0]+200.*[1])+expo(0)",min,max);
  f_qcd2->SetParameters(func_qcd->GetParameter(0),func_qcd->GetParameter(1));//,func_qcd->GetParameter(2));
  f_sub2->SetParameters(func_sub->GetParameter(0),func_sub->GetParameter(1));//,func_sub->GetParameter(2));
  TF1 *f_qcd_model = new TF1("f_qcd_model","expo(0)+[2]",min,max);
  f_qcd_model->SetParameters(func_qcd_model->GetParameter(0),func_qcd_model->GetParameter(1),func_qcd_model->GetParameter(2));

  for(size_t i = 0; i < fSamples.size(); ++i){
   
    Long64_t nentries =  fSamples[i].tree->GetEntries();
    Long64_t nbytes = 0, nb = 0;
    
    Double_t weight = fSamples[i].xsection * fSamples[i].kfact * fSamples[i].lumi / (fSamples[i].nevents);
    if(fVerbose>2) cout << "ABCD: looping over " << fSamples[i].name << endl;
    if(fVerbose>2) cout << "      sample has weight " << weight << endl; 

    TString weights   = fSamples[i].type!="data" ?  TString::Format("(%.15f*pileUp.Weight)",weight) : TString::Format("(%.15f)",weight);
    TString selection = TString::Format("(%s) * (%s)"      ,weights.Data(),basecut.Data());
    TString sel_up    = TString::Format("(%s) * (%s && %s)",weights.Data(),basecut.Data(),upper_cut.Data());
    TString sel_lo    = TString::Format("(%s) * (%s && %s)",weights.Data(),basecut.Data(),lower_cut.Data());
    if(fVerbose>2) cout << "      sel_lo = " << weights << endl; 
//     TString selection = TString::Format("(%f) * (%s)"      ,weight,basecut.Data());
//     TString sel_up    = TString::Format("(%f) * (%s && %s)",weight,basecut.Data(),upper_cut.Data());
//     TString sel_lo    = TString::Format("(%f) * (%s && %s)",weight,basecut.Data(),lower_cut.Data());

    TString variable;
    TString var2 = "misc.MT2";
    if (fSamples[i].type == "data"){
      variable  = TString::Format("%s>>+%s",var2.Data(),h_pred_lower_y_band_data->GetName());
      int nevL = fSamples[i].tree->Draw(variable,sel_lo,"goff");
      variable  = TString::Format("%s>>+%s",var2.Data(),h_pred_upper_y_band_data->GetName());
      int nevU = fSamples[i].tree->Draw(variable,sel_up,"goff");
    }
    else if (fSamples[i].type == "mc" && fSamples[i].name == "QCD_Pt_120to170"){   // WARNING: checking statistical fluctuation in this pt bin!!!
      TString sel_u    = TString::Format("(%s) * (%s && misc.MT2<80&& %s)",weights.Data(),basecut.Data(),upper_cut.Data());
      TString sel_l    = TString::Format("(%s) * (%s && misc.MT2<80&& %s)",weights.Data(),basecut.Data(),lower_cut.Data());
      variable  = TString::Format("%s>>+%s",var2.Data(),h_pred_lower_y_band_qcd->GetName());
      int nevL = fSamples[i].tree->Draw(variable,sel_l,"goff");
      variable  = TString::Format("%s>>+%s",var2.Data(),h_pred_upper_y_band_qcd->GetName());
      int nevU = fSamples[i].tree->Draw(variable,sel_u,"goff");
    }
    else if (fSamples[i].type == "mc" && fSamples[i].name == "QCD_Pt_170to300"){   // WARNING: checking statistical fluctuation in this pt bin!!!
      TString sel_u    = TString::Format("(%s) * (%s && misc.MT2<140&& %s)",weights.Data(),basecut.Data(),upper_cut.Data());
      TString sel_l    = TString::Format("(%s) * (%s && misc.MT2<170&& %s)",weights.Data(),basecut.Data(),lower_cut.Data());
      variable  = TString::Format("%s>>+%s",var2.Data(),h_pred_lower_y_band_qcd->GetName());
      int nevL = fSamples[i].tree->Draw(variable,sel_l,"goff");
      variable  = TString::Format("%s>>+%s",var2.Data(),h_pred_upper_y_band_qcd->GetName());
      int nevU = fSamples[i].tree->Draw(variable,sel_u,"goff");
    }
    else if (fSamples[i].type == "mc" && fSamples[i].name == "QCD_Pt_300to470"){   // WARNING: checking statistical fluctuation in this pt bin!!!
      TString sel_u    = TString::Format("(%s) * (%s && misc.MT2<200&& %s)",weights.Data(),basecut.Data(),upper_cut.Data());
      TString sel_l    = TString::Format("(%s) * (%s && %s)",weights.Data(),basecut.Data(),lower_cut.Data());
      variable  = TString::Format("%s>>+%s",var2.Data(),h_pred_lower_y_band_qcd->GetName());
      int nevL = fSamples[i].tree->Draw(variable,sel_l,"goff");
      variable  = TString::Format("%s>>+%s",var2.Data(),h_pred_upper_y_band_qcd->GetName());
      int nevU = fSamples[i].tree->Draw(variable,sel_u,"goff");
    }
    else if (fSamples[i].type == "mc" && fSamples[i].sname == "QCD"){
      variable  = TString::Format("%s>>+%s",var2.Data(),h_pred_lower_y_band_qcd->GetName());
      int nevL = fSamples[i].tree->Draw(variable,sel_lo,"goff");
      variable  = TString::Format("%s>>+%s",var2.Data(),h_pred_upper_y_band_qcd->GetName());
      int nevU = fSamples[i].tree->Draw(variable,sel_up,"goff");
    }
    else if (fSamples[i].type == "mc"){
      variable  = TString::Format("%s>>+%s",var2.Data(),h_pred_lower_y_band_mc->GetName());
      int nevL = fSamples[i].tree->Draw(variable,sel_lo,"goff");
      variable  = TString::Format("%s>>+%s",var2.Data(),h_pred_upper_y_band_mc->GetName());
      int nevU = fSamples[i].tree->Draw(variable,sel_up,"goff");
    }
    else if (fSamples[i].type == "susy"){
      variable  = TString::Format("%s>>+%s",var2.Data(),h_pred_lower_y_band_susy->GetName());
      int nevL = fSamples[i].tree->Draw(variable,sel_lo,"goff");
      variable  = TString::Format("%s>>+%s",var2.Data(),h_pred_upper_y_band_susy->GetName());
      int nevU = fSamples[i].tree->Draw(variable,sel_up,"goff");
    }
  }
  
  
  h_pred_lower_y_band_sub ->Add(h_pred_lower_y_band_data,h_pred_lower_y_band_mc,1,-1);
  h_pred_upper_y_band_sub ->Add(h_pred_upper_y_band_data,h_pred_upper_y_band_mc,1,-1);
  h_pred_lower_y_band_susy->Add(h_pred_lower_y_band_sub); // adding susy contamination to subtracted data

  for (int i= 0; i<=nbins+1; i++){
    if (h_pred_lower_y_band_sub->GetBinContent(i)<0)  h_pred_lower_y_band_sub->SetBinContent(i, 0.);
    if (h_pred_upper_y_band_sub->GetBinContent(i)<0)  h_pred_upper_y_band_sub->SetBinContent(i, 0.);
  }

  TH1D *h_pred_lower_y_band_data_2   = (TH1D*)h_pred_lower_y_band_data->Clone("h_pred_lower_y_band_data_2");         h_pred_lower_y_band_data_2   ->Sumw2();
  TH1D *h_pred_lower_y_band_qcd_2    = (TH1D*)h_pred_lower_y_band_qcd ->Clone("h_pred_lower_y_band_qcd_2" );         h_pred_lower_y_band_qcd_2    ->Sumw2();
  TH1D *h_pred_lower_y_band_sub_2    = (TH1D*)h_pred_lower_y_band_sub ->Clone("h_pred_lower_y_band_sub_2" );         h_pred_lower_y_band_sub_2    ->Sumw2();
  TH1D *h_pred_lower_y_band_data_c   = (TH1D*)h_pred_lower_y_band_data->Clone("h_pred_lower_y_band_data_c");         h_pred_lower_y_band_data_c   ->Sumw2();
  TH1D *h_pred_lower_y_band_qcd_c    = (TH1D*)h_pred_lower_y_band_qcd ->Clone("h_pred_lower_y_band_qcd_c" );         h_pred_lower_y_band_qcd_c    ->Sumw2();
  TH1D *h_pred_lower_y_band_sub_c    = (TH1D*)h_pred_lower_y_band_sub ->Clone("h_pred_lower_y_band_sub_c" );         h_pred_lower_y_band_sub_c    ->Sumw2();
  TH1D *h_pred_lower_y_band_susy_c   = (TH1D*)h_pred_lower_y_band_susy->Clone("h_pred_lower_y_band_susy_c" );        h_pred_lower_y_band_susy_c  ->Sumw2();
  TH1D *h_pred_lower_y_band_data_2_c = (TH1D*)h_pred_lower_y_band_data->Clone("h_pred_lower_y_band_data_2_c");       h_pred_lower_y_band_data_2_c ->Sumw2();
  TH1D *h_pred_lower_y_band_qcd_2_c  = (TH1D*)h_pred_lower_y_band_qcd ->Clone("h_pred_lower_y_band_qcd_2_c" );       h_pred_lower_y_band_qcd_2_c  ->Sumw2();
  TH1D *h_pred_lower_y_band_sub_2_c  = (TH1D*)h_pred_lower_y_band_sub ->Clone("h_pred_lower_y_band_sub_2_c" );       h_pred_lower_y_band_sub_2_c  ->Sumw2();
  TH1D *h_pred_lower_y_band_qcd_model= (TH1D*)h_pred_lower_y_band_qcd ->Clone("h_pred_lower_y_band_qcd_model");      h_pred_lower_y_band_qcd_model->Sumw2();

  // from qcd fit (exp+const) in whole range (50-500)
  h_pred_lower_y_band_qcd_model->Multiply(f_qcd_model);
  // from qcd fit

  h_pred_lower_y_band_data  ->Multiply(f_qcd);
  h_pred_lower_y_band_qcd   ->Multiply(f_qcd);
  h_pred_lower_y_band_sub   ->Multiply(f_qcd);
  // from data fit
  h_pred_lower_y_band_data_2->Multiply(f_sub);
  h_pred_lower_y_band_qcd_2 ->Multiply(f_sub);
  h_pred_lower_y_band_sub_2 ->Multiply(f_sub);
  // from qcd fit (exp+const)
  h_pred_lower_y_band_data_c  ->Multiply(f_qcd2);
  h_pred_lower_y_band_qcd_c   ->Multiply(f_qcd2);
  h_pred_lower_y_band_sub_c   ->Multiply(f_qcd2);
  // from data fit (exp+const)
  h_pred_lower_y_band_data_2_c->Multiply(f_sub2);
  h_pred_lower_y_band_qcd_2_c ->Multiply(f_sub2);
  h_pred_lower_y_band_sub_2_c ->Multiply(f_sub2);
  // with susy contamination from data fit (exp+const)
  h_pred_lower_y_band_susy  ->Multiply(f_sub);
  h_pred_lower_y_band_susy_c->Multiply(f_sub2);

  cout << "SUSY prediction:" << endl;
  printEstimation(h_pred_upper_y_band_susy,h_pred_upper_y_band_susy, nbins, min, max);
  cout << "QCD prediction:" << endl;
  printEstimation(h_pred_upper_y_band_qcd,h_pred_upper_y_band_qcd, nbins, min, max);
  cout << "------------------------------------" << endl;
  cout << "QCD estimation from QCD model (exp+const):" << endl;
  printEstimation(h_pred_lower_y_band_qcd_model,h_pred_lower_y_band_qcd_model, nbins, min, max);
  cout << "------------------------------------" << endl;
  cout << "QCD estimation from QCD (fit to QCD):" << endl;
  printEstimation(h_pred_lower_y_band_qcd,h_pred_lower_y_band_qcd_c, nbins, min, max);
  cout << "------------------------------------" << endl;
  cout << "QCD estimation from QCD (fit to data):" << endl;
  printEstimation(h_pred_lower_y_band_qcd_2,h_pred_lower_y_band_qcd_2_c, nbins, min, max);
  cout << "------------------------------------" << endl;
  cout << "QCD estimation from data (fit to QCD):" << endl;
  printEstimation(h_pred_lower_y_band_data,h_pred_lower_y_band_data_c, nbins, min, max);
  cout << "------------------------------------" << endl;
  cout << "QCD estimation from data (fit to data):" << endl;
  printEstimation(h_pred_lower_y_band_data_2,h_pred_lower_y_band_data_2_c, nbins, min, max);
  cout << "------------------------------------" << endl;
  cout << "QCD estimation from EWK subtracted data (fit to QCD):" << endl;
  printEstimation(h_pred_lower_y_band_sub,h_pred_lower_y_band_sub_c, nbins, min, max);
  cout << "------------------------------------" << endl;
  cout << "QCD estimation from EWK subtracted data (fit to data):" << endl;
  printEstimation(h_pred_lower_y_band_sub_2,h_pred_lower_y_band_sub_2_c, nbins, min, max);
  cout << "------------------------------------" << endl;
  cout << "QCD estimation from QCD with susy contamination (fit to data):" << endl;
  printEstimation(h_pred_lower_y_band_susy,h_pred_lower_y_band_susy_c, nbins, min, max);
  cout << "------------------------------------" << endl;

}

void MassPlotter::printEstimation(TH1D* h_pred, TH1D* h_pred_c, int nbins, float min, float max){
  float delta = (max-min)/((float)nbins);
  double yield[5][2], error[5][2];
  yield[0][0] = Util::IntegralAndError(h_pred,  (int)((100.-min)/delta)+1, nbins+1, error[0][0]);
  cout << "\tMT2>100 = " << yield[0][0]  << " +/- " << error[0][0];
  yield[0][1] = Util::IntegralAndError(h_pred_c,(int)((100.-min)/delta)+1, nbins+1, error[0][1]);
  cout << "\tMT2>100 = " << yield[0][1]  << " +/- " << error[0][1] << endl;
  yield[1][0] = Util::IntegralAndError(h_pred,  (int)((150.-min)/delta)+1, nbins+1, error[1][0]);
  cout << "\tMT2>150 = " << yield[1][0]  << " +/- " << error[1][0];
  yield[1][1] = Util::IntegralAndError(h_pred_c,(int)((150.-min)/delta)+1, nbins+1, error[1][1]);
  cout << "\tMT2>150 = " << yield[1][1]  << " +/- " << error[1][1] << endl;
  yield[2][0] = Util::IntegralAndError(h_pred,  (int)((200.-min)/delta)+1, nbins+1, error[2][0]);
  cout << "\tMT2>200 = " << yield[2][0]  << " +/- " << error[2][0];
  yield[2][1] = Util::IntegralAndError(h_pred_c,(int)((200.-min)/delta)+1, nbins+1, error[2][1]);
  cout << "\tMT2>200 = " << yield[2][1]  << " +/- " << error[2][1] << endl;
  yield[3][0] = Util::IntegralAndError(h_pred,  (int)((250.-min)/delta)+1, nbins+1, error[3][0]);
  cout << "\tMT2>250 = " << yield[3][0]  << " +/- " << error[3][0];
  yield[3][1] = Util::IntegralAndError(h_pred_c,(int)((250.-min)/delta)+1, nbins+1, error[3][1]);
  cout << "\tMT2>250 = " << yield[3][1]  << " +/- " << error[3][1] << endl;
  yield[4][0] = Util::IntegralAndError(h_pred,  (int)((300.-min)/delta)+1, nbins+1, error[4][0]);
  cout << "\tMT2>300 = " << yield[4][0]  << " +/- " << error[4][0];
  yield[4][1] = Util::IntegralAndError(h_pred_c,(int)((300.-min)/delta)+1, nbins+1, error[4][1]);
  cout << "\tMT2>300 = " << yield[4][1]  << " +/- " << error[4][1] << endl;

  for (int i=0; i<2; i++){
    cout << "yield_" << (i==0 ? "opt" : "pes") << " = {" << endl;
    for (int j=0; j<5; j++){
      cout << yield[j][i] << (j<4 ? ", " : "\n}\n");
    }
    cout << "error_" << (i==0 ? "opt" : "pes") << " = {" << endl;
    for (int j=0; j<5; j++){
      cout << error[j][i] << (j<4 ? ", " : "\n}\n");
    }
  }

}

//_________________________________________________________________________________
void MassPlotter::plotRatioStack(THStack* hstack, TH1* h1_orig, TH1* h2_orig, bool logflag, bool normalize, TString name, TLegend* leg, TString xtitle, TString ytitle,int njets, int nleps, TString saveMacro){
	// define canvas and pads 
	TH1D *h1 = (TH1D*)h1_orig->Clone("h1_copy");
	TH1D *h2 = (TH1D*)h2_orig->Clone("h2_copy");

	h1->SetStats(0);	
	h2->SetStats(0);	
	
	TCanvas* c1 = new TCanvas(name,"", 20,100,1000,700);
	c1 -> cd();
	
	float border = 0.3;
 	float scale = (1-border)/border;
 
 	TPad *p_plot  = new TPad(name+"_plotpad",  "Pad containing the overlay plot", 0.00, border, 1.00, 1.00, 0, 0);
 	p_plot->SetBottomMargin(0);
 	p_plot->Draw();
 	TPad *p_ratio = new TPad(name+"_ratiopad", "Pad containing the ratio",        0.00, 0.00, 1.00, border, 0, 0);
 	p_ratio->SetTopMargin(0);
 	p_ratio->SetBottomMargin(0.35);
 	p_ratio->Draw();
 
	// draw overlay plot
 	p_plot ->cd();

	if(logflag) gPad->SetLogy(1);
	gPad->SetFillStyle(0);
		
	// Scaling
	if(normalize){
		h1->Scale(1.0/h1->Integral());
		h2->Scale(1.0/h2->Integral());
	}
	
	// Determine plotting range
	double max1 = h1->GetMaximum();
	double max2 = h2->GetMaximum();
	double max  = (max1>max2)?max1:max2;
	if(logflag) max = 2*max;
	else max = 1.05*max;
	h1->SetMaximum(max);
	h2->SetMaximum(max);
	hstack->SetMaximum(max);
//	h1->SetMinimum(0.000000001);
//	h2->SetMinimum(0.000000001);

	hstack->SetMinimum(0.02);
	hstack ->Draw("hist");
	//h2     ->Draw("sameEX0");
        h2     ->Draw("sameE1"); //LEO MOD

	// title
	TLatex lat;
	lat.SetNDC(1);
	lat.SetTextColor(4);
	lat.SetTextSize(0.1);

	// y axis title 
	lat.SetTextAlign(33); 
	lat.SetTextColor(1);
	lat.SetTextAngle(90);
	lat.DrawLatex(0.035, 0.9, ytitle);
		
	TLatex TitleBox;
	TitleBox.SetNDC();
	TitleBox.SetTextSize(0.05);
	TString text;
	if (njets>=10)
	  text = TString::Format("%d-%d Jets",njets/10,njets%10);
	else
	  text = njets < 0 ? TString::Format("#geq %d Jets",abs(njets)) : TString::Format("%d Jets",abs(njets));
	text += nleps == 1 ? ", 1 Lepton" : "";
	TitleBox.DrawLatex(0.18,0.943,text.Data());

 	p_plot ->Draw();
	gPad->RedrawAxis();

	if(leg != NULL ){
		leg -> SetFillColor(0);
		leg -> SetBorderSize(0);
		leg -> Draw();
	} 
	
	// draw the ratio plot
 	p_ratio ->cd();
	//gPad->SetLogy();
 
	TH1D *h_ratio = (TH1D*)h2_orig->Clone("h2_copy");	

	h_ratio->GetXaxis()->SetTitleSize(scale * h1->GetXaxis()->GetTitleSize());
	h_ratio->GetYaxis()->SetTitleSize(scale * h1->GetYaxis()->GetTitleSize());
	h_ratio->GetXaxis()->SetLabelSize(scale * h1->GetXaxis()->GetLabelSize());
	h_ratio->GetYaxis()->SetLabelSize(scale * h1->GetYaxis()->GetLabelSize());
	h_ratio->GetXaxis()->SetTickLength(scale * h1->GetXaxis()->GetTickLength());
	h_ratio->GetYaxis()->SetTickLength(h1->GetYaxis()->GetTickLength());
	h_ratio->SetLineWidth(2);
	//h_ratio->SetFillColor(kBlue);//leo
	//h_ratio->SetLineColor(kBlue);
	h_ratio ->SetStats(0);
	h_ratio ->SetMarkerStyle(20);
	//h_ratio ->SetMarkerSize(0.1);

 	h_ratio ->Divide(h2, h1);
 	h_ratio ->SetMinimum(0.4);
 	h_ratio ->SetMaximum(3.0);
	h_ratio ->GetYaxis()->SetTitleOffset(h1->GetYaxis()->GetTitleOffset());

	h_ratio ->DrawCopy("E"/*2"*/);
 
	TLine *l3 = new TLine(h1->GetXaxis()->GetXmin(), 1.00, h1->GetXaxis()->GetXmax(), 1.00);
	l3->SetLineWidth(2);
	l3->SetLineStyle(7);
	l3->Draw();
	// ratio title
	lat.SetTextAlign(33);
	lat.SetTextColor(1);
	lat.SetTextSize(0.15);
	lat.SetNDC(1);
	lat.SetTextAngle(90);
	lat.DrawLatex(0.035, 1.0, "data / MC");
	
	// x axis title
	lat.SetTextSize(0.2);
	lat.SetTextAngle(0);
	float ypos = xtitle.Contains("slash") || xtitle.Contains("_") || xtitle.Contains("^") ? 0.23 : 0.16;
	lat.DrawLatex(0.9, ypos, xtitle);
	gPad->RedrawAxis();
 	p_ratio ->Draw();
 	c1->Update();

	TString save=name+"_ratio";
	if(fSave) Util::Print(c1, save, fOutputDir, fOutputFile);
	if(saveMacro == "gif")
		c1->SaveAs(save+".gif");

// 	delete h1;
// 	delete h2;
// 	delete h_ratio;
// 	delete p_plot;
// 	delete p_ratio;
// 	delete c1;

}
//_________________________________________________________________________________
void MassPlotter::plotRatioStack(THStack* hstack, TH1* h1_orig, TH1* h2_orig, TH1* h3, bool logflag, bool normalize, TString name, TLegend* leg, TString xtitle, TString ytitle,int njets,int nbjets, int nleps, float overlayScale, TString saveMacro){
  //LEO TRUE USE THIS

	// define canvas and pads 
	TH1D *h1 = (TH1D*)h1_orig->Clone("h1_copy");
	TH1D *h2 = (TH1D*)h2_orig->Clone("h2_copy");

	h1->SetStats(0);	
	h2->SetStats(0);	
	h1->SetMarkerStyle(20);
	h2->SetMarkerStyle(20);

	
	//TCanvas* c1 = new TCanvas(name,"", 20,100,1000,700);
	TCanvas* c1 = new TCanvas(name+"c_ratio","",0,0,600,600 /*37, 60,636,670*/);
	c1->SetFrameLineWidth(1);
	c1 -> cd();
	
	float border = 0.2;
 	float scale = (1-border)/border;
 
 	TPad *p_plot  = new TPad(name+"_plotpad",  "Pad containing the overlay plot", 0,0.211838,1,1 /*0.00, border, 1.00, 1.00, 0, 0*/);
 	//p_plot->SetBottomMargin(0.05);
	//p_plot->SetTopMargin(0.09);
	//p_plot->SetLeftMargin(0.1669107);
      	//p_plot->SetRightMargin(0.02);
	p_plot->SetLeftMargin(0.131579);
	p_plot->SetRightMargin(0.08);
	p_plot->SetTopMargin(0.06895515);
	p_plot->SetBottomMargin(0.07206074);
 	p_plot->Draw();
 	TPad *p_ratio = new TPad(name+"_ratiopad", "Pad containing the ratio",   0,0.01863354,0.9967105,0.2189441/*     0.00, 0.05, 1.00, border, 0, 0*/);
 	//p_ratio->SetTopMargin(0.03);
 	//p_ratio->SetBottomMargin(0.05/*5*/);
	//p_ratio->SetRightMargin(0.02);
	p_ratio->SetLeftMargin(0.1336634);	
	p_ratio->SetRightMargin(0.075);
	p_ratio->SetTopMargin(0.06976745);
	p_ratio->SetBottomMargin(0.2790698);

	p_ratio->Draw();
 
	// draw overlay plot
 	p_plot ->cd();

	if(logflag) gPad->SetLogy(1);
	gPad->SetFillStyle(0);
		
	// Scaling
	if(normalize){
		h1->Scale(1.0/h1->Integral());
		h2->Scale(1.0/h2->Integral());
	}
	
	// Determine plotting range
	double max1 = h1->GetMaximum();
	double max2 = h2->GetMaximum();
	double max  = (max1>max2)?max1:max2;
	if(logflag) max = 2.5*max;
	else max = 1.5*max;
	h1->SetMaximum(max);
	h2->SetMaximum(max);
	hstack->SetMaximum(max);
//	h1->SetMinimum(0.000000001);
//	h2->SetMinimum(0.000000001);
	stringstream yTitle;
	if(fEventsPerGeV){
		if(fabs(h1_orig->GetBinWidth(1) -h1_orig->GetBinWidth(h1_orig->GetNbinsX()-1))<0.01){
		double binwidth = h1_orig->GetBinWidth(1);
		yTitle.precision(3);
		yTitle << ytitle.Data();
		yTitle << " / ";
		yTitle << binwidth;
		yTitle << " GeV";
		} else{
		cout << h1_orig->GetBinWidth(1) << " " << h1_orig->GetBinWidth(h1_orig->GetNbinsX()-1) << endl;
		}
	}else{
		yTitle << ytitle.Data();
	}
	hstack->GetYaxis()->SetTitle(yTitle.str().c_str());
	hstack->GetYaxis()->SetLabelSize(0.05);
	hstack->GetYaxis()->SetTitleSize(0.05);
	hstack->GetYaxis()->SetTitleOffset(1.3);

	
	//MT2_bSel[0]->SetTitleSize(0.03);
	///MT2_bSel[0]->SetTitleOffset(1.);
	hstack->SetMinimum(0.02);
	hstack->Draw("hist");
	h2    ->Draw("sameE");
	h3->Scale(overlayScale ? overlayScale : h2->Integral() / h3->Integral());
	h3->SetFillColor(0);
	h3->SetLineStyle(kDotted);
	h3->SetLineWidth(4);
	h3->Draw("samehist");

	TLatex TitleBox;
	TitleBox.SetNDC();
	TitleBox.SetTextSize(0.03);

	TString text;
	if (njets>=10)
	  text = TString::Format("%d-%d jets",njets/10,njets%10);
	else
	  text = njets < 0 ? TString::Format("#geq %d jets",abs(njets)) : TString::Format("%d jets",abs(njets));
	text += nbjets==-10 ? "" : nbjets < 0 ? TString::Format(", #geq %d b-tag",abs(nbjets)) : TString::Format(", %d b-tag",abs(nbjets));
	text += nleps == 1 ? ", 1 lepton" : "";


// 	TString text ="";
// 	text = fMT2Analysis?  "M_{T2} Analysis                                          ":"";
// 	text +=fMT2bAnalysis? "M_{T2}b Analysis                                         ":"";
// 	TString lumi = TString::Format("%1.2f",fSamples[0].lumi/1000.);
// 	text +="CMS Preliminary, #sqrt{s} = 8 TeV, L = "+lumi+" fb^{-1}";
	TitleBox.DrawLatex(0.13,0.943,text.Data());
	TLatex LumiBox;
	LumiBox.SetNDC();
	LumiBox.SetTextSize(0.0305);
	TString lumi = TString::Format("%1.2f",fSamples[0].lumi/1000.);
	LumiBox.DrawLatex(0.68,0.943,"#sqrt{s} = 8 TeV, L = "+lumi+" fb^{-1}");//standard
	//LumiBox.DrawLatex(0.49,0.943,"CMS Preliminary, #sqrt{s} = 8 TeV, L = "+lumi+" fb^{-1}");//for CMS Preliminary
	//LumiBox.DrawLatex(0.62,0.943,"CMS, #sqrt{s} = 8 TeV, L = "+lumi+" fb^{-1}");//for CMS

 	p_plot ->Draw();
	gPad->RedrawAxis();

	if(leg != NULL ){
		leg -> SetFillColor(0);
		leg -> SetBorderSize(0);
		leg -> Draw();
	} 
	
	// draw the ratio plot
 	p_ratio ->cd();
	//gPad->SetLogy();
 
	TH1D *h_ratio = (TH1D*)h2_orig->Clone("h2_copy");	
	h_ratio ->SetStats(0);
	h_ratio ->SetMarkerStyle(20);

 	h_ratio ->Divide(h2, h1);
 	h_ratio ->SetMinimum(0.4);
 	h_ratio ->SetMaximum(3.0);
	h_ratio ->GetYaxis()->SetTitleOffset(h1->GetYaxis()->GetTitleOffset());

	//MC with errors
	TH1D*h_ratio_mc = (TH1D*)h1_orig->Clone("h1_copy");
	h_ratio_mc->Divide(h1);
	h_ratio_mc->GetYaxis()->SetRangeUser(0,2);
	h_ratio_mc->GetXaxis()->SetLabelSize( 0.);
	h_ratio_mc->GetYaxis()->SetTitle("Data / MC");
	h_ratio_mc->GetXaxis()->SetTitle(xtitle);
	h_ratio_mc->GetXaxis()->SetTitleSize(0.2);
	h_ratio_mc->GetXaxis()->SetTitleOffset(0.5);
	h_ratio_mc->GetYaxis()->SetLabelSize(0.19);
	h_ratio_mc->GetXaxis()->SetTickLength(0.09);
	h_ratio_mc	->GetYaxis()->SetTitleSize(0.18);
	h_ratio_mc->GetYaxis()->SetTitleOffset(0.36);
	h_ratio_mc->GetYaxis()->SetNdivisions(509);

	
	h_ratio_mc->SetFillStyle(3001);
	h_ratio_mc->Draw("E2");
	h_ratio ->DrawCopy("Esame");//LEO MOD
 
	TLine *l3 = new TLine(h1->GetXaxis()->GetXmin(), 1.00, h1->GetXaxis()->GetXmax(), 1.00);
	l3->SetLineWidth(2);
	l3->SetLineStyle(7);
	l3->Draw();
	
	gPad->RedrawAxis();
 	p_ratio ->Draw();
 	c1->Update();

	TString save=name+"_ratio";
	if(fSave)Util::Print(c1, save, fOutputDir, fOutputFile);
    if(saveMacro = "gif")	
		c1->SaveAs(save + ".gif");


}
//_________________________________________________________________________________
void MassPlotter::plotRatio(TH1* h1_orig, TH1* h2_orig, bool logflag, bool normalize, TString name, TLegend* leg, TString xtitle, TString ytitle){
	// define canvas and pads
	TH1D *h1 = (TH1D*)h1_orig->Clone("h1_copy");
	TH1D *h2 = (TH1D*)h2_orig->Clone("h2_copy");

	h1->SetStats(0);	
	h2->SetStats(0);	
	
	TCanvas* c1 = new TCanvas("c1","", 20,100,1000,700);
	c1 -> cd();
	
	float border = 0.3;
 	float scale = (1-border)/border;
 
 	TPad *p_plot  = new TPad("plotpad",  "Pad containing the overlay plot", 0.00, border, 1.00, 1.00, 0, 0);
 	p_plot->SetBottomMargin(0);
 	p_plot->Draw();
 	TPad *p_ratio = new TPad("ratiopad", "Pad containing the ratio",        0.00, 0.00, 1.00, border, 0, 0);
 	p_ratio->SetTopMargin(0);
 	p_ratio->SetBottomMargin(0.35);
 	p_ratio->Draw();
 
	gPad->SetFillStyle(0);
	// draw overlay plot
 	p_plot ->cd();

	if(logflag) gPad->SetLogy(1);
	gPad->SetFillStyle(0);
		
	// Scaling
	if(normalize){
		h1->Scale(1.0/h1->Integral());
		h2->Scale(1.0/h2->Integral());
	}
	
	// Determine plotting range
	double max1 = h1->GetMaximum();
	double max2 = h2->GetMaximum();
	double max  = (max1>max2)?max1:max2;
	if(logflag) max = 2*max;
	else max = 1.05*max;
	h1->SetMaximum(max);
	h2->SetMaximum(max);
//	h1->SetMinimum(0.000000001);
//	h2->SetMinimum(0.000000001);
	h2     ->Draw();
	h1     ->Draw("same");
	
	// title
	TLatex lat;
	lat.SetNDC(1);
	lat.SetTextColor(4);
	lat.SetTextSize(0.07);

	// y axis title 
	lat.SetTextAlign(33); 
	lat.SetTextColor(1);
	lat.SetTextAngle(90);
	lat.DrawLatex(0.035, 0.9, ytitle);
		
 	p_plot ->Draw();
	gPad->RedrawAxis();

	if(leg != NULL ){
		leg -> SetFillColor(0);
		leg -> SetBorderSize(0);
		leg -> Draw();
	} 
	
	// draw the ratio plot
 	p_ratio ->cd();
 
	TH1D *h_ratio = (TH1D*)h2_orig->Clone("h2_copy");	

	h_ratio->GetXaxis()->SetTitleSize(scale * h1->GetXaxis()->GetTitleSize());
	h_ratio->GetYaxis()->SetTitleSize(scale * h1->GetYaxis()->GetTitleSize());
	h_ratio->GetXaxis()->SetLabelSize(scale * h1->GetXaxis()->GetLabelSize());
	h_ratio->GetYaxis()->SetLabelSize(scale * h1->GetYaxis()->GetLabelSize());
	h_ratio->GetXaxis()->SetTickLength(scale * h1->GetXaxis()->GetTickLength());
	h_ratio->GetYaxis()->SetTickLength(h1->GetYaxis()->GetTickLength());
	h_ratio->SetLineWidth(2);
	h_ratio->SetFillColor(kBlue);
	h_ratio->SetLineColor(kBlue);
	h_ratio ->SetStats(0);
	h_ratio ->SetMarkerStyle(20);
	h_ratio ->SetMarkerSize(0.1);

 	h_ratio ->Divide(h1, h2);
 	h_ratio ->SetMinimum(0.4);
 	h_ratio ->SetMaximum(3.0);
	h_ratio ->GetYaxis()->SetTitleOffset(h1->GetYaxis()->GetTitleOffset());

	gPad->SetFillStyle(0);
	h_ratio ->DrawCopy("E2");
 
	TLine *l3 = new TLine(h1->GetXaxis()->GetXmin(), 1.00, h1->GetXaxis()->GetXmax(), 1.00);
	l3->SetLineWidth(2);
	l3->SetLineStyle(7);
	l3->Draw();
	// ratio title
	lat.SetTextAlign(33);
	lat.SetTextColor(1);
	lat.SetTextSize(0.07);
	lat.SetNDC(1);
	lat.SetTextAngle(90);
	lat.DrawLatex(0.035, 1.0, "ratio");
	
	// x axis title
	lat.SetTextAngle(0);
	float ypos = xtitle.Contains("slash") || xtitle.Contains("_") || xtitle.Contains("^") ? 0.23 : 0.16;
	lat.DrawLatex(0.9, ypos, xtitle);
	//gPad->SetLogy(1);
	gPad->RedrawAxis();
 	p_ratio ->Draw();
 	c1->Update();

	TString save=name+"_ratio";
	if(fSave)Util::PrintEPS(c1, save, fOutputDir);	

// 	delete h1;
//	delete h2;
//	delete h_ratio;
// 	delete p_plot;
// 	delete p_ratio;
// 	delete c1;

}
//___________________________________________________________________________
void MassPlotter::printHisto(TH1* h_orig, TString canvname,  Option_t *drawopt, bool logflag){
	TH1D* h = (TH1D*)h_orig->Clone("h1_copy");
	TCanvas *col = new TCanvas(canvname, "", 0, 0, 900, 700);
	col->SetFillStyle(0);
	col->SetFrameFillStyle(0);
	col->cd();
	gPad->SetFillStyle(0);
	gStyle->SetPalette(1);
	if(logflag) gPad->SetLogy(1);
	h->SetMinimum(0.01);
	gPad->SetRightMargin(0.15);
	gPad->SetLogz();
	h->SetStats(false);
	h->DrawCopy(drawopt);
	gPad->RedrawAxis();
	if(fSave)Util::PrintNoEPS(col, canvname, fOutputDir, fOutputFile);
	h->Write();
	delete col;
	delete h;

}

//____________________________________________________________________________
void MassPlotter::printHisto(THStack* h, TString canvname, Option_t *drawopt, bool logflag){
	
	TCanvas *col = new TCanvas(canvname, "", 0, 0, 900, 700);
	col->SetFillStyle(0);
	col->SetFrameFillStyle(0);
	col->cd();
	gPad->SetFillStyle(0);
	if(logflag) gPad->SetLogy(1);
	h->Draw(drawopt);
	gPad->RedrawAxis();
	col ->Update();
	gPad->RedrawAxis();
	if(fSave)Util::PrintNoEPS(col, canvname, fOutputDir, fOutputFile);
	delete col;

}
//____________________________________________________________________________
void MassPlotter::printHisto(THStack* h, TH1* h_data, TH1* h_mc_sum, TLegend* leg,  TString canvname, Option_t *drawopt, bool logflag, TString xtitle, TString ytitle,int njets, int nleps, bool stacked){

	TCanvas *col = new TCanvas(canvname, "", 0, 0, 900, 700);
	col->SetFillStyle(0);
	col->SetFrameFillStyle(0);
	col->cd();
	gPad->SetFillStyle(0);
	if(logflag) {
		gPad->SetLogy(1);
		h        -> SetMinimum(stacked ? 0.05 : 0.0001);
		h_mc_sum -> SetMinimum(0.05);
		h_data   -> SetMinimum(0.05);
	}else{
		h->SetMinimum(0);
	}

	// Determine plotting range
	double max1 = h_data->GetMaximum();
	double max2 = h_mc_sum->GetMaximum();
	double max  = (max1>max2)?max1:max2;
	if(logflag) max = 2*max;
	else max = 1.05*max;
	if(stacked){
	  h_data  ->SetMaximum(max);
	  h_mc_sum->SetMaximum(max);
	  h       ->SetMaximum(max);
	}

	h->Draw(stacked ? drawopt : "histnostack");
	//h_mc_sum -> Draw("same, E2");
	if(h_data->Integral()>0 && stacked) {
		h_data       ->Draw("same");
	}
	if(leg != NULL ){
		leg -> SetFillColor(0);
		leg -> SetBorderSize(0);
		leg -> Draw();
	} 
	gPad->RedrawAxis();
	col ->Update();

	TLatex TitleBox;
	TitleBox.SetNDC();
	TitleBox.SetTextSize(0.05);
	TString text;
	if (njets>=10)
	  text = TString::Format("%d-%d Jets",njets/10,njets%10);
	else
	  text = njets < 0 ? TString::Format("#geq %d Jets",abs(njets)) : TString::Format("%d Jets",abs(njets));
	text += nleps == 1 ? ", 1 Lepton" : "";
	TitleBox.DrawLatex(0.18,0.943,text.Data());

	// title
	TLatex lat;
	lat.SetNDC(1);
	lat.SetTextColor(4);
	lat.SetTextSize(0.07);

	// y axis title 
	lat.SetTextAlign(33); 
	lat.SetTextColor(1);
	lat.SetTextAngle(90);
	lat.DrawLatex(0.03, 0.9, ytitle);
	
	// x axis title
	lat.SetTextAngle(0);
	float ypos = xtitle.Contains("slash") || xtitle.Contains("_") || xtitle.Contains("^") ? 0.1 : 0.06;
	lat.DrawLatex(0.9, ypos, xtitle);
	gPad->RedrawAxis();
	if(fSave)Util::PrintNoEPS(col, canvname, fOutputDir, fOutputFile);
	if(fSave)Util::PrintEPS(col, canvname, fOutputDir);
// 	delete col;

}
//____________________________________________________________________________
void MassPlotter::printHisto(THStack* h, TH1* h_data, TH1* h_mc_sum, TH1* h_susy, TLegend* legend,  TString canvname, Option_t *drawopt, bool logflag, TString xtitle, TString ytitle,int njets, int nbjets, int nleps, float overlayScale){

	TCanvas *col = new TCanvas(canvname, "", 0, 0, 600, 600);
	col->SetRightMargin(0.08);
	col->SetLeftMargin(0.18);
	col->SetTopMargin(0.07);
	col->SetBottomMargin(0.17);

	TLegend *leg = (TLegend*) legend->Clone("leg");

	col->SetFillStyle(0);
	col->SetFrameFillStyle(0);
	col->cd();
	gPad->SetFillStyle(0);
	if(logflag) {
		gPad->SetLogy(1);
		h        -> SetMinimum(0.05);
		h_mc_sum -> SetMinimum(0.05);
		h_data   -> SetMinimum(0.05);
		h_susy   -> SetMinimum(0.05);
	}else{
		h->SetMinimum(0);
	}

	// Determine plotting range
	double max1 = h_data->GetMaximum();
	double max2 = h_mc_sum->GetMaximum();
	double max  = (max1>max2)?max1:max2;
	if(logflag) max = 2*max;
	else max = 1.05*max;
	h_data  ->SetMaximum(max);
	h_mc_sum->SetMaximum(max);
	h       ->SetMaximum(max);

	h->Draw(drawopt);
	//h_mc_sum -> Draw("same, E2");
	if(h_data->Integral()>0) {
		h_data       ->Draw("sameE");
	}
	h_susy->Scale(overlayScale ? overlayScale : h_data->Integral() / h_susy->Integral());
	h_susy->SetLineStyle(kDotted);
	h_susy->SetFillColor(0);
	h_susy->Draw("samehist");
	if(leg != NULL ){
		leg -> SetY1NDC(0.68);
		leg -> SetFillColor(0);
		leg -> SetBorderSize(0);
		leg -> Draw();
	} 
	gPad->RedrawAxis();
	col ->Update();

	TLatex TitleBox;
	TitleBox.SetNDC();
	TitleBox.SetTextSize(0.0305);
	//TString text = fMT2Analysis? "M_{T2} Analysis        ":"";
	//text +=fMT2bAnalysis? "M_{T2}b Analysis     ":"";
	TString text;
	if (njets>=10)
	  text = TString::Format("%d-%d jets",njets/10,njets%10);
	else
	  text = njets < 0 ? TString::Format("#geq %d jets",abs(njets)) : TString::Format("%d jets",abs(njets));
	text += nbjets==-10 ? "" : nbjets < 0 ? TString::Format(", #geq %d b-tag",abs(nbjets)) : TString::Format(", %d b-tag",abs(nbjets));
	text += nleps == 1 ? ", 1 lepton" : "";
	TString lumi = TString::Format("%1.2f",fSamples[0].lumi/1000.);
	text +="CMS Preliminary, #sqrt{s} = 8 TeV, L = "+lumi+" fb^{-1}";
	TitleBox.DrawLatex(0.18,0.943,text.Data());

// 	TLatex TitleBox;
// 	TitleBox.SetNDC();
// 	TitleBox.SetTextSize(0.0305);
// 	TString text = fMT2Analysis? "M_{T2} Analysis        ":"";
// 	text +=fMT2bAnalysis? "M_{T2}b Analysis     ":"";
// 	TString lumi = TString::Format("%1.2f",fSamples[0].lumi/1000.);
// 	text +="CMS Preliminary, #sqrt{s} = 8 TeV, L = "+lumi+" fb^{-1}";
// 	TitleBox.DrawLatex(0.18,0.943,text.Data());


	h->GetXaxis()->SetTitle(xtitle);
	h->GetXaxis()->SetLabelSize(0.05);
	h->GetXaxis()->SetTitleSize(0.05);
	h->GetXaxis()->SetTitleOffset(1.1);
	
	stringstream yTitle;
	if(fEventsPerGeV){
		if(fabs(h_data->GetBinWidth(1) -h_data->GetBinWidth(h_data->GetNbinsX()-1))<0.01){
		double binwidth = h_data->GetBinWidth(1);
		yTitle.precision(3);
		yTitle << ytitle.Data();
		yTitle << " / ";
		yTitle << binwidth;
		yTitle << " GeV";
		} else{
		cout << h_data->GetBinWidth(1) << " " << h_data->GetBinWidth(h_data->GetNbinsX()-1) << endl;
		}
	}else{
		yTitle << ytitle.Data();
	}
	h->GetYaxis()->SetTitle(yTitle.str().c_str());
	h->GetYaxis()->SetLabelSize(0.05);
	h->GetYaxis()->SetTitleSize(0.05);
	h->GetYaxis()->SetTitleOffset(1.47);

	gPad->RedrawAxis();
	if(fSave)Util::PrintNoEPS(col, canvname, fOutputDir, fOutputFile);
	if(fSave)Util::PrintEPS(col, canvname, fOutputDir);
// 	delete col;

}
//____________________________________________________________________________
void MassPlotter::printHisto(THStack* h, TH1* h_data, TH1* h_mc_sum, vector<TH1D*> h_sig, int nsig, TLegend* legend,  TString canvname, Option_t *drawopt, bool logflag, TString xtitle, TString ytitle,int njets,int nbjets, int nleps, float overlayScale){

	TCanvas *col = new TCanvas(canvname, "", 0, 0, 600, 600);
	col->SetRightMargin(0.08);
	col->SetLeftMargin(0.18);
	col->SetTopMargin(0.07);
	col->SetBottomMargin(0.17);

	TLegend *leg = (TLegend*) legend->Clone("leg");

	col->SetFillStyle(0);
	col->SetFrameFillStyle(0);
	col->cd();
	gPad->SetFillStyle(0);
	if(logflag) {
		gPad->SetLogy(1);
		h        -> SetMinimum(0.05);
		h_mc_sum -> SetMinimum(0.05);
		h_data   -> SetMinimum(0.05);
		for (int i=0; i<nsig; i++)  h_sig[i]-> SetMinimum(0.05);
	}else{
		h->SetMinimum(0);
	}

	// Determine plotting range
	double max1 = h_data->GetMaximum();
	double max2 = h_mc_sum->GetMaximum();
	double max  = (max1>max2)?max1:max2;
	if(logflag) max = 2*max;
	else max = 1.05*max;
	h_data  ->SetMaximum(max);
	h_mc_sum->SetMaximum(max);
	h       ->SetMaximum(max);

	h->Draw(drawopt);
	//h_mc_sum -> Draw("same, E2");
	if(h_data->Integral()>0) {
		h_data       ->Draw("sameE");
	}
	for (int i=0; i<nsig; i++) {
	  h_sig[i]->Scale(overlayScale ? overlayScale : h_data->Integral() / h_sig[i]->Integral());
	  h_sig[i]->SetLineStyle(2);
	  //h_sig[i]->SetLineWidth(2);
	  h_sig[i]->SetFillColor(0);
	  h_sig[i]->Draw("samehist");
	}
	if(leg != NULL ){
		leg -> SetY1NDC(0.68);
		leg -> SetFillColor(0);
		leg -> SetBorderSize(0);
		leg -> Draw();
	} 
	gPad->RedrawAxis();
	col ->Update();

	TLatex TitleBox;
	TitleBox.SetNDC();
	TitleBox.SetTextSize(0.0305);
	TString text;
	if (njets>=10)
	  text = TString::Format("%d-%d jets",njets/10,njets%10);
	else
	  text = njets < 0 ? TString::Format("#geq %d jets",abs(njets)) : TString::Format("%d jets",abs(njets));
	text += nbjets==-10 ? "" : nbjets < 0 ? TString::Format(", #geq %d b-tag",abs(nbjets)) : TString::Format(", %d b-tag",abs(nbjets));
	text += nleps == 1 ? ", 1 lepton" : "";
	TitleBox.DrawLatex(0.18,0.943,text.Data());//standard
//	TitleBox.DrawLatex(0.247,0.895,text.Data());//for CMS Preliminary
//	TitleBox.DrawLatex(0.257,0.895,text.Data());//for CMS, and 1 lepton
//	TitleBox.DrawLatex(0.18,0.943,text.Data());//for CMS
	TLatex LumiBox;
	LumiBox.SetNDC();
	LumiBox.SetTextSize(0.0305);
	TString lumi = TString::Format("%1.2f",fSamples[0].lumi/1000.);
	LumiBox.DrawLatex(0.615,0.943,"#sqrt{s} = 8 TeV, L = "+lumi+" fb^{-1}");//standard
//	LumiBox.DrawLatex(0.18,0.943,"CMS Preliminary, #sqrt{s} = 8 TeV, L = "+lumi+" fb^{-1}";//for CMS Preliminary
//	LumiBox.DrawLatex(0.18,0.943,"CMS, #sqrt{s} = 8 TeV, L = "+lumi+" fb^{-1}");//for CMS, and 1 lepton
//	LumiBox.DrawLatex(0.55,0.943,"CMS, #sqrt{s} = 8 TeV, L = "+lumi+" fb^{-1}");//for CMS


// 	TLatex TitleBox;
// 	TitleBox.SetNDC();
// 	TitleBox.SetTextSize(0.0305);
// 	TString text = fMT2Analysis? "M_{T2} Analysis        ":"";
// 	text +=fMT2bAnalysis? "M_{T2}b Analysis     ":"";
// 	TString lumi = TString::Format("%1.2f",fSamples[0].lumi/1000.);
// 	text +="CMS Preliminary, #sqrt{s} = 8 TeV, L = "+lumi+" fb^{-1}";
// 	TitleBox.DrawLatex(0.18,0.943,text.Data());

	h->GetXaxis()->SetTitle(xtitle);
	h->GetXaxis()->SetLabelSize(0.05);
	h->GetXaxis()->SetTitleSize(0.05);
	h->GetXaxis()->SetTitleOffset(1.1);
	
	stringstream yTitle;
	if(fEventsPerGeV){
		if(fabs(h_data->GetBinWidth(1) -h_data->GetBinWidth(h_data->GetNbinsX()-1))<0.01){
		double binwidth = h_data->GetBinWidth(1);
		yTitle.precision(3);
		yTitle << ytitle.Data();
		yTitle << " / ";
		yTitle << binwidth;
		yTitle << " GeV";
		} else{
		cout << h_data->GetBinWidth(1) << " " << h_data->GetBinWidth(h_data->GetNbinsX()-1) << endl;
		}
	}else{
		yTitle << ytitle.Data();
	}
	h->GetYaxis()->SetTitle(yTitle.str().c_str());
	h->GetYaxis()->SetLabelSize(0.05);
	h->GetYaxis()->SetTitleSize(0.05);
	h->GetYaxis()->SetTitleOffset(1.47);

	gPad->RedrawAxis();
	if(fSave)Util::PrintNoEPS(col, canvname, fOutputDir, fOutputFile);
	if(fSave)Util::PrintEPS(col, canvname, fOutputDir);
// 	delete col;

}
//____________________________________________________________________________
void MassPlotter::printHisto(THStack* h, TH1* h_data, TLegend* leg,  TString canvname, Option_t *drawopt, bool logflag , TString xtitle, TString ytitle){

	TCanvas *col = new TCanvas(canvname, "", 0, 0, 900, 700);
	col->SetFillStyle(0);
	col->SetFrameFillStyle(0);
	col->cd();
	gPad->SetFillStyle(0);
	if(logflag) {
		gPad->SetLogy(1);
		h->SetMinimum(0.01);
		h_data->SetMinimum(0.01);
	}else{
		h->SetMinimum(0);
	}

	// Determine plotting range
	double max1 = h_data->GetMaximum();
	double max2 = h     ->GetMaximum();
	double max  = (max1>max2)?max1:max2;
	if(logflag) max = 2*max;
	else max = 1.05*max;
	h_data  ->SetMaximum(max);
	h       ->SetMaximum(max);

	h->Draw(drawopt);
	if(h_data->Integral()>0) {
		h_data  ->Draw("same, E1");
	}
	if(leg != NULL ){
		leg -> SetFillColor(0);
		leg -> SetBorderSize(0);
		leg -> Draw();
	} 
	gPad->RedrawAxis();
	col ->Update();

	// title
	TLatex lat;
	lat.SetNDC(1);
	lat.SetTextColor(4);
	lat.SetTextSize(0.07);

	// y axis title 
	lat.SetTextAlign(33); 
	lat.SetTextColor(1);
	lat.SetTextAngle(90);
	lat.DrawLatex(0.03, 0.9, ytitle);
	
	// x axis title
	lat.SetTextAngle(0);
	float ypos = xtitle.Contains("slash") || xtitle.Contains("_") || xtitle.Contains("^") ? 0.1 : 0.06;
	lat.DrawLatex(0.9, ypos, xtitle);
	gPad->RedrawAxis();
	col ->Update();
	if(fSave) Util::PrintNoEPS(col, canvname, fOutputDir, fOutputFile);
	delete col;

}
//____________________________________________________________________________
void MassPlotter::loadSamples(const char* filename){
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
			cout <<"My path: " <<fPath << endl;
			
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

			IN.getline(buffer, 400, '\n');
			sscanf(buffer, "File\t%s", StringValue);
			TString file =fPath+StringValue;
			cout<<"my file: "<<file<<endl;
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

			if(s.type!="data"){
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
			} else{
			  s.PU_avg_weight =1;
			  s.nevents       =1;
			}

			// DON'T DO THIS AT HOME !!!!
			if ( s.name == "T1tttt_mGlu-1000_mLSP-400" ) s.nevents = 20000;
			if ( s.name.Contains("T2bb") ) s.nevents = 10000;
			if ( s.name.Contains("T2tt") ) s.nevents = 50000;
			// DON'T DO THAT AT HOME !!!!

			if ( s.type == "susy" && s.name.Contains("Tau")){
			  TH2F * h_SMSEvents = (TH2F*) s.file->Get("h_SMSEvents");
			  int binNumber = h_SMSEvents->FindBin(150.0, 150.0);//350,50//saeid
			  
			  s.nevents = h_SMSEvents->GetBinContent(binNumber); //saeid
			  s.nevents = 1000;//350,50//saeid
			  s.PU_avg_weight =1;		//saeid	  
			}//saeid
			else{

			if ( s.type == "susy" && !s.name.Contains("TStau")){//saeid
			  TH2F * h_SMSEvents = (TH2F*) s.file->Get("h_SMSEvents");
			  int binNumber = h_SMSEvents->FindBin(350.0, 50.0);//350,50//saeid
			  
			  s.nevents = h_SMSEvents->GetBinContent(binNumber); //saeid
			  //			  s.nevents = 128333;//350,50//saeid
			  s.PU_avg_weight =1;		//saeid	  
			}//saeid
			else
			  if ( s.type == "susy" && s.name.Contains("TStau")){
			    s.PU_avg_weight =1;	
			    s.nevents = 10000;
			}
			}
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

/*
//_________________________________________________________________________
void MassPlotter::TauContamination(int sample_index, Long64_t nevents, int flag){
  float sys = 0.05; //systematic uncertanty on the fake rate and efficiency
  float IsoCut = 0.5;//cut on the isolation 0.5:VLoose, 1.5:Loose...
  float maxDPhiTauMET = 11.5;//To cut on the DeltaPhi between the tau and MET

  setFlags(flag); 

  TH1::SetDefaultSumw2();
  TString  cnames[NumberOfSamples] = {"QCD", "Wjets", "Zjets", "Top", "MC", "data"};
  int      ccolor[NumberOfSamples] = {  401,     417,     419,   600,  500,  632};
  TString varname = "MT2";
  for (int i=0; i<NumberOfSamples; i++){
    MT2[i] = new TH1F(varname+"_"+cnames[i], "", 150, 0, 750);
    MT2[i] -> SetFillColor  (ccolor[i]);
    MT2[i] -> SetLineColor  (ccolor[i]);
    MT2[i] -> SetLineWidth  (2);
    MT2[i] -> SetMarkerColor(ccolor[i]);
    MT2[i] -> SetStats(false);
  }

  MT2[5] -> SetMarkerStyle(20);
  MT2[5] -> SetMarkerColor(kBlack);
  MT2[5] -> SetLineColor(kBlack);
  
  MT2[4] -> SetFillStyle(3004);
  MT2[4] -> SetFillColor(kBlack);

  TH2F* hEtaPtJet           = new TH2F("hEtaPtJet",           "hEtaPtJet",          1,0,1000.0,1,0,2.5);//tauable jets in MC in QCD dominated region (MT2 < 70)
  TH2F* hEtaPtTauFake       = new TH2F("hEtaPtTauFake",       "hEtaPtTauFake",      1,0,1000.0,1,0,2.5);//taus in MC in QCD dominated region  (MT2 < 70)
  TH2F* hEtaPtJetData       = new TH2F("hEtaPtJetData",       "hEtaPtJetData",      1,0,1000.0,1,0,2.5);//tauable jets in data in QCD dominated region (MT2 < 70)
  TH2F* hEtaPtTauFakeData   = new TH2F("hEtaPtTauFakeData",   "hEtaPtTauFakeData",  1,0,1000.0,1,0,2.5);//taus in data in QCD dominated region  (MT2 < 70)

  TH1F* hMT2TauFound        = new TH1F("hMT2TauFound",   "hMT2TauFound",   NumberOfMT2Bins - 1, xMT2bin);//Rec/Id taus
  TH1F* hMT2JetSig          = new TH1F("hMT2JetSig",     "hMT2JetSig",     NumberOfMT2Bins - 1, xMT2bin);//tauable jets in signal region (MC)
  TH1F* hMT2TauSig          = new TH1F("hMT2TauSig",     "hMT2TauSig",     NumberOfMT2Bins - 1, xMT2bin);//taus in the signal region (MC)
  TH1F* hMT2JetSigData      = new TH1F("hMT2JetSigData", "hMT2JetSigData", NumberOfMT2Bins - 1, xMT2bin);//tauable jets in signal region (data)
  TH1F* hMT2TauSigData      = new TH1F("hMT2TauSigData", "hMT2TauSigData", NumberOfMT2Bins - 1, xMT2bin);//taus in the signal region (data)
  TH1F* hRatioTight         = new TH1F("hRatioTight",       "hRatioTight",       NumberOfMT2Bins - 1, xMT2bin);//Used for the final rescaling the estimation from relaxed to appplied cut
  TH1F* hRatioRelaxed       = new TH1F("hRatioRelaxed",       "hRatioRelaxed",       NumberOfMT2Bins - 1, xMT2bin);//Used for the final rescaling the estimation from relaxed to appplied cut

  TH1F* hNAllTaus           = new TH1F("hNAllTaus",      "hNAllTaus",      NumberOfMT2Bins - 1, xMT2bin);//All of the generated hadronic taus
  TH1F* hNAllTauJets        = new TH1F("hNAllTauJets",   "hNAllTauJets",   NumberOfMT2Bins - 1, xMT2bin);//All of the generated hadronic taus matched with a rec jet


  for(unsigned int i = 0; i < fSamples.size(); i++){

    if(sample_index != -1 && i != sample_index)
      continue;

    int data = 0;
    sample Sample = fSamples[i];
    
    if(Sample.type == "data")
      data = 1;

    if(Sample.type == "susy") 
      continue;
     
    fMT2tree = new MT2tree();
    Sample.tree->SetBranchAddress("MT2tree", &fMT2tree);
    Long64_t nentries =  Sample.tree->GetEntries();

    float Weight = Sample.xsection * Sample.kfact * Sample.lumi / (Sample.nevents * Sample.PU_avg_weight );

    std::cout << setfill('=') << std::setw(70) << "" << std::endl;
    cout << "looping over :     " <<endl;	
    cout << "   Name:           " << Sample.name << endl;
    cout << "   File:           " << (Sample.file)->GetName() << endl;
    cout << "   Events:         " << Sample.nevents  << endl;
    cout << "   Events in tree: " << Sample.tree->GetEntries() << endl; 
    cout << "   Xsection:       " << Sample.xsection << endl;
    cout << "   kfactor:        " << Sample.kfact << endl;
    cout << "   avg PU weight:  " << Sample.PU_avg_weight << endl;
    cout << "   Weight:         " << Weight <<endl;
    std::cout << setfill('-') << std::setw(70) << "" << std::endl;
   
    for (Long64_t jentry=0; jentry<min(nentries, nevents);jentry++) {
      Sample.tree->GetEntry(jentry); 
      if ( fVerbose>2 && jentry % 100000 == 0 ){  cout << "+++ Proccessing event " << jentry << endl;}
 
      float weight = Weight;
      if(data == 0)
	weight *= fMT2tree->pileUp.Weight;
      else
	weight = 1.0;

      //For 2011 needs to be updated for 2012
// 	  if (data == 1){
// 	    if(!((fMT2tree->trigger.HLT_HT440_v2 ==1 && fMT2tree->misc.Run<161216)                                                               ||
// 		 (fMT2tree->trigger.HLT_HT450_v2 ==1 && fMT2tree->misc.Run>=161216 && fMT2tree->misc.Run< 163269)                                ||
// 		 (fMT2tree->trigger.HLT_HT500_v3 ==1 && fMT2tree->misc.Run>=163269 && fMT2tree->misc.Run<=163869)                                ||
// 		 (fMT2tree->trigger.HLT_HT500_v4 ==1 && fMT2tree->misc.Run>=165088 && fMT2tree->misc.Run< 165970)                                ||
// 		 (fMT2tree->trigger.HLT_HT550_v5 ==1 && fMT2tree->misc.Run>=165970 && fMT2tree->misc.Run<=167043 && fMT2tree->misc.Run!=166346)  ||
// 		 (fMT2tree->trigger.HLT_HT550_v6 ==1 && fMT2tree->misc.Run==166346)                                                              ||
// 		 (fMT2tree->trigger.HLT_HT550_v7 ==1 && fMT2tree->misc.Run>=167078 && fMT2tree->misc.Run< 170249)                                ||
// 		 (fMT2tree->trigger.HLT_HT550_v8 ==1 && fMT2tree->misc.Run>=170249 && fMT2tree->misc.Run< 173236)                                ||
// 		 (fMT2tree->trigger.HLT_HT600_v1 ==1 && fMT2tree->misc.Run>=173236 && fMT2tree->misc.Run< 178420)                                || 
// 		 (fMT2tree->trigger.HLT_HT650_v4 ==1 && fMT2tree->misc.Run>=178420 && fMT2tree->misc.Run< 179959)))
// 	      continue;
// 	  }


      if (!(fMT2tree->misc.MET             >  30                                   &&
	    fMT2tree->misc.HT              >  750                                  &&
	    (fMT2tree->NEles==0  ||  fMT2tree->ele[0].lv.Pt()<10)                  &&
	    (fMT2tree->NMuons==0 ||  fMT2tree->muo[0].lv.Pt()<10)                  &&
	    fMT2tree->misc.Jet0Pass        == 1                                    &&
	    fMT2tree->misc.Jet1Pass        == 1                                    &&
	    fMT2tree->misc.SecondJPt       >  100                                  &&
	    fMT2tree->misc.PassJetID       == 1                                    &&
	    fMT2tree->misc.Vectorsumpt     <  70                                   &&
//            fMT2tree->misc.HBHENoiseFlagIso == 0                                   &&
//	    fMT2tree->misc.CSCTightHaloID   == 0                                   &&
            fMT2tree->misc.CrazyHCAL        == 0))
	  continue;

 	  if(High && fMT2tree->misc.HT             <  950)
 	    continue;

 	  if(Low && fMT2tree->misc.HT              >  950)
	    continue;

	  if(mT2){
	    if(!( fMT2tree->misc.LeadingJPt      >  100                                  &&
//  		  fMT2tree->misc.MinMetJetDPhi   >  0.3                                  &&
		  fMT2tree->NJetsIDLoose40       >= 3 
		  ))
	      continue;
	  }else{
	    if(!( fMT2tree->NJetsIDLoose40       >= 4                                    &&
// 		  fMT2tree->NBJets               >  0                                    &&
//  		  fMT2tree->NBJetsHE             >  0                                    &&
//  		  fMT2tree->misc.MinMetJetDPhi4  >  0.3                                  &&
		  fMT2tree->misc.LeadingJPt      >  150                                  
 		  ))
 	    continue;
	    
	    if(!((BTagRelaxed == 0 && fMT2tree->NBJets40CSVM > 0) || (BTagRelaxed == 1 && fMT2tree->NBJets40CSVM > 0)))
	      continue;	
	  }

	
	   
	  float METPhi = fMT2tree->misc.METPhi; //To be used if a cone close to met is used to look for taus
	  
   for(int t = 0; t < fMT2tree->NTaus; t++){
 	if(fMT2tree->tau[t].MT > 100 || fMT2tree->tau[t].Isolation < IsoCut)
 	  continue;

	float minDeltaR = 100.0;

 	float DeltaPhi = Util::DeltaPhi(fMT2tree->tau[t].lv.Phi(), fMT2tree->misc.METPhi);

	if(DeltaPhi > maxDPhiTauMET)
	  continue;

	  if(data == 1){
	    
	    MT2[5]->Fill(fMT2tree->misc.MT2, weight);//data
	    
	  }else{
	    MT2[4]->Fill(fMT2tree->misc.MT2, weight);
	    
	    if(Sample.sname == "Top")  
	      MT2[3]->Fill(fMT2tree->misc.MT2, weight);
	    else
	      if(Sample.sname == "DY")	
		MT2[2]->Fill(fMT2tree->misc.MT2, weight);
	      else
		if(Sample.name == "WJetsToLNu")
		  MT2[1]->Fill(fMT2tree->misc.MT2, weight);
		else
		  if(Sample.sname == "QCD")
		    MT2[0]->Fill(fMT2tree->misc.MT2, weight);
	  }}

   for(int j = 0; j < fMT2tree->NGenLepts; j++){
     if(abs(fMT2tree->genlept[j].ID) != 15)
       continue;
	
	    int HadTau = 1;
	    
	    for(int iii = 0; iii < fMT2tree->NGenLepts; iii++){
	      if((abs(fMT2tree->genlept[iii].ID) == 11 && fMT2tree->genlept[iii].MID == fMT2tree->genlept[j].ID) || 
		 (abs(fMT2tree->genlept[iii].ID) == 13 && fMT2tree->genlept[iii].MID == fMT2tree->genlept[j].ID)){
		HadTau = 0;
	      }
	    }

	    if(HadTau == 0)
	      continue;
        
	if((mT2 && fMT2tree->misc.MinMetJetDPhi   >  0.3) || (mT2b && fMT2tree->NBJets40CSVM >  0) || (mT2b && fMT2tree->misc.MinMetJetDPhi4   >  0.3 && fMT2tree->NBJets40CSVM >  0))
	  hRatioTight->Fill(fMT2tree->misc.MT2, weight);
	else
	  hRatioRelaxed->Fill(fMT2tree->misc.MT2, weight);	  
	
	hNAllTaus->Fill(fMT2tree->misc.MT2, weight);
	

	//To find the jet matched to this gen had tau
	int jetIndex = -1;
	float minDR = 100.0;

	float genleptEta = fMT2tree->genlept[j].lv.Eta();

	float genleptPhi = fMT2tree->genlept[j].lv.Phi();

	for(int k = 0; k < fMT2tree->NJets; k++){
	  float deltaR = Util::GetDeltaR(genleptEta, fMT2tree->jet[k].lv.Eta(), genleptPhi, fMT2tree->jet[k].lv.Phi());
	  if(deltaR < minDR){
	    minDR = deltaR;
	    jetIndex = k;
	  }
	}	

	if(minDR < 0.1){
	  hNAllTauJets->Fill(fMT2tree->misc.MT2, weight);
	}

	jetIndex = -1;
	minDR = 100.0;

	//To find the rec/id tau matched to this gen had tau
	for(int k = 0; k < fMT2tree->NTaus; k++){
	  float deltaR = Util::GetDeltaR(genleptEta, fMT2tree->tau[k].lv.Eta(), genleptPhi, fMT2tree->tau[k].lv.Phi());
	  if(deltaR < minDR){
	    minDR = deltaR;
	    jetIndex = k;
	  }
	}

	float DeltaPhi = Util::DeltaPhi(fMT2tree->tau[jetIndex].lv.Phi(), METPhi);

	if(DeltaPhi > maxDPhiTauMET)
	  continue;
	
	if(minDR < 0.1){
	  if(fMT2tree->tau[jetIndex].Isolation > IsoCut && fMT2tree->tau[jetIndex].MT < 100)
	    hMT2TauFound->Fill(fMT2tree->misc.MT2, weight);
	}
      }
      
      //QCD dominated region, use these events to find the fake rate
      if(fMT2tree->misc.MT2 < 70){
	if(data == 1){
	  for(int j = 0; j < fMT2tree->NJets; j++){
	    
	    if(fMT2tree->jet[j].isTauMatch < 0 )// || fMT2tree->tau[fMT2tree->jet[j].isTauMatch].MT > 100)
	      continue;

	    float DeltaPhi = Util::DeltaPhi(fMT2tree->jet[j].lv.Phi(), METPhi);

	    if(DeltaPhi > maxDPhiTauMET)
	      continue;
	
	    hEtaPtJetData->Fill(fMT2tree->jet[j].lv.Pt(), fabs(fMT2tree->jet[j].lv.Eta()), weight);

	    if(fMT2tree->tau[fMT2tree->jet[j].isTauMatch].Isolation > IsoCut && fMT2tree->tau[fMT2tree->jet[j].isTauMatch].MT < 100){  
	      hEtaPtTauFakeData->Fill(fMT2tree->jet[j].lv.Pt(), fabs(fMT2tree->jet[j].lv.Eta()), weight);
	    }
	  }
	}else{
	  for(int j = 0; j < fMT2tree->NJets; j++){

	    if(fMT2tree->jet[j].isTauMatch < 0 )// || fMT2tree->tau[fMT2tree->jet[j].isTauMatch].MT > 100)
	      continue;	    

	    float DeltaPhi = Util::DeltaPhi(fMT2tree->jet[j].lv.Phi(), METPhi);

	    if(DeltaPhi > maxDPhiTauMET)
	      continue;

	    hEtaPtJet->Fill(fMT2tree->jet[j].lv.Pt(), fabs(fMT2tree->jet[j].lv.Eta()), weight);
	    
	    if(fMT2tree->tau[fMT2tree->jet[j].isTauMatch].Isolation > IsoCut && fMT2tree->tau[fMT2tree->jet[j].isTauMatch].MT < 100){ 
	      
	      hEtaPtTauFake->Fill(fMT2tree->jet[j].lv.Pt(), fabs(fMT2tree->jet[j].lv.Eta()), weight);
	    }
	  }
	}
      }
      //Fill The number of Rec/Id taus and tauable jets in the signal region
      if(fMT2tree->misc.MT2 >= MT2Min){
	if(data == 1){
	  for(int j = 0; j < fMT2tree->NJets; j++){
	    if(fMT2tree->jet[j].isTauMatch < 0)// || fMT2tree->tau[fMT2tree->jet[j].isTauMatch].MT > 100)
	      continue;

	    float DeltaPhi = Util::DeltaPhi(fMT2tree->jet[j].lv.Phi(), METPhi);

	    if(DeltaPhi > maxDPhiTauMET)
	      continue;
	    hMT2JetSigData->Fill(fMT2tree->misc.MT2, weight);
	  }
	
	  for(int j = 0; j < fMT2tree->NTaus; j++){

	    float DeltaPhi = Util::DeltaPhi(fMT2tree->tau[j].lv.Phi(), METPhi);

	    if(DeltaPhi > maxDPhiTauMET)
	      continue;
	    
	    if(fMT2tree->tau[j].MT < 100 && fMT2tree->tau[j].Isolation > IsoCut){

	      hMT2TauSigData->Fill(fMT2tree->misc.MT2, weight);
	    }
	  }
	}else{
	  
	  for(int j = 0; j < fMT2tree->NJets; j++){
	    if(fMT2tree->jet[j].isTauMatch < 0)// || fMT2tree->tau[fMT2tree->jet[j].isTauMatch].MT > 100) 
	      continue;

	    float DeltaPhi = Util::DeltaPhi(fMT2tree->jet[j].lv.Phi(), METPhi);

	    if(DeltaPhi > maxDPhiTauMET)
	      continue;
	
	    hMT2JetSig->Fill(fMT2tree->misc.MT2, weight);
	  }
	
	  for(int j = 0; j < fMT2tree->NTaus; j++){

	    float DeltaPhi = Util::DeltaPhi(fMT2tree->tau[j].lv.Phi(), METPhi);

	    if(DeltaPhi > maxDPhiTauMET)
	      continue;
	
	    if(fMT2tree->tau[j].MT < 100 && fMT2tree->tau[j].Isolation > IsoCut){
	      hMT2TauSig->Fill(fMT2tree->misc.MT2, weight);

	      float minDR = 100.0;
	      
	      float tauEta = fMT2tree->tau[j].lv.Eta();

	      float tauPhi = fMT2tree->tau[j].lv.Phi();
	    }
	  }
	}
      }
    }//for (Long64_t jentry=0; jentry<min(nentries, nevents);jentry++) 
  }// for(int i = 0; i < fSamples.size(); i++){


 for(int j = 1; j < NumberOfMT2Bins; j++){
   float id = hMT2TauFound  ->GetBinContent(j);
   float jet= hNAllTauJets->GetBinContent(j);
   float gen= hNAllTaus   ->GetBinContent(j);
   float EffTotal = -1.0;
   if(gen != 0) EffTotal = id/gen;
   float Eff = -1.0;
   if(jet != 0) Eff = id/jet;
   float Acc = -1.0;
   if(gen != 0) Acc = jet/gen;

   cout<<"$"<<xMT2bin[j-1]<<"-"<<xMT2bin[j]<<" Eff = "<<EffTotal<<", Eff * Acc = "<<Eff<<" * "<<Acc<<endl;
 }

 int p = 0;
 hMT2TauFound->Divide(hNAllTaus);
 Cout(p++, hMT2TauFound);

 //To add the syatematic uncertainty
 for(int j = 1; j < NumberOfMT2Bins; j++){
   float val = hMT2TauFound->GetBinContent(j);
   float Err = hMT2TauFound->GetBinError(j);
   hMT2TauFound->SetBinError(j, sqrt(Err * Err + sys * sys * val * val));
 }

 Cout(p++, hEtaPtTauFake);
 Cout(p++, hEtaPtJet);
 float TauFake  = hEtaPtTauFake->Integral();
  float Jet      = hEtaPtJet    ->Integral();
  float Fake     = TauFake/Jet;
  cout<<" Fake "<<Fake<<endl;
  hEtaPtTauFake->Divide(hEtaPtJet);
  Cout(p++, hEtaPtTauFake);

  Cout(p++, hMT2JetSig);

  hMT2JetSig->Scale(Fake);

  Cout(p++, hMT2JetSig);
  Cout(p++, hMT2TauSig);
  hMT2JetSig->Add(hMT2TauSig, -1);
  hMT2JetSig->Scale(-1.0);
  Cout(p++, hMT2JetSig);

  hMT2JetSig->Divide(hMT2TauFound);
  Cout(p++, hMT2JetSig);

  float TauFakeData  = hEtaPtTauFakeData->Integral();
  float JetData      = hEtaPtJetData    ->Integral();
  float FakeData     = TauFakeData/JetData;
  cout<<" FakeData "<<FakeData<<endl;
  hEtaPtTauFakeData->Divide(hEtaPtJetData);
  
  Cout(p++, hMT2JetSigData);
  hMT2JetSigData->Scale(FakeData);
  
  Cout(p++, hMT2JetSigData);
  Cout(p++, hMT2TauSigData);
  hMT2JetSigData->Add(hMT2TauSigData, -1);
  hMT2JetSigData->Scale(-1.0);
  Cout(p++, hMT2JetSigData);
    

  hMT2JetSigData->Divide(hMT2TauFound);
  
  Cout(p++, hMT2JetSigData);
  
  int k = 1;

  int Chi2  = 0;

  float Sim = 0;

  float Pre = 0;

  float PreData    = 0;

  float SimErr     = 0;

  float PreErr     = 0;

  float PreErrData = 0;

  Cout(p++, hRatioRelaxed);
  Cout(p++, hRatioTight);

  hRatioRelaxed->Add(hRatioTight);
  
  hRatioTight->Divide(hRatioRelaxed);

  Cout(p++, hRatioTight);

  hMT2JetSig->Multiply(hRatioTight);

  hMT2JetSigData->Multiply(hRatioTight);

  hNAllTaus->Multiply(hRatioTight);
 
  for(int j = 1; j < NumberOfMT2Bins; j++){
    
    float TotalErrSim = hMT2JetSig    ->GetBinError(j);
    float StatErrSim  = (hMT2TauSig    ->GetBinError(j)/hMT2TauFound->GetBinContent(j))*hRatioTight->GetBinContent(j);
    float SysErrSim   = sqrt(fabs(TotalErrSim * TotalErrSim - StatErrSim * StatErrSim));

    float TotalErrData = hMT2JetSigData    ->GetBinError(j);
    float StatErrData  = (hMT2TauSigData    ->GetBinError(j)/hMT2TauFound->GetBinContent(j))*hRatioTight->GetBinContent(j);
    float SysErrData   = sqrt(fabs(TotalErrData * TotalErrData - StatErrData * StatErrData));
    
    cout<<"$"<<xMT2bin[j-1]<<"-"<<xMT2bin[j]<<"$: Sim  "<<hNAllTaus     ->GetBinContent(j)<<" +/- "<<hNAllTaus        ->GetBinError(j)*hRatioTight->GetBinContent(j)<<
                                                 " Pre  "<<hMT2JetSig    ->GetBinContent(j)<<" +/- "<<StatErrSim <<" +/- "<< SysErrSim <<                              
                                                 " Data "<<hMT2JetSigData->GetBinContent(j)<<" +/- "<<StatErrData<<" +/- "<< SysErrData<<endl;
  }

  for(int j = 0; j < NumberOfSamples; j++){
     AddOverAndUnderFlow(MT2[j], true, true);
  }
  printYield();

THStack* h_stack     = new THStack(varname, "");
  for(int j = 0; j < NumberOfSamples; j++){
    MT2[j]->Rebin(3);
    TH1F* mt2 = (TH1F*)MT2[j]->Clone();
    mt2->SetName("mt2");
    if(j < (NumberOfSamples - 2))
      h_stack  -> Add(MT2[j]);
    delete mt2;
  }  

  TLegend* Legend1 = new TLegend(.71,.54,.91,.92);
  Legend1->AddEntry(MT2[0], "QCD", "f");
  Legend1->AddEntry(MT2[1], "W+jets", "f");
  Legend1->AddEntry(MT2[2], "Z+jets", "f");
  Legend1->AddEntry(MT2[3], "Top", "f");
  Legend1->AddEntry(MT2[5], "data", "l");
  
  //  TLegend *Legend1;
  printHisto(h_stack, MT2[5], MT2[4], Legend1 , "MTC", "hist", true, "MT2", "Events", -4, 0);

  plotRatioStack(h_stack,  MT2[4], MT2[5], true, false, "MT2_ratio", Legend1, "MT2", "Events", -4, 0);

  TCanvas *CanvasMT2 = new TCanvas("myMT2", "myMT2");
  hNAllTaus      ->SetLineColor(kBlue);
  hMT2JetSig     ->SetLineColor(kRed);
  hMT2JetSigData ->SetLineColor(kGreen);
  hNAllTaus      ->Draw();
  hMT2JetSig     ->Draw("sames");
  hMT2JetSigData ->Draw("sames");
}*/


void MassPlotter::TauContamination(int sample_index, Long64_t nevents, int flag, TString cuts, TString trigger, int signalContamination){
  float sys = 0.05; //systematic uncertanty on the fake rate and efficiency
  float IsoCut = 1.5;//cut on the isolation 0.5:VLoose, 1.5:Loose...
  float IsoCutLoose = IsoCut - 1.0;
  float maxDPhiTauMET = 11.5;//To cut on the DeltaPhi between the tau and MET

  setFlags(flag); 
  TH1::SetDefaultSumw2();

  TString  cnames[NumberOfSamples] = {"QCD", "Wjets", "Zjets", "Top", "MC", "data"};
  int      ccolor[NumberOfSamples] = {  401,     417,     419,   600,  500,  632};
  TString varname = "MT2";
  for (int i=0; i<NumberOfSamples; i++){
    MT2[i] = new TH1F(varname+"_"+cnames[i], "", 150, 0, 750);
    MT2[i] -> SetFillColor  (ccolor[i]);
    MT2[i] -> SetLineColor  (ccolor[i]);
    MT2[i] -> SetLineWidth  (2);
    MT2[i] -> SetMarkerColor(ccolor[i]);
    MT2[i] -> SetStats(false);
  }

  MT2[5] -> SetMarkerStyle(20);
  MT2[5] -> SetMarkerColor(kBlack);
  MT2[5] -> SetLineColor(kBlack);
  
  MT2[4] -> SetFillStyle(3004);
  MT2[4] -> SetFillColor(kBlack);

  TH2F* hEtaPtJet           = new TH2F("hEtaPtJet",           "hEtaPtJet",          1,0,1000.0,1,0,2.5);//tauable jets in MC in QCD dominated region (MT2 < 70)
  TH2F* hEtaPtTauFake       = new TH2F("hEtaPtTauFake",       "hEtaPtTauFake",      1,0,1000.0,1,0,2.5);//taus in MC in QCD dominated region  (MT2 < 70)
  TH2F* hEtaPtJetData       = new TH2F("hEtaPtJetData",       "hEtaPtJetData",      1,0,1000.0,1,0,2.5);//tauable jets in data in QCD dominated region (MT2 < 70)
  TH2F* hEtaPtTauFakeData   = new TH2F("hEtaPtTauFakeData",   "hEtaPtTauFakeData",  1,0,1000.0,1,0,2.5);//taus in data in QCD dominated region  (MT2 < 70)

  TH1F* hMT2TauFound        = new TH1F("hMT2TauFound",   "hMT2TauFound",   NumberOfMT2Bins - 1, xMT2bin);//Rec/Id taus
  TH1F* hMT2JetSig          = new TH1F("hMT2JetSig",     "hMT2JetSig",     NumberOfMT2Bins - 1, xMT2bin);//tauable jets in signal region (MC)
  TH1F* hMT2TauSig          = new TH1F("hMT2TauSig",     "hMT2TauSig",     NumberOfMT2Bins - 1, xMT2bin);//taus in the signal region (MC)
  TH1F* hMT2JetSigData      = new TH1F("hMT2JetSigData", "hMT2JetSigData", NumberOfMT2Bins - 1, xMT2bin);//tauable jets in signal region (data)
  TH1F* hMT2TauSigData      = new TH1F("hMT2TauSigData", "hMT2TauSigData", NumberOfMT2Bins - 1, xMT2bin);//taus in the signal region (data)
  TH1F* hRatioTight         = new TH1F("hRatioTight",       "hRatioTight",       NumberOfMT2Bins - 1, xMT2bin);//Used for the final rescaling the estimation from relaxed to appplied cut
  TH1F* hRatioRelaxed       = new TH1F("hRatioRelaxed",       "hRatioRelaxed",       NumberOfMT2Bins - 1, xMT2bin);//Used for the final rescaling the estimation from relaxed to appplied cut

  TH1F* hNAllTaus           = new TH1F("hNAllTaus",      "hNAllTaus",      NumberOfMT2Bins - 1, xMT2bin);//All of the generated hadronic taus
  TH1F* hNAllTauJets        = new TH1F("hNAllTauJets",   "hNAllTauJets",   NumberOfMT2Bins - 1, xMT2bin);//All of the generated hadronic taus matched with a rec jet

  cout<<" trigger "<<trigger<<endl;
  cout<<" cuts "<<cuts<<endl;

  for(int ii = 0; ii < fSamples.size(); ii++){
    float nAllTausBeforeThisSample = hNAllTaus->Integral();
    if(sample_index != -1 && ii != sample_index)
      continue;
    
    TString myCuts = cuts;
 
    int data = 0;
    sample Sample = fSamples[ii];
    
    if(Sample.type == "data"){
      data = 1;
      myCuts += " && " + trigger;
    }

    if(Sample.type == "susy" && signalContamination == 0)// || Sample.sname == "QCD") 
      continue;
     
    fMT2tree = new MT2tree();
    Sample.tree->SetBranchAddress("MT2tree", &fMT2tree);

    float Weight = Sample.xsection * Sample.kfact * Sample.lumi / (Sample.nevents);

    std::cout << setfill('=') << std::setw(70) << "" << std::endl;
    cout << "looping over :     " <<endl;	
    cout << "   Name:           " << Sample.name << endl;
    cout << "   File:           " << (Sample.file)->GetName() << endl;
    cout << "   Events:         " << Sample.nevents  << endl;
    cout << "   Events in tree: " << Sample.tree->GetEntries() << endl; 
    cout << "   Xsection:       " << Sample.xsection << endl;
    cout << "   kfactor:        " << Sample.kfact << endl;
    cout << "   avg PU weight:  " << Sample.PU_avg_weight << endl;
    cout << "   Weight:         " << Weight <<endl;
    std::cout << setfill('-') << std::setw(70) << "" << std::endl;
   

    Sample.tree->Draw(">>selList", myCuts);
    TEventList *myEvtList = (TEventList*)gDirectory->Get("selList");
    Sample.tree->SetEventList(myEvtList);

    Long64_t nentries =  myEvtList->GetN();//Sample.tree->GetEntries();


    for (Long64_t jentry=0; jentry<min(nentries, nevents);jentry++) {
      //Sample.tree->GetEntry(jentry); 
      Sample.tree->GetEntry(myEvtList->GetEntry(jentry));

      if ( fVerbose>2 && jentry % 100000 == 0 ){  
	fprintf(stdout, "\rProcessed events: %6d of %6d ", jentry + 1, nentries);
	fflush(stdout);
      }
 
      float weight = Weight;
      if(data == 1)
 	weight = 1.0;
      else
	if(Sample.type != "susy")
	  weight *= (fMT2tree->pileUp.Weight * fMT2tree->SFWeight.BTagCSV40ge1 * fMT2tree->SFWeight.TauTagge1/Sample.PU_avg_weight);//

      int Filled = 0;

      float METPhi = fMT2tree->misc.METPhi; //To be used if a cone close to met is used to look for taus
      if(Sample.type != "susy"){
      for(int t = 0; t < fMT2tree->NTaus; t++){

 	if(fMT2tree->tau[t].MT > 100 || fMT2tree->tau[t].Isolation3Hits < IsoCut || Filled == 1)
 	  continue;

	Filled = 1;

	float minDeltaR = 100.0;

 	float DeltaPhi = Util::DeltaPhi(fMT2tree->tau[t].lv.Phi(), fMT2tree->misc.METPhi);

	if(DeltaPhi > maxDPhiTauMET)
	  continue;

	if(data == 1){
	    
	  MT2[5]->Fill(fMT2tree->misc.MT2, weight);//data
	    
	}else{
	  MT2[4]->Fill(fMT2tree->misc.MT2, weight);
	    
	  if(Sample.sname == "Top")  
	    MT2[3]->Fill(fMT2tree->misc.MT2, weight);
	  else
	    if(Sample.sname == "DY")	
	      MT2[2]->Fill(fMT2tree->misc.MT2, weight);
	    else
	      if(Sample.sname == "Wtolnu")
		MT2[1]->Fill(fMT2tree->misc.MT2, weight);
	      else
		if(Sample.sname == "QCD")
		  MT2[0]->Fill(fMT2tree->misc.MT2, weight);
	  //	  if(fMT2tree->misc.MT2 >= MT2Min)
	  //cout<<"Tau:: "<<fMT2tree->misc.MT2<<endl;
	}}}

      Filled = 0;
      int GenTau = 0;
      if(Sample.type != "susy"){
       if(fMT2tree->misc.TopDecayMode != 48){//skip di_tau_had  events

      for(int j = 0; j < fMT2tree->NGenLepts; j++){
	if(abs(fMT2tree->genlept[j].ID) != 15 || abs(fMT2tree->genlept[j].MID) != 24)
	  continue;
		

	int HadTau = 1;
	    
	for(int i = 0; i < fMT2tree->NGenLepts; i++){
	  if((abs(fMT2tree->genlept[i].ID) == 11 && fMT2tree->genlept[i].MID == fMT2tree->genlept[j].ID && abs(fMT2tree->genlept[i].GMID) == 24) || 
	     (abs(fMT2tree->genlept[i].ID) == 13 && fMT2tree->genlept[i].MID == fMT2tree->genlept[j].ID && abs(fMT2tree->genlept[i].GMID) == 24)){
	    HadTau = 0;
	  }
	}

 	if(HadTau == 0 || Filled != 0){
// 	  if(Sample.sname == "Top") cout<<"ele "<<fMT2tree->TopDecayModeResult(11)<<endl;
// 	  if(Sample.sname == "Top") cout<<"mu "<<fMT2tree->TopDecayModeResult(13)<<endl;
 	  continue;
	}

	GenTau++;

	Filled = 1;
	if((mT2 && fMT2tree->misc.MinMetJetDPhi   >  0.3) || (mT2b && fMT2tree->misc.MinMetJetDPhi4   >  0.3))
	  hRatioTight->Fill(fMT2tree->misc.MT2, weight);
	else
	  hRatioRelaxed->Fill(fMT2tree->misc.MT2, weight);	  
	  
	
	hNAllTaus->Fill(fMT2tree->misc.MT2, weight);
	
	//To find the jet matched to this gen had tau
	int jetIndex = -1;
	float minDR = 100.0;

	float genleptEta = fMT2tree->genlept[j].lv.Eta();

	float genleptPhi = fMT2tree->genlept[j].lv.Phi();

	for(int k = 0; k < fMT2tree->NJets; k++){
	  float deltaR = Util::GetDeltaR(genleptEta, fMT2tree->jet[k].lv.Eta(), genleptPhi, fMT2tree->jet[k].lv.Phi());
	  if(deltaR < minDR){
	    minDR = deltaR;
	    jetIndex = k;
	  }
	}	

	if(minDR < 0.1){
	  hNAllTauJets->Fill(fMT2tree->misc.MT2, weight);
	}

	jetIndex = -1;
	minDR = 100.0;

	//To find the rec/id tau matched to this gen had tau
	for(int k = 0; k < fMT2tree->NTaus; k++){
	  float deltaR = Util::GetDeltaR(genleptEta, fMT2tree->tau[k].lv.Eta(), genleptPhi, fMT2tree->tau[k].lv.Phi());
	  if(deltaR < minDR){
	    minDR = deltaR;
	    jetIndex = k;
	  }
	}

	float DeltaPhi = Util::DeltaPhi(fMT2tree->tau[jetIndex].lv.Phi(), METPhi);

	if(DeltaPhi > maxDPhiTauMET)
	  continue;
	
	if(minDR < 0.1){
	  if(fMT2tree->tau[jetIndex].Isolation3Hits > IsoCut && fMT2tree->tau[jetIndex].MT < 100)
	    hMT2TauFound->Fill(fMT2tree->misc.MT2, weight);
	}
      }}//TopDeayMode
      Filled = 0;
      if(GenTau == 0 && data != 1){
	for(int j = 0; j < fMT2tree->NTaus; j++){
	  if(fMT2tree->tau[j].Isolation3Hits > IsoCut && fMT2tree->tau[j].MT < 100 && Filled == 0){
	    Filled = 1;
	    hMT2JetSig->Fill(fMT2tree->misc.MT2, weight);
	    hMT2JetSigData->Fill(fMT2tree->misc.MT2, weight);
	  }}}
      }//!=susy


      /*
      //QCD dominated region, use these events to find the fake rate
      if(fMT2tree->misc.MT2 < 70){
	if(data == 1){
	  for(int j = 0; j < fMT2tree->NJets; j++){
	    
	    if(fMT2tree->jet[j].isTauMatch < 0  || fMT2tree->tau[fMT2tree->jet[j].isTauMatch].MT > 100)
	      continue;

	    float DeltaPhi = Util::DeltaPhi(fMT2tree->jet[j].lv.Phi(), METPhi);

	    if(DeltaPhi > maxDPhiTauMET)
	      continue;
	    if(fMT2tree->tau[fMT2tree->jet[j].isTauMatch].Isolation3Hits > IsoCutLoose)
	      hEtaPtJetData->Fill(fMT2tree->jet[j].lv.Pt(), fabs(fMT2tree->jet[j].lv.Eta()), weight);

	    if(fMT2tree->tau[fMT2tree->jet[j].isTauMatch].Isolation3Hits > IsoCut && fMT2tree->tau[fMT2tree->jet[j].isTauMatch].MT < 100){  
	      hEtaPtTauFakeData->Fill(fMT2tree->jet[j].lv.Pt(), fabs(fMT2tree->jet[j].lv.Eta()), weight);
	    }
	  }
	}else{
	  for(int j = 0; j < fMT2tree->NJets; j++){

	    if(fMT2tree->jet[j].isTauMatch < 0  || fMT2tree->tau[fMT2tree->jet[j].isTauMatch].MT > 100)
	      continue;	    

	    float DeltaPhi = Util::DeltaPhi(fMT2tree->jet[j].lv.Phi(), METPhi);

	    if(DeltaPhi > maxDPhiTauMET)
	      continue;

	    if(fMT2tree->tau[fMT2tree->jet[j].isTauMatch].Isolation3Hits > IsoCutLoose)
	      hEtaPtJet->Fill(fMT2tree->jet[j].lv.Pt(), fabs(fMT2tree->jet[j].lv.Eta()), weight);
	    
	    if(fMT2tree->tau[fMT2tree->jet[j].isTauMatch].Isolation3Hits > IsoCut && fMT2tree->tau[fMT2tree->jet[j].isTauMatch].MT < 100){ 
	      
	      hEtaPtTauFake->Fill(fMT2tree->jet[j].lv.Pt(), fabs(fMT2tree->jet[j].lv.Eta()), weight);
	    }
	  }
	}
      }*/
      
      Filled = 0;

      //Fill The number of Rec/Id taus and tauable jets in the signal region
      if(fMT2tree->misc.MT2 >= MT2Min){
	if(data == 1){
	  for(int j = 0; j < fMT2tree->NJets; j++){
	    if(fMT2tree->jet[j].isTauMatch < 0 || fMT2tree->tau[fMT2tree->jet[j].isTauMatch].MT > 100)
	      continue;

	    float DeltaPhi = Util::DeltaPhi(fMT2tree->jet[j].lv.Phi(), METPhi);

	    if(DeltaPhi > maxDPhiTauMET)
	      continue;
 
	    
// 	    if(fMT2tree->tau[fMT2tree->jet[j].isTauMatch].Isolation3Hits > IsoCutLoose){
// 	      hMT2JetSigData->Fill(fMT2tree->misc.MT2, weight);
// 	      cout<<"JetSigData:: "<<fMT2tree->misc.MT2<<endl;
// 	    }
	  }
	
	  for(int j = 0; j < fMT2tree->NTaus; j++){

	    float DeltaPhi = Util::DeltaPhi(fMT2tree->tau[j].lv.Phi(), METPhi);

	    if(DeltaPhi > maxDPhiTauMET)
	      continue;
	    
	    if(fMT2tree->tau[j].MT < 100 && fMT2tree->tau[j].Isolation3Hits > IsoCut &&  Filled == 0){
	      Filled = 1;
	      hMT2TauSigData->Fill(fMT2tree->misc.MT2, weight);
	      //	      cout<<"TauSigData:: "<<fMT2tree->misc.MT2<<endl;
	    }
	  }
	}else{
	  
	  for(int j = 0; j < fMT2tree->NJets; j++){
	    if(fMT2tree->jet[j].isTauMatch < 0 || fMT2tree->tau[fMT2tree->jet[j].isTauMatch].MT > 100) 
	      continue;

	    float DeltaPhi = Util::DeltaPhi(fMT2tree->jet[j].lv.Phi(), METPhi);

	    if(DeltaPhi > maxDPhiTauMET)
	      continue;
// 	    if(fMT2tree->tau[fMT2tree->jet[j].isTauMatch].Isolation3Hits > IsoCutLoose){
// 	      hMT2JetSig->Fill(fMT2tree->misc.MT2, weight);
// 	      cout<<"JetSig:: "<<fMT2tree->misc.MT2<<endl;
// 	    }
	  }
	
	  for(int j = 0; j < fMT2tree->NTaus; j++){

	    float DeltaPhi = Util::DeltaPhi(fMT2tree->tau[j].lv.Phi(), METPhi);

	    if(DeltaPhi > maxDPhiTauMET)
	      continue;
	
	    if(fMT2tree->tau[j].MT < 100 && fMT2tree->tau[j].Isolation3Hits > IsoCut &&  Filled == 0){
	      Filled = 1;
	      hMT2TauSig->Fill(fMT2tree->misc.MT2, weight);
	      //	      cout<<"TauSig:: "<<fMT2tree->misc.MT2<<endl;
	      float minDR = 100.0;
	      
	      float tauEta = fMT2tree->tau[j].lv.Eta();

	      float tauPhi = fMT2tree->tau[j].lv.Phi();
	    }
	  }
	}
      }
    }//for (Long64_t jentry=0; jentry<min(nentries, nevents);jentry++) 
    cout<<" hNAllTaus->Integral() "<<hNAllTaus->Integral()<<endl;
    cout<<" nAllTaus In ThisSample "<<(hNAllTaus->Integral() - nAllTausBeforeThisSample)<<endl;
  }// for(int i = 0; i < fSamples.size(); i++){

 for(int j = 1; j < NumberOfMT2Bins; j++){
   float id = hMT2TauFound  ->GetBinContent(j);
   float jet= hNAllTauJets->GetBinContent(j);
   float gen= hNAllTaus   ->GetBinContent(j);
   float EffTotal = -1.0;
   if(gen != 0) EffTotal = id/gen;
   float Eff = -1.0;
   if(jet != 0) Eff = id/jet;
   float Acc = -1.0;
   if(gen != 0) Acc = jet/gen;

   cout<<"$"<<xMT2bin[j-1]<<"-"<<xMT2bin[j]<<" Eff = "<<EffTotal<<", Eff * Acc = "<<Eff<<" * "<<Acc<<endl;
 }

 int p = 0;
 hMT2TauFound->Divide(hNAllTaus);
 Cout(p++, hMT2TauFound);

 //To add the syatematic uncertainty
 for(int j = 1; j < NumberOfMT2Bins; j++){
   float val = hMT2TauFound->GetBinContent(j);
   float Err = hMT2TauFound->GetBinError(j);
   hMT2TauFound->SetBinError(j, sqrt(Err * Err + sys * sys * val * val));
 }

 Cout(p++, hEtaPtTauFake);
 Cout(p++, hEtaPtJet);
 float TauFake  = hEtaPtTauFake->Integral();
  float Jet      = hEtaPtJet    ->Integral();
  float Fake     = TauFake/Jet;
  cout<<" Fake "<<Fake<<endl;
  hEtaPtTauFake->Divide(hEtaPtJet);
  Cout(p++, hEtaPtTauFake);

  Cout(p++, hMT2JetSig);

  //  hMT2JetSig->Scale(Fake);

  Cout(p++, hMT2JetSig);
  Cout(p++, hMT2TauSig);
  hMT2JetSig->Add(hMT2TauSig, -1);
  hMT2JetSig->Scale(-1.0);
  Cout(p++, hMT2JetSig);

  hMT2JetSig->Divide(hMT2TauFound);
  Cout(p++, hMT2JetSig);

  float TauFakeData  = hEtaPtTauFakeData->Integral();
  float JetData      = hEtaPtJetData    ->Integral();
  float FakeData     = TauFakeData/JetData;
  cout<<" FakeData "<<FakeData<<endl;
  hEtaPtTauFakeData->Divide(hEtaPtJetData);
  
  Cout(p++, hMT2JetSigData);
  //   hMT2JetSigData->Scale(FakeData);
  
  Cout(p++, hMT2JetSigData);
  Cout(p++, hMT2TauSigData);
  hMT2JetSigData->Add(hMT2TauSigData, -1);
  hMT2JetSigData->Scale(-1.0);
  Cout(p++, hMT2JetSigData);
    

  hMT2JetSigData->Divide(hMT2TauFound);
  
  Cout(p++, hMT2JetSigData);
  

  Cout(p++, hRatioRelaxed);
  Cout(p++, hRatioTight);

  hRatioRelaxed->Add(hRatioTight);
  
  hRatioTight->Divide(hRatioRelaxed);

  Cout(p++, hRatioTight);

//   hMT2JetSig->Multiply(hRatioTight);

//   hMT2JetSigData->Multiply(hRatioTight);

//   hNAllTaus->Multiply(hRatioTight);
 
  for(int j = 1; j < NumberOfMT2Bins; j++){
    
    float TotalErrSim = hMT2JetSig    ->GetBinError(j);
    float StatErrSim  = (hMT2TauSig    ->GetBinError(j)/hMT2TauFound->GetBinContent(j))*hRatioTight->GetBinContent(j);
    float SysErrSim   = sqrt(fabs(TotalErrSim * TotalErrSim - StatErrSim * StatErrSim));

    float TotalErrData = hMT2JetSigData    ->GetBinError(j);
    float StatErrData  = (hMT2TauSigData    ->GetBinError(j)/hMT2TauFound->GetBinContent(j))*hRatioTight->GetBinContent(j);
    float SysErrData   = sqrt(fabs(TotalErrData * TotalErrData - StatErrData * StatErrData));
    
    cout<<"$"<<xMT2bin[j-1]<<"-"<<xMT2bin[j]<<"$: Sim  "<<hNAllTaus     ->GetBinContent(j)<<" +/- "<<hNAllTaus        ->GetBinError(j)*hRatioTight->GetBinContent(j)<<
                                                 " Pre  "<<hMT2JetSig    ->GetBinContent(j)<<" +/- "<<StatErrSim <<" +/- "<< SysErrSim <<                              
                                                 " Data "<<hMT2JetSigData->GetBinContent(j)<<" +/- "<<StatErrData<<" +/- "<< SysErrData<<endl;
  }

  for(int j = 0; j < NumberOfSamples; j++){
     AddOverAndUnderFlow(MT2[j]);
  }
  printYield();

THStack* h_stack     = new THStack(varname, "");
 double xbin[7] = {0.0, 35.0, 70.0, 110.0, 150.0, 250.0, 2000};
 TH1F *hnew[NumberOfSamples];
  for(int j = 0; j < NumberOfSamples; j++){
    hnew[j] = (TH1F*)MT2[j]->Rebin(6,"hnew",xbin);
    TH1F* mt2 = (TH1F*)MT2[j]->Clone();
    mt2->SetName("mt2");
    if(j < (NumberOfSamples - 2))
      //h_stack  -> Add(MT2[j]);
      h_stack  -> Add(hnew[j]);
    delete mt2;
  }  

  TLegend* Legend1 = new TLegend(.71,.54,.91,.92);
  Legend1->AddEntry(MT2[0], "QCD", "f");
  Legend1->AddEntry(MT2[1], "W+jets", "f");
  Legend1->AddEntry(MT2[2], "Z+jets", "f");
  Legend1->AddEntry(MT2[3], "Top", "f");
  Legend1->AddEntry(MT2[5], "data", "l");
  
  //  TLegend *Legend1;
  printHisto(h_stack, MT2[5]->Rebin(6,"hnew1",xbin), MT2[4]->Rebin(6,"hnew4",xbin), Legend1 , "MTC", "hist", true, "MT2", "Events", -4, 0);

  plotRatioStack(h_stack,  MT2[4]->Rebin(6,"hnew2",xbin), MT2[5]->Rebin(6,"hnew3",xbin), true, false, "MT2_ratio", Legend1, "MT2", "Events", -4, 0);

  TCanvas *CanvasMT2 = new TCanvas("myMT2", "myMT2");
  hNAllTaus      ->SetLineColor(kBlue);
  hMT2JetSig     ->SetLineColor(kRed);
  hMT2JetSigData ->SetLineColor(kGreen);
  hNAllTaus      ->Draw();
  hMT2JetSig     ->Draw("sames");
  hMT2JetSigData ->Draw("sames");
}


void MassPlotter::NewTauContamination(int sample_index, Long64_t nevents, int flag){
  float sys = 0.00; //systematic uncertanty on the fake rate and efficiency
  float IsoCut = 1.5;//cut on the isolation 0.5:VLoose, 1.5:Loose...
  float IsoCutLoose = IsoCut - 1.0;
  float maxDPhiTauMET = 11.5;//To cut on the DeltaPhi between the tau and MET

  setFlags(flag); 

  TH1::SetDefaultSumw2();
  TString  cnames[NumberOfSamples] = {"QCD", "Wjets", "Zjets", "Top", "MC", "data"};
  int      ccolor[NumberOfSamples] = {  401,     417,     419,   600,  500,  632};
  TString varname = "MT2";
  for (int i=0; i<NumberOfSamples; i++){
    MT2[i] = new TH1F(varname+"_"+cnames[i], "", 150, 0, 750);
    MT2[i] -> SetFillColor  (ccolor[i]);
    MT2[i] -> SetLineColor  (ccolor[i]);
    MT2[i] -> SetLineWidth  (2);
    MT2[i] -> SetMarkerColor(ccolor[i]);
    MT2[i] -> SetStats(false);
  }

  MT2[5] -> SetMarkerStyle(20);
  MT2[5] -> SetMarkerColor(kBlack);
  MT2[5] -> SetLineColor(kBlack);
  
  MT2[4] -> SetFillStyle(3004);
  MT2[4] -> SetFillColor(kBlack);

  TH2F* hEtaPtJet           = new TH2F("hEtaPtJet",           "hEtaPtJet",          1,0,1000.0,1,0,2.5);//tauable jets in MC in QCD dominated region (MT2 < 70)
  TH2F* hEtaPtTauFake       = new TH2F("hEtaPtTauFake",       "hEtaPtTauFake",      1,0,1000.0,1,0,2.5);//taus in MC in QCD dominated region  (MT2 < 70)
  TH2F* hEtaPtJetData       = new TH2F("hEtaPtJetData",       "hEtaPtJetData",      1,0,1000.0,1,0,2.5);//tauable jets in data in QCD dominated region (MT2 < 70)
  TH2F* hEtaPtTauFakeData   = new TH2F("hEtaPtTauFakeData",   "hEtaPtTauFakeData",  1,0,1000.0,1,0,2.5);//taus in data in QCD dominated region  (MT2 < 70)

  TH1F* hMT2TauTightFound        = new TH1F("hMT2TauTightFound",   "hMT2TauTightFound",   NumberOfMT2Bins - 1, xMT2bin);//Rec/IdT taus
  TH1F* hMT2TauLooseFound        = new TH1F("hMT2TauLooseFound",   "hMT2TauLooseFound",   NumberOfMT2Bins - 1, xMT2bin);//Rec/IdL taus

  TH1F* hMT2Tau0Sig         = new TH1F("hMT2Tau0Sig",     "hMT2Tau0Sig",     NumberOfMT2Bins - 1, xMT2bin);//events without tau(MC)
  TH1F* hMT2Tau1Sig         = new TH1F("hMT2Tau1Sig",     "hMT2Tau1Sig",     NumberOfMT2Bins - 1, xMT2bin);//taus in the signal region (MC)
  TH1F* hMT2Tau0SigData     = new TH1F("hMT2Tau0SigData", "hMT2Tau0SigData", NumberOfMT2Bins - 1, xMT2bin);//events without tau (data)
  TH1F* hMT2Tau1SigData     = new TH1F("hMT2Tau1SigData", "hMT2Tau1SigData", NumberOfMT2Bins - 1, xMT2bin);//taus in the signal region (data)
  TH1F* hRatioTight         = new TH1F("hRatioTight",       "hRatioTight",       NumberOfMT2Bins - 1, xMT2bin);//Used for the final rescaling the estimation from relaxed to appplied cut
  TH1F* hRatioRelaxed       = new TH1F("hRatioRelaxed",       "hRatioRelaxed",       NumberOfMT2Bins - 1, xMT2bin);//Used for the final rescaling the estimation from relaxed to appplied cut

  TH1F* hNAllTaus           = new TH1F("hNAllTaus",      "hNAllTaus",      NumberOfMT2Bins - 1, xMT2bin);//All of the generated hadronic taus
  TH1F* hNAllTauJets        = new TH1F("hNAllTauJets",   "hNAllTauJets",   NumberOfMT2Bins - 1, xMT2bin);//All of the generated hadronic taus matched with a rec jet

  TH1F* hMT2TauFake         = new TH1F("hMT2TauFake",   "hMT2TauFake",   NumberOfMT2Bins - 1, xMT2bin);//Fake rate (MC)
  TH1F* hMT2TauFakeData     = new TH1F("hMT2TauFakeData", "hMT2TauFakeData",   NumberOfMT2Bins - 1, xMT2bin);//Fake rate (data)
  TH1F* hMT2Unit            = new TH1F("hMT2Unit",      "hMT2Unit",   NumberOfMT2Bins - 1, xMT2bin);//Fake rate (MC)
 
  for(int ii = 0; ii < fSamples.size(); ii++){

    if(sample_index != -1 && ii != sample_index)
      continue;

    int data = 0;
    sample Sample = fSamples[ii];
    
    if(Sample.type == "data")
      data = 1;

    if(Sample.type == "susy")// || Sample.sname == "QCD") 
      continue;

    fMT2tree = new MT2tree();
    Sample.tree->SetBranchAddress("MT2tree", &fMT2tree);
    Long64_t nentries =  Sample.tree->GetEntries();

    float Weight = Sample.xsection * Sample.kfact * Sample.lumi / (Sample.nevents * Sample.PU_avg_weight );

    std::cout << setfill('=') << std::setw(70) << "" << std::endl;
    cout << "looping over :     " <<endl;	
    cout << "   Name:           " << Sample.name << endl;
    cout << "   File:           " << (Sample.file)->GetName() << endl;
    cout << "   Events:         " << Sample.nevents  << endl;
    cout << "   Events in tree: " << Sample.tree->GetEntries() << endl; 
    cout << "   Xsection:       " << Sample.xsection << endl;
    cout << "   kfactor:        " << Sample.kfact << endl;
    cout << "   avg PU weight:  " << Sample.PU_avg_weight << endl;
    cout << "   Weight:         " << Weight <<endl;
    std::cout << setfill('-') << std::setw(70) << "" << std::endl;
   
    for (Long64_t jentry=0; jentry<min(nentries, nevents);jentry++) {
      Sample.tree->GetEntry(jentry); 
      if ( fVerbose>2 && jentry % 100000 == 0 ){  
	fprintf(stdout, "\rProcessed events: %6d of %6d ", jentry + 1, nentries);
	fflush(stdout);
      }
 
      float weight = Weight;
      if(data == 0)
	weight *= fMT2tree->pileUp.Weight;
      else
	weight = 1.0;

      if (!(fMT2tree->misc.MET             >  30                                  &&
	    fMT2tree->NBJetsCSVT           > 0                                    &&
            fMT2tree->NJetsIDLoose40       >= 4                                   && 
           (fMT2tree->NEles==0  ||  fMT2tree->ele[0].lv.Pt()<5)                   &&
 	   (fMT2tree->NMuons==0 ||  fMT2tree->muo[0].lv.Pt()<5)                   &&
 	    fMT2tree->misc.Jet0Pass        == 1                                   &&
 	    fMT2tree->misc.Jet1Pass        == 1                                   &&
 	    fMT2tree->misc.PassJetID       == 1                                   &&
	    fMT2tree->misc.Vectorsumpt     <  70                                  &&
	    //          (fMT2tree->misc.MinMetJetDPhi >0.3||fMT2tree->misc.MinMetJetDPhiIndex>3) &&
	    fMT2tree->misc.MinMetJetDPhi4  >  0.3                                  &&
	    //Six Or Quad
	    ((fMT2tree->NJetsIDLoose40  >=4 && fMT2tree->jet[3].lv.Pt() >= 100.0)  || 
	     (fMT2tree->NJetsIDLoose40  >=5 && fMT2tree->jet[4].lv.Pt() >= 85.0)   || 
	     (fMT2tree->NJetsIDLoose40  >=6 && fMT2tree->jet[5].lv.Pt() >= 65.0)   || 
	     (fMT2tree->NJetsIDLoose40  >=7 && fMT2tree->jet[6].lv.Pt() >= 55.0))  &&
	    
	    (fMT2tree->misc.ProcessID!=0 || 
	     ((fMT2tree->misc.CrazyHCAL==0          && fMT2tree->misc.NegativeJEC==0   && 
	       fMT2tree->misc.CSCTightHaloIDFlag==0 && fMT2tree->misc.HBHENoiseFlag==0 && 
	       fMT2tree->misc.hcalLaserEventFlag==0 && fMT2tree->misc.trackingFailureFlag==0 && 
	       fMT2tree->misc.eeBadScFlag==0 && fMT2tree->misc.EcalDeadCellTriggerPrimitiveFlag==0 )&&
	      //Six Or Quad
	      ((((fMT2tree->trigger.HLT_QuadJet80_v1 ==1)||(fMT2tree->trigger.HLT_QuadJet80_v2 ==1)||
		 (fMT2tree->trigger.HLT_QuadJet80_v3 ==1)) && 
		((fMT2tree->NJetsIDLoose40  >=4 && fMT2tree->jet[3].lv.Pt() >= 100.0 ) || 
		 (fMT2tree->NJetsIDLoose40  >=5 && fMT2tree->jet[4].lv.Pt() >= 85.0 )))|| 

	       (((fMT2tree->trigger.HLT_SixJet45_v1 ==1)||(fMT2tree->trigger.HLT_SixJet45_v2 ==1)||
		 (fMT2tree->trigger.HLT_SixJet45_v3 ==1)) && 
		((fMT2tree->NJetsIDLoose40  >=6 && fMT2tree->jet[5].lv.Pt() >= 65.0 ) || 
		 (fMT2tree->NJetsIDLoose40  >=7 && fMT2tree->jet[6].lv.Pt() >= 55.0 ))))))))
	continue;

      float METPhi = fMT2tree->misc.METPhi; //To be used if a cone close to met is used to look for taus
	  
      for(int t = 0; t < fMT2tree->NTaus; t++){

 	if(fMT2tree->tau[t].MT > 100 || fMT2tree->tau[t].Isolation < IsoCut)
 	  continue;

	float minDeltaR = 100.0;

 	float DeltaPhi = Util::DeltaPhi(fMT2tree->tau[t].lv.Phi(), fMT2tree->misc.METPhi);

	if(DeltaPhi > maxDPhiTauMET)
	  continue;

	if(data == 1){
	    
	  MT2[5]->Fill(fMT2tree->misc.MT2, weight);//data
	    
	}else{
	  MT2[4]->Fill(fMT2tree->misc.MT2, weight);
	    
	  if(Sample.sname == "Top")  
	    MT2[3]->Fill(fMT2tree->misc.MT2, weight);
	  else
	    if(Sample.sname == "DY")	
	      MT2[2]->Fill(fMT2tree->misc.MT2, weight);
	    else
	      if(Sample.sname == "Wtolnu")
		MT2[1]->Fill(fMT2tree->misc.MT2, weight);
	      else
		if(Sample.sname == "QCD")
		  MT2[0]->Fill(fMT2tree->misc.MT2, weight);
	  if(fMT2tree->misc.MT2 >= MT2Min)
	    cout<<"Tau:: "<<fMT2tree->misc.MT2<<endl;
	}}

      int Filled = 0;

      for(int j = 0; j < fMT2tree->NGenLepts; j++){
	if(abs(fMT2tree->genlept[j].ID) != 15)
	  continue;
	
	int HadTau = 1;
	    
	for(int i = 0; i < fMT2tree->NGenLepts; i++){
	  if((abs(fMT2tree->genlept[i].ID) == 11 && fMT2tree->genlept[i].MID == fMT2tree->genlept[j].ID) || 
	     (abs(fMT2tree->genlept[i].ID) == 13 && fMT2tree->genlept[i].MID == fMT2tree->genlept[j].ID)){
	    HadTau = 0;
	  }
	}

	if(HadTau == 0)
	  continue;
        
	if(Filled == 0){
	  //Filled = 1;
	  if((mT2 && fMT2tree->misc.MinMetJetDPhi   >  0.3) || (mT2b && fMT2tree->misc.MinMetJetDPhi4   >  0.3))
	    hRatioTight->Fill(fMT2tree->misc.MT2, weight);
	  else
	    hRatioRelaxed->Fill(fMT2tree->misc.MT2, weight);	  
	}
	
	hNAllTaus->Fill(fMT2tree->misc.MT2, weight);
	

	//To find the jet matched to this gen had tau
	int jetIndex = -1;
	float minDR = 100.0;

	float genleptEta = fMT2tree->genlept[j].lv.Eta();

	float genleptPhi = fMT2tree->genlept[j].lv.Phi();

	for(int k = 0; k < fMT2tree->NJets; k++){
	  float deltaR = Util::GetDeltaR(genleptEta, fMT2tree->jet[k].lv.Eta(), genleptPhi, fMT2tree->jet[k].lv.Phi());
	  if(deltaR < minDR){
	    minDR = deltaR;
	    jetIndex = k;
	  }
	}	

	if(minDR < 0.1){
	  hNAllTauJets->Fill(fMT2tree->misc.MT2, weight);
	}

	jetIndex = -1;
	minDR = 100.0;

	//To find the rec/id tau matched to this gen had tau
	for(int k = 0; k < fMT2tree->NTaus; k++){
	  float deltaR = Util::GetDeltaR(genleptEta, fMT2tree->tau[k].lv.Eta(), genleptPhi, fMT2tree->tau[k].lv.Phi());
	  if(deltaR < minDR){
	    minDR = deltaR;
	    jetIndex = k;
	  }
	}

	float DeltaPhi = Util::DeltaPhi(fMT2tree->tau[jetIndex].lv.Phi(), METPhi);

	if(DeltaPhi > maxDPhiTauMET)
	  continue;
	
	if(minDR < 0.1){
	  if(fMT2tree->tau[jetIndex].Isolation > IsoCutLoose && fMT2tree->tau[jetIndex].MT < 100)
	    hMT2TauLooseFound->Fill(fMT2tree->misc.MT2, weight);

	  if(fMT2tree->tau[jetIndex].Isolation > IsoCut && fMT2tree->tau[jetIndex].MT < 100)
	    hMT2TauTightFound->Fill(fMT2tree->misc.MT2, weight);
	}
      }
      
      //QCD dominated region, use these events to find the fake rate
      if(fMT2tree->misc.MT2 < 70){
	if(data == 1){
	  for(int j = 0; j < fMT2tree->NJets; j++){
	    
	    if(fMT2tree->jet[j].isTauMatch < 0  || fMT2tree->tau[fMT2tree->jet[j].isTauMatch].MT > 100)
	      continue;

	    float DeltaPhi = Util::DeltaPhi(fMT2tree->jet[j].lv.Phi(), METPhi);

	    if(DeltaPhi > maxDPhiTauMET)
	      continue;
	    if(fMT2tree->tau[fMT2tree->jet[j].isTauMatch].Isolation > IsoCutLoose){
	      hEtaPtJetData->Fill(fMT2tree->jet[j].lv.Pt(), fabs(fMT2tree->jet[j].lv.Eta()), weight);
	      
	      if(fMT2tree->tau[fMT2tree->jet[j].isTauMatch].Isolation > IsoCut && fMT2tree->tau[fMT2tree->jet[j].isTauMatch].MT < 100){  
		hEtaPtTauFakeData->Fill(fMT2tree->jet[j].lv.Pt(), fabs(fMT2tree->jet[j].lv.Eta()), weight);
	      }}}}else{
	  for(int j = 0; j < fMT2tree->NJets; j++){
	    
	    if(fMT2tree->jet[j].isTauMatch < 0  || fMT2tree->tau[fMT2tree->jet[j].isTauMatch].MT > 100)
	      continue;	    

	    float DeltaPhi = Util::DeltaPhi(fMT2tree->jet[j].lv.Phi(), METPhi);

	    if(DeltaPhi > maxDPhiTauMET)
	      continue;

	    if(fMT2tree->tau[fMT2tree->jet[j].isTauMatch].Isolation > IsoCutLoose){
	      hEtaPtJet->Fill(fMT2tree->jet[j].lv.Pt(), fabs(fMT2tree->jet[j].lv.Eta()), weight);
	      
	      if(fMT2tree->tau[fMT2tree->jet[j].isTauMatch].Isolation > IsoCut && fMT2tree->tau[fMT2tree->jet[j].isTauMatch].MT < 100){ 
		
		hEtaPtTauFake->Fill(fMT2tree->jet[j].lv.Pt(), fabs(fMT2tree->jet[j].lv.Eta()), weight);
	      }}}}}
      
      //Fill The number of Rec/Id taus and tauable jets in the signal region
      if(fMT2tree->misc.MT2 >= MT2Min){
	if(data == 1){
	  int nTausLoose = 0;
	  int nTaus = 0;
	  for(int j = 0; j < fMT2tree->NTaus; j++){

	    float DeltaPhi = Util::DeltaPhi(fMT2tree->tau[j].lv.Phi(), METPhi);

	    if(DeltaPhi > maxDPhiTauMET)
	      continue;
	    
	    if(fMT2tree->tau[j].MT < 100){
	      if(fMT2tree->tau[j].Isolation > IsoCutLoose){
		nTausLoose++;
		cout<<"TauSigLooseData:: "<<fMT2tree->misc.MT2<<endl;
		if(fMT2tree->tau[j].Isolation > IsoCut){
		  nTaus++;
		  cout<<"TauSigData:: "<<fMT2tree->misc.MT2<<endl;}}}}
	    
	    if(nTausLoose > 0){
	      if(nTaus == 0)
		//	nTaus=0;
		hMT2Tau0SigData->Fill(fMT2tree->misc.MT2, weight);
	      else
		hMT2Tau1SigData->Fill(fMT2tree->misc.MT2, weight);}	    
	  }else{
	  int nTausLoose = 0;
	  int nTaus = 0;
	  
	  for(int j = 0; j < fMT2tree->NTaus; j++){

	    float DeltaPhi = Util::DeltaPhi(fMT2tree->tau[j].lv.Phi(), METPhi);

	    if(DeltaPhi > maxDPhiTauMET)
	      continue;
	
	    if(fMT2tree->tau[j].MT < 100){
	      if(fMT2tree->tau[j].Isolation > IsoCutLoose){
		nTausLoose++;
		cout<<"TauSigLoose:: "<<fMT2tree->misc.MT2<<endl;
		if(fMT2tree->tau[j].Isolation > IsoCut){
		  nTaus++;
		  cout<<"TauSig:: "<<fMT2tree->misc.MT2<<endl;}}}}
	    
	  if(nTausLoose > 0){
	    if(nTaus == 0){
	      hMT2Tau0Sig->Fill(fMT2tree->misc.MT2, weight);}
	    else
	      hMT2Tau1Sig->Fill(fMT2tree->misc.MT2, weight);}}}
    }//for (Long64_t jentry=0; jentry<min(nentries, nevents);jentry++) 
  }// for(int i = 0; i < fSamples.size(); i++){


 for(int j = 1; j < NumberOfMT2Bins; j++){
   float idT = hMT2TauTightFound  ->GetBinContent(j);
   float idL = hMT2TauLooseFound  ->GetBinContent(j);
   float jet= hNAllTauJets->GetBinContent(j);
   float gen= hNAllTaus   ->GetBinContent(j);
   float EffTotal = -1.0;
   if(gen != 0) EffTotal = idT/gen;
   float EffLtoT = -1.0;
   if(idL != 0) EffLtoT = idT/idL;
   float Eff = -1.0;
   if(jet != 0) Eff = idL/jet;
   float Acc = -1.0;
   if(gen != 0) Acc = jet/gen;

   cout<<"$"<<xMT2bin[j-1]<<"-"<<xMT2bin[j]<<" Eff = "<<EffTotal<<", LtoT * Eff * Acc = "<<EffLtoT<<" * "<<Eff<<" * "<<Acc<<endl;
 }

 int p = 0;
 Cout(p++, hMT2TauTightFound);
 Cout(p++, hMT2TauLooseFound);
 Cout(p++, hNAllTauJets);
 Cout(p++, hNAllTaus);
 hMT2TauTightFound->Divide(hMT2TauLooseFound);
 cout<<"hMT2TauTightFound->Divide(hMT2TauLooseFound);"<<endl;
 Cout(p++, hMT2TauTightFound);

 //To add the syatematic uncertainty
 for(int j = 1; j < NumberOfMT2Bins; j++){
   float val = hMT2TauTightFound->GetBinContent(j);
   float Err = hMT2TauTightFound->GetBinError(j);
   hMT2TauTightFound->SetBinError(j, sqrt(Err * Err + sys * sys * val * val));
 }

 Cout(p++, hEtaPtTauFake);
 Cout(p++, hEtaPtJet);
 float TauFake  = hEtaPtTauFake->Integral();
 float Jet      = hEtaPtJet    ->Integral();
 float Fake     = TauFake/Jet;
 cout<<" Fake "<<Fake<<endl;
 hEtaPtTauFake->Divide(hEtaPtJet);
 Cout(p++, hEtaPtTauFake);
 for(int j = 1; j < NumberOfMT2Bins; j++){
   hMT2TauFake->SetBinContent(j,hEtaPtTauFake->GetBinContent(4));
   hMT2TauFake->SetBinError(j,hEtaPtTauFake->GetBinError(4));
   hMT2Unit->SetBinContent(j,1.0);
 }

  Cout(p++, hMT2TauFake);
  Cout(p++, hMT2Unit);

  Cout(p++, hMT2Tau0Sig);

  hMT2Tau0Sig->Multiply(hMT2TauFake);
  cout<<"hMT2Tau0Sig->Multiply(hMT2TauFake);"<<endl;
  Cout(p++, hMT2Tau0Sig);

  hMT2Unit->Add(hMT2TauFake, -1.0);
  cout<<"hMT2Unit->Add(hMT2TauFake, -1.0);"<<endl;
  Cout(p++, hMT2Unit);

  Cout(p++, hMT2Tau1Sig);

  hMT2Tau1Sig->Multiply(hMT2Unit);
  cout<<"hMT2Tau1Sig->Multiply(hMT2Unit);"<<endl;
  Cout(p++, hMT2Tau1Sig);

  hMT2Tau0Sig->Add(hMT2Tau1Sig, -1.0);
  cout<<"hMT2Tau0Sig->Add(hMT2Tau1Sig, -1.0);"<<endl;
  Cout(p++, hMT2Tau0Sig);

  hMT2TauFake->Add(hMT2TauTightFound, -1.0);
  cout<<"hMT2TauFake->Add(hMT2TauTightFound, -1.0);"<<endl;
  Cout(p++, hMT2TauFake);

  hMT2Tau0Sig->Divide(hMT2TauFake);
  cout<<"hMT2Tau0Sig->Divide(hMT2TauFake);"<<endl;
  Cout(p++, hMT2Tau0Sig);

  Cout(p++, hNAllTaus);

  hMT2TauLooseFound->Divide(hNAllTaus);
  cout<<"hMT2TauLooseFound->Divide(hNAllTaus);"<<endl;
  Cout(p++, hMT2TauLooseFound);
  
  hMT2Tau0Sig->Divide(hMT2TauLooseFound);
  cout<<"hMT2Tau0Sig->Divide(hMT2TauLooseFound);"<<endl;
  Cout(p++, hMT2Tau0Sig);



  float TauFakeData  = hEtaPtTauFakeData->Integral();
  float JetData      = hEtaPtJetData    ->Integral();
  float FakeData     = TauFakeData/JetData;
  Cout(p++, hEtaPtTauFakeData);
  Cout(p++, hEtaPtJetData);

  cout<<" FakeData "<<FakeData<<endl;
  hEtaPtTauFakeData->Divide(hEtaPtJetData);

  for(int j = 1; j < NumberOfMT2Bins; j++){
    hMT2TauFakeData->SetBinContent(j,hEtaPtTauFakeData->GetBinContent(4));
    hMT2TauFakeData->SetBinError(j,hEtaPtTauFakeData->GetBinError(4));
    hMT2Unit->SetBinContent(j,1.0);
    hMT2Unit->SetBinError(j,0.0);
  }
  Cout(p++, hMT2TauFakeData);
  Cout(p++, hMT2Unit);

  Cout(p++, hMT2Tau0SigData);

  hMT2Tau0SigData->Multiply(hMT2TauFakeData);
  cout<<"hMT2Tau0SigData->Multiply(hMT2TauFakeData);"<<endl;
  Cout(p++, hMT2Tau0SigData);

  hMT2Unit->Add(hMT2TauFakeData, -1.0);
  cout<<"hMT2Unit->Add(hMT2TauFakeData, -1.0);"<<endl;
  Cout(p++, hMT2Unit);

  Cout(p++, hMT2Tau1SigData);
  
  hMT2Tau1SigData->Multiply(hMT2Unit);
  cout<<"hMT2Tau1SigData->Multiply(hMT2Unit);"<<endl;
  Cout(p++, hMT2Tau1SigData);

  hMT2Tau0SigData->Add(hMT2Tau1SigData, -1.0);
  cout<<"hMT2Tau0SigData->Add(hMT2Tau1SigData, -1.0);"<<endl;
  Cout(p++, hMT2Tau0SigData);

  hMT2TauFakeData->Add(hMT2TauTightFound, -1.0);
  cout<<"hMT2TauFakeData->Add(hMT2TauTightFound, -1.0);"<<endl;
  Cout(p++, hMT2TauFakeData);

  hMT2Tau0SigData->Divide(hMT2TauFakeData);
  cout<<"hMT2Tau0SigData->Divide(hMT2TauFakeData);"<<endl;
  Cout(p++, hMT2Tau0SigData);

  Cout(p++, hMT2TauLooseFound);
  
  hMT2Tau0SigData->Divide(hMT2TauLooseFound);
  cout<<"hMT2Tau0SigData->Divide(hMT2TauLooseFound);"<<endl;
  Cout(p++, hMT2Tau0SigData);


  Cout(p++, hRatioRelaxed);
  Cout(p++, hRatioTight);

  hRatioRelaxed->Add(hRatioTight);
  
  hRatioTight->Divide(hRatioRelaxed);

  Cout(p++, hRatioTight);

//   hMT2Tau0Sig->Multiply(hRatioTight);

//   hMT2Tau0SigData->Multiply(hRatioTight);

//   hNAllTaus->Multiply(hRatioTight);
 
  for(int j = 1; j < NumberOfMT2Bins; j++){
    
    float TotalErrSim = hMT2Tau0Sig    ->GetBinError(j);
//     float StatErrSim  = (hMT2Tau1Sig    ->GetBinError(j)/hMT2TauFound->GetBinContent(j))*hRatioTight->GetBinContent(j);
//     float SysErrSim   = sqrt(fabs(TotalErrSim * TotalErrSim - StatErrSim * StatErrSim));

     float TotalErrData = hMT2Tau0SigData    ->GetBinError(j);
//     float StatErrData  = (hMT2Tau1SigData    ->GetBinError(j)/hMT2TauFound->GetBinContent(j))*hRatioTight->GetBinContent(j);
//     float SysErrData   = sqrt(fabs(TotalErrData * TotalErrData - StatErrData * StatErrData));
    
     cout<<"$"<<xMT2bin[j-1]<<"-"<<xMT2bin[j]<<"$: Sim  "<<hNAllTaus     ->GetBinContent(j)<<" +/- "<<hNAllTaus        ->GetBinError(j)<<//*hRatioTight->GetBinContent(j)<<
       " Pre  "<<hMT2Tau0Sig    ->GetBinContent(j)<<" +/- "<<TotalErrSim<<//" +/- "<<StatErrSim <<" +/- "<< SysErrSim <<                              
      " Data "<<hMT2Tau0SigData->GetBinContent(j)<<" +/- "<<TotalErrData<<//" +/- "<<StatErrData<<" +/- "<< SysErrData<<
      endl;
  }

  for(int j = 0; j < NumberOfSamples; j++){
     AddOverAndUnderFlow(MT2[j]);
  }
  printYield();

THStack* h_stack     = new THStack(varname, "");
 double xbin[5] = {0.0, 35.0, 70.0, 110.0, 750.0};
 TH1F *hnew[NumberOfSamples];
  for(int j = 0; j < NumberOfSamples; j++){
    hnew[j] = (TH1F*)MT2[j]->Rebin(4,"hnew",xbin);
    TH1F* mt2 = (TH1F*)MT2[j]->Clone();
    mt2->SetName("mt2");
    if(j < (NumberOfSamples - 2))
      //h_stack  -> Add(MT2[j]);
      h_stack  -> Add(hnew[j]);
    delete mt2;
  }  

  TLegend* Legend1 = new TLegend(.71,.54,.91,.92);
  Legend1->AddEntry(MT2[0], "QCD", "f");
  Legend1->AddEntry(MT2[1], "W+jets", "f");
  Legend1->AddEntry(MT2[2], "Z+jets", "f");
  Legend1->AddEntry(MT2[3], "Top", "f");
  Legend1->AddEntry(MT2[5], "data", "l");
  
  //  TLegend *Legend1;
  printHisto(h_stack, MT2[5]->Rebin(4,"hnew1",xbin), MT2[4]->Rebin(4,"hnew4",xbin), Legend1 , "MTC", "hist", true, "MT2", "Events", -4, 0);

  plotRatioStack(h_stack,  MT2[4]->Rebin(4,"hnew2",xbin), MT2[5]->Rebin(4,"hnew3",xbin), true, false, "MT2_ratio", Legend1, "MT2", "Events", -4, 0);

  TCanvas *CanvasMT2 = new TCanvas("myMT2", "myMT2");
  hNAllTaus       ->SetLineColor(kBlue);
  hMT2Tau0Sig     ->SetLineColor(kRed);
  hMT2Tau0SigData ->SetLineColor(kGreen);
  hNAllTaus       ->Draw();
  hMT2Tau0Sig     ->Draw("sames");
  hMT2Tau0SigData ->Draw("sames");
}

void MassPlotter::setFlags(int flag){
    mT2  = 0;
    mT2b = 0;
    High = 0;
    Low  = 0;
    BTagRelaxed = 0;

    if(flag < 10)  mT2 = 1;
    if(flag > 10)  mT2b= 1;
    if(flag == 5 || flag == 15)  Low = 1;
    if(flag == 7 || flag == 17)  High = 1;

    if(flag == 19){
      High = 1;
      BTagRelaxed = 1;}

    if(mT2){
      MT2Min = 150;

      xMT2bin[1] = 200.0;
      xMT2bin[2] = 275.0;
      xMT2bin[3] = 375.0;
      xMT2bin[4] = 500.0;
      xMT2bin[5] = 2000.0;
      
      newxMT2bin[0] = 0.0;
      newxMT2bin[1] = 30.0;
      newxMT2bin[2] = 90.0;
      newxMT2bin[3] = 120.0;
      newxMT2bin[4] = 150.0;
      newxMT2bin[5] = 200.0;
      newxMT2bin[6] = 275.0;
      newxMT2bin[7] = 375.0;
      newxMT2bin[8] = 500.0;
      newxMT2bin[9] = 750.0;
    }

    if(mT2b){
      MT2Min = 125;

      xMT2bin[1] = 150.0;
      xMT2bin[2] = 200.0;
      xMT2bin[3] = 300.0;
      xMT2bin[4] = 2000.0;
      xMT2bin[5] = 3000.0;
  
      newxMT2bin[0] = 0.0;
      newxMT2bin[1] = 25.0;
      newxMT2bin[2] = 50.0;
      newxMT2bin[3] = 75.0;
      newxMT2bin[4] = 100.0;
      newxMT2bin[5] = 125.0;
      newxMT2bin[6] = 150.0;
      newxMT2bin[7] = 200.0;
      newxMT2bin[8] = 300.0;
      newxMT2bin[9] = 750.0;
    }

    if(mT2b && High){
      MT2Min = 125;

      xMT2bin[1] = 150.0;
      xMT2bin[2] = 180.0;
      xMT2bin[3] = 260.0;
      xMT2bin[4] = 2000.0;
      xMT2bin[5] = 3000.0;

      newxMT2bin[0] = 0.0;
      newxMT2bin[1] = 25.0;
      newxMT2bin[2] = 50.0;
      newxMT2bin[3] = 75.0;
      newxMT2bin[4] = 100.0;
      newxMT2bin[5] = 125.0;
      newxMT2bin[6] = 150.0;
      newxMT2bin[7] = 180.0;
      newxMT2bin[8] = 260.0;
      newxMT2bin[9] = 750.0;
    }

    if(flag == 27){
      High = 1;
      MT2Min = 125;
      
      xMT2bin[1] = 1500.0;
      xMT2bin[2] = 2000.0;
      xMT2bin[3] = 2500.0;
      xMT2bin[4] = 20000.0;
      xMT2bin[5] = 30000.0;

      newxMT2bin[0] = 0.0;
      newxMT2bin[1] = 25.0;
      newxMT2bin[2] = 50.0;
      newxMT2bin[3] = 75.0;
      newxMT2bin[4] = 125.0;
      newxMT2bin[5] = 1500.0;
      newxMT2bin[6] = 2000.0;
      newxMT2bin[7] = 2500.0;
      newxMT2bin[8] = 20000.0;
      newxMT2bin[9] = 30000.0;
    }

    xMT2bin[0] = MT2Min;

}






void MassPlotter::AddOverAndUnderFlow(TH1 * Histo, bool overflow, bool underflow){
  // Add underflow & overflow bins
  // This failed for older ROOT version when the first(last) bin is empty
  // and there are underflow (overflow) events --- must check whether this 
  // is still the case
  if(underflow){
  Histo->SetBinContent(1,
		       Histo->GetBinContent(0) + Histo->GetBinContent(1));
  Histo->SetBinError(1,
		     sqrt(Histo->GetBinError(0)*Histo->GetBinError(0)+
			  Histo->GetBinError(1)*Histo->GetBinError(1) ));
  Histo->SetBinContent(0, 0.0);
  } if(overflow){
  Histo->SetBinContent(Histo->GetNbinsX(),
		       Histo->GetBinContent(Histo->GetNbinsX()  )+ 
		       Histo->GetBinContent(Histo->GetNbinsX()+1) );
  Histo->SetBinError(Histo->GetNbinsX(),
		     sqrt(Histo->GetBinError(Histo->GetNbinsX()  )*
			  Histo->GetBinError(Histo->GetNbinsX()  )+
			  Histo->GetBinError(Histo->GetNbinsX()+1)*
			  Histo->GetBinError(Histo->GetNbinsX()+1)  ));
  Histo->SetBinContent(Histo->GetNbinsX() + 1, 0.0);
  }
}

void MassPlotter::Cout(int k, TH1F * Histo){
  cout<<k<<" "<<Histo->GetName()<<endl;
  for(int j = 1; j <= Histo->GetNbinsX(); j++)
    cout<<"Bin :: "<<j<<" "<<Histo->GetBinContent(j)<<" +/- "<<Histo->GetBinError(j)<<endl;
}

void MassPlotter::Cout(int k, TH2F * Histo){
  cout<<k<<" "<<Histo->GetName()<<endl;
  int nYbins = Histo->GetNbinsY();
  int nXbins = Histo->GetNbinsX();

  for(int i = 1; i <= nYbins; i++){
    for(int j = 1; j <= nXbins; j++){
      cout<<Histo->GetBinContent((nXbins + 2) * i + j)<<" +/- "<<Histo->GetBinError((nXbins + 2) * i + j)<<endl;
    }
  }
}

void MassPlotter::printYield(){

    cout<<"           & "<<MT2[0]->GetName()     << " & "<<MT2[2]->GetName()      <<" & "<<MT2[1]->GetName()      <<" & "<<MT2[3]->GetName()     <<" & "<<MT2[4]->GetName()<<"& "<<MT2[5]->GetName()     <<"\\"<<endl;
  
    int tmpxMT2bin[10];
  
    for(int i = 0; i < NumberOfMT2Bins; i++)
      tmpxMT2bin[i] = xMT2bin[i]/5;

    double err;

    MT2[4]->IntegralAndError(0,tmpxMT2bin[NumberOfMT2Bins],err);
    cout<<" full selection & "<<MT2[0]->Integral()     << " & "<<MT2[2]->Integral()      <<" & "<<MT2[1]->Integral()      <<" & "<<MT2[3]->Integral()     <<" & "<<MT2[4]->Integral(0,tmpxMT2bin[NumberOfMT2Bins])<<" $pm$ "<<err<<" & "<<MT2[5]->Integral()     <<"\\"<<endl;


    for(int i = 0; i < NumberOfMT2Bins-1; i++){
      MT2[4]->IntegralAndError(tmpxMT2bin[i]+1,tmpxMT2bin[i+1],err);
      cout<<"$"<<xMT2bin[i]<<"-"<<xMT2bin[i+1]<<"$ &      "<<MT2[0]->Integral(tmpxMT2bin[i]+1,tmpxMT2bin[i+1])<< " & "<<MT2[2]->Integral(tmpxMT2bin[i]+1,tmpxMT2bin[i+1]) <<" & "<<MT2[1]->Integral(tmpxMT2bin[i]+1,tmpxMT2bin[i+1]) <<" & "<<MT2[3]->Integral(tmpxMT2bin[i]+1,tmpxMT2bin[i+1])<<" & "<<MT2[4]->Integral(tmpxMT2bin[i]+1,tmpxMT2bin[i+1])<<" $pm$ "<<err<<" & "<<MT2[5]->Integral(tmpxMT2bin[i]+1,tmpxMT2bin[i+1])<<"\\"<<endl;
    }
  }




void MassPlotter::TriggerEfficiency(int sixjet, int quadjet, int plusNJets, Long64_t nevents){

  int SixJet = sixjet;
  int QuadJet = quadjet;
    
 TH1F *hQuadJet50_v5 = new TH1F("HLT_QuadJet50_v5", "HLT_QuadJet50_v5", 500,0,250);
  TH1F *hQuadJet50_v5_PFHT350_v4567 = new TH1F("HLT_QuadJet50_v5_PFHT350_v4567", "HLT_QuadJet50_v5_PFHT350_v4567", 500,0,250);


  TH1F *hSixJet35_v1 = new TH1F("HLT_SixJet35_v1", "HLT_SixJet35_v1", 250,0,250);
  TH1F *hSixJet45_v1 = new TH1F("HLT_SixJet45_v1", "HLT_SixJet45_v1", 250,0,250);

  TH1F *hSixJet35_v2 = new TH1F("HLT_SixJet35_v2", "HLT_SixJet35_v2", 250,0,250);
  TH1F *hSixJet45_v2 = new TH1F("HLT_SixJet45_v2", "HLT_SixJet45_v2", 250,0,250);

  TH1F *hSixJet35_v3 = new TH1F("HLT_SixJet35_v3", "HLT_SixJet35_v3", 250,0,250);
  TH1F *hSixJet45_v3 = new TH1F("HLT_SixJet45_v3", "HLT_SixJet45_v3", 250,0,250);

  TH1F *hQuadJet70_v1 = new TH1F("HLT_QuadJet70_v1", "HLT_QuadJet70_v1", 250,0,250);
  TH1F *hQuadJet80_v1 = new TH1F("HLT_QuadJet80_v1", "HLT_QuadJet80_v1", 250,0,250);

  TH1F *hQuadJet70_v2 = new TH1F("HLT_QuadJet70_v2", "HLT_QuadJet70_v2", 250,0,250);
  TH1F *hQuadJet80_v2 = new TH1F("HLT_QuadJet80_v2", "HLT_QuadJet80_v2", 250,0,250);

  TH1F *hQuadJet70_v3 = new TH1F("HLT_QuadJet70_v3", "HLT_QuadJet70_v3", 250,0,250);
  TH1F *hQuadJet80_v3 = new TH1F("HLT_QuadJet80_v3", "HLT_QuadJet80_v3", 250,0,250);

  TH1F *hHT350V3 = new TH1F("HLT_PFHT350_v3", "HLT_PFHT350_v3", 1000,0,1000);
  TH1F *hHT650V5 = new TH1F("HLT_PFHT650_v5", "HLT_PFHT650_v5", 1000,0,1000);
 
  TH1F *hHT350V4 = new TH1F("HLT_PFHT350_v4", "HLT_PFHT350_v4", 1000,0,1000);
  TH1F *hHT650V6 = new TH1F("HLT_PFHT650_v6", "HLT_PFHT650_v6", 1000,0,1000);
 
  TH1F *hHT350V5 = new TH1F("HLT_PFHT350_v5", "HLT_PFHT350_v5", 1000,0,1000);
  TH1F *hHT650V7 = new TH1F("HLT_PFHT650_v7", "HLT_PFHT650_v7", 1000,0,1000);
 
  TH1F *hHT350V6 = new TH1F("HLT_PFHT350_v6", "HLT_PFHT350_v6", 1000,0,1000);
  TH1F *hHT650V8 = new TH1F("HLT_PFHT650_v8", "HLT_PFHT650_v8", 1000,0,1000);

  TH1F *hHT350V7 = new TH1F("HLT_PFHT350_v7", "HLT_PFHT350_v7", 1000,0,1000);
  TH1F *hHT650V9 = new TH1F("HLT_PFHT650_v9", "HLT_PFHT650_v9", 1000,0,1000);

  TH1F *hSixJet45V1_35V1 = new TH1F("hSixJet45V1_35V1","", 250,0,250);
  TH1F *hSixJet45V1_35V2 = new TH1F("hSixJet45V1_35V2","", 250,0,250);
  TH1F *hSixJet45V1_35V3 = new TH1F("hSixJet45V1_35V3","", 250,0,250);

  TH1F *hSixJet45V2_35V1 = new TH1F("hSixJet45V2_35V1","", 250,0,250);
  TH1F *hSixJet45V2_35V2 = new TH1F("hSixJet45V2_35V2","", 250,0,250);
  TH1F *hSixJet45V2_35V3 = new TH1F("hSixJet45V2_35V3","", 250,0,250);

  TH1F *hSixJet45V3_35V1 = new TH1F("hSixJet45V3_35V1","", 250,0,250);
  TH1F *hSixJet45V3_35V2 = new TH1F("hSixJet45V3_35V2","", 250,0,250);
  TH1F *hSixJet45V3_35V3 = new TH1F("hSixJet45V3_35V3","", 250,0,250);

  TH1F *hQuadJet80V1_70V1 = new TH1F("hQuadJet80V1_70V1","", 250,0,250);
  TH1F *hQuadJet80V1_70V2 = new TH1F("hQuadJet80V1_70V2","", 250,0,250);
  TH1F *hQuadJet80V1_70V3 = new TH1F("hQuadJet80V1_70V3","", 250,0,250);

  TH1F *hQuadJet80V2_70V1 = new TH1F("hQuadJet80V2_70V1","", 250,0,250);
  TH1F *hQuadJet80V2_70V2 = new TH1F("hQuadJet80V2_70V2","", 250,0,250);
  TH1F *hQuadJet80V2_70V3 = new TH1F("hQuadJet80V2_70V3","", 250,0,250);

  TH1F *hQuadJet80V3_70V1 = new TH1F("hQuadJet80V3_70V1","", 250,0,250);
  TH1F *hQuadJet80V3_70V2 = new TH1F("hQuadJet80V3_70V2","", 250,0,250);
  TH1F *hQuadJet80V3_70V3 = new TH1F("hQuadJet80V3_70V3","", 250,0,250);

  TH1F *h650V5_350V3 = new TH1F("HLT_PFHT650_v5_HLT_PFHT350_v3","", 1000,0,1000);
  TH1F *h650V5_350V4 = new TH1F("HLT_PFHT650_v5_HLT_PFHT350_v4","", 1000,0,1000);
  TH1F *h650V5_350V5 = new TH1F("HLT_PFHT650_v5_HLT_PFHT350_v5","", 1000,0,1000);
  TH1F *h650V5_350V6 = new TH1F("HLT_PFHT650_v5_HLT_PFHT350_v6","", 1000,0,1000);
  TH1F *h650V5_350V7 = new TH1F("HLT_PFHT650_v5_HLT_PFHT350_v7","", 1000,0,1000);
  TH1F *h650V6_350V3 = new TH1F("HLT_PFHT650_v6_HLT_PFHT350_v3","", 1000,0,1000);
  TH1F *h650V6_350V4 = new TH1F("HLT_PFHT650_v6_HLT_PFHT350_v4","", 1000,0,1000);
  TH1F *h650V6_350V5 = new TH1F("HLT_PFHT650_v6_HLT_PFHT350_v5","", 1000,0,1000);
  TH1F *h650V6_350V6 = new TH1F("HLT_PFHT650_v6_HLT_PFHT350_v6","", 1000,0,1000);
  TH1F *h650V6_350V7 = new TH1F("HLT_PFHT650_v6_HLT_PFHT350_v7","", 1000,0,1000);
  TH1F *h650V7_350V3 = new TH1F("HLT_PFHT650_v7_HLT_PFHT350_v3","", 1000,0,1000);
  TH1F *h650V7_350V4 = new TH1F("HLT_PFHT650_v7_HLT_PFHT350_v4","", 1000,0,1000);
  TH1F *h650V7_350V5 = new TH1F("HLT_PFHT650_v7_HLT_PFHT350_v5","", 1000,0,1000);
  TH1F *h650V7_350V6 = new TH1F("HLT_PFHT650_v7_HLT_PFHT350_v6","", 1000,0,1000);
  TH1F *h650V7_350V7 = new TH1F("HLT_PFHT650_v7_HLT_PFHT350_v7","", 1000,0,1000);
  TH1F *h650V8_350V3 = new TH1F("HLT_PFHT650_v8_HLT_PFHT350_v3","", 1000,0,1000);
  TH1F *h650V8_350V4 = new TH1F("HLT_PFHT650_v8_HLT_PFHT350_v4","", 1000,0,1000);
  TH1F *h650V8_350V5 = new TH1F("HLT_PFHT650_v8_HLT_PFHT350_v5","", 1000,0,1000);
  TH1F *h650V8_350V6 = new TH1F("HLT_PFHT650_v8_HLT_PFHT350_v6","", 1000,0,1000);
  TH1F *h650V8_350V7 = new TH1F("HLT_PFHT650_v8_HLT_PFHT350_v7","", 1000,0,1000);
  TH1F *h650V9_350V3 = new TH1F("HLT_PFHT650_v9_HLT_PFHT350_v3","", 1000,0,1000);
  TH1F *h650V9_350V4 = new TH1F("HLT_PFHT650_v9_HLT_PFHT350_v4","", 1000,0,1000);
  TH1F *h650V9_350V5 = new TH1F("HLT_PFHT650_v9_HLT_PFHT350_v5","", 1000,0,1000);
  TH1F *h650V9_350V6 = new TH1F("HLT_PFHT650_v9_HLT_PFHT350_v6","", 1000,0,1000);
  TH1F *h650V9_350V7 = new TH1F("HLT_PFHT650_v9_HLT_PFHT350_v7","", 1000,0,1000);

  TH1F *hRunSixJet45V1_35V1 = new TH1F("RunHLT_SixJet45_v1_HLT_SixJet35_v1","", 10000,190000,200000);
  TH1F *hRunSixJet45V2_35V2 = new TH1F("RunHLT_SixJet45_v2_HLT_SixJet35_v2","", 10000,190000,200000);
  TH1F *hRunSixJet45V3_35V3 = new TH1F("RunHLT_SixJet45_v3_HLT_SixJet35_v3","", 10000,190000,200000);

  TH1F *hRunQuadJet80V1_70V1 = new TH1F("RunHLT_QuadJet80_v1_HLT_QuadJet70_v1","", 10000,190000,200000);
  TH1F *hRunQuadJet80V2_70V2 = new TH1F("RunHLT_QuadJet80_v2_HLT_QuadJet70_v2","", 10000,190000,200000);
  TH1F *hRunQuadJet80V3_70V3 = new TH1F("RunHLT_QuadJet80_v3_HLT_QuadJet70_v3","", 10000,190000,200000);

  TH1F *hRun650V5_350V3 = new TH1F("RunHLT_PFHT650_v5_HLT_PFHT350_v3","", 10000,190000,200000);
  TH1F *hRun650V6_350V4 = new TH1F("RunHLT_PFHT650_v6_HLT_PFHT350_v4","", 10000,190000,200000);
  TH1F *hRun650V7_350V5 = new TH1F("RunHLT_PFHT650_v7_HLT_PFHT350_v5","", 10000,190000,200000);
  TH1F *hRun650V8_350V6 = new TH1F("RunHLT_PFHT650_v8_HLT_PFHT350_v6","", 10000,190000,200000);
  TH1F *hRun650V9_350V7 = new TH1F("RunHLT_PFHT650_v9_HLT_PFHT350_v7","", 10000,190000,200000);

  TH1F *hQuadJetPt80 = new TH1F("QuadJetPt80","", 75, 0.0, 1.5);
  TH1F *hQuadJetPt85 = new TH1F("QuadJetPt85","", 75, 0.0, 1.5);
  TH1F *hQuadJetPt90 = new TH1F("QuadJetPt90","", 75, 0.0, 1.5);
  TH1F *hQuadJetPt95 = new TH1F("QuadJetPt95","", 75, 0.0, 1.5);
  TH1F *hQuadJetPt100 = new TH1F("QuadJetPt100","", 75, 0.0, 1.5);
  TH1F *hQuadJetPt105 = new TH1F("QuadJetPt105","", 75, 0.0, 1.5);
  TH1F *hQuadJetPt110 = new TH1F("QuadJetPt110","", 75, 0.0, 1.5);
  TH1F *hQuadJetPt115 = new TH1F("QuadJetPt115","", 75, 0.0, 1.5);
    
  TH1F *hSixJetPt50 = new TH1F("SixJetPt50","", 75, 0.0, 1.5);
  TH1F *hSixJetPt55 = new TH1F("SixJetPt55","", 75, 0.0, 1.5);
  TH1F *hSixJetPt60 = new TH1F("SixJetPt60","", 75, 0.0, 1.5);
  TH1F *hSixJetPt65 = new TH1F("SixJetPt65","", 75, 0.0, 1.5);
  TH1F *hSixJetPt70 = new TH1F("SixJetPt70","", 75, 0.0, 1.5);
  TH1F *hSixJetPt75 = new TH1F("SixJetPt75","", 75, 0.0, 1.5);
  TH1F *hSixJetPt80 = new TH1F("SixJetPt80","", 75, 0.0, 1.5);
  TH1F *hSixJetPt85 = new TH1F("SixJetPt85","", 75, 0.0, 1.5);
 
  TH1F *hHT650 = new TH1F("HT650","", 75, 0.0, 1.5);
  TH1F *hHT675 = new TH1F("HT675","", 75, 0.0, 1.5);
  TH1F *hHT700 = new TH1F("HT700","", 75, 0.0, 1.5);
  TH1F *hHT725 = new TH1F("HT725","", 75, 0.0, 1.5);
  TH1F *hHT750 = new TH1F("HT750","", 75, 0.0, 1.5);
  TH1F *hHT775 = new TH1F("HT775","", 75, 0.0, 1.5);
  TH1F *hHT800 = new TH1F("HT800","", 75, 0.0, 1.5);
  TH1F *hHT825 = new TH1F("HT825","", 75, 0.0, 1.5);

  for(int i = 0; i < fSamples.size(); i++){
    sample Sample = fSamples[i];
    // if(Sample.name != "Run2012A")
    if(Sample.type != "data" || Sample.name != "Run2012D1F")
      continue;
 
    fMT2tree = new MT2tree();

    Sample.tree->SetBranchAddress("MT2tree", &fMT2tree);
    Long64_t nentries =  Sample.tree->GetEntries();

    std::cout << setfill('=') << std::setw(70) << "" << std::endl;
    cout << "looping over :     " <<endl;	
    cout << "   Name:           " << Sample.name << endl;
    cout << "   File:           " << (Sample.file)->GetName() << endl;
    cout << "   Events:         " << Sample.nevents  << endl;
    cout << "   Events in tree: " << Sample.tree->GetEntries() << endl; 
    cout << "   Xsection:       " << Sample.xsection << endl;
    cout << "   kfactor:        " << Sample.kfact << endl;
    cout << "   avg PU weight:  " << Sample.PU_avg_weight << endl;
    std::cout << setfill('-') << std::setw(70) << "" << std::endl;

 
    for (Long64_t jentry=0; jentry<min(nentries, nevents);jentry++) {
      Sample.tree->GetEntry(jentry); 
      
      if ( fVerbose>2 && jentry % 100000 == 0 ){  
	fprintf(stdout, "\rProcessed events: %6d of %6d ", jentry + 1, nentries);
	fflush(stdout);}

      if(fMT2tree->NJets > (SixJet+plusNJets)){
	if(fMT2tree->trigger.HLT_SixJet35_v1 && fMT2tree->trigger.HLT_SixJet45_v1)
	  hSixJet45V1_35V1->Fill(fMT2tree->jet[SixJet].lv.Pt());
	if(fMT2tree->trigger.HLT_SixJet35_v2 && fMT2tree->trigger.HLT_SixJet45_v1)
	  hSixJet45V1_35V2->Fill(fMT2tree->jet[SixJet].lv.Pt());
	if(fMT2tree->trigger.HLT_SixJet35_v3 && fMT2tree->trigger.HLT_SixJet45_v1)
	  hSixJet45V1_35V3->Fill(fMT2tree->jet[SixJet].lv.Pt());
	if(fMT2tree->trigger.HLT_SixJet35_v1 && fMT2tree->trigger.HLT_SixJet45_v2)
	  hSixJet45V2_35V1->Fill(fMT2tree->jet[SixJet].lv.Pt());
	if(fMT2tree->trigger.HLT_SixJet35_v2 && fMT2tree->trigger.HLT_SixJet45_v2)
	  hSixJet45V2_35V2->Fill(fMT2tree->jet[SixJet].lv.Pt());
	if(fMT2tree->trigger.HLT_SixJet35_v3 && fMT2tree->trigger.HLT_SixJet45_v2)
	  hSixJet45V2_35V3->Fill(fMT2tree->jet[SixJet].lv.Pt());
	if(fMT2tree->trigger.HLT_SixJet35_v1 && fMT2tree->trigger.HLT_SixJet45_v3)
	  hSixJet45V3_35V1->Fill(fMT2tree->jet[SixJet].lv.Pt());
	if(fMT2tree->trigger.HLT_SixJet35_v2 && fMT2tree->trigger.HLT_SixJet45_v3)
	  hSixJet45V3_35V2->Fill(fMT2tree->jet[SixJet].lv.Pt());
	if(fMT2tree->trigger.HLT_SixJet35_v3 && fMT2tree->trigger.HLT_SixJet45_v3)
	  hSixJet45V3_35V3->Fill(fMT2tree->jet[SixJet].lv.Pt());
      }

      if(fMT2tree->NJets > (QuadJet+plusNJets)){
	if(fMT2tree->trigger.HLT_QuadJet70_v1 && fMT2tree->trigger.HLT_QuadJet80_v1)
	  hQuadJet80V1_70V1->Fill(fMT2tree->jet[QuadJet].lv.Pt());
	if(fMT2tree->trigger.HLT_QuadJet70_v2 && fMT2tree->trigger.HLT_QuadJet80_v1)
	  hQuadJet80V1_70V2->Fill(fMT2tree->jet[QuadJet].lv.Pt());
	if(fMT2tree->trigger.HLT_QuadJet70_v3 && fMT2tree->trigger.HLT_QuadJet80_v1)
	  hQuadJet80V1_70V3->Fill(fMT2tree->jet[QuadJet].lv.Pt());
	if(fMT2tree->trigger.HLT_QuadJet70_v1 && fMT2tree->trigger.HLT_QuadJet80_v2)
	  hQuadJet80V2_70V1->Fill(fMT2tree->jet[QuadJet].lv.Pt());
	if(fMT2tree->trigger.HLT_QuadJet70_v2 && fMT2tree->trigger.HLT_QuadJet80_v2)
	  hQuadJet80V2_70V2->Fill(fMT2tree->jet[QuadJet].lv.Pt());
	if(fMT2tree->trigger.HLT_QuadJet70_v3 && fMT2tree->trigger.HLT_QuadJet80_v2)
	  hQuadJet80V2_70V3->Fill(fMT2tree->jet[QuadJet].lv.Pt());
	if(fMT2tree->trigger.HLT_QuadJet70_v1 && fMT2tree->trigger.HLT_QuadJet80_v3)
	  hQuadJet80V3_70V1->Fill(fMT2tree->jet[QuadJet].lv.Pt());
	if(fMT2tree->trigger.HLT_QuadJet70_v2 && fMT2tree->trigger.HLT_QuadJet80_v3)
	  hQuadJet80V3_70V2->Fill(fMT2tree->jet[QuadJet].lv.Pt());
	if(fMT2tree->trigger.HLT_QuadJet70_v3 && fMT2tree->trigger.HLT_QuadJet80_v3)
	  hQuadJet80V3_70V3->Fill(fMT2tree->jet[QuadJet].lv.Pt());
      }


      if(fMT2tree->trigger.HLT_PFHT350_v3 && fMT2tree->trigger.HLT_PFHT650_v5)
     	h650V5_350V3->Fill(fMT2tree->misc.HT);
      if(fMT2tree->trigger.HLT_PFHT350_v4 && fMT2tree->trigger.HLT_PFHT650_v5)
     	h650V5_350V4->Fill(fMT2tree->misc.HT);
      if(fMT2tree->trigger.HLT_PFHT350_v5 && fMT2tree->trigger.HLT_PFHT650_v5)
     	h650V5_350V5->Fill(fMT2tree->misc.HT);
      if(fMT2tree->trigger.HLT_PFHT350_v6 && fMT2tree->trigger.HLT_PFHT650_v5)
     	h650V5_350V6->Fill(fMT2tree->misc.HT);
      if(fMT2tree->trigger.HLT_PFHT350_v7 && fMT2tree->trigger.HLT_PFHT650_v5)
     	h650V5_350V7->Fill(fMT2tree->misc.HT);
   
      if(fMT2tree->trigger.HLT_PFHT350_v3 && fMT2tree->trigger.HLT_PFHT650_v6)
     	h650V6_350V3->Fill(fMT2tree->misc.HT);
      if(fMT2tree->trigger.HLT_PFHT350_v4 && fMT2tree->trigger.HLT_PFHT650_v6)
     	h650V6_350V4->Fill(fMT2tree->misc.HT);
      if(fMT2tree->trigger.HLT_PFHT350_v5 && fMT2tree->trigger.HLT_PFHT650_v6)
     	h650V6_350V5->Fill(fMT2tree->misc.HT);
      if(fMT2tree->trigger.HLT_PFHT350_v6 && fMT2tree->trigger.HLT_PFHT650_v6)
     	h650V6_350V6->Fill(fMT2tree->misc.HT);
      if(fMT2tree->trigger.HLT_PFHT350_v7 && fMT2tree->trigger.HLT_PFHT650_v6)
     	h650V6_350V7->Fill(fMT2tree->misc.HT);
   
     if(fMT2tree->trigger.HLT_PFHT350_v3 && fMT2tree->trigger.HLT_PFHT650_v7)
     	h650V7_350V3->Fill(fMT2tree->misc.HT);
      if(fMT2tree->trigger.HLT_PFHT350_v4 && fMT2tree->trigger.HLT_PFHT650_v7)
     	h650V7_350V4->Fill(fMT2tree->misc.HT);
      if(fMT2tree->trigger.HLT_PFHT350_v5 && fMT2tree->trigger.HLT_PFHT650_v7)
     	h650V7_350V5->Fill(fMT2tree->misc.HT);
      if(fMT2tree->trigger.HLT_PFHT350_v6 && fMT2tree->trigger.HLT_PFHT650_v7)
     	h650V7_350V6->Fill(fMT2tree->misc.HT);
      if(fMT2tree->trigger.HLT_PFHT350_v7 && fMT2tree->trigger.HLT_PFHT650_v7)
     	h650V7_350V7->Fill(fMT2tree->misc.HT);
 
      if(fMT2tree->trigger.HLT_PFHT350_v3 && fMT2tree->trigger.HLT_PFHT650_v8)
     	h650V8_350V3->Fill(fMT2tree->misc.HT);
      if(fMT2tree->trigger.HLT_PFHT350_v4 && fMT2tree->trigger.HLT_PFHT650_v8)
     	h650V8_350V4->Fill(fMT2tree->misc.HT);
      if(fMT2tree->trigger.HLT_PFHT350_v5 && fMT2tree->trigger.HLT_PFHT650_v8)
     	h650V8_350V5->Fill(fMT2tree->misc.HT);
      if(fMT2tree->trigger.HLT_PFHT350_v6 && fMT2tree->trigger.HLT_PFHT650_v8)
      	h650V8_350V6->Fill(fMT2tree->misc.HT);
      if(fMT2tree->trigger.HLT_PFHT350_v7 && fMT2tree->trigger.HLT_PFHT650_v8)
      	h650V8_350V7->Fill(fMT2tree->misc.HT);
    
      if(fMT2tree->trigger.HLT_PFHT350_v3 && fMT2tree->trigger.HLT_PFHT650_v9)
     	h650V9_350V3->Fill(fMT2tree->misc.HT);
      if(fMT2tree->trigger.HLT_PFHT350_v4 && fMT2tree->trigger.HLT_PFHT650_v9)
     	h650V9_350V4->Fill(fMT2tree->misc.HT);
      if(fMT2tree->trigger.HLT_PFHT350_v5 && fMT2tree->trigger.HLT_PFHT650_v9)
     	h650V9_350V5->Fill(fMT2tree->misc.HT);
      if(fMT2tree->trigger.HLT_PFHT350_v6 && fMT2tree->trigger.HLT_PFHT650_v9)
      	h650V9_350V6->Fill(fMT2tree->misc.HT);
      if(fMT2tree->trigger.HLT_PFHT350_v7 && fMT2tree->trigger.HLT_PFHT650_v9)
      	h650V9_350V7->Fill(fMT2tree->misc.HT);
    
  
      if(fMT2tree->NJets > (QuadJet+plusNJets)){
        //if(fMT2tree->trigger.HLT_PFHT350_v4 || fMT2tree->trigger.HLT_PFHT350_v5 || fMT2tree->trigger.HLT_PFHT350_v6 || fMT2tree->trigger.HLT_PFHT350_v7){
        //if(fMT2tree->trigger.HLT_PFHT650_v6 || fMT2tree->trigger.HLT_PFHT650_v7 || fMT2tree->trigger.HLT_PFHT650_v8 || fMT2tree->trigger.HLT_PFHT650_v9){
	//if(fMT2tree->trigger.HLT_Mu40_HT200_v1 || fMT2tree->trigger.HLT_Mu40_HT200_v2){
	//if(fMT2tree->trigger.HLT_QuadJet45_v1 && fMT2tree->misc.Run > 203777){
	//if(fMT2tree->trigger.HLT_DiPFJetAve40_v6 && fMT2tree->misc.Run > 194270){
	//if(fMT2tree->trigger.HLT_IsoMu24_v17){
	//if(fMT2tree->trigger.HLT_IsoMu24_eta2p1_v15){

	    hQuadJet50_v5_PFHT350_v4567->Fill(fMT2tree->jet[QuadJet].lv.Pt());
	    
	    if(fMT2tree->trigger.HLT_QuadJet50_v1 || fMT2tree->trigger.HLT_QuadJet50_v2  || fMT2tree->trigger.HLT_QuadJet50_v3|| fMT2tree->trigger.HLT_QuadJet50_v5)
	      hQuadJet50_v5->Fill(fMT2tree->jet[QuadJet].lv.Pt());
      }//}

      int Ref = 0;
      int Tri = 0;
 
      if(fMT2tree->trigger.HLT_PFHT350_v3){
  	Ref = 1;
 	if(fMT2tree->trigger.HLT_PFHT650_v5){
  	  Tri = 1;
	  hRun650V5_350V3->Fill(fMT2tree->misc.Run);
	}
      }
      hHT350V3->Fill(fMT2tree->misc.HT, Ref);
      hHT650V5->Fill(fMT2tree->misc.HT, Tri);
      
      Ref = 0;
      Tri = 0;
      
      if((fMT2tree->trigger.HLT_PFHT350_v4)){
  	Ref = 1;
 	if((fMT2tree->trigger.HLT_PFHT650_v6)){
  	  Tri = 1;
	  hRun650V6_350V4->Fill(fMT2tree->misc.Run);
	}
      }
      hHT350V4->Fill(fMT2tree->misc.HT, Ref);
      hHT650V6->Fill(fMT2tree->misc.HT, Tri);

      Ref = 0;
      Tri = 0;
      
      if((fMT2tree->trigger.HLT_PFHT350_v5)){
  	Ref = 1;
 	if((fMT2tree->trigger.HLT_PFHT650_v7)){
  	  Tri = 1;
	  hRun650V7_350V5->Fill(fMT2tree->misc.Run);
	}
      }
      
      hHT350V5->Fill(fMT2tree->misc.HT, Ref);
      hHT650V7->Fill(fMT2tree->misc.HT, Tri);

      Ref = 0;
      Tri = 0;

      if((fMT2tree->trigger.HLT_PFHT350_v6)){
  	Ref = 1;
 	if((fMT2tree->trigger.HLT_PFHT650_v8)){
  	  Tri = 1;
	  hRun650V8_350V6->Fill(fMT2tree->misc.Run);
	}
      }
      hHT350V6->Fill(fMT2tree->misc.HT, Ref);
      hHT650V8->Fill(fMT2tree->misc.HT, Tri);

      Ref = 0;
      Tri = 0;

      if((fMT2tree->trigger.HLT_PFHT350_v7)){
  	Ref = 1;
 	if((fMT2tree->trigger.HLT_PFHT650_v9)){
  	  Tri = 1;
	  hRun650V9_350V7->Fill(fMT2tree->misc.Run);
	}
      }
      hHT350V7->Fill(fMT2tree->misc.HT, Ref);
      hHT650V9->Fill(fMT2tree->misc.HT, Tri);

      if(fMT2tree->NJets > (SixJet+plusNJets)){
	Ref = 0;
	Tri = 0;
	
	if((fMT2tree->trigger.HLT_SixJet35_v1)){
	  Ref = 1;
	  if((fMT2tree->trigger.HLT_SixJet45_v1)){
	    Tri = 1;
	    hRunSixJet45V1_35V1->Fill(fMT2tree->misc.Run);
	  }
	}
	hSixJet35_v1->Fill(fMT2tree->jet[SixJet].lv.Pt(), Ref);
	hSixJet45_v1->Fill(fMT2tree->jet[SixJet].lv.Pt(), Tri);

	Ref = 0;
	Tri = 0;
	
	if((fMT2tree->trigger.HLT_SixJet35_v2)){
	  Ref = 1;
	  if((fMT2tree->trigger.HLT_SixJet45_v2)){
	    Tri = 1;
	    hRunSixJet45V2_35V2->Fill(fMT2tree->misc.Run);
	  }
	}
	hSixJet35_v2->Fill(fMT2tree->jet[SixJet].lv.Pt(), Ref);
	hSixJet45_v2->Fill(fMT2tree->jet[SixJet].lv.Pt(), Tri);

	Ref = 0;
	Tri = 0;
	
	if((fMT2tree->trigger.HLT_SixJet35_v3)){
	  Ref = 1;
	  if((fMT2tree->trigger.HLT_SixJet45_v3)){
	    Tri = 1;
	    hRunSixJet45V3_35V3->Fill(fMT2tree->misc.Run);
	  }
	}
	hSixJet35_v3->Fill(fMT2tree->jet[SixJet].lv.Pt(), Ref);
	hSixJet45_v3->Fill(fMT2tree->jet[SixJet].lv.Pt(), Tri);
      }

      if(fMT2tree->NJets > (QuadJet+plusNJets)){
	Ref = 0;
	Tri = 0;
	
	if((fMT2tree->trigger.HLT_QuadJet70_v1)){
	  Ref = 1;
	  if((fMT2tree->trigger.HLT_QuadJet80_v1)){
	    Tri = 1;
	    hRunQuadJet80V1_70V1->Fill(fMT2tree->misc.Run);
	  }
	}
	hQuadJet70_v1->Fill(fMT2tree->jet[QuadJet].lv.Pt(), Ref);
	hQuadJet80_v1->Fill(fMT2tree->jet[QuadJet].lv.Pt(), Tri);

	Ref = 0;
	Tri = 0;
	
	if((fMT2tree->trigger.HLT_QuadJet70_v2)){
	  Ref = 1;
	  if((fMT2tree->trigger.HLT_QuadJet80_v2)){
	    Tri = 1;
	    hRunQuadJet80V2_70V2->Fill(fMT2tree->misc.Run);
	  }
	}
	hQuadJet70_v2->Fill(fMT2tree->jet[QuadJet].lv.Pt(), Ref);
	hQuadJet80_v2->Fill(fMT2tree->jet[QuadJet].lv.Pt(), Tri);

	Ref = 0;
	Tri = 0;
	
	if((fMT2tree->trigger.HLT_QuadJet70_v3)){
	  Ref = 1;
	  if((fMT2tree->trigger.HLT_QuadJet80_v3)){
	    Tri = 1;
	    hRunQuadJet80V3_70V3->Fill(fMT2tree->misc.Run);
	  }
	}
	hQuadJet70_v3->Fill(fMT2tree->jet[QuadJet].lv.Pt(), Ref);
	hQuadJet80_v3->Fill(fMT2tree->jet[QuadJet].lv.Pt(), Tri);
      }
    }
  }  
     TCanvas *MyHT = new TCanvas("MyHT","MyHT");
     MyHT->Divide(5,5);
     MyHT->cd(1);
     h650V5_350V3->Draw();
     MyHT->cd(2);
     h650V5_350V4->Draw();
     MyHT->cd(3);
     h650V5_350V5->Draw();
     MyHT->cd(4);
     h650V5_350V6->Draw();
     MyHT->cd(5);
     h650V5_350V7->Draw();
     MyHT->cd(6);
     h650V6_350V3->Draw();
     MyHT->cd(7);
     h650V6_350V4->Draw();
     MyHT->cd(8);
     h650V6_350V5->Draw();
     MyHT->cd(9);
     h650V6_350V6->Draw();
     MyHT->cd(10);
     h650V6_350V7->Draw();
     MyHT->cd(11);
     h650V7_350V3->Draw();
     MyHT->cd(12);
     h650V7_350V4->Draw();
     MyHT->cd(13);
     h650V7_350V5->Draw();
     MyHT->cd(14);
     h650V7_350V6->Draw();
     MyHT->cd(15);
     h650V7_350V7->Draw();
     MyHT->cd(16);
     h650V8_350V3->Draw();
     MyHT->cd(17);
     h650V8_350V4->Draw();
     MyHT->cd(18);
     h650V8_350V5->Draw();
     MyHT->cd(19);
     h650V8_350V6->Draw();
     MyHT->cd(20);
     h650V8_350V7->Draw();
     MyHT->cd(21);
     h650V9_350V3->Draw();
     MyHT->cd(22);
     h650V9_350V4->Draw();
     MyHT->cd(23);
     h650V9_350V5->Draw();
     MyHT->cd(24);
     h650V9_350V6->Draw();
     MyHT->cd(25);
     h650V9_350V7->Draw();


     TCanvas *MyQuadJet = new TCanvas("MyQuadJet","MyQuadJet");
     MyQuadJet->Divide(3,3);
     MyQuadJet->cd(1);
     hQuadJet80V1_70V1->Draw();
     MyQuadJet->cd(2);
     hQuadJet80V1_70V2->Draw();
     MyQuadJet->cd(3);
     hQuadJet80V1_70V3->Draw();
     MyQuadJet->cd(4);
     hQuadJet80V2_70V1->Draw();
     MyQuadJet->cd(5);
     hQuadJet80V2_70V2->Draw();
     MyQuadJet->cd(6);
     hQuadJet80V2_70V3->Draw();
     MyQuadJet->cd(7);
     hQuadJet80V3_70V1->Draw();
     MyQuadJet->cd(8);
     hQuadJet80V3_70V2->Draw();
     MyQuadJet->cd(9);
     hQuadJet80V3_70V3->Draw();

     TCanvas *MySixJet = new TCanvas("MySixJet","MySixJet");
     MySixJet->Divide(3,3);
     MySixJet->cd(1);
     hSixJet45V1_35V1->Draw();
     MySixJet->cd(2);
     hSixJet45V1_35V2->Draw();
     MySixJet->cd(3);
     hSixJet45V1_35V3->Draw();
     MySixJet->cd(4);
     hSixJet45V2_35V1->Draw();
     MySixJet->cd(5);
     hSixJet45V2_35V2->Draw();
     MySixJet->cd(6);
     hSixJet45V2_35V3->Draw();
     MySixJet->cd(7);
     hSixJet45V3_35V1->Draw();
     MySixJet->cd(8);
     hSixJet45V3_35V2->Draw();
     MySixJet->cd(9);
     hSixJet45V3_35V3->Draw();


     TCanvas *MyEff = new TCanvas("MyEff","MyEff");
     MyEff->Divide(3,5);
     MyEff->cd(1);
     hHT350V3->Draw();
     MyEff->cd(2);
     hHT650V5->Draw();
     MyEff->cd(3);
     TGraphAsymmErrors *gEff650V5_350V3 = new TGraphAsymmErrors(hHT650V5, hHT350V3);
     gEff650V5_350V3->Draw("AP");

     MyEff->cd(4);
     hHT350V4->Draw();
     MyEff->cd(5);
     hHT650V6->Draw();
     MyEff->cd(6);
     TGraphAsymmErrors *gEff650V6_350V4 = new TGraphAsymmErrors(hHT650V6, hHT350V4);
     gEff650V6_350V4->Draw("AP");

     MyEff->cd(7);
     hHT350V5->Draw();
     MyEff->cd(8);
     hHT650V7->Draw();
     MyEff->cd(9);
     TGraphAsymmErrors *gEff650V7_350V5 = new TGraphAsymmErrors(hHT650V7, hHT350V5);
     gEff650V7_350V5->Draw("AP");

     MyEff->cd(10);
     hHT350V6->Draw();
     MyEff->cd(11);
     hHT650V8->Draw();
     MyEff->cd(12);
     TGraphAsymmErrors *gEff650V8_350V6 = new TGraphAsymmErrors(hHT650V8, hHT350V6);
     gEff650V8_350V6->Draw("AP");

     MyEff->cd(13);
     hHT350V7->Draw();
     MyEff->cd(14);
     hHT650V9->Draw();
     MyEff->cd(15);
     TGraphAsymmErrors *gEff650V9_350V7 = new TGraphAsymmErrors(hHT650V9, hHT350V7);
     gEff650V9_350V7->Draw("AP");

     for(int i = 651; i < 1001; i++){
       float V5_V3 = hHT650V5->GetBinContent(i)/hHT350V3->GetBinContent(i);
       float V5_V3W = hHT650V5->GetBinContent(i);

       float V6_V4 = hHT650V6->GetBinContent(i)/hHT350V4->GetBinContent(i);
       float V6_V4W = hHT650V6->GetBinContent(i);

       float V7_V5 = hHT650V7->GetBinContent(i)/hHT350V5->GetBinContent(i);
       float V7_V5W = hHT650V7->GetBinContent(i);

       float V8_V6 = hHT650V8->GetBinContent(i)/hHT350V6->GetBinContent(i);
       float V8_V6W = hHT650V8->GetBinContent(i);

       float V9_V7 = hHT650V9->GetBinContent(i)/hHT350V7->GetBinContent(i);
       float V9_V7W = hHT650V9->GetBinContent(i);

       hHT650->Fill(V5_V3, V5_V3W);
       hHT650->Fill(V6_V4, V6_V4W);
       hHT650->Fill(V7_V5, V7_V5W);
       hHT650->Fill(V8_V6, V8_V6W);
       hHT650->Fill(V9_V7, V9_V7W);
       
       if(i > 675){
	 hHT675->Fill(V5_V3, V5_V3W);
	 hHT675->Fill(V6_V4, V6_V4W);
	 hHT675->Fill(V7_V5, V7_V5W);
	 hHT675->Fill(V8_V6, V8_V6W);
	 hHT675->Fill(V9_V7, V9_V7W);
       }

       if(i > 700){
	 hHT700->Fill(V5_V3, V5_V3W);
	 hHT700->Fill(V6_V4, V6_V4W);
	 hHT700->Fill(V7_V5, V7_V5W);
	 hHT700->Fill(V8_V6, V8_V6W);
	 hHT700->Fill(V9_V7, V9_V7W);
       }
       
       if(i > 725){
	 hHT725->Fill(V5_V3, V5_V3W);
	 hHT725->Fill(V6_V4, V6_V4W);
	 hHT725->Fill(V7_V5, V7_V5W);
	 hHT725->Fill(V8_V6, V8_V6W);
	 hHT725->Fill(V9_V7, V9_V7W);
       }
       
       if(i > 750){
	 hHT750->Fill(V5_V3, V5_V3W);
	 hHT750->Fill(V6_V4, V6_V4W);
	 hHT750->Fill(V7_V5, V7_V5W);
	 hHT750->Fill(V8_V6, V8_V6W);
	 hHT750->Fill(V9_V7, V9_V7W);
       }
       
       if(i > 775){
	 hHT775->Fill(V5_V3, V5_V3W);
	 hHT775->Fill(V6_V4, V6_V4W);
	 hHT775->Fill(V7_V5, V7_V5W);
	 hHT775->Fill(V8_V6, V8_V6W);
	 hHT775->Fill(V9_V7, V9_V7W);
       }
       
       if(i > 800){
	 hHT800->Fill(V5_V3, V5_V3W);
	 hHT800->Fill(V6_V4, V6_V4W);
	 hHT800->Fill(V7_V5, V7_V5W);
	 hHT800->Fill(V8_V6, V8_V6W);
	 hHT800->Fill(V9_V7, V9_V7W);
       }
       
       if(i > 825){
	 hHT825->Fill(V5_V3, V5_V3W);
	 hHT825->Fill(V6_V4, V6_V4W);
	 hHT825->Fill(V7_V5, V7_V5W);
	 hHT825->Fill(V8_V6, V8_V6W);
	 hHT825->Fill(V9_V7, V9_V7W);
       }
     }
       
     TCanvas *MyEffHT = new TCanvas("MyEffHT","MyEffHT");
     MyEffHT->Divide(4,2);
     MyEffHT->cd(1);
     hHT650->Draw();
     MyEffHT->cd(2);
     hHT675->Draw();
     MyEffHT->cd(3);
     hHT700->Draw();
     MyEffHT->cd(4);
     hHT725->Draw();
     MyEffHT->cd(5);
     hHT750->Draw();
     MyEffHT->cd(6);
     hHT775->Draw();
     MyEffHT->cd(7);
     hHT800->Draw();
     MyEffHT->cd(8);
     hHT825->Draw();

     TCanvas *MyEff2 = new TCanvas("MyEff2","MyEff2");
     MyEff2->Divide(3,3);
     MyEff2->cd(1);
     hQuadJet70_v1->Draw();
     MyEff2->cd(2);
     hQuadJet80_v1->Draw();
     MyEff2->cd(3);
     TGraphAsymmErrors *gEffQuadJet80V1_70V1 = new TGraphAsymmErrors(hQuadJet80_v1, hQuadJet70_v1);
     gEffQuadJet80V1_70V1->Draw("AP");


     MyEff2->cd(4);
     hQuadJet70_v2->Draw();
     MyEff2->cd(5);
     hQuadJet80_v2->Draw();
     MyEff2->cd(6);
     TGraphAsymmErrors *gEffQuadJet80V2_70V2 = new TGraphAsymmErrors(hQuadJet80_v2, hQuadJet70_v2);
     gEffQuadJet80V2_70V2->Draw("AP");

     MyEff2->cd(7);
     hQuadJet70_v3->Draw();
     MyEff2->cd(8);
     hQuadJet80_v3->Draw();
     MyEff2->cd(9);
     TGraphAsymmErrors *gEffQuadJet80V3_70V3 = new TGraphAsymmErrors(hQuadJet80_v3, hQuadJet70_v3);
     gEffQuadJet80V3_70V3->Draw("AP");

     for(int i = 81; i < 251; i++){
       float V1_V1 = hQuadJet80_v1->GetBinContent(i)/hQuadJet70_v1->GetBinContent(i);
       float V1_V1W = hQuadJet80_v1->GetBinContent(i);

       float V2_V2 = hQuadJet80_v2->GetBinContent(i)/hQuadJet70_v2->GetBinContent(i);
       float V2_V2W = hQuadJet80_v2->GetBinContent(i);

       float V3_V3 = hQuadJet80_v3->GetBinContent(i)/hQuadJet70_v3->GetBinContent(i);
       float V3_V3W = hQuadJet80_v3->GetBinContent(i);

       hQuadJetPt80->Fill(V1_V1, V1_V1W);
       hQuadJetPt80->Fill(V2_V2, V2_V2W);
       hQuadJetPt80->Fill(V3_V3, V3_V3W);

       if(i > 85){
	 hQuadJetPt85->Fill(V1_V1, V1_V1W);
	 hQuadJetPt85->Fill(V2_V2, V2_V2W);
	 hQuadJetPt85->Fill(V3_V3, V3_V3W);
       }

       if(i > 90){
	 hQuadJetPt90->Fill(V1_V1, V1_V1W);
	 hQuadJetPt90->Fill(V2_V2, V2_V2W);
	 hQuadJetPt90->Fill(V3_V3, V3_V3W);
       }

       if(i > 95){
	 hQuadJetPt95->Fill(V1_V1, V1_V1W);
	 hQuadJetPt95->Fill(V2_V2, V2_V2W);
	 hQuadJetPt95->Fill(V3_V3, V3_V3W);
       }
       
       if(i > 100){
	 hQuadJetPt100->Fill(V1_V1, V1_V1W);
	 hQuadJetPt100->Fill(V2_V2, V2_V2W);
	 hQuadJetPt100->Fill(V3_V3, V3_V3W);
       }

       if(i > 105){
	 hQuadJetPt105->Fill(V1_V1, V1_V1W);
	 hQuadJetPt105->Fill(V2_V2, V2_V2W);
	 hQuadJetPt105->Fill(V3_V3, V3_V3W);
       }
      
       if(i > 110){
	 hQuadJetPt110->Fill(V1_V1, V1_V1W);
	 hQuadJetPt110->Fill(V2_V2, V2_V2W);
	 hQuadJetPt110->Fill(V3_V3, V3_V3W);
       }

       if(i > 115){
	 hQuadJetPt115->Fill(V1_V1, V1_V1W);
	 hQuadJetPt115->Fill(V2_V2, V2_V2W);
	 hQuadJetPt115->Fill(V3_V3, V3_V3W);
       }
     }

     TCanvas *MyEffQuad = new TCanvas("MyEffQuad","MyEffQuad");
     MyEffQuad->Divide(4,2);
     MyEffQuad->cd(1);
     hQuadJetPt80->Draw();
     MyEffQuad->cd(2);
     hQuadJetPt85->Draw();
     MyEffQuad->cd(3);
     hQuadJetPt90->Draw();
     MyEffQuad->cd(4);
     hQuadJetPt95->Draw();
     MyEffQuad->cd(5);
     hQuadJetPt100->Draw();
     MyEffQuad->cd(6);
     hQuadJetPt105->Draw();
     MyEffQuad->cd(7);
     hQuadJetPt110->Draw();
     MyEffQuad->cd(8);
     hQuadJetPt115->Draw();

     TCanvas *MyEff3 = new TCanvas("MyEff3","MyEff3");
     MyEff3->Divide(3,3);
     MyEff3->cd(1);
     hSixJet35_v1->Draw();
     MyEff3->cd(2);
     hSixJet45_v1->Draw();
     MyEff3->cd(3);
     TGraphAsymmErrors *gEffSixJet45V1_35V1 = new TGraphAsymmErrors(hSixJet45_v1, hSixJet35_v1);
     gEffSixJet45V1_35V1->Draw("AP");
   
     MyEff3->cd(4);
     hSixJet35_v2->Draw();
     MyEff3->cd(5);
     hSixJet45_v2->Draw();
     MyEff3->cd(6);
     TGraphAsymmErrors *gEffSixJet45V2_35V2 = new TGraphAsymmErrors(hSixJet45_v2, hSixJet35_v2);
     gEffSixJet45V2_35V2->Draw("AP");

     MyEff3->cd(7);
     hSixJet35_v3->Draw();
     MyEff3->cd(8);
     hSixJet45_v3->Draw();
     MyEff3->cd(9);
     TGraphAsymmErrors *gEffSixJet45V3_35V3 = new TGraphAsymmErrors(hSixJet45_v3, hSixJet35_v3);
     gEffSixJet45V3_35V3->Draw("AP");



     for(int i = 51; i < 251; i++){
       float V1_V1 = hSixJet45_v1->GetBinContent(i)/hSixJet35_v1->GetBinContent(i);
       float V1_V1W = hSixJet45_v1->GetBinContent(i);

       float V2_V2 = hSixJet45_v2->GetBinContent(i)/hSixJet35_v2->GetBinContent(i);
       float V2_V2W = hSixJet45_v2->GetBinContent(i);

       float V3_V3 = hSixJet45_v3->GetBinContent(i)/hSixJet35_v3->GetBinContent(i);
       float V3_V3W = hSixJet45_v3->GetBinContent(i);

       hSixJetPt50->Fill(V1_V1, V1_V1W);
       hSixJetPt50->Fill(V2_V2, V2_V2W);
       hSixJetPt50->Fill(V3_V3, V3_V3W);

       if(i > 55){
	 hSixJetPt55->Fill(V1_V1, V1_V1W);
	 hSixJetPt55->Fill(V2_V2, V2_V2W);
	 hSixJetPt55->Fill(V3_V3, V3_V3W);
       }
       
       if(i > 60){
	 hSixJetPt60->Fill(V1_V1, V1_V1W);
	 hSixJetPt60->Fill(V2_V2, V2_V2W);
	 hSixJetPt60->Fill(V3_V3, V3_V3W);
       }

       if(i > 65){
	 hSixJetPt65->Fill(V1_V1, V1_V1W);
	 hSixJetPt65->Fill(V2_V2, V2_V2W);
	 hSixJetPt65->Fill(V3_V3, V3_V3W);
       }
       
       if(i > 70){
	 hSixJetPt70->Fill(V1_V1, V1_V1W);
	 hSixJetPt70->Fill(V2_V2, V2_V2W);
	 hSixJetPt70->Fill(V3_V3, V3_V3W);
       }

      if(i > 75){
	 hSixJetPt75->Fill(V1_V1, V1_V1W);
	 hSixJetPt75->Fill(V2_V2, V2_V2W);
	 hSixJetPt75->Fill(V3_V3, V3_V3W);
       }

       if(i > 80){
	 hSixJetPt80->Fill(V1_V1, V1_V1W);
	 hSixJetPt80->Fill(V2_V2, V2_V2W);
	 hSixJetPt80->Fill(V3_V3, V3_V3W);
       }

      if(i > 85){
	 hSixJetPt85->Fill(V1_V1, V1_V1W);
	 hSixJetPt85->Fill(V2_V2, V2_V2W);
	 hSixJetPt85->Fill(V3_V3, V3_V3W);
       }
     }

     TCanvas *MyEffSix = new TCanvas("MyEffSix","MyEffSix");
     MyEffSix->Divide(4,2);
     MyEffSix->cd(1);
     hSixJetPt50->Draw();
     MyEffSix->cd(2);
     hSixJetPt55->Draw();
     MyEffSix->cd(3);
     hSixJetPt60->Draw();
     MyEffSix->cd(4);
     hSixJetPt65->Draw();
     MyEffSix->cd(5);
     hSixJetPt70->Draw();
     MyEffSix->cd(6);
     hSixJetPt75->Draw();
     MyEffSix->cd(7);
     hSixJetPt80->Draw();
     MyEffSix->cd(8);
     hSixJetPt85->Draw();

     TCanvas *MyCRun = new TCanvas("MyCRun","MyCRun");
     MyCRun->Divide(2,3);
     MyCRun->cd(1);
     hRun650V5_350V3->Draw();
     MyCRun->cd(2);
     hRun650V6_350V4->Draw();
     MyCRun->cd(3);
     hRun650V7_350V5->Draw();
     MyCRun->cd(4);
     hRun650V8_350V6->Draw();
     MyCRun->cd(5);
     hRun650V9_350V7->Draw();

     TCanvas *MyCRun2 = new TCanvas("MyCRun2","MyCRun2");
     MyCRun2->Divide(1,3);
     MyCRun2->cd(1);
     hRunQuadJet80V1_70V1->Draw();
     MyCRun2->cd(2);
     hRunQuadJet80V2_70V2->Draw();
     MyCRun2->cd(3);
     hRunQuadJet80V3_70V3->Draw();

     TCanvas *MyCRun3 = new TCanvas("MyCRun3","MyCRun3");
     MyCRun3->Divide(1,3);
     MyCRun3->cd(1);
     hRunSixJet45V1_35V1->Draw();
     MyCRun3->cd(2);
     hRunSixJet45V2_35V2->Draw();
     MyCRun3->cd(3);
     hRunSixJet45V3_35V3->Draw();


    TCanvas *MyCQuadHT = new TCanvas("MyCQuadHT","MyCQuadHT");
     MyCQuadHT->Divide(4,2);
     MyCQuadHT->cd(1);
     hQuadJet50_v5_PFHT350_v4567->SetLineColor(4);
     TH1F *hQuadJet50_v5_PFHT350_v4567_2 = (TH1F*)hQuadJet50_v5_PFHT350_v4567->Rebin(2, "hQuadJet50_v5_PFHT350_v4567_2");
     TH1F *hQuadJet50_v5_2 = (TH1F*)hQuadJet50_v5->Rebin(2, "hQuadJet50_v5_2");
     hQuadJet50_v5_PFHT350_v4567_2->Draw();
     hQuadJet50_v5_2->Draw("same");

     MyCQuadHT->cd(5);
     TGraphAsymmErrors *gQuadJet50_v5_PFHT350_v4567_2 = new TGraphAsymmErrors(hQuadJet50_v5_2, hQuadJet50_v5_PFHT350_v4567_2);
     gQuadJet50_v5_PFHT350_v4567_2->Draw("AP");

     MyCQuadHT->cd(2);
     TH1F *hQuadJet50_v5_PFHT350_v4567_5 = (TH1F*)hQuadJet50_v5_PFHT350_v4567->Rebin(5, "hQuadJet50_v5_PFHT350_v4567_5");
     TH1F *hQuadJet50_v5_5 = (TH1F*)hQuadJet50_v5->Rebin(5, "hQuadJet50_v5_5");
     hQuadJet50_v5_PFHT350_v4567_5->Draw();
     hQuadJet50_v5_5->Draw("same");

     MyCQuadHT->cd(6);
     TGraphAsymmErrors *gQuadJet50_v5_PFHT350_v4567_5 = new TGraphAsymmErrors(hQuadJet50_v5_5, hQuadJet50_v5_PFHT350_v4567_5);
     gQuadJet50_v5_PFHT350_v4567_5->Draw("AP");

     MyCQuadHT->cd(3);
     TH1F *hQuadJet50_v5_PFHT350_v4567_10 = (TH1F*)hQuadJet50_v5_PFHT350_v4567->Rebin(10, "hQuadJet50_v5_PFHT350_v4567_10");
     TH1F *hQuadJet50_v5_10 = (TH1F*)hQuadJet50_v5->Rebin(10, "hQuadJet50_v5_10");
     hQuadJet50_v5_PFHT350_v4567_10->Draw();
     hQuadJet50_v5_10->Draw("same");

     MyCQuadHT->cd(7);
     TGraphAsymmErrors *gQuadJet50_v5_PFHT350_v4567_10 = new TGraphAsymmErrors(hQuadJet50_v5_10, hQuadJet50_v5_PFHT350_v4567_10);
     gQuadJet50_v5_PFHT350_v4567_10->Draw("AP");


     MyCQuadHT->cd(4);
     TH1F *hQuadJet50_v5_PFHT350_v4567_20 = (TH1F*)hQuadJet50_v5_PFHT350_v4567->Rebin(20, "hQuadJet50_v5_PFHT350_v4567_20");
     TH1F *hQuadJet50_v5_20 = (TH1F*)hQuadJet50_v5->Rebin(20, "hQuadJet50_v5_20");
     hQuadJet50_v5_PFHT350_v4567_20->Draw();
     hQuadJet50_v5_20->Draw("same");

     MyCQuadHT->cd(8);
     TGraphAsymmErrors *gQuadJet50_v5_PFHT350_v4567_20 = new TGraphAsymmErrors(hQuadJet50_v5_20, hQuadJet50_v5_PFHT350_v4567_20);
     gQuadJet50_v5_PFHT350_v4567_20->Draw("AP");
     MyCQuadHT->SaveAs("MyC.C");
}

void MassPlotter::GetGenInfo(MT2GenParticle *T1, MT2GenParticle *T2, MT2GenParticle *W1, MT2GenParticle *W2, MT2GenParticle *B1, MT2GenParticle *B2, MT2GenParticle *U1, MT2GenParticle *U2, MT2GenParticle *D1, MT2GenParticle *D2){
  T1->Reset();
  T2->Reset();
  W1->Reset();
  W2->Reset();
  B1->Reset();
  B2->Reset();
  U1->Reset();
  U2->Reset();
  D1->Reset();
  D2->Reset();
  
  MT2GenParticle *Mo;

  for(int i = 0; i < fMT2tree->NGenParticles; i++){
    // cout<<i<<" Id "<<fMT2tree->genparticle[i].ID<<" Mo "<<fMT2tree->genparticle[i].MIndex<<" GMo "<<fMT2tree->genparticle[i].GMIndex<<endl;  

    Mo  = fMT2tree->genparticle[i].GetMother(fMT2tree->genparticle, fMT2tree->NGenParticles);

    if(abs(fMT2tree->genparticle[i].ID) == 6){
      if(T1->Index == -9)
	*T1 = fMT2tree->genparticle[i];
      else
	*T2 = fMT2tree->genparticle[i];
    }

    if(abs(fMT2tree->genparticle[i].ID) == 24){
      if(Mo->Index == T1->Index)
	*W1 = fMT2tree->genparticle[i];
      else
	*W2 = fMT2tree->genparticle[i];
    }

    if(abs(fMT2tree->genparticle[i].ID) == 5 && abs(Mo->ID) == 6){
      if(Mo->Index == T1->Index)
	*B1 = fMT2tree->genparticle[i];
      else
	*B2 = fMT2tree->genparticle[i];
    }


    if(abs(fMT2tree->genparticle[i].ID) < 6 && abs(Mo->ID) == 24){
      if(Mo->Index == W1->Index){
	if(U1->Index == -9)
	  *U1 =  fMT2tree->genparticle[i];
	else 
	  *D1 =  fMT2tree->genparticle[i];
      }else{
	if(U2->Index == -9)
	  *U2 =  fMT2tree->genparticle[i];
	else 
	  *D2 =  fMT2tree->genparticle[i];
      }}
  }
  // cout<<" T1->ID "<<T1->ID<<endl;
  // cout<<" T2->ID "<<T2->ID<<endl;
  // cout<<" W1->ID "<<W1->ID<<endl;
  // cout<<" W2->ID "<<W2->ID<<endl;
  // cout<<" B1->ID "<<B1->ID<<endl;
  // cout<<" B2->ID "<<B2->ID<<endl;
  // cout<<" U1->ID "<<U1->ID<<endl;
  // cout<<" U2->ID "<<U2->ID<<endl;
  // cout<<" D1->ID "<<D1->ID<<endl;
  // cout<<" D2->ID "<<D2->ID<<endl;
}

float MassPlotter::MatchedJets(MT2GenParticle a, int &b){
  float minDR = 100.0;
  for(int i=0; i<fMT2tree->NJets; ++i){
    if(fMT2tree->jet[i].IsGoodPFJet(20.0, 2.4, 1) ==false) continue;
    
    if(a.lv.DeltaR(fMT2tree->jet[i].lv) < minDR){
      minDR = a.lv.DeltaR(fMT2tree->jet[i].lv);
      b = i;
    }
  }
    return minDR;
}

void MassPlotter::InitializeTopSearch(TopSearch * myTopSearch, int nPart, TLorentzVector lv, int bjet){
	myTopSearch->SetEvt(nPart, 0, 1); 
	myTopSearch->SetEvt(nPart, 1, lv.Px()); 
	myTopSearch->SetEvt(nPart, 2, lv.Py()); 
	myTopSearch->SetEvt(nPart, 3, lv.Pz()); 
	myTopSearch->SetEvt(nPart, 4, lv.E()); 
	myTopSearch->SetEvt(nPart, 5, lv.M()); 
	myTopSearch->SetEvt(nPart, 6, lv.Pt()); 
	myTopSearch->SetEvt(nPart, 7, lv.Eta()); 
	myTopSearch->SetEvt(nPart, 8, lv.Phi()); 
	myTopSearch->SetEvt(nPart, 9, bjet);
}

void MassPlotter::InitializeHemisphereFinder(vector<float> &px, vector<float> &py, vector<float> &pz, vector<float> &E, TLorentzVector lv){
  	px.push_back(lv.Px());
	py.push_back(lv.Py());
	pz.push_back(lv.Pz());
	E .push_back(lv.E());
}

void MassPlotter::TopStudy(TString mySample, Long64_t nevents){

  TH1F* hWMass = new TH1F("W","W",100,0,300);
  TH1F* hTMass = new TH1F("T","T",50,0,300);
  TH2F* hTWMass = new TH2F("TW","TW",50,0,150,50,0,300);
  TH2F* hHemi12 = new TH2F("Hemi12","Hemi12",10,-0.5,9.5,10,-0.5,9.5);
  TH1F* hChi2T = new TH1F("Chi2T", "Chi2T",15,0,5);

  TH1F* hMinDR = new TH1F("MinDR","MinDR",1020,-1,101);
  TH1F* hDPtTop = new TH1F("DPtTop","DPtTop", 10000, -500, 500);

  double x[11] = {0.0, 30.0, 60.0, 90.0, 130.0, 170.0, 250.0, 350.0, 450.0, 1000.0, 1300.0}; 
  double xEta[6] = {0, 0.5, 1.0, 1.5, 2.5, 4.5};
 
  TH1F* hPtTop = new TH1F("PtTop","PtTop", 10, x);
  TH1F* hMT2Top = new TH1F("MT2Top","MT2Top", 10, x);
  TH1F* hEtaTop = new TH1F("EtaTop","EtaTop", 5, xEta);
  TH1F* hPtTopTogether = new TH1F("PtTopTogether","PtTopTogether", 10, x);
  TH1F* hPtTopEff = new TH1F("PtTopEff","PtTopEff", 10, x);
  TH1F* hEtaTopEff = new TH1F("EtaTopEff","EtaTopEff", 5, xEta);
  TH1F* hMT2TopEff = new TH1F("MT2TopEff","MT2TopEff", 10, x);

  TH1F* hPtW = new TH1F("PtW","PtW", 10, x);
  TH1F* hPtWTogether = new TH1F("PtWTogether","PtWTogether", 10, x);

  int NTopsW1b = 0;
  int NTopsW1j = 0;
  int NTopsRightB = 0;
  int NTopsRecTop = 0;


  int NTops0 = 0;

  int Together1 = 0;
  int SplittedW1 = 0;
  int SplittedT1 = 0;
  int NTops1 = 0;

  int Together2 = 0;
  int SplittedW2 = 0;
  int SplittedT2 = 0;
  int NTops2 = 0;

  int together = 0;
  int splittedW = 0;
  int splittedT = 0;
  int nTops = 0;

  int nGenFullHad = 0;
  int TogetherGenFullHad = 0;
  int SplittedBGenFullHad = 0;
  int SplittedWGenFullHad = 0;
  int CorrectTopGenFullHad = 0;

  int TogetherGenTBW = 0;

  int TogetherGenTBQQ = 0;
  int SplittedBGenTBQQ = 0;
  int SplittedWGenTBQQ = 0;

  MT2GenParticle T1, T2, W1, W2, B1, B2, U1, U2, D1, D2;

  for(int ii = 0; ii <(int) fSamples.size(); ii++){
    sample Sample = fSamples[ii];
    // if(Sample.name != "Run2012A")
    if(Sample.name != mySample)
    //   //  if(Sample.sname != mySample)
      continue;

    // debug level for TopSearch:
    // = 0 if no printout at all
    // = 1 if only final statistics
    // = 2 if also short event per event output
    // = 3 if also longer event per event output
    const int debug = 3;
   
    TopSearch *myTopSearch = new TopSearch(); 
    myTopSearch->SetDebug(debug);
    myTopSearch->SetWithbtag(3); // = 1: btag preferred over chisq, = 2: correct b-jet requested for top
    myTopSearch->SetChisqTmax(2.0);
    myTopSearch->SetChisqWmax(2.0);
    myTopSearch->SetChisqTWlepmax(0.0);
    myTopSearch->SetChisqW1b(1.0);//chisqW1bpenlty
    myTopSearch->SetChisqWlept(1.0);//chisqWleptpenlty
    myTopSearch->SetRescaleW(1);// = 1: rescale W 4-vector to correct W mass, = 0: do not rescale

    //Changed from the default values
    myTopSearch->SetCkWoverlap(0);

    fMT2tree = new MT2tree();

    Sample.tree->SetBranchAddress("MT2tree", &fMT2tree);
    Long64_t nentries =  Sample.tree->GetEntries();

    Long64_t minNentriesNevents = min(nentries, nevents);

    std::cout << setfill('=') << std::setw(70) << "" << std::endl;
    cout << "looping over :     " <<endl;	
    cout << "   Name:           " << Sample.name << endl;
    cout << "   File:           " << (Sample.file)->GetName() << endl;
    cout << "   Events:         " << Sample.nevents  << endl;
    cout << "   Events in tree: " << nentries << endl; 
    cout << "   Events to read: " << minNentriesNevents << endl; 
    cout << "   Xsection:       " << Sample.xsection << endl;
    cout << "   kfactor:        " << Sample.kfact << endl;
    cout << "   avg PU weight:  " << Sample.PU_avg_weight << endl;
    std::cout << setfill('-') << std::setw(70) << "" << std::endl;

    for (Long64_t jentry=0; jentry<minNentriesNevents;jentry++) {
     
      Sample.tree->GetEntry(jentry); 
      together = 0;
      splittedW = 0;
      splittedT = 0;
      nTops = 0;
 
      if ( fVerbose>2 && jentry % 100000 == 0 ){ 
	fprintf(stdout, "\rProcessed events: %6d of %6d ", jentry + 1, nentries);
	fflush(stdout);
      }
      // if(fMT2tree->misc.TopDecayMode != 0)
      // 	continue;
      if((fMT2tree->NEles!=0  && fMT2tree->ele[0].lv.Pt() >= 10) || (fMT2tree->NMuons!=0 && fMT2tree->muo[0].lv.Pt() >= 10))
	continue;
      // if(fMT2tree->NJetsIDLoose < 3 )
      //  	continue;
  
     if (!(fMT2tree->misc.MET            >  30                                   &&
          fMT2tree->NBJets               > 0                                     &&
         (fMT2tree->NEles==0  ||  fMT2tree->ele[0].lv.Pt()<10)                   &&
         (fMT2tree->NMuons==0 ||  fMT2tree->muo[0].lv.Pt()<10)                   &&
	  fMT2tree->misc.Jet0Pass        == 1                                    &&
	  fMT2tree->misc.Jet1Pass        == 1                                    &&
	  fMT2tree->misc.LeadingJPt      >  150                                  &&
	  fMT2tree->misc.SecondJPt       >  100                                  &&
	  fMT2tree->misc.PassJetID       == 1                                    &&
	  fMT2tree->misc.Vectorsumpt     <  70                                   &&
         (fMT2tree->misc.MinMetJetDPhi >0.3||fMT2tree->misc.MinMetJetDPhiIndex>3) &&
	 ((fMT2tree->NJetsIDLoose40  >=4 && fMT2tree->jet[3].lv.Pt() >= 120.0) || (fMT2tree->NJetsIDLoose40  >=6 && fMT2tree->jet[5].lv.Pt() >= 80.0))
	   ))
	  continue;


     // if(fMT2tree->misc.MT2 < 170.0)
     //   continue;

     TLorentzVector w1, t1, w2, t2;

     // fill Pseudojets with selected objects
     vector<float> px, py, pz, E;
 
     int nPart = 0;
 
     // struct to keep track of indices in hemi association
     struct HemiObj {vector<TString> type; vector<int> index; vector<int> hemisphere;
     } hemiobjs;
      

     GetGenInfo(&T1, &T2, &W1, &W2, &B1, &B2, &U1, &U2, &D1, &D2);

     if(B1.Index == B2.Index){
       //Wjets
       for(int i=0; i<fMT2tree->NJets; ++i){
       if(fMT2tree->jet[i].IsGoodPFJet(20.0, 2.4, 1) ==false) continue;
       //      cout<<nPart<<" jet: "<<i<<" ";
       //      jet[i].lv.Print();

       InitializeTopSearch(myTopSearch, nPart, fMT2tree->jet[i].lv, fMT2tree->jet[i].IsBJet(3));
       nPart++;

       InitializeHemisphereFinder(px, py, pz, E, fMT2tree->jet[i].lv);

       hemiobjs.index.push_back(i); 
       hemiobjs.type.push_back("jet"), hemiobjs.hemisphere.push_back(0);
     }
     myTopSearch->SetNpart(nPart);
   
     myTopSearch->GoSearch();

     NTopsRecTop += myTopSearch->NpT();

     for(int i = 0; i < myTopSearch->NpT(); i++){
       TLorentzVector Top = myTopSearch->PT(i);

       hChi2T->Fill(myTopSearch->ChisqT(i));
       hPtTop->Fill(Top.Pt());
       hEtaTop->Fill(fabs(Top.Eta()));
       hMT2Top->Fill(fMT2tree->misc.MT2);

       int Wnumber = myTopSearch->IpT(i,0);
       int W_1 = myTopSearch->IpW(Wnumber, 1);

       if(W_1 == -990)
	 NTopsW1b++;
       else{
	 if(W_1 == -999)
	   NTopsW1j++;
	   
	 if(myTopSearch->heavyJet(i) > 0)
	   NTopsRightB++;
       }
     }
     }//Wjets

     if(U1.Index == D2.Index)      
       continue;
    
     MT2GenParticle TopSemiLep;
     int TopSemiLepIndex = -1;
     if(abs(U1.Index + D1.Index + U2.Index + D2.Index) > 100){
       if(U1.Index == D1.Index){
	 TopSemiLep = T2;
	 if(T1.lv.Pt() > T2.lv.Pt())
	   TopSemiLepIndex = 1;
	 else
	   TopSemiLepIndex = 0;
       }
       else{
	 TopSemiLep = T1;
	 if(T1.lv.Pt() > T2.lv.Pt())
	   TopSemiLepIndex = 0;
	 else
	   TopSemiLepIndex = 1;
	
       }
       // cout<<" SemiLep MT2 :: "<<fMT2tree->misc.MT2<<" TopEta "<<TopSemiLep.lv.Eta()<<endl;
       // TopSemiLep.lv.Print();
       hMT2Top->Fill(fMT2tree->misc.MT2);
       hPtTop->Fill(TopSemiLep.lv.Pt());
       hEtaTop->Fill(fabs(TopSemiLep.lv.Eta()));


       hemiobjs.index.clear(); 
       hemiobjs.type.clear();
       hemiobjs.hemisphere.clear();
       nPart = 0;
       myTopSearch->SetNpart(nPart);   
       px.clear();
       py.clear();
       pz.clear();
       E.clear();

       for(int i=0; i<fMT2tree->NJets; ++i){
	 if(fMT2tree->jet[i].IsGoodPFJet(20.0, 2.4, 1) ==false) continue;
	 //      cout<<nPart<<" jet: "<<i<<" ";
	 //      jet[i].lv.Print();

	 InitializeTopSearch(myTopSearch, nPart, fMT2tree->jet[i].lv, fMT2tree->jet[i].IsBJet(3));
	 nPart++;

	 InitializeHemisphereFinder(px, py, pz, E, fMT2tree->jet[i].lv);

	 hemiobjs.index.push_back(i); 
	 hemiobjs.type.push_back("jet"), hemiobjs.hemisphere.push_back(0);
       }

       myTopSearch->SetNpart(nPart);

       myTopSearch->GoSearch();
       
       float minDR = 100.0;
       int MatchedRecTopIndex = -1;

       NTopsRecTop += myTopSearch->NpT();

       for(int i = 0; i < myTopSearch->NpT(); i++){
	 TLorentzVector Top = myTopSearch->PT(i);
	 
	 if(TopSemiLep.lv.DeltaR(Top) < minDR){
	   minDR = TopSemiLep.lv.DeltaR(Top);
	   MatchedRecTopIndex = i;
	 }
       }
       
       hMinDR->Fill(minDR);
       if(minDR < 1.0){
	 int Wnumber = myTopSearch->IpT(MatchedRecTopIndex,0);
	 int W_1 = myTopSearch->IpW(Wnumber, 1);

	 if(W_1 == -990)
	   NTopsW1b++;
	 else{
	   if(W_1 == -999)
	     NTopsW1j++;
	   
	   if(myTopSearch->heavyJet(MatchedRecTopIndex) > 0)
	     NTopsRightB++;
	 }

	 // cout<<"Rec SemiLep MT2 :: "<<fMT2tree->misc.MT2<<" DR "<<minDR<<" MatchedRecTopIndex "<<MatchedRecTopIndex<<endl;
	 // myTopSearch->PT(MatchedRecTopIndex).Print();
	 hDPtTop->Fill((TopSemiLep.lv.Pt() - myTopSearch->PT(MatchedRecTopIndex).Pt())/TopSemiLep.lv.Pt());
	 hMT2TopEff->Fill(fMT2tree->misc.MT2);
	 hPtTopEff->Fill(TopSemiLep.lv.Pt());
	 hEtaTopEff->Fill(fabs(TopSemiLep.lv.Eta()));
       }
     }//SemiLept
     
     int b1, b2, u1, u2, d1, d2;

     float MaxDR = 0.0;
     
     if(B1.Index != -9 && B2.Index != -9 && U1.Index != -9 && U2.Index != -9 && D1.Index != -9 && D2.Index != -9){
       // hMinDR->Fill(MatchedJets(B1, b1));
       // hMinDR->Fill(MatchedJets(B2, b2));
       // hMinDR->Fill(MatchedJets(U1, u1));
       // hMinDR->Fill(MatchedJets(U2, u2));
       // hMinDR->Fill(MatchedJets(D1, d1));
       // hMinDR->Fill(MatchedJets(D2, d2));

       if(MatchedJets(B1, b1) > MaxDR)
	 MaxDR = MatchedJets(B1, b1);
       if(MatchedJets(B2, b2) > MaxDR)
	 MaxDR = MatchedJets(B2, b2);
       if(MatchedJets(U1, u1) > MaxDR)
	 MaxDR = MatchedJets(U1, u1);
       if(MatchedJets(U2, u2) > MaxDR)
	 MaxDR = MatchedJets(U2, u2);
       if(MatchedJets(D1, d1) > MaxDR)
	 MaxDR = MatchedJets(D1, d1);
       if(MatchedJets(D2, d2) > MaxDR)
	 MaxDR = MatchedJets(D2, d2);
     }else{
       // cout<<" B1 "<<B1.Index<<" B2 "<<B2.Index<<" D1 "<<D1.Index<<" D2 "<<D2.Index<<" U1 "<<U1.Index<<" U2 "<<U2.Index<<endl;
       continue;
     }

     // if(fabs(T1.lv.Eta()) > 2.4 || T1.lv.Pt() < 100.0 || fabs(T2.lv.Eta()) > 2.4 || T2.lv.Pt() < 100.0)
     // 	continue;

     // if(!(fabs(B1.lv.Eta()) < 2.4 && fabs(U1.lv.Eta()) < 2.4 && fabs(D1.lv.Eta()) < 2.4  && fabs(B2.lv.Eta()) < 2.4 && fabs(U2.lv.Eta()) < 2.4 && fabs(D2.lv.Eta()) < 2.4 && B1.lv.Pt() > 20.0 && U1.lv.Pt() > 20.0  && D1.lv.Pt() > 20.0  && B2.lv.Pt() > 20.0  && U2.lv.Pt() > 20.0  && D2.lv.Pt() > 20.0))
     //   continue;
     vector<MT2GenParticle> GenTops;

     if(T1.lv.Pt() > T2.lv.Pt() ){
       GenTops.push_back(T1);
       GenTops.push_back(T2);
     }else{
       GenTops.push_back(T2);
       GenTops.push_back(T1);
     }


     nGenFullHad++;
     hMT2Top->Fill(fMT2tree->misc.MT2);
     hMT2Top->Fill(fMT2tree->misc.MT2);
     // cout<<" Full Had MT2 :: "<<fMT2tree->misc.MT2<<" TopEta "<<T1.lv.Eta()<<endl;
     // T1.lv.Print();
     // cout<<" Full Had MT2 :: "<<fMT2tree->misc.MT2<<" TopEta "<<T2.lv.Eta()<<endl;
     // T2.lv.Print();
     hPtTop->Fill(T1.lv.Pt());
     hPtTop->Fill(T2.lv.Pt());
     hEtaTop->Fill(fabs(T1.lv.Eta()));
     hEtaTop->Fill(fabs(T2.lv.Eta()));
     hPtW->Fill(W1.lv.Pt());
     hPtW->Fill(W2.lv.Pt());

     if(T1.ID < T2.ID)
       for(int i = 0; i < fMT2tree->NGenParticles; i++)
	 cout<<i<<" Id "<<fMT2tree->genparticle[i].ID<<" Mo "<<fMT2tree->genparticle[i].MIndex<<" GMo "<<fMT2tree->genparticle[i].GMIndex<<endl; 

     int Overlap = 0;

     if(b1 == b2 || u1 == u2 || d1 == d2 || u1 == d2 || d1 == u2 || u1 == b2 || d1 == b2 || b1 == u2 || b1 == d2)
       Overlap += 1;

     if(b1 == u1 || b1 == d1 || u1 == d1)
       Overlap += 2;

     if(b2 == u2 || b2 == d2 || u2 == d2)
       Overlap += 4;

  
     nPart = 0;
     myTopSearch->SetNpart(nPart);   
  
     InitializeTopSearch(myTopSearch, nPart, B1.lv, 1);
     InitializeHemisphereFinder(px, py, pz, E, B1.lv);
     hemiobjs.index.push_back(nPart); 
     hemiobjs.type.push_back("jet"), hemiobjs.hemisphere.push_back(0);
     nPart++;
  
     InitializeTopSearch(myTopSearch, nPart, U1.lv, 0);
     InitializeHemisphereFinder(px, py, pz, E, U1.lv);
     hemiobjs.index.push_back(nPart); 
     hemiobjs.type.push_back("jet"), hemiobjs.hemisphere.push_back(0);
     nPart++;
      
     InitializeTopSearch(myTopSearch, nPart, D1.lv, 0);
     InitializeHemisphereFinder(px, py, pz, E, D1.lv);
     hemiobjs.index.push_back(nPart); 
     hemiobjs.type.push_back("jet"), hemiobjs.hemisphere.push_back(0);
     nPart++;

     InitializeTopSearch(myTopSearch, nPart, B2.lv, 1);
     InitializeHemisphereFinder(px, py, pz, E, B2.lv);
     hemiobjs.index.push_back(nPart); 
     hemiobjs.type.push_back("jet"), hemiobjs.hemisphere.push_back(0);
     nPart++;

     InitializeTopSearch(myTopSearch, nPart, U2.lv, 0);
     InitializeHemisphereFinder(px, py, pz, E, U2.lv);
     hemiobjs.index.push_back(nPart); 
     hemiobjs.type.push_back("jet"), hemiobjs.hemisphere.push_back(0);
     nPart++;
      
     InitializeTopSearch(myTopSearch, nPart, D2.lv, 0);
     InitializeHemisphereFinder(px, py, pz, E, D2.lv);
     hemiobjs.index.push_back(nPart); 
     hemiobjs.type.push_back("jet"), hemiobjs.hemisphere.push_back(0);
     nPart++;

     myTopSearch->SetNpart(nPart);

     myTopSearch->GoSearch();

     // get hemispheres (seed 2: max inv mass, association method: default 3 = minimal lund distance)
     Hemisphere* hemisp = new Hemisphere(px, py, pz, E, 2, 3);
     //      2,1  0.38
     //      1,3  0.33
     vector<int> grouping = hemisp->getGrouping();

     //    TLorentzVector pseudojet1(0.,0.,0.,0.);
     //     TLorentzVector pseudojet2(0.,0.,0.,0.);
	
     //     for(int i=0; i<px.size(); ++i){
     // 	if(grouping[i]==1){
     // 	  pseudojet1.SetPx(pseudojet1.Px() + px[i]);
     // 	  pseudojet1.SetPy(pseudojet1.Py() + py[i]);
     // 	  pseudojet1.SetPz(pseudojet1.Pz() + pz[i]);
     // 	  pseudojet1.SetE( pseudojet1.E()  + E[i]);	
     // 	}else if(grouping[i] == 2){
     // 	  pseudojet2.SetPx(pseudojet2.Px() + px[i]);
     // 	  pseudojet2.SetPy(pseudojet2.Py() + py[i]);
     // 	  pseudojet2.SetPz(pseudojet2.Pz() + pz[i]);
     // 	  pseudojet2.SetE( pseudojet2.E()  + E[i]);
     // 	}
     //     }

     // // define MET
     // TLorentzVector MET = B1.lv + B2.lv + U1.lv + U2.lv + D1.lv + D2.lv;

     // Float_t MT2 = fMT2tree->CalcMT2(0,0, pseudojet1, pseudojet2, MET);
     // hTMass->Fill(MT2);

 
     for(int i = 0; i < myTopSearch->NpT(); i++){
       // if(myTopSearch->NpT() == 1)
       //   continue;

       TLorentzVector Top = myTopSearch->PT(i);
       int b = myTopSearch->IpT(i,1);
       int Wnumber = myTopSearch->IpT(i,0);
       int W_0 = myTopSearch->IpW(Wnumber, 0);
       int W_1 = myTopSearch->IpW(Wnumber, 1);
       // cout<<" myTopSearch->NpT() "<<myTopSearch->NpT()<<" b "<<b<<" W_0 "<<W_0<<" W_1 "<<W_1<<endl;
       // cout<<" W "<<endl;
       // myTopSearch->PW(i).Print();
       // cout<<" W1.lv "<<endl;
       // W1.lv.Print();
       // cout<<" T2.lv "<<endl;
       // W2.lv.Print();
	
       if((b == 0 && fabs((myTopSearch->PW(i) - W1.lv).M()) < 1 && fabs((myTopSearch->PW(i) - W1.lv).E()) < 1 ) || 
	  (b == 3 && fabs((myTopSearch->PW(i) - W2.lv).M()) < 1 && fabs((myTopSearch->PW(i) - W2.lv).E()) < 1 ))
	 CorrectTopGenFullHad++;
       else
	 if((b + W_0 + W_1) == 3 || (b + W_0 + W_1) == 12 )
	   CorrectTopGenFullHad++;
	   

       //   else{
       //     //	    if(!((Top - T1.lv).M() < 0.1 || (Top - T2.lv).M() < 0.1)){
       //   if(!((b + W_0 + W_1) == 3 || (b + W_0 + W_1) == 12 )){
       //     cout<<" myTopSearch->NpT() "<<myTopSearch->NpT()<<" b "<<b<<" W_0 "<<W_0<<" W_1 "<<W_1<<endl;
       //     cout<<" Top "<<Top.M()<<endl;
       //     Top.Print();
       //     cout<<" T1.lv "<<T1.lv.M()<<endl;
       //     T1.lv.Print();
       //     cout<<" T2.lv "<<T2.lv.M()<<endl;
       //     T2.lv.Print();
       //     cout<<" W "<<myTopSearch->PW(i).M()<<endl;
       //     myTopSearch->PW(i).Print();
       //     cout<<" W1.lv "<<W1.lv.M()<<endl;
       //     W1.lv.Print();
       //     cout<<" W2.lv "<<W2.lv.M()<<endl;
       //     W2.lv.Print();
       //    }
       // }
     }

     delete hemisp;
     int B1Hemi = grouping[0];
     int B2Hemi = grouping[3];
     if(B1Hemi == grouping[1] && B1Hemi == grouping[2]){
       TogetherGenFullHad++;
       // hPtTopTogether->Fill(T1.lv.Pt());
       // hPtWTogether->Fill(W1.lv.Pt());
     }


     if(B2Hemi == grouping[4] && B2Hemi == grouping[5]){
       TogetherGenFullHad++;
       // hPtTopTogether->Fill(T2.lv.Pt());
       // hPtWTogether->Fill(W2.lv.Pt());
     }

     if(B1Hemi != grouping[1] && B1Hemi != grouping[2])
       SplittedBGenFullHad++;

     if(B2Hemi != grouping[4] && B2Hemi != grouping[5])
       SplittedBGenFullHad++;

     if(grouping[1] != grouping[2])
       SplittedWGenFullHad++;

     if(grouping[4] != grouping[5])
       SplittedWGenFullHad++;

     if(Overlap == 0 && MaxDR < 0.5){
       //if((b_0 == b1 && ((W_0 == u1 && W_1 == d1) || (W_0 == d1 && W_1 == u1))) || (b0 == b2 && ((W_0 == u2 && W_1 == d2) || (W_0 == d2 && W_1 == u2))))

       w1 = fMT2tree->jet[u1].lv + fMT2tree->jet[d1].lv;
       t1 = fMT2tree->jet[b1].lv + w1;

       w2 = fMT2tree->jet[u2].lv + fMT2tree->jet[d2].lv;
       t2 = fMT2tree->jet[b2].lv + w2;

       // hWMass->Fill(w1.M());
       // hWMass->Fill(w2.M());
       // hTMass->Fill(t1.M());
       // hTMass->Fill(t2.M());
       // hTWMass->Fill(w1.M(), t1.M());
       // hTWMass->Fill(w2.M(), t2.M());
     }

     // struct to keep track of indices in hemi association
     // struct HemiObj {vector<TString> type; vector<int> index; vector<int> hemisphere;
     //} hemiobjs;
     hemiobjs.index.clear(); 
     hemiobjs.type.clear();
     hemiobjs.hemisphere.clear();
     grouping.clear();
     px.clear();
     py.clear();
     pz.clear();
     E.clear();

     nPart = 0;
     InitializeHemisphereFinder(px, py, pz, E, B2.lv);
     hemiobjs.index.push_back(nPart); 
     hemiobjs.type.push_back("jet"), hemiobjs.hemisphere.push_back(0);
     nPart++;
  
     InitializeHemisphereFinder(px, py, pz, E, W2.lv);
     hemiobjs.index.push_back(nPart); 
     hemiobjs.type.push_back("jet"), hemiobjs.hemisphere.push_back(0);
     nPart++;
   
     InitializeHemisphereFinder(px, py, pz, E, T1.lv);
     hemiobjs.index.push_back(nPart); 
     hemiobjs.type.push_back("Top"), hemiobjs.hemisphere.push_back(0);
     nPart++;

     hemisp = new Hemisphere(px, py, pz, E, 2, 3);
     grouping = hemisp->getGrouping();
     delete hemisp;

     B1Hemi = grouping[0];
     if(B1Hemi == grouping[1]){
       TogetherGenTBW++;
       // hPtTopTogether->Fill(T2.lv.Pt());
       // hPtWTogether->Fill(W2.lv.Pt());
     }

     if(B1Hemi == grouping[1] && B1Hemi == grouping[2])
       cout<<"TogetherGenTBW::Everything in one hemi !!!!"<<endl;

     hemiobjs.index.clear(); 
     hemiobjs.type.clear();
     hemiobjs.hemisphere.clear();
     grouping.clear();
     px.clear();
     py.clear();
     pz.clear();
     E.clear();

     nPart = 0;
     InitializeHemisphereFinder(px, py, pz, E, B1.lv);
     hemiobjs.index.push_back(nPart); 
     hemiobjs.type.push_back("jet"), hemiobjs.hemisphere.push_back(0);
     nPart++;
  
     InitializeHemisphereFinder(px, py, pz, E, W1.lv);
     hemiobjs.index.push_back(nPart); 
     hemiobjs.type.push_back("jet"), hemiobjs.hemisphere.push_back(0);
     nPart++;
   
     InitializeHemisphereFinder(px, py, pz, E, T2.lv);
     hemiobjs.index.push_back(nPart); 
     hemiobjs.type.push_back("Top"), hemiobjs.hemisphere.push_back(0);
     nPart++;

     hemisp = new Hemisphere(px, py, pz, E, 2, 3);
     grouping = hemisp->getGrouping();
     delete hemisp;

     B1Hemi = grouping[0];
     if(B1Hemi == grouping[1]){
       TogetherGenTBW++;
       // hPtTopTogether->Fill(T1.lv.Pt());
       // hPtWTogether->Fill(W1.lv.Pt());
     }

     if(B1Hemi == grouping[1] && B1Hemi == grouping[2])
       cout<<"TogetherGenTBW::Everything in one hemi !!!!"<<endl;

     hemiobjs.index.clear(); 
     hemiobjs.type.clear();
     hemiobjs.hemisphere.clear();
     grouping.clear();
     px.clear();
     py.clear();
     pz.clear();
     E.clear();

     nPart = 0;
     InitializeHemisphereFinder(px, py, pz, E, B2.lv);
     hemiobjs.index.push_back(nPart); 
     hemiobjs.type.push_back("jet"), hemiobjs.hemisphere.push_back(0);
     nPart++;
  
     InitializeHemisphereFinder(px, py, pz, E, U2.lv);
     hemiobjs.index.push_back(nPart); 
     hemiobjs.type.push_back("jet"), hemiobjs.hemisphere.push_back(0);
     nPart++;
   
     InitializeHemisphereFinder(px, py, pz, E, D2.lv);
     hemiobjs.index.push_back(nPart); 
     hemiobjs.type.push_back("jet"), hemiobjs.hemisphere.push_back(0);
     nPart++;
   
     InitializeHemisphereFinder(px, py, pz, E, T1.lv);
     hemiobjs.index.push_back(nPart); 
     hemiobjs.type.push_back("Top"), hemiobjs.hemisphere.push_back(0);
     nPart++;

     hemisp = new Hemisphere(px, py, pz, E, 2, 3);
     grouping = hemisp->getGrouping();
     delete hemisp;

     B1Hemi = grouping[0];
     if(B1Hemi == grouping[1] && B1Hemi == grouping[2]){
       TogetherGenTBQQ++;
       hPtTopTogether->Fill(T2.lv.Pt());
       hPtWTogether->Fill(W2.lv.Pt());
     }else
       if(B1Hemi != grouping[1] && B1Hemi != grouping[2])
	 SplittedBGenTBQQ++;
       else
	 if(grouping[1] != grouping[2])
	   SplittedWGenTBQQ++;
	 else
	   cout<<" B1Hemi  = "<<B1Hemi<<" grouping[1] "<<grouping[1]<<" grouping[2] "<<grouping[2]<<" grouping[3] "<<grouping[3]<<endl;

     if(B1Hemi == grouping[1] && B1Hemi == grouping[2] && B1Hemi == grouping[3])
       cout<<"TogetherGenTBQQ::Everything in one hemi !!!!"<<endl;


     hemiobjs.index.clear(); 
     hemiobjs.type.clear();
     hemiobjs.hemisphere.clear();
     grouping.clear();
     px.clear();
     py.clear();
     pz.clear();
     E.clear();

     nPart = 0;
     InitializeHemisphereFinder(px, py, pz, E, B1.lv);
     hemiobjs.index.push_back(nPart); 
     hemiobjs.type.push_back("jet"), hemiobjs.hemisphere.push_back(0);
     nPart++;
  
     InitializeHemisphereFinder(px, py, pz, E, U1.lv);
     hemiobjs.index.push_back(nPart); 
     hemiobjs.type.push_back("jet"), hemiobjs.hemisphere.push_back(0);
     nPart++;
   
     InitializeHemisphereFinder(px, py, pz, E, D1.lv);
     hemiobjs.index.push_back(nPart); 
     hemiobjs.type.push_back("jet"), hemiobjs.hemisphere.push_back(0);
     nPart++;
   
     InitializeHemisphereFinder(px, py, pz, E, T2.lv);
     hemiobjs.index.push_back(nPart); 
     hemiobjs.type.push_back("Top"), hemiobjs.hemisphere.push_back(0);
     nPart++;

     hemisp = new Hemisphere(px, py, pz, E, 2, 3);
     grouping = hemisp->getGrouping();
     delete hemisp;

     B1Hemi = grouping[0];
     if(B1Hemi == grouping[1] && B1Hemi == grouping[2]){
       TogetherGenTBQQ++;
       hPtTopTogether->Fill(T1.lv.Pt());
       hPtWTogether->Fill(W1.lv.Pt());
     }else
       if(B1Hemi != grouping[1] && B1Hemi != grouping[2])
	 SplittedBGenTBQQ++;
       else
	 if(grouping[1] != grouping[2])
	   SplittedWGenTBQQ++;
	 else
	   cout<<" B1Hemi  = "<<B1Hemi<<" grouping[1] "<<grouping[1]<<" grouping[2] "<<grouping[2]<<" grouping[3] "<<grouping[3]<<endl;

     if(B1Hemi == grouping[1] && B1Hemi == grouping[2] && B1Hemi == grouping[3])
       cout<<"TogetherGenTBQQ::Everything in one hemi !!!!"<<endl;


     hemiobjs.index.clear(); 
     hemiobjs.type.clear();
     hemiobjs.hemisphere.clear();
     grouping.clear();
     nPart = 0;
     myTopSearch->SetNpart(nPart);   
     px.clear();
     py.clear();
     pz.clear();
     E.clear();

     for(int i=0; i<fMT2tree->NJets; ++i){
       if(fMT2tree->jet[i].IsGoodPFJet(20.0, 2.4, 1) ==false) continue;
       //      cout<<nPart<<" jet: "<<i<<" ";
       //      jet[i].lv.Print();

       InitializeTopSearch(myTopSearch, nPart, fMT2tree->jet[i].lv, fMT2tree->jet[i].IsBJet(3));
       nPart++;

       InitializeHemisphereFinder(px, py, pz, E, fMT2tree->jet[i].lv);

       hemiobjs.index.push_back(i); 
       hemiobjs.type.push_back("jet"), hemiobjs.hemisphere.push_back(0);
     }

     myTopSearch->SetNpart(nPart);

     myTopSearch->GoSearch();

     // get hemispheres (seed 2: max inv mass, association method: default 3 = minimal lund distance)
     hemisp = new Hemisphere(px, py, pz, E, 2, 3);
     grouping = hemisp->getGrouping();
     delete hemisp;
   
     int hemi1 = 0;
     int hemi2 = 0;
     for(int i = 0; i <(int) px.size(); i++){
       if(grouping[i] == 1)
	 hemi1++;
       if(grouping[i] == 2)
	 hemi2++;
     }
     if((hemi1 + hemi2) != px.size())
       cout<<" hemi1 "<<hemi1<<" hemi2 "<<hemi2<<" px.size() "<<px.size()<<endl;
      
     hHemi12->Fill(hemi1, hemi2);

     for(int i = 0; i < myTopSearch->NpT(); i++){
       TLorentzVector Top = myTopSearch->PT(i);


       // if(abs(myTopSearch->heavyJet(i)) > 200 || nTops > 2)
       // if(nTops > 2)
       //   continue;

       int b_0 = myTopSearch->IpT(i,1);
       int WNumber = myTopSearch->IpT(i,0);
       int W_0 = myTopSearch->IpW(WNumber, 0);
       int W_1 = myTopSearch->IpW(WNumber, 1);

       // if(fabs(myTopSearch->Mw(WNumber) - 80.4) > 0.5 || fabs((myTopSearch->PT(i)).M() - 172.0) > 1.0)
       //   continue;

       nTops++;
       // cout<<" grouping size "<<grouping.size()<<" E.size() "<<E.size()<<" fMT2tree->NJets "<<fMT2tree->NJets<<endl;
       int bHemi = grouping[b_0];
       // cout<<" nTops "<<nTops<<endl;
       // cout<<" heavyJet(i) "<<myTopSearch->heavyJet(i)<<" mW "<<myTopSearch->Mw(WNumber)<<" mT "<<(myTopSearch->PT(i)).M()<<endl;
       hWMass->Fill(myTopSearch->Mw(WNumber));
       hTMass->Fill((myTopSearch->PT(i)).M());
       hTWMass->Fill(myTopSearch->Mw(WNumber),(myTopSearch->PT(i)).M()); 
       hChi2T->Fill(myTopSearch->ChisqT(i));

       // cout<<" b_0, hemi "<<b_0<<" ,"<<bHemi<<endl;
       // cout<<" W_0, hemi "<<W_0<<" ,"<<grouping[W_0]<<endl;
       // cout<<" W_1, hemi "<<W_1<<" ,"<<grouping[W_1]<<endl;
	
       if(abs(myTopSearch->heavyJet(i)) < 200){
	 if(grouping[W_0]==bHemi && grouping[W_1]==bHemi)
	   together++;
	 else
	   if(grouping[W_0]!=bHemi && grouping[W_1]!=bHemi)
	     splittedT++;
	   else
	     if(grouping[W_0]!= grouping[W_1])
	       splittedW++;
	     else{
	       cout<<"What is this? "<<endl;
	       cout<<"b_0, hemi "<<b_0<<", "<<bHemi<<endl;
	       cout<<"W_0, hemi "<<W_0<<", "<<grouping[W_0]<<endl;
	       cout<<"W_1, hemi "<<W_1<<", "<<grouping[W_1]<<endl;
	     }
       }else{
	 if(grouping[W_0]==bHemi)
	   together++;
	 else
	   splittedT++;
       }
     }
     // for(int i = 0; i < myTopSearch->NpT(); i++){
     //   cout<<"Rec Full Had  "<<myTopSearch->PT(i).Eta()<<endl;}
     NTopsRecTop += myTopSearch->NpT();
     vector<int> MatchedGenTop;
     for(int j = 0; j < 2; j++){
       float minDR = 100.0;
       int MatchedRecTopIndex = -1;

       for(int i = 0; i < myTopSearch->NpT(); i++){
	 
	 if(MatchedGenTop.size() > 0 && i == MatchedGenTop[0])
	   continue;

	 TLorentzVector Top = myTopSearch->PT(i);
	 
	 if(GenTops[j].lv.DeltaR(Top) < minDR){
	   minDR = GenTops[j].lv.DeltaR(Top);
	   MatchedRecTopIndex = i;
	 }
       }

       hMinDR->Fill(minDR);
       if(minDR < 1.0){
	 MatchedGenTop.push_back(MatchedRecTopIndex);
	 int Wnumber = myTopSearch->IpT(MatchedRecTopIndex,0);
	 int W_1 = myTopSearch->IpW(Wnumber, 1);
	 if(W_1 == -990)
	   NTopsW1b++;
	 else{
	   if(W_1 == -999)
	     NTopsW1j++;
	   
	   if(myTopSearch->heavyJet(MatchedRecTopIndex) > 0)
	     NTopsRightB++;
	 }

	 hDPtTop->Fill((GenTops[j].lv.Pt() - myTopSearch->PT(MatchedGenTop[j]).Pt())/GenTops[j].lv.Pt());
	 hMT2TopEff->Fill(fMT2tree->misc.MT2);
	 // cout<<"Rec Full Had MT2 :: "<<fMT2tree->misc.MT2<<" DR "<<minDR<<" MatchedRecTopIndex "<<MatchedRecTopIndex<<endl;
	 // myTopSearch->PT(MatchedGenTop[j]).Print();
	 hPtTopEff->Fill(GenTops[j].lv.Pt());
	 hEtaTopEff->Fill(fabs(GenTops[j].lv.Eta()));
       }
     }
     
     if(nTops == 0)
       NTops0++;

     if(nTops == 1){
       NTops1 += nTops;
       Together1  += together;
       SplittedW1 += splittedW;
       SplittedT1 += splittedT;
     }
   
     if(nTops >= 2){
       NTops2 += nTops;
       Together2  += together;
       SplittedW2 += splittedW;
       SplittedT2 += splittedT;
     }

    }

  }



  cout<<" Total number of events with GenFullHad = "<<nGenFullHad<<endl;
  cout<<"TogetherGenTBW = "<<TogetherGenTBW<<", "<<(TogetherGenTBW/(2.0 * nGenFullHad))<<endl;
  cout<<"==================="<<endl;

  cout<<"TogetherGenTBQQ = "<<TogetherGenTBQQ<<", "<<(TogetherGenTBQQ/(2.0 * nGenFullHad))<<endl;
  cout<<"SplittedTGenTBQQ = "<<SplittedBGenTBQQ<<", "<<(SplittedBGenTBQQ/(2.0 * nGenFullHad))<<endl;
  cout<<"SplittedWGenTBQQ = "<<SplittedWGenTBQQ<<", "<<(SplittedWGenTBQQ/(2.0 * nGenFullHad))<<endl;
  cout<<"==================="<<endl;
 
  cout<<"CorrectTopGenFullHad = "<<CorrectTopGenFullHad<<", "<<(CorrectTopGenFullHad/(2.0 * nGenFullHad))<<endl;

  cout<<"==================="<<endl;
  cout<<"TogetherGenFullHad = "<<TogetherGenFullHad<<", "<<(TogetherGenFullHad/(2.0 * nGenFullHad))<<endl;
  cout<<"SplittedTGenFullHad = "<<SplittedBGenFullHad<<", "<<(SplittedBGenFullHad/(2.0 * nGenFullHad))<<endl;
  cout<<"SplittedWGenFullHad = "<<SplittedWGenFullHad<<", "<<(SplittedWGenFullHad/(2.0 * nGenFullHad))<<endl;
  cout<<"==================="<<endl;

  // cout<<"Total number of events = "<<(NTops0 + NTops1 + NTops2/2)<<endl;  

  cout<<"Events with 0 Top = "<<NTops0<<endl;
  cout<<"==================="<<endl;
  cout<<"Events with 1 Top = "<<NTops1<<endl;
  cout<<"Together = "<<Together1<<", "<<(Together1/(1.0 * NTops1))<<endl;
  cout<<"SplittedT = "<<SplittedT1<<", "<<(SplittedT1/(1.0 * NTops1))<<endl;
  cout<<"SplittedW = "<<SplittedW1<<", "<<(SplittedW1/(1.0 * NTops1))<<endl;
  cout<<"==================="<<endl;
  cout<<"Tops in events with >= 2 Top = "<<NTops2<<endl;
  cout<<"Together = "<<Together2<<", "<<(Together2/(1.0 * NTops2))<<endl;
  cout<<"SplittedT = "<<SplittedT2<<", "<<(SplittedT2/(1.0 * NTops2))<<endl;
  cout<<"SplittedW = "<<SplittedW2<<", "<<(SplittedW2/(1.0 * NTops2))<<endl;


  cout<< " W1b "<<NTopsW1b<< " W1j "<<NTopsW1j<<" RightB "<<(NTopsRightB + NTopsW1b)<<" NTopsRecTop "<<NTopsRecTop<<endl;


  TCanvas *MyC = new TCanvas("MyC","MyC");
  MyC->Divide(4,4);
  MyC->cd(1);
  hWMass->GetXaxis()->SetRangeUser(0.0,300.0);
  hWMass->Draw();
  hTMass->Draw("sames");
  MyC->cd(2);
  hTWMass->Draw();
  MyC->cd(3);
  hHemi12->Draw();
  MyC->cd(4);
  hChi2T->GetXaxis()->SetRangeUser(0.0,5.0);
  hChi2T->Draw();

  MyC->cd(5);
  hPtTop->Draw();
  MyC->cd(6);
  //hPtTopTogether->Divide(hPtTop);
  TGraphAsymmErrors *gPtTopTogether = new TGraphAsymmErrors(hPtTopTogether, hPtTop);
  gPtTopTogether->GetXaxis()->SetTitle("Top Pt (GeV)");
  gPtTopTogether->GetYaxis()->SetTitle("Together Fraction");
  gPtTopTogether->Draw("AP");
  MyC->cd(7);
  //hPtWTogether->Divide(hPtW);
  TGraphAsymmErrors *gPtWTogether = new TGraphAsymmErrors(hPtWTogether, hPtW);
  gPtWTogether->GetXaxis()->SetTitle("W Pt (GeV)");
  gPtWTogether->GetYaxis()->SetTitle("Together Fraction");
  gPtWTogether->Draw("AP");
  
  MyC->cd(8);

  cout<<"Number of Tops = "<<hPtTop->Integral()<<endl;
  cout<<"Number of Rec Tops = "<<hPtTopEff->Integral()<<endl;
  cout<<"Eff = "<<hPtTopEff->Integral()/hPtTop->Integral()<<endl;
  
  //hPtTopEff->Divide(hPtTop);
  TGraphAsymmErrors *gPtTopEff = new TGraphAsymmErrors(hPtTopEff, hPtTop);
  gPtTopEff->GetXaxis()->SetTitle("Top Pt (GeV)");
  gPtTopEff->GetYaxis()->SetTitle("Top Finding Eff");
  gPtTopEff->Draw("AP");

  MyC->cd(9);
  hEtaTop->Draw();
  MyC->cd(10);
  //hEtaTopEff->Divide(hEtaTop);
  TGraphAsymmErrors *gEtaTopEff = new TGraphAsymmErrors(hEtaTopEff, hEtaTop);
  gEtaTopEff->GetXaxis()->SetTitle("Top Eta");
  gEtaTopEff->GetYaxis()->SetTitle("Top Finding Eff");
  gEtaTopEff->Draw("AP");

  MyC->cd(11);
  hMT2Top->Draw();
  MyC->cd(12);
  cout<<"Number of Tops (MT2) = "<<hMT2Top->Integral()<<endl;
  cout<<"Number of Rec Tops = "<<hMT2TopEff->Integral()<<endl;
  cout<<"Eff = "<<hMT2TopEff->Integral()/hMT2Top->Integral()<<endl;

  // hMT2TopEff->Divide(hMT2Top);
  TGraphAsymmErrors *gMT2TopEff = new TGraphAsymmErrors(hMT2TopEff, hMT2Top);
  gMT2TopEff->GetXaxis()->SetTitle("Event MT2");
  gMT2TopEff->GetYaxis()->SetTitle("Top Finding Eff");
  gMT2TopEff->Draw("AP");

  MyC->cd(13);
  hDPtTop->Draw();
  MyC->cd(14);
  hMinDR->Draw();
}



void MassPlotter::MyMakePlot(Long64_t nevents){
  static const int NumberOfSamples = 7;
  TH1F* myHisto[NumberOfSamples];


  TH1::SetDefaultSumw2();
  TString  cnames[NumberOfSamples] = {"QCD", "Wjets", "Zjets", "Top", "MC", "susy", "data"};
  int      ccolor[NumberOfSamples] = {  401,     417,     419,   600,  500,      0,  632};
  TString varname = "MT2";
  for (int i=0; i<NumberOfSamples; i++){
    myHisto[i] = new TH1F(varname+"_"+cnames[i], "", 35, 0, 3.5);
    myHisto[i] -> SetFillColor  (ccolor[i]);
    myHisto[i] -> SetLineColor  (ccolor[i]);
    myHisto[i] -> SetLineWidth  (2);
    myHisto[i] -> SetMarkerColor(ccolor[i]);
    myHisto[i] -> SetStats(false);
  }

  myHisto[6] -> SetMarkerStyle(20);
  myHisto[6] -> SetMarkerColor(kBlack);
  myHisto[6] -> SetLineColor(kBlack);
  
  myHisto[4] -> SetFillStyle(3004);
  myHisto[4] -> SetFillColor(kBlack);
  
  myHisto[5] -> SetLineColor(kBlack);
  myHisto[5] -> SetLineStyle(kDotted);


  for(int ii = 0; ii <(int) fSamples.size(); ii++){
    sample Sample = fSamples[ii];
 
    fMT2tree = new MT2tree();
    int counter = 0;
    Sample.tree->SetBranchAddress("MT2tree", &fMT2tree);
    Long64_t nentries =  Sample.tree->GetEntries();
    float weight = Sample.xsection * Sample.kfact * Sample.lumi / (Sample.nevents * Sample.PU_avg_weight );
    Long64_t minNentriesNevents = min(nentries, nevents);

    std::cout << setfill('=') << std::setw(70) << "" << std::endl;
    cout << "looping over :     " <<endl;	
    cout << "   Name:           " << Sample.name << endl;
    cout << "   File:           " << (Sample.file)->GetName() << endl;
    cout << "   Events:         " << Sample.nevents  << endl;
    cout << "   Events in tree: " << nentries << endl; 
    cout << "   Xsection:       " << Sample.xsection << endl;
    cout << "   kfactor:        " << Sample.kfact << endl;
    cout << "   avg PU weight:  " << Sample.PU_avg_weight << endl;
    cout << "   Weight:         " << weight <<endl;
    std::cout << setfill('-') << std::setw(70) << "" << std::endl;

    for (Long64_t jentry=0; jentry<minNentriesNevents;jentry++) {
     
      Sample.tree->GetEntry(jentry); 
      if ( fVerbose>2 && jentry % 100000 == 0 ){
	fprintf(stdout, "\rProcessed events: %6d of %6d ", jentry + 1, nentries);
	fflush(stdout);
      }

    //   if(jentry ==64){
// 	//	if(((fMT2tree->NJetsIDLoose40  >=4 && fMT2tree->jet[3].lv.Pt() >= 100.0) || (fMT2tree->NJetsIDLoose40  >=5 && fMT2tree->jet[4].lv.Pt() >= 85.0) || (fMT2tree->NJetsIDLoose40  >=6 && fMT2tree->jet[5].lv.Pt() >= 65.0) || (fMT2tree->NJetsIDLoose40  >=7 && fMT2tree->jet[6].lv.Pt() >= 55.0) && (fMT2tree->NEles==0  ||  fMT2tree->ele[0].lv.Pt()<5) ))

// 	cout<<"NEle "<<fMT2tree->NEles<<" pt0 "<<fMT2tree->ele[0].lv.Pt()<<endl;
// 	if(fMT2tree->NEles==0  ||  fMT2tree->ele[0].lv.Pt()<5)                
// 	  cout<<jentry<<endl;
// 	cout<<"nJ "<<fMT2tree->NJetsIDLoose40<<" j4pt "<<fMT2tree->jet[3].lv.Pt()<<" j5pt "<<fMT2tree->jet[4].lv.Pt()<<" j6pt "<<fMT2tree->jet[5].lv.Pt()<<" j7pt "<<fMT2tree->jet[6].lv.Pt()<<endl;

//       }

      if (!(fMT2tree->misc.MET             >  30                                  &&
	    fMT2tree->misc.MT2             >  100                                 &&
	    fMT2tree->NBJetsCSVT           > 0                                     &&
            fMT2tree->NJetsIDLoose40       >= 4                                   && 
           (fMT2tree->NEles==0  ||  fMT2tree->ele[0].lv.Pt()<5)                   &&
 	   (fMT2tree->NMuons==0 ||  fMT2tree->muo[0].lv.Pt()<5)                   &&
 	    fMT2tree->misc.Jet0Pass        == 1                                    &&
 	    fMT2tree->misc.Jet1Pass        == 1                                    &&
// 	    // 	  fMT2tree->misc.SecondJPt       >  100                                  &&
// 	    // 	  fMT2tree->misc.LeadingJPt      >  150                                  &&
 	    fMT2tree->misc.PassJetID       == 1                                    &&
	    fMT2tree->misc.Vectorsumpt     <  70                                   &&
	    //          (fMT2tree->misc.MinMetJetDPhi >0.3||fMT2tree->misc.MinMetJetDPhiIndex>3) &&
	    fMT2tree->misc.MinMetJetDPhi4  >  0.3                                  &&
	    //Six Or Quad
	   ((fMT2tree->NJetsIDLoose40  >=4 && fMT2tree->jet[3].lv.Pt() >= 100.0) || (fMT2tree->NJetsIDLoose40  >=5 && fMT2tree->jet[4].lv.Pt() >= 85.0) || (fMT2tree->NJetsIDLoose40  >=6 && fMT2tree->jet[5].lv.Pt() >= 65.0) || (fMT2tree->NJetsIDLoose40  >=7 && fMT2tree->jet[6].lv.Pt() >= 55.0))                                  &&
	    //HT Trigger
	    // 	   fMT2tree->misc.HT             > 750                                   && 

	    (fMT2tree->misc.ProcessID!=0 || ((fMT2tree->misc.CrazyHCAL==0 &&     fMT2tree->misc.NegativeJEC==0 && fMT2tree->misc.CSCTightHaloIDFlag==0 && fMT2tree->misc.HBHENoiseFlag==0 && fMT2tree->misc.hcalLaserEventFlag==0 &&  fMT2tree->misc.trackingFailureFlag==0 && fMT2tree->misc.eeBadScFlag==0 && fMT2tree->misc.EcalDeadCellTriggerPrimitiveFlag==0 )&&
	     //Six Or Quad
	     ((((fMT2tree->trigger.HLT_QuadJet80_v1 ==1)||(fMT2tree->trigger.HLT_QuadJet80_v2 ==1)||(fMT2tree->trigger.HLT_QuadJet80_v3 ==1)) && ((fMT2tree->NJetsIDLoose40  >=4 && fMT2tree->jet[3].lv.Pt() >= 100.0 ) || (fMT2tree->NJetsIDLoose40  >=5 && fMT2tree->jet[4].lv.Pt() >= 85.0 )))|| (((fMT2tree->trigger.HLT_SixJet45_v1 ==1)||(fMT2tree->trigger.HLT_SixJet45_v2 ==1)||(fMT2tree->trigger.HLT_SixJet45_v3 ==1)) && ((fMT2tree->NJetsIDLoose40  >=6 && fMT2tree->jet[5].lv.Pt() >= 65.0 ) || (fMT2tree->NJetsIDLoose40  >=7 && fMT2tree->jet[6].lv.Pt() >= 55.0 )))))
	     //HT Trigger
	     //    ((fMT2tree->trigger.HLT_PFHT650_v5 ==1) || (fMT2tree->trigger.HLT_PFHT650_v6 ==1) || (fMT2tree->trigger.HLT_PFHT650_v7 ==1) || (fMT2tree->trigger.HLT_PFHT650_v8 ==1) || (fMT2tree->trigger.HLT_PFHT650_v9 ==1))
	        )
))
	continue;
      cout<<++counter<<" : "<<jentry<<endl;
      double myVar = fMT2tree->misc.MinMetBJetDPhi;

      if(Sample.type == "data"){
	    
	myHisto[6]->Fill(myVar, weight);//data
	    
      }else{
	if(Sample.sname == "SMS")
	  myHisto[5]->Fill(myVar, weight);
	else{
	  myHisto[4]->Fill(myVar, weight);
	    
	  if(Sample.sname == "Top")  
	    myHisto[3]->Fill(myVar, weight);
	  else
	    if(Sample.sname == "DY")	
	      myHisto[2]->Fill(myVar, weight);
	    else
	      if(Sample.name == "WJetsToLNu")
		myHisto[1]->Fill(myVar, weight);
	      else
		if(Sample.sname == "QCD")
		  myHisto[0]->Fill(myVar, weight);
	}
      }
    }}
  for(int j = 0; j < NumberOfSamples; j++){
    AddOverAndUnderFlow(myHisto[j]);
  }
  THStack* h_stack     = new THStack(varname, "");
  for(int j = 0; j < NumberOfSamples; j++){

    TH1F* mt2 = (TH1F*)myHisto[j]->Clone();
    mt2->SetName("mt2");
    if(j < (NumberOfSamples - 3))
      h_stack  -> Add(myHisto[j]);
    delete mt2;
  }
    
  TLegend* Legend1 = new TLegend(.71,.54,.91,.92);
  Legend1->AddEntry(myHisto[0], "QCD", "f");
  Legend1->AddEntry(myHisto[1], "W+jets", "f");
  Legend1->AddEntry(myHisto[2], "Z+jets", "f");
  Legend1->AddEntry(myHisto[3], "Top", "f");
  Legend1->AddEntry(myHisto[5], "SUSY", "l");
  Legend1->AddEntry(myHisto[6], "data", "l");
  
  //  TLegend *Legend1;
  printHisto(h_stack, myHisto[6], myHisto[4], Legend1 , "MTC", "hist", true, "MT2", "Events", -4, 0, 1);

  plotRatioStack(h_stack,  myHisto[4], myHisto[6], true, false, "_ratio", Legend1, "MT2", "Events", -4, 0);
}
  


void MassPlotter::MySmallMakePlot(Long64_t nevents){
  static const int NumberOfSamples = 7;
  TH1F* myHisto[NumberOfSamples];


  TH1::SetDefaultSumw2();
  TString  cnames[NumberOfSamples] = {"QCD", "Wjets", "Zjets", "Top", "MC", "susy", "data"};
  int      ccolor[NumberOfSamples] = {  401,     417,     419,   600,  500,      0,  632};
  TString varname = "MT2";
  for (int i=0; i<NumberOfSamples; i++){
    myHisto[i] = new TH1F(varname+"_"+cnames[i], "", 35, 0, 3.5);
    myHisto[i] -> SetFillColor  (ccolor[i]);
    myHisto[i] -> SetLineColor  (ccolor[i]);
    myHisto[i] -> SetLineWidth  (2);
    myHisto[i] -> SetMarkerColor(ccolor[i]);
    myHisto[i] -> SetStats(false);
  }

  myHisto[6] -> SetMarkerStyle(20);
  myHisto[6] -> SetMarkerColor(kBlack);
  myHisto[6] -> SetLineColor(kBlack);
  
  myHisto[4] -> SetFillStyle(3004);
  myHisto[4] -> SetFillColor(kBlack);
  
  myHisto[5] -> SetLineColor(kBlack);
  myHisto[5] -> SetLineStyle(kDotted);


  float MinMetJetDPhi4, MET, VectorSumPt, MT2, MT2Top, MT2MassiveTop, MT2MassiveOnlyTop, HemiDPhi, MinMetBCSVLJetDPhi, MinMetBCSVMJetDPhi, MinMetBCSVTJetDPhi;
  int NJetsIDLoose40, NBJets, NBJetsHE, NBJetsCSVL, NBJetsCSVM, NBJetsCSVT, NTops;
  bool HTTrigger, Trigger4OR6Jets;

  double MET4[4];
  double BestTop4[4];
  double SecondTop4[4];
  double ClosestRBtoMET4[4];
  double LeadingJet4[4];
  double SeconsJet4[4]; 

  for(int ii = 0; ii <(int) fSamples.size(); ii++){
    sample Sample = fSamples[ii];
  
    //    fMT2tree = new MT2tree();
    //    Sample.tree->SetBranchAddress("MT2tree", &fMT2tree);
   
    Sample.tree->SetBranchAddress("MinMetJetDPhi4", &MinMetJetDPhi4);
    Sample.tree->SetBranchAddress("MET", &MET);     
    Sample.tree->SetBranchAddress("VectorSumPt", &VectorSumPt);     
    Sample.tree->SetBranchAddress("MT2", &MT2);     
    Sample.tree->SetBranchAddress("MT2Top", &MT2Top);     
    Sample.tree->SetBranchAddress("MT2MassiveTop", &MT2MassiveTop);     
    Sample.tree->SetBranchAddress("MT2MassiveOnlyTop", &MT2MassiveOnlyTop);
    Sample.tree->SetBranchAddress("NJetsIDLoose40", &NJetsIDLoose40); 
    Sample.tree->SetBranchAddress("NBJets", &NBJets); 
    Sample.tree->SetBranchAddress("NBJetsHE", &NBJetsHE); 
    Sample.tree->SetBranchAddress("NBJetsCSVL", &NBJetsCSVL); 
    Sample.tree->SetBranchAddress("NBJetsCSVM", &NBJetsCSVM); 
    Sample.tree->SetBranchAddress("NBJetsCSVT", &NBJetsCSVT); 
    Sample.tree->SetBranchAddress("NTops", &NTops); 
    Sample.tree->SetBranchAddress("HTTrigger", &HTTrigger); 
    Sample.tree->SetBranchAddress("Trigger4OR6Jets", &Trigger4OR6Jets);
    Sample.tree->SetBranchAddress("MET4", &MET4); 
    Sample.tree->SetBranchAddress("BestTop4", &BestTop4);
    Sample.tree->SetBranchAddress("SecondTop4", &SecondTop4);
    Sample.tree->SetBranchAddress("ClosestRBtoMET4", &ClosestRBtoMET4);
    Sample.tree->SetBranchAddress("LeadingJet4", &LeadingJet4);
    Sample.tree->SetBranchAddress("SeconsJet4", &SeconsJet4);
    Sample.tree->SetBranchAddress("HemiDPhi", &HemiDPhi);
    Sample.tree->SetBranchAddress("MinMetBCSVLJetDPhi", &MinMetBCSVLJetDPhi);
    Sample.tree->SetBranchAddress("MinMetBCSVMJetDPhi", &MinMetBCSVMJetDPhi);
    Sample.tree->SetBranchAddress("MinMetBCSVTJetDPhi", &MinMetBCSVTJetDPhi);

    Long64_t nentries =  Sample.tree->GetEntries();
    float weight = Sample.xsection * Sample.kfact * Sample.lumi / (Sample.nevents * Sample.PU_avg_weight );
    Long64_t minNentriesNevents = min(nentries, nevents);

    std::cout << setfill('=') << std::setw(70) << "" << std::endl;
    cout << "looping over :     " <<endl;	
    cout << "   Name:           " << Sample.name << endl;
    cout << "   File:           " << (Sample.file)->GetName() << endl;
    cout << "   Events:         " << Sample.nevents  << endl;
    cout << "   Events in tree: " << nentries << endl; 
    cout << "   Xsection:       " << Sample.xsection << endl;
    cout << "   kfactor:        " << Sample.kfact << endl;
    cout << "   avg PU weight:  " << Sample.PU_avg_weight << endl;
    cout << "   Weight:         " << weight <<endl;
    std::cout << setfill('-') << std::setw(70) << "" << std::endl;

    for (Long64_t jentry=0; jentry<minNentriesNevents;jentry++) {
     
      Sample.tree->GetEntry(jentry); 
      if ( fVerbose>2 && jentry % 100000 == 0 ){
	fprintf(stdout, "\rProcessed events: %6d of %6d ", jentry + 1, nentries);
	fflush(stdout);
     }

         if (!(MET               >  30                                   &&
 	      MT2                >  100                                  &&
	      NBJetsCSVL         >  0                                    &&
 	      VectorSumPt        <  70                                   &&
// 	      //          (MinMetJetDPhi >0.3||MinMetJetDPhiIndex>3) &&
	      MinMetJetDPhi4     >  0.3                                  &&
	      HemiDPhi           <  2.25                                 &&
	    //Six Or Quad
	      Trigger4OR6Jets              
	       ))
	continue;

	 //	 TLorentzVector Jet1(LeadingJet4[0], LeadingJet4[1], LeadingJet4[2], LeadingJet4[3]);
	 //	 TLorentzVector Jet2(SeconsJet4[0], SeconsJet4[1], SeconsJet4[2], SeconsJet4[3]);
	 
	 double myVar = HemiDPhi;//Util::DeltaPhi(Jet1.Phi(), Jet2.Phi());//Jet1.DeltaPhi(Jet2);

      if(Sample.type == "data"){
	    
	myHisto[6]->Fill(myVar, weight);//data
	    
      }else{
	if(Sample.sname == "SMS")
	  myHisto[5]->Fill(myVar, weight);
	else{
	  myHisto[4]->Fill(myVar, weight);
	    
	  if(Sample.sname == "Top")  
	    myHisto[3]->Fill(myVar, weight);
	  else
	    if(Sample.sname == "DY")	
	      myHisto[2]->Fill(myVar, weight);
	    else
	      if(Sample.name == "WJetsToLNu")
		myHisto[1]->Fill(myVar, weight);
	      else
		if(Sample.sname == "QCD")
		  myHisto[0]->Fill(myVar, weight);
	}
      }
    }}
  for(int j = 0; j < NumberOfSamples; j++){
    AddOverAndUnderFlow(myHisto[j]);
  }
  THStack* h_stack     = new THStack(varname, "");
  for(int j = 0; j < NumberOfSamples; j++){

    TH1F* mt2 = (TH1F*)myHisto[j]->Clone();
    mt2->SetName("mt2");
    if(j < (NumberOfSamples - 3))
      h_stack  -> Add(myHisto[j]);
    delete mt2;
  }
      cout << "------------------------------------"                << endl
	   << "QCD Integral:      " << myHisto[0]->Integral()  << endl
	   << "W+jets Integral:   " << myHisto[1]->Integral()  << endl
	   << "Z+jets Integral:   " << myHisto[2]->Integral()  << endl
	   << "Top Integral:      " << myHisto[3]->Integral()  << endl
	   << "TOTAL BG:          " << myHisto[4]->Integral() <<endl
	   << "SUSY:              " << myHisto[5]->Integral()  << endl
	   << "Data:              " << myHisto[6]->Integral()  << endl;
		
  TLegend* Legend1 = new TLegend(.71,.54,.91,.92);
  Legend1->AddEntry(myHisto[0], "QCD", "f");
  Legend1->AddEntry(myHisto[1], "W+jets", "f");
  Legend1->AddEntry(myHisto[2], "Z+jets", "f");
  Legend1->AddEntry(myHisto[3], "Top", "f");
  Legend1->AddEntry(myHisto[5], "SUSY", "l");
  Legend1->AddEntry(myHisto[6], "data", "l");
  
  printHisto(h_stack, myHisto[6], myHisto[4], myHisto[5], Legend1 , "MTC", "hist", true, "MT2", "Events", -4, 0);

  plotRatioStack(h_stack,  myHisto[4], myHisto[5], true, false, "_ratio", Legend1, "MT2", "Events", -4, 0);
}
  






void MassPlotter::SortHighMT2(float MT2cut, Long64_t nevents){

  vector<int> RunNo;
  vector<int> LSNo;
  vector<int> EventNo;
  vector<int> counterNo;
  vector<float> MT2;
  int counter = 0;
  for(int ii = 0; ii <(int) fSamples.size(); ii++){
    sample Sample = fSamples[ii];
    if(Sample.type != "data")
      continue;

    fMT2tree = new MT2tree();

    Sample.tree->SetBranchAddress("MT2tree", &fMT2tree);
    Long64_t nentries =  Sample.tree->GetEntries();

    Long64_t minNentriesNevents = min(nentries, nevents);
 
    std::cout << setfill('=') << std::setw(70) << "" << std::endl;
    cout << "looping over :     " <<endl;	
    cout << "   Name:           " << Sample.name << endl;
    cout << "   File:           " << (Sample.file)->GetName() << endl;
    cout << "   Events:         " << Sample.nevents  << endl;
    cout << "   Events in tree: " << nentries << endl; 
    cout << "   Events to read: " << minNentriesNevents << endl; 
    cout << "   Xsection:       " << Sample.xsection << endl;
    cout << "   kfactor:        " << Sample.kfact << endl;
    cout << "   avg PU weight:  " << Sample.PU_avg_weight << endl;
    std::cout << setfill('-') << std::setw(70) << "" << std::endl;

    for (Long64_t jentry=0; jentry<minNentriesNevents;jentry++) {
     
      Sample.tree->GetEntry(jentry); 
      if ( fVerbose>2 && jentry % 100000 == 0 ){ 
	fprintf(stdout, "\rProcessed events: %6d of %6d ", jentry + 1, nentries);
	fflush(stdout);
      }

      if (!(fMT2tree->misc.MET            >  30                                   &&
	    fMT2tree->NBJets               > 0                                     &&
	    (fMT2tree->NEles==0  ||  fMT2tree->ele[0].lv.Pt()<10)                   &&
	    (fMT2tree->NMuons==0 ||  fMT2tree->muo[0].lv.Pt()<10)                   &&
	    fMT2tree->misc.Jet0Pass        == 1                                    &&
	    fMT2tree->misc.Jet1Pass        == 1                                    &&
	    fMT2tree->misc.SecondJPt       >  100                                  &&
	    fMT2tree->misc.LeadingJPt      >  150                                  &&
	    fMT2tree->misc.PassJetID       == 1                                    &&
	    fMT2tree->misc.Vectorsumpt     <  70                                   &&
	    (fMT2tree->misc.MinMetJetDPhi >0.3||fMT2tree->misc.MinMetJetDPhiIndex>3) &&
	    //Six Or Quad
	    //((fMT2tree->NJetsIDLoose40  >=4 && fMT2tree->jet[3].lv.Pt() >= 120.0) || (fMT2tree->NJetsIDLoose40  >=6 && fMT2tree->jet[5].lv.Pt() >= 80.0))                                  &&
	    //HT Trigger
	    fMT2tree->misc.HT             > 750                                   && 

	    ((  fMT2tree->misc.ProcessID==10 || (fMT2tree->misc.CrazyHCAL==0 &&     fMT2tree->misc.NegativeJEC==0 && fMT2tree->misc.CSCTightHaloIDFlag==0 && fMT2tree->misc.HBHENoiseFlag==0 && fMT2tree->misc.hcalLaserEventFlag==0 &&  fMT2tree->misc.trackingFailureFlag==0 && fMT2tree->misc.eeBadScFlag==0 && fMT2tree->misc.EcalDeadCellTriggerPrimitiveFlag==0 ))&&
	     //Six Or Quad
	     //	   ((((fMT2tree->trigger.HLT_QuadJet80_v1 ==1)||(fMT2tree->trigger.HLT_QuadJet80_v2 ==1)||(fMT2tree->trigger.HLT_QuadJet80_v3 ==1)) && fMT2tree->NJetsIDLoose40  >=4 && fMT2tree->jet[3].lv.Pt() >= 120.0 ) || (((fMT2tree->trigger.HLT_SixJet45_v1 ==1)||(fMT2tree->trigger.HLT_SixJet45_v2 ==1)||(fMT2tree->trigger.HLT_SixJet45_v3 ==1)) && fMT2tree->NJetsIDLoose40  >=6 && fMT2tree->jet[5].lv.Pt() >= 80.0 ))
	     //HT Trigger
	     ((fMT2tree->trigger.HLT_PFHT650_v5 ==1) || (fMT2tree->trigger.HLT_PFHT650_v6 ==1) || (fMT2tree->trigger.HLT_PFHT650_v7 ==1) || (fMT2tree->trigger.HLT_PFHT650_v8 ==1) || (fMT2tree->trigger.HLT_PFHT650_v9 ==1))
	

	     )
	    ))
	continue;

      if(fMT2tree->misc.MT2 > MT2cut){
	RunNo.push_back(fMT2tree->misc.Run);
	LSNo.push_back(fMT2tree->misc.LumiSection);
	EventNo.push_back(fMT2tree->misc.Event);
	counterNo.push_back(counter);
	MT2.push_back(fMT2tree->misc.MT2);
	counter++;
      }
    }}
  
  // cout<<" MT2 > "<< MT2cut<<endl;
  // for(int i = 0; i < counter; i++){
  //   cout<<counterNo[i]<<" MT2 = "<< MT2[counterNo[i]]<<" Run::LS::Event "<<RunNo[counterNo[i]]<<"::"<<LSNo[counterNo[i]]<<"::"<<EventNo[counterNo[i]]<<endl;
  // }
    
  counterNo = Util::VSort(counterNo, MT2, true);
  cout<<" MT2 > "<< MT2cut<<endl;
  for(int i = 0; i < counter; i++){
    cout<<counterNo[i]<<" MT2 = "<< MT2[counterNo[i]]<<" Run::LS::Event "<<RunNo[counterNo[i]]<<"::"<<LSNo[counterNo[i]]<<"::"<<EventNo[counterNo[i]]<<endl;
  }
}


void MassPlotter::Efficiency(TString mySample){
  
  TString s = ";m_{Glu};m_{LSP}";
  int nBins = 50;
  float min = 0;
  float max = 2500;
  
  if(mySample == "SMST2"){
    s = ";m_{Stop};m_{LSP}";
    nBins = 100;
  }

  TH2F *SMSEventsAllCuts4or6 = new TH2F("SMSEventsAllCuts4or6", "AllCuts4or6" + s, nBins, min, max, nBins, min, max);
  TH2F *SMSEventsAllCuts4or6MT2_100DelPhi = new TH2F("SMSEventsAllCuts4or6MT2_100DelPhi", "4or6MT2_100DelPhi" + s, nBins, min, max, nBins, min, max);
  TH2F *SMSEventsAllCutsHT = new TH2F("SMSEventsAllCutsHT", "AllCutsHT" + s, nBins, min, max, nBins, min, max);
  TH2F *SMSEventsAllCutsHTMT2_100DelPhi = new TH2F("SMSEventsAllCutsHTMT2_100DelPhi", "HTMT2_100DelPhi" + s, nBins, min, max, nBins, min, max);

  TH2F *SMSEvents = new TH2F("SMSEvents", "PreSelection" + s, nBins, min, max, nBins, min, max);
  TH2F *SMSEventsHT = new TH2F("SMSEventsHT", "HT" + s, nBins, min, max, nBins, min, max);
  TH2F *SMSEvents4or6 = new TH2F("SMSEvents4or6", "4 or 6 Jets" + s, nBins, min, max, nBins, min, max);
  TH2F *SMSEventsMT2 = new TH2F("SMSEventsMT2", "MT2" + s, nBins, min, max, nBins, min, max);
  TH2F *SMSEventsMT2100 = new TH2F("SMSEventsMT2100", "MT2" + s, nBins, min, max, nBins, min, max);
  TH2F *SMSEventsMT2Top = new TH2F("SMSEventsMT2Top", "MT2Top" + s, nBins, min, max, nBins, min, max);
  TH2F *SMSEventsNTop1 = new TH2F("SMSEventsNTop1", "NTop1" + s, nBins, min, max, nBins, min, max);
  TH2F *SMSEventsNTop2 = new TH2F("SMSEventsNTop2", "NTop2" + s, nBins, min, max, nBins, min, max);
  TH2F *SMSEventsLeptonVeto = new TH2F("SMSEventsLeptonVeto", "LeptonVeto" + s, nBins, min, max, nBins, min, max);
  TH2F *SMSEventsJetPass = new TH2F("SMSEventsJetPass", "JetPass" + s, nBins, min, max, nBins, min, max);
  TH2F *SMSEventsSecondJet = new TH2F("SMSEventsSecondJet", "SecondJet" + s, nBins, min, max, nBins, min, max);
  TH2F *SMSEventsLeadingJet = new TH2F("SMSEventsLeadingJet", "LeadingJet" + s, nBins, min, max, nBins, min, max);
  TH2F *SMSEventsMetJetDPhi = new TH2F("SMSEventsMetJetDPhi", "MetJetDPhi" + s, nBins, min, max, nBins, min, max);
  TH2F *SMSEventsMetJetDPhi5 = new TH2F("SMSEventsMetJetDPhi5", "MetJetDPhi 0.5 " + s, nBins, min, max, nBins, min, max);
  TH2F *SMSEventsMetJetDPhi10 = new TH2F("SMSEventsMetJetDPhi10", "MetJetDPhi 1.0 " + s, nBins, min, max, nBins, min, max);
  TH1F *hDeltaPhiMetTop = new TH1F("DeltaPhiMetTop", "DeltaPhiMetTop", 32, 0, 3.14);
  TH1F *hDeltaPhiMetTopTTBar = new TH1F("DeltaPhiMetTopTTBar", "DeltaPhiMetTopTTBar", 32, 0, 3.14);

  TH2F *h_SMSEvents;
  for(int ii = 0; ii <(int) fSamples.size(); ii++){
    sample Sample = fSamples[ii];

    if(Sample.name != mySample && Sample.sname != "Top")
      continue;

    Long64_t nentries = Sample.tree->GetEntries();
 
    std::cout << setfill('=') << std::setw(70) << "" << std::endl;
    cout << "looping over :     " <<endl;	
    cout << "   Name:           " << Sample.name << endl;
    cout << "   File:           " << (Sample.file)->GetName() << endl;
    cout << "   Events:         " << Sample.nevents  << endl;
    cout << "   Events in tree: " << nentries << endl; 
    cout << "   Xsection:       " << Sample.xsection << endl;
    cout << "   kfactor:        " << Sample.kfact << endl;
    cout << "   avg PU weight:  " << Sample.PU_avg_weight << endl;
    std::cout << setfill('-') << std::setw(70) << "" << std::endl;
 

    fMT2tree = new MT2tree();
    Sample.tree->SetBranchAddress("MT2tree", &fMT2tree);

    if(Sample.sname == "Top"){
      for (Long64_t jentry=0; jentry<nentries;jentry++) {
	Sample.tree->GetEntry(jentry); 
	if ( jentry % 100000 == 0 ){ 
	fprintf(stdout, "\rProcessed events: %6d of %6d ", jentry + 1, nentries);
	fflush(stdout);
	}
	if (!(fMT2tree->misc.MET                > 30                                    &&
	      fMT2tree->misc.MT2                > 130                                   &&
	      fMT2tree->NBJetsCSVM              > 0                                     &&
	      (fMT2tree->NEles==0  ||  fMT2tree->ele[0].lv.Pt()<5)                   &&
	      (fMT2tree->NMuons==0 ||  fMT2tree->muo[0].lv.Pt()<5)                   &&
	      fMT2tree->misc.Jet0Pass        == 1                                    &&
	      fMT2tree->misc.Jet1Pass        == 1                                    &&
	      // fMT2tree->misc.SecondJPt       >  100                                  &&
	      // fMT2tree->misc.LeadingJPt      >  150                                  &&
	      fMT2tree->misc.PassJetID       == 1                                      &&
	      //   	     fMT2tree->misc.Vectorsumpt        < 70                                    &&
	      // fMT2tree->misc.MinMetJetDPhi4            >  0.6                                  &&
	      ((fMT2tree->NJetsIDLoose40  >=4 && fMT2tree->jet[3].lv.Pt() >= 100.0) || (fMT2tree->NJetsIDLoose40  >=6 && fMT2tree->jet[5].lv.Pt() >= 65.0) ||
	       (fMT2tree->NJetsIDLoose40  >=5 && fMT2tree->jet[4].lv.Pt() >= 85.0) || (fMT2tree->NJetsIDLoose40  >=7 && fMT2tree->jet[6].lv.Pt() >= 55.0)) &&
	      fMT2tree->NJetsIDLoose40          >=4))
	  continue;

	if(fMT2tree->NTops                   >  0){
	  hDeltaPhiMetTopTTBar->Fill(Util::DeltaPhi(fMT2tree->fitTop[0].lv.Phi(), fMT2tree->misc.METPhi));
	}}
      continue;
    }
    h_SMSEvents    = (TH2F*) Sample.file->Get("h_SMSEvents");
	
    for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Sample.tree->GetEntry(jentry); 

      if ( jentry % 100000 == 0 ){ 
	fprintf(stdout, "\rProcessed events: %6d of %6d ", jentry + 1, nentries);
	fflush(stdout);}

      if (!(fMT2tree->misc.MET                > 30                                    &&
	    fMT2tree->NBJets                  > 0                                     &&
	    fMT2tree->misc.Vectorsumpt        < 70                                    &&
	    fMT2tree->NJetsIDLoose40          >=4))
	continue;
      

      if ((fMT2tree->NEles==0  ||  fMT2tree->ele[0].lv.Pt()<10)                   &&
	  (fMT2tree->NMuons==0 ||  fMT2tree->muo[0].lv.Pt()<10)                   &&
	  fMT2tree->misc.Jet0Pass        == 1                                    &&
	  fMT2tree->misc.Jet1Pass        == 1                                    &&
	  fMT2tree->misc.SecondJPt       >  100                                  &&
	  fMT2tree->misc.LeadingJPt      >  150                                  &&
	  fMT2tree->misc.PassJetID       == 1                                    
	  ){
	//Six Or Quad
	if((fMT2tree->NJetsIDLoose40  >=4 && fMT2tree->jet[3].lv.Pt() >= 120.0) || (fMT2tree->NJetsIDLoose40  >=6 && fMT2tree->jet[5].lv.Pt() >= 80.0)){ 
	  if((fMT2tree->misc.MinMetJetDPhi >0.3||fMT2tree->misc.MinMetJetDPhiIndex>3) && fMT2tree->misc.MT2  >  125)
	    SMSEventsAllCuts4or6->Fill(fMT2tree->Susy.MassGlu, fMT2tree->Susy.MassLSP);
      
	  if((fMT2tree->misc.MinMetJetDPhi >0.5||fMT2tree->misc.MinMetJetDPhiIndex>3) && fMT2tree->misc.MT2  >  100)
	    SMSEventsAllCuts4or6MT2_100DelPhi->Fill(fMT2tree->Susy.MassGlu, fMT2tree->Susy.MassLSP);
	}
	//HT Trigger
	if(fMT2tree->misc.HT             > 750){
	  if((fMT2tree->misc.MinMetJetDPhi >0.3||fMT2tree->misc.MinMetJetDPhiIndex>3) && fMT2tree->misc.MT2  >  125)
	    SMSEventsAllCutsHT->Fill(fMT2tree->Susy.MassGlu, fMT2tree->Susy.MassLSP);
      
	  if((fMT2tree->misc.MinMetJetDPhi >0.5||fMT2tree->misc.MinMetJetDPhiIndex>3) && fMT2tree->misc.MT2  >  100)
	    SMSEventsAllCutsHTMT2_100DelPhi->Fill(fMT2tree->Susy.MassGlu, fMT2tree->Susy.MassLSP);
	}
      }

      SMSEvents->Fill(fMT2tree->Susy.MassGlu, fMT2tree->Susy.MassLSP);
  
      if(fMT2tree->misc.HT              > 750)
	SMSEventsHT->Fill(fMT2tree->Susy.MassGlu, fMT2tree->Susy.MassLSP);
  
      if(fMT2tree->misc.MT2             > 125)
	SMSEventsMT2->Fill(fMT2tree->Susy.MassGlu, fMT2tree->Susy.MassLSP);

      if(fMT2tree->misc.MT2             > 100)
	SMSEventsMT2100->Fill(fMT2tree->Susy.MassGlu, fMT2tree->Susy.MassLSP);
  
      if((fMT2tree->NEles==0  ||  fMT2tree->ele[0].lv.Pt()<10) && (fMT2tree->NMuons==0 ||  fMT2tree->muo[0].lv.Pt()<10))
	SMSEventsLeptonVeto->Fill(fMT2tree->Susy.MassGlu, fMT2tree->Susy.MassLSP);
  
      if((fMT2tree->misc.Jet0Pass        == 1) && (fMT2tree->misc.Jet1Pass        == 1) && (fMT2tree->misc.PassJetID       == 1 ))
	SMSEventsJetPass->Fill(fMT2tree->Susy.MassGlu, fMT2tree->Susy.MassLSP);
  
      if(fMT2tree->misc.SecondJPt       >  100) 
	SMSEventsSecondJet->Fill(fMT2tree->Susy.MassGlu, fMT2tree->Susy.MassLSP);
  
      if(fMT2tree->misc.LeadingJPt      >  150)
	SMSEventsLeadingJet->Fill(fMT2tree->Susy.MassGlu, fMT2tree->Susy.MassLSP);
  
      if(fMT2tree->misc.MinMetJetDPhi >0.3||fMT2tree->misc.MinMetJetDPhiIndex>3)
	SMSEventsMetJetDPhi->Fill(fMT2tree->Susy.MassGlu, fMT2tree->Susy.MassLSP);

      if(fMT2tree->misc.MinMetJetDPhi >0.5||fMT2tree->misc.MinMetJetDPhiIndex>3)
	SMSEventsMetJetDPhi5->Fill(fMT2tree->Susy.MassGlu, fMT2tree->Susy.MassLSP);

      if(fMT2tree->misc.MinMetJetDPhi >1.0||fMT2tree->misc.MinMetJetDPhiIndex>3)
	SMSEventsMetJetDPhi10->Fill(fMT2tree->Susy.MassGlu, fMT2tree->Susy.MassLSP);

      if((fMT2tree->NJetsIDLoose40  >=4 && fMT2tree->jet[3].lv.Pt() >= 120.0) || (fMT2tree->NJetsIDLoose40  >=6 && fMT2tree->jet[5].lv.Pt() >= 80.0))
	SMSEvents4or6->Fill(fMT2tree->Susy.MassGlu, fMT2tree->Susy.MassLSP);
  
      if(fMT2tree->misc.MT2Top             > 125)
	SMSEventsMT2Top->Fill(fMT2tree->Susy.MassGlu, fMT2tree->Susy.MassLSP);

      if(fMT2tree->NTops                   >  0){
	SMSEventsNTop1->Fill(fMT2tree->Susy.MassGlu, fMT2tree->Susy.MassLSP);
    
      }
  
      if(fMT2tree->NTops                   >  1)
	SMSEventsNTop2->Fill(fMT2tree->Susy.MassGlu, fMT2tree->Susy.MassLSP);
    

      if (!(fMT2tree->misc.MET                > 30                                    &&
	    fMT2tree->misc.MT2                > 130                                   &&
	    fMT2tree->NBJetsCSVM              > 0                                     &&
	    (fMT2tree->NEles==0  ||  fMT2tree->ele[0].lv.Pt()<5)                   &&
	    (fMT2tree->NMuons==0 ||  fMT2tree->muo[0].lv.Pt()<5)                   &&
	    fMT2tree->misc.Jet0Pass        == 1                                    &&
	    fMT2tree->misc.Jet1Pass        == 1                                    &&
	    // fMT2tree->misc.SecondJPt       >  100                                  &&
	    // fMT2tree->misc.LeadingJPt      >  150                                  &&
	    fMT2tree->misc.PassJetID       == 1                                      &&
	    //   	     fMT2tree->misc.Vectorsumpt        < 70                                    &&
	    // fMT2tree->misc.MinMetJetDPhi4            >  0.6                                  &&
	    ((fMT2tree->NJetsIDLoose40  >=4 && fMT2tree->jet[3].lv.Pt() >= 100.0) || (fMT2tree->NJetsIDLoose40  >=6 && fMT2tree->jet[5].lv.Pt() >= 65.0) ||
	     (fMT2tree->NJetsIDLoose40  >=5 && fMT2tree->jet[4].lv.Pt() >= 85.0) || (fMT2tree->NJetsIDLoose40  >=7 && fMT2tree->jet[6].lv.Pt() >= 55.0)) &&
	    fMT2tree->NJetsIDLoose40          >=4))
	continue;
  
      if(fMT2tree->NTops                   >  0){
	hDeltaPhiMetTop->Fill(Util::DeltaPhi(fMT2tree->fitTop[0].lv.Phi(), fMT2tree->misc.METPhi));
      }
    }}
  
  if(mySample == "SMST2")
    h_SMSEvents->Rebin2D(5,5);
  else
    h_SMSEvents->Rebin2D(10,10);

  TH2F* PreSelected = (TH2F*)SMSEvents->Clone();
  PreSelected->SetName("PreSelected");
  
  TCanvas *MyC = new TCanvas("MyC","MyC");
  MyC->Divide(4,3);
  MyC->cd(1);
  SMSEvents->Divide(h_SMSEvents);
  SMSEvents->Draw("COLZ");
  MyC->cd(2);
  SMSEventsHT->Divide(PreSelected);
  SMSEventsHT->Draw("COLZ");
  MyC->cd(3);
  SMSEvents4or6->Divide(PreSelected);
  SMSEvents4or6->Draw("COLZ");
  MyC->cd(4);
  SMSEventsMT2->Divide(PreSelected);
  SMSEventsMT2->Draw("COLZ");
  MyC->cd(5);
  // SMSEventsMT2Top->Divide(PreSelected);
  // SMSEventsMT2Top->Draw("COLZ");
  SMSEventsMT2100->Divide(PreSelected);
  SMSEventsMT2100->Draw("COLZ");
  MyC->cd(6);
  SMSEventsSecondJet->Divide(PreSelected);
  SMSEventsSecondJet->Draw("COLZ");
  MyC->cd(7);
  SMSEventsLeadingJet->Divide(PreSelected);
  SMSEventsLeadingJet->Draw("COLZ");
  MyC->cd(8);
  SMSEventsMetJetDPhi->Divide(PreSelected);
  SMSEventsMetJetDPhi->Draw("COLZ");
  MyC->cd(9);
  SMSEventsMetJetDPhi5->Divide(PreSelected);
  SMSEventsMetJetDPhi5->Draw("COLZ");
  MyC->cd(10);
  SMSEventsMetJetDPhi10->Divide(PreSelected);
  SMSEventsMetJetDPhi10->Draw("COLZ");
  // SMSEventsLeptonVeto->Divide(PreSelected);
  // SMSEventsLeptonVeto->Draw("COLZ");
  // SMSEventsJetPass->Divide(PreSelected);
  // SMSEventsJetPass->Draw("COLZ");
  MyC->cd(11);
  SMSEventsNTop1->Divide(PreSelected);
  SMSEventsNTop1->Draw("COLZ");
  MyC->cd(12);
  SMSEventsNTop2->Divide(PreSelected);
  SMSEventsNTop2->Draw("COLZ");

  TCanvas *MyC3 = new TCanvas("MyC3","MyC3");
  hDeltaPhiMetTopTTBar->Scale(1.0/hDeltaPhiMetTopTTBar->Integral());
  hDeltaPhiMetTop->Scale(1.0/hDeltaPhiMetTop->Integral());
  hDeltaPhiMetTopTTBar->SetLineColor(4);
  hDeltaPhiMetTop->Draw();
  hDeltaPhiMetTopTTBar->Draw("sames");

  TCanvas *MyC2 = new TCanvas("MyC2","MyC2");
  MyC2->Divide(2,2);

  MyC2->cd(1);
  SMSEventsAllCuts4or6->Divide(h_SMSEvents);
  SMSEventsAllCuts4or6->Draw("COLZ");

  MyC2->cd(2);
  SMSEventsAllCuts4or6MT2_100DelPhi->Divide(h_SMSEvents);
  SMSEventsAllCuts4or6MT2_100DelPhi->Draw("COLZ");

  MyC2->cd(3);
  SMSEventsAllCutsHT->Divide(h_SMSEvents);
  SMSEventsAllCutsHT->Draw("COLZ");

  MyC2->cd(4);
  SMSEventsAllCutsHTMT2_100DelPhi->Divide(h_SMSEvents);
  SMSEventsAllCutsHTMT2_100DelPhi->Draw("COLZ");
}


void MassPlotter::TopStudy2(TString mySample, Long64_t nevents){
  double x[11] = {0.0, 30.0, 60.0, 90.0, 130.0, 170.0, 250.0, 350.0, 450.0, 1000.0, 1300.0}; 
  TH1F* hPtTop = new TH1F("PtTop","PtTop", 10, x);
  TH1F* hPtTopEff = new TH1F("PtTopEff","PtTopEff", 10, x);

  MT2GenParticle T1, T2, W1, W2, B1, B2, U1, U2, D1, D2;

  for(int ii = 0; ii <(int) fSamples.size(); ii++){
    sample Sample = fSamples[ii];
    if(Sample.name != mySample)
      //   //  if(Sample.sname != mySample)
      continue;

    // debug level for TopSearch:
    // = 0 if no printout at all
    // = 1 if only final statistics
    // = 2 if also short event per event output
    // = 3 if also longer event per event output
    const int debug = 0;
   
    TopSearch *myTopSearch = new TopSearch(); 
    myTopSearch->SetDebug(debug);
    myTopSearch->SetWithbtag(1); // = 1: btag preferred over chisq, = 2: correct b-jet requested for top
    myTopSearch->SetChisqTmax(2.0);
    myTopSearch->SetChisqWmax(2.0);
    myTopSearch->SetChisqTWlepmax(0.0);
    myTopSearch->SetChisqW1b(1.0);//chisqW1bpenlty
    myTopSearch->SetChisqWlept(1.0);//chisqWleptpenlty
    myTopSearch->SetRescaleW(1);// = 1: rescale W 4-vector to correct W mass, = 0: do not rescale

    fMT2tree = new MT2tree();

    Sample.tree->SetBranchAddress("MT2tree", &fMT2tree);
    Long64_t nentries =  Sample.tree->GetEntries();

    Long64_t minNentriesNevents = min(nentries, nevents);

    std::cout << setfill('=') << std::setw(70) << "" << std::endl;
    cout << "looping over :     " <<endl;	
    cout << "   Name:           " << Sample.name << endl;
    cout << "   File:           " << (Sample.file)->GetName() << endl;
    cout << "   Events:         " << Sample.nevents  << endl;
    cout << "   Events in tree: " << nentries << endl; 
    cout << "   Events to read: " << minNentriesNevents << endl; 
    cout << "   Xsection:       " << Sample.xsection << endl;
    cout << "   kfactor:        " << Sample.kfact << endl;
    cout << "   avg PU weight:  " << Sample.PU_avg_weight << endl;
    std::cout << setfill('-') << std::setw(70) << "" << std::endl;

    for (Long64_t jentry=0; jentry<minNentriesNevents;jentry++) {
     
      Sample.tree->GetEntry(jentry); 

      if ( fVerbose>2 && jentry % 100000 == 0 ){  
	fprintf(stdout, "\rProcessed events: %6d of %6d ", jentry + 1, nentries);
	fflush(stdout);}

      if (!(fMT2tree->misc.MET            >  30                                   &&
	    fMT2tree->NBJets               > 0                                     &&
	    (fMT2tree->NEles==0  ||  fMT2tree->ele[0].lv.Pt()<10)                   &&
	    (fMT2tree->NMuons==0 ||  fMT2tree->muo[0].lv.Pt()<10)                   &&
	    fMT2tree->misc.Jet0Pass        == 1                                    &&
	    fMT2tree->misc.Jet1Pass        == 1                                    &&
	    fMT2tree->misc.LeadingJPt      >  150                                  &&
	    fMT2tree->misc.SecondJPt       >  100                                  &&
	    fMT2tree->misc.PassJetID       == 1                                    &&
	    fMT2tree->misc.Vectorsumpt     <  70                                   &&
	    (fMT2tree->misc.MinMetJetDPhi >0.3||fMT2tree->misc.MinMetJetDPhiIndex>3) &&
	    ((fMT2tree->NJetsIDLoose40  >=4 && fMT2tree->jet[3].lv.Pt() >= 120.0) || (fMT2tree->NJetsIDLoose40  >=6 && fMT2tree->jet[5].lv.Pt() >= 80.0))
	    ))
	continue;
      int nPart = 0;
      //Read thegen info and find the particle indices 
      GetGenInfo(&T1, &T2, &W1, &W2, &B1, &B2, &U1, &U2, &D1, &D2);

      //Reject Fully Hadronic (quarks do not have mother)
      if(U1.Index == D2.Index)      
	continue;

      MT2GenParticle TopSemiLep;
      //To select semileptonic ttbar 
      if(abs(U1.Index + D1.Index + U2.Index + D2.Index) > 100){
	if(U1.Index == D1.Index){
	  TopSemiLep = T2;
	}else{
	  TopSemiLep = T1;}

	hPtTop->Fill(TopSemiLep.lv.Pt());

	//clear the memory
	nPart = 0;
	myTopSearch->SetNpart(nPart);   


	for(int i=0; i<fMT2tree->NJets; ++i){
	  if(fMT2tree->jet[i].IsGoodPFJet(20.0, 2.4, 1) ==false) continue;
	  InitializeTopSearch(myTopSearch, nPart, fMT2tree->jet[i].lv, fMT2tree->jet[i].IsBJet(3));
	  nPart++;
	}

	myTopSearch->SetNpart(nPart);

	myTopSearch->GoSearch();
       
	float minDR = 100.0;
	int MatchedRecTopIndex = -1;
  
	for(int i = 0; i < myTopSearch->NpT(); i++){
	  TLorentzVector Top = myTopSearch->PT(i);
	 
	  if(TopSemiLep.lv.DeltaR(Top) < minDR){
	    minDR = TopSemiLep.lv.DeltaR(Top);
	    MatchedRecTopIndex = i;
	  }
	}

	if(minDR < 1.0){
	  hPtTopEff->Fill(TopSemiLep.lv.Pt());
	}
      }}}

  TCanvas *MyC = new TCanvas("MyC", "MyC");
  // hPtTopEff->Divide(hPtTop);
  //hPtTopEff->Draw();

  TGraphAsymmErrors *gPtTopEff = new TGraphAsymmErrors(hPtTopEff, hPtTop);
  gPtTopEff->GetXaxis()->SetTitle("Top Pt (GeV)");
  gPtTopEff->GetYaxis()->SetTitle("Top Finding Eff");
  gPtTopEff->Draw("AP");
}


void MassPlotter::QCD(){
  
  TH1::SetDefaultSumw2();
  double x[14] = {0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0, 125.0, 500.0, 2000.0}; 
  TH1F *MT2HighDPhi = new TH1F("MT2HighDPhi", "MT2HighDPhi", 13, x);//40,   0, 1000);//
  TH1F *MT2LowDPhi = new TH1F("MT2LowDPhi", "MT2LowDPhi", 13, x);//40,   0, 1000);//


    for(int ii = 0; ii <(int) fSamples.size(); ii++){
    sample Sample = fSamples[ii];
//    if(Sample.sname == "QCD" || Sample.sname == "SMS")
//       continue;
    if(Sample.type != "data")
      continue;

    fMT2tree = new MT2tree();

    Sample.tree->SetBranchAddress("MT2tree", &fMT2tree);
    Long64_t nentries =  Sample.tree->GetEntries();
    float Weight = Sample.xsection * Sample.kfact * Sample.lumi / (Sample.nevents * Sample.PU_avg_weight );
 
    if(Sample.type == "mc")
      Weight *= -1.0;
    
    std::cout << setfill('=') << std::setw(70) << "" << std::endl;
    cout << "looping over :     " <<endl;	
    cout << "   Name:           " << Sample.name << endl;
    cout << "   File:           " << (Sample.file)->GetName() << endl;
    cout << "   Events:         " << Sample.nevents  << endl;
    cout << "   Events in tree: " << nentries << endl; 
    cout << "   Xsection:       " << Sample.xsection << endl;
    cout << "   kfactor:        " << Sample.kfact << endl;
    cout << "   avg PU weight:  " << Sample.PU_avg_weight << endl;
    cout << "   Weight:         " << Weight <<endl;
    std::cout << setfill('-') << std::setw(70) << "" << std::endl;

    for (Long64_t jentry=0; jentry<nentries;jentry++) {
     
      Sample.tree->GetEntry(jentry); 
      if ( fVerbose>2 && jentry % 100000 == 0 ){ 
	fprintf(stdout, "\rProcessed events: %6d of %6d ", jentry + 1, nentries);
	fflush(stdout);
      }
      
      if (!((fMT2tree->NEles==0  ||  fMT2tree->ele[0].lv.Pt()<5)                   &&
	    (fMT2tree->NMuons==0 ||  fMT2tree->muo[0].lv.Pt()<5)))
	continue;

      if(fMT2tree->misc.MinMetJetDPhi4 > 0.3)
	MT2HighDPhi->Fill(fMT2tree->misc.MT2, Weight);
      else
	//if(fMT2tree->misc.MinMetJetDPhi4 <= 0.2)
	MT2LowDPhi->Fill(fMT2tree->misc.MT2, Weight);
    }}

    TCanvas *MyC = new TCanvas("MyC", "MyC");
    MyC->Divide(2,1);

    MyC->cd(1);
    cout<<"Low  "<<MT2LowDPhi->GetBinContent(14)<<" "<<MT2LowDPhi->GetBinError(14)<<endl;
    MT2LowDPhi->Draw();
    TH1F* hMT2HighDPhi = (TH1F*) MT2HighDPhi->Clone("hMT2HighDPhi");
    hMT2HighDPhi->SetName("hMT2HighDPhi");
    hMT2HighDPhi->SetLineColor(4);
    hMT2HighDPhi->Draw("sames");

    MyC->cd(2);
    cout<<"High "<<MT2HighDPhi->GetBinContent(14)<<" "<<MT2HighDPhi->GetBinError(14)<<endl;
    MT2HighDPhi->Divide(MT2LowDPhi);
    cout<<"Divi "<<MT2HighDPhi->GetBinContent(14)<<" "<<MT2HighDPhi->GetBinError(14)<<endl;
    MT2HighDPhi->Draw();
}



void MassPlotter::vs(){ 
 
         TH2F *METvsVSPT = new TH2F("METvsVSPT", "METvsVSPT", 100, 0,1000, 100,0,1000);

         for(int i = 0; i <(int) fSamples.size(); i++){
         sample Sample = fSamples[i];

    fMT2tree = new MT2tree();
 
    Sample.tree->SetBranchAddress("MT2tree", &fMT2tree);
    Long64_t n =  Sample.tree->GetEntries();
    for (Long64_t j=0; j<10000;j++) {
     
      Sample.tree->GetEntry(j); 
               METvsVSPT->Fill(fMT2tree->misc.MET, fMT2tree->misc.Vectorsumpt);
    }
}

    TCanvas *MyC = new TCanvas("MyC", "MyC");

    MyC->cd(1);
   METvsVSPT->Draw();
}

void MassPlotter::SpecialMakePlot(int nevents, TString cuts, TString trigger){
  TH1::SetDefaultSumw2();

  TH1F *NBJets[NumberOfSamples+1];
  TString  cnames[NumberOfSamples+1] = {"QCD", "Wjets", "Zjets", "Top", "MC", "data", "SMS"};
  int      ccolor[NumberOfSamples+1] = {  401,     417,     419,   855,  500,    632,     0};
  TString varname = "NBJetsCSVT";
  for (int i=0; i<(NumberOfSamples+1); i++){
    NBJets[i] = new TH1F(varname+"_"+cnames[i], "", 4, 0, 4);
    NBJets[i] -> SetFillColor  (ccolor[i]);
    NBJets[i] -> SetLineColor  (ccolor[i]);
    NBJets[i] -> SetLineWidth  (2);
    NBJets[i] -> SetMarkerColor(ccolor[i]);
    NBJets[i] -> SetStats(false);
  }

  NBJets[5] -> SetMarkerStyle(20);
  NBJets[5] -> SetMarkerColor(kBlack);
  NBJets[5] -> SetLineColor(kBlack);
  
  NBJets[4] -> SetFillStyle(3004);
  NBJets[4] -> SetFillColor(kBlack);

  NBJets[6] -> SetLineColor(kBlack);
  NBJets[6] -> SetLineStyle(kDotted);

  cout<<" trigger "<<trigger<<endl;
  cout<<" cuts "<<cuts<<endl;

  for(int ii = 0; ii < fSamples.size(); ii++){

//     if(fSamples[ii].type != "data" || fSamples[ii].name != "Run2012D1F")
//       continue;
    float Weight = fSamples[ii].xsection * fSamples[ii].kfact * fSamples[ii].lumi / (fSamples[ii].nevents);

     std::cout << setfill('=') << std::setw(70) << "" << std::endl;
    cout << "looping over :     " <<endl;	
    cout << "   Name:           " << fSamples[ii].name << endl;
    cout << "   File:           " << (fSamples[ii].file)->GetName() << endl;
    cout << "   Events:         " << fSamples[ii].nevents  << endl;
    cout << "   Events in tree: " << fSamples[ii].tree->GetEntries() << endl; 
    cout << "   Xsection:       " << fSamples[ii].xsection << endl;
    cout << "   kfactor:        " << fSamples[ii].kfact << endl;
    cout << "   avg PU weight:  " << fSamples[ii].PU_avg_weight << endl;
    cout << "   Weight:         " << Weight <<endl;
    std::cout << setfill('-') << std::setw(70) << "" << std::endl;

    fMT2tree = new MT2tree();
    fSamples[ii].tree->SetBranchAddress("MT2tree", &fMT2tree);

    TString myCuts = cuts;

    if(fSamples[ii].type == "data"){
      myCuts += " && " + trigger;
    }


    fSamples[ii].tree->Draw(">>selList", myCuts);
    TEventList *myEvtList = (TEventList*)gDirectory->Get("selList");
    fSamples[ii].tree->SetEventList(myEvtList);


    Int_t nentries =  myEvtList->GetN();


    for (Int_t i=0;i<min(nentries, nevents); i++) {
      if(i > nevents) continue;

      if ( i % 100000 == 0 ){ 	
	fprintf(stdout, "\rProcessed events: %6d of %6d ", i + 1, nentries);
	fflush(stdout);
      }
 
      fSamples[ii].tree->GetEntry(myEvtList->GetEntry(i));
 
     float weight = Weight;
     if(fSamples[ii].type == "data")
       weight = 1.0;
     else{
       if(fSamples[ii].type == "susy")
	 weight *= 1.0;
       else{
	 weight *= 1.0;//(fMT2tree->pileUp.Weight/fSamples[ii].PU_avg_weight);//
	 if(fMT2tree->NBJetsCSVT == 0)
	   weight *=fMT2tree->SFWeight.BTagCSV40eq0;
	 if(fMT2tree->NBJetsCSVT == 1)
	   weight *=fMT2tree->SFWeight.BTagCSV40eq1;
	 if(fMT2tree->NBJetsCSVT == 2)
	   weight *=fMT2tree->SFWeight.BTagCSV40eq2;
	 if(fMT2tree->NBJetsCSVT >= 3)
	   weight *=fMT2tree->SFWeight.BTagCSV40ge3;
       }}
   
     if(fSamples[ii].type == "data"){
       NBJets[5]->Fill(fMT2tree->NBJetsCSVT, weight);//data
     }else{
       if(fSamples[ii].type == "susy")
	 NBJets[6]->Fill(fMT2tree->NBJetsCSVT, weight);
       else
	 NBJets[4]->Fill(fMT2tree->NBJetsCSVT, weight);
	    
       if(fSamples[ii].sname == "Top")  
	 NBJets[3]->Fill(fMT2tree->NBJetsCSVT, weight);
       else
	 if(fSamples[ii].sname == "DY")	
	   NBJets[2]->Fill(fMT2tree->NBJetsCSVT, weight);
	 else
	   if(fSamples[ii].sname == "Wtolnu")
	     NBJets[1]->Fill(fMT2tree->NBJetsCSVT, weight);
	   else
	     if(fSamples[ii].sname == "QCD")
	       NBJets[0]->Fill(fMT2tree->NBJetsCSVT, weight);
     }
     
   
    }
    cout << endl;
  }
  
  THStack* h_stack     = new THStack(varname, "");
  for(int j = 0; j < (NumberOfSamples+ 1); j++){
    AddOverAndUnderFlow(NBJets[j]);

    if(j < (NumberOfSamples - 2))
      h_stack  -> Add(NBJets[j]);
  }  

  TLegend* Legend1 = new TLegend(.71,.54,.91,.92);
  Legend1->AddEntry(NBJets[0], "QCD", "f");
  Legend1->AddEntry(NBJets[1], "W+jets", "f");
  Legend1->AddEntry(NBJets[2], "Z+jets", "f");
  Legend1->AddEntry(NBJets[3], "Top", "f");
  Legend1->AddEntry(NBJets[5], "data", "l");
  Legend1->AddEntry(NBJets[6], "SMS", "l");

  cout << "------------------------------------"                << endl
       << "QCD Integral:      " << NBJets[0]->Integral()  << endl
       << "W+jets Integral:   " << NBJets[1]->Integral()  << endl
       << "Z+jets Integral:   " << NBJets[2]->Integral()  << endl
       << "Top Integral:      " << NBJets[3]->Integral()  << endl
       << "TOTAL BG:          " << NBJets[4]->Integral()  << endl
       << "SUSY:              " << NBJets[6]->Integral()  << endl
       << "Data:              " << NBJets[5]->Integral()  << endl;
        

  cout<< "SUSY1:              " << NBJets[6]->Integral()  << endl;

  //  TLegend *Legend1;
  printHisto(h_stack, NBJets[5], NBJets[4], NBJets[6],Legend1 , "MTC", "hist", true, varname, "Events", -4, 0, 1,1);
  cout<< "SUSY2:              " << NBJets[6]->Integral()  << endl;

  plotRatioStack(h_stack,  NBJets[4], NBJets[5], NBJets[6],true, false, "MT2_ratio", Legend1, varname, "Events", -4, -10, 0, 1);
  cout<< "SUSY3:              " << NBJets[6]->Integral()  << endl;

}


void MassPlotter::MakeCutFlowTable( std::vector<std::string> all_cuts ){

  TH1::SetDefaultSumw2();
  std::vector< TH1* > CutFlowHistos;
  string TotalBkgtitle = "Total Bkg";
  TH1 *hTotalBkg = new TH1D(TotalBkgtitle.c_str(), TotalBkgtitle.c_str(), all_cuts.size() , 0 , all_cuts.size() );
  for(int ii = 0; ii < fSamples.size(); ii++){

    string hName = "h" + string(fSamples[ii].name.Data()) + "_cutflow" ;
    TH1* hCutFlowSample = new TH1D( hName.c_str() , fSamples[ii].sname  , all_cuts.size() , 0 , all_cuts.size() );

    Double_t weight=0;
    if(fPUReweight) weight = fSamples[ii].xsection * fSamples[ii].kfact * fSamples[ii].lumi / (fSamples[ii].nevents*fSamples[ii].PU_avg_weight);
    else            weight = fSamples[ii].xsection * fSamples[ii].kfact * fSamples[ii].lumi / (fSamples[ii].nevents);
    
    TTree* theTree = fSamples[ii].tree ;
    cout << "|" << endl << fSamples[ii].name ;
    TEventList* latest_event_list = NULL;

    string full_cut = "1 ";
    int cut_number = 0 ; 

    for( vector<string>::const_iterator cut = all_cuts.begin(); cut != all_cuts.end() ; cut++){

      full_cut += ("&&" + *cut);

      TString btagweight = "1.00"; //stored btag weights up to >=3, default is weight==1 to avoid non-existing weights
      if(full_cut.find("NBJetsCSVT") != string::npos ) btagweight = TString::Format("SFWeight.BTagCSV40ge%d",abs(1));

      TString selection;
      if(     fSamples[ii].type!="data" && fPUReweight && fbSFReWeight) selection = TString::Format("(%.15f*pileUp.Weight*%s)",weight, btagweight.Data());
      else if(fSamples[ii].type!="data" && fPUReweight                ) selection = TString::Format("(%.15f*pileUp.Weight)",   weight);
      else if(fSamples[ii].type!="data" &&                fbSFReWeight) selection = TString::Format("(%.15f*%s)",              weight, btagweight.Data());
      else                                                            selection = TString::Format("(%.15f)",                 weight); 
    
      if(     fSamples[ii].type=="susy")  selection = TString::Format("(%.15f)",weight);
      string full_weight = selection.Data();

      if( cut_number == 0 )
	theTree->Draw( ">>event_list" ,   full_cut.c_str() );
      else{
	theTree->SetEventList( latest_event_list );
	theTree->Draw( ">>event_list" ,  full_cut.c_str() );
      }
      latest_event_list = (TEventList*)( gDirectory->Get("event_list")->Clone("tmp_event_list") );
      theTree->SetEventList( latest_event_list );
      theTree->Draw( "1>>temp_hist" , full_weight.c_str() );
      TH1* hTemp = (TH1*)( gDirectory->Get("temp_hist") );
      double err_yields;
      double yields = hTemp->IntegralAndError( 0 , hTemp->GetNbinsX() , err_yields );
      hCutFlowSample->SetBinContent( cut_number+1 , yields );
      hCutFlowSample->SetBinError( cut_number+1 , err_yields );
      cout << "|" << yields << "\\pm" << err_yields << flush;
      
      cut_number++;
    }

 
    bool was_found_in_snames = false;
    for( std::vector< TH1* >::iterator sname_histo = CutFlowHistos.begin() ; sname_histo != CutFlowHistos.end() ; sname_histo++){
      if( string((*sname_histo)->GetTitle()) == string(fSamples[ii].sname.Data()) ){
	(*sname_histo)->Add( hCutFlowSample );
	was_found_in_snames = true;
      }
    }
    if( !was_found_in_snames ){
      string hName = "h" + string(fSamples[ii].sname.Data()) + "_cutflow" ;
      gROOT->cd();
      CutFlowHistos.push_back( (TH1*)(hCutFlowSample->Clone( hName.c_str() )) );
    }
    if(!(fSamples[ii].sname.Contains("ata") || fSamples[ii].sname.Contains("SUSY")))
      hTotalBkg->Add(hCutFlowSample);

    delete hCutFlowSample;
  }
  std::vector< TH1* >::iterator TotalBkgPlace = CutFlowHistos.end();
  TotalBkgPlace--;

  CutFlowHistos.insert(TotalBkgPlace, hTotalBkg);

  for( std::vector< TH1* >::iterator sname_histo = CutFlowHistos.begin() ; sname_histo != CutFlowHistos.end() ; sname_histo++){
    cout << endl << (*sname_histo)->GetTitle() << "&" ;
    for( int i = 1 ; i <= (*sname_histo)->GetNbinsX() ; i++){
      cout << (*sname_histo)->GetBinContent(i) << "\\pm" << (*sname_histo)->GetBinError(i) << "&" ;
    }
    cout << endl;
  }

  for( int i = 1 ; i <=  all_cuts.size(); i++){
    if(i == 1)
      cout <<  "          " << "&" ;
    else   
      cout << (all_cuts[i-1]) << " &" ;
    for( std::vector< TH1* >::iterator sname_histo = CutFlowHistos.begin() ; sname_histo != CutFlowHistos.end() ; sname_histo++){
      if(i == 1)
	cout <<  (*sname_histo)->GetTitle() << "&" ;
      else{
	if((TString) (*sname_histo)->GetTitle() != (TString) TotalBkgtitle.c_str())
	  cout << (*sname_histo)->GetBinContent(i) << "&" ;
	else
	  cout << (*sname_histo)->GetBinContent(i) << "\\pm" << (*sname_histo)->GetBinError(i) << "&" ;
      }}
    cout << endl;
  }
}










	double MassPlotter::DeltaPhi(double phi1, double phi2){
		// From cmssw reco::deltaPhi()
		double result = phi1 - phi2;
		while( result >   TMath::Pi() ) result -= TMath::TwoPi();
		while( result <= -TMath::Pi() ) result += TMath::TwoPi();
		return TMath::Abs(result);
	}
 
