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
#include <vector>
#include "TROOT.h"
//#include "../MT2Code/include/MT2Shapes.hh"
using namespace std;
using namespace RooFit ;

// options ---------------
Int_t    fVerbose                   = 6;
Bool_t   fSaveResults               = true;
Bool_t   fDraw                      = true;
TString  fOutDir                    = "../Results/GammaJetsPredictionFromZll/20110117";
TString  fRootFile                  = "Plots.root";
Bool_t   fWriteToFile               = false; // writes couts to file
// Cfg ----------------------------------------------------------------------------
TString fSamplesPhotonPt           ="../RootMacros/samples/samples_RemovedPhotons.dat";
TString fSamplesPhotonPtMConly     ="../RootMacros/samples/samples_RemovedPhotonsMC.dat";
TString fSamplesZllPt              ="../RootMacros/samples/samples_BoostedZll.dat";
TString fSamplesZllPtMConly        ="../RootMacros/samples/samples_BoostedZllMC.dat";
TString fSamplesZnunu              ="../RootMacros/samples/samples_Znunu.dat";
// stearing ------------
Bool_t fDoData                       = false;
Bool_t fDoDataZllToPhotonRatio       = false;     // take Zll to Photon ratio from data (or MC)
Bool_t fDoPrediction                 = false;
Bool_t fDoFits                       = false;    // take ratio from fits rather than ratio of integrals
Bool_t fDoPtSpectraComparison        = true;     // produce fancy pt comparison plot
Bool_t fDoPtSpectraComparisonAccCorr = false;     // acceptance efficiency correction for Zll in Zll/photon pt plot
Bool_t fDoPtSpectraComparisonScaling = true; // scale normalization of MC pt spectra tp match data
Bool_t fDoMT2SpectraCompaison        = false;
Bool_t fDoMT2SpectraCompaisonScaling = false;
Bool_t fDoZnunuGammaGenPtSpectraComparison = false; // gen-level Znunu to Photon Pt spectra comparison
// stearing of steps
Bool_t fDoPhotonPtShape          = true;
Bool_t fDoZllPtShape             = true;
Bool_t fDoGenZllShape            = false;
Bool_t fDoGenAccZllShape         = false;
Bool_t fDoGenAccZllRecoShape     = false;
Bool_t fDoPhotonMT2Shape         = false;
Bool_t fDoZnunuMT2Shape          = false;
Bool_t fDoPhotonRecoEff          = false;
// 
Float_t fZllToGammaFitMin        = 350;
Float_t fZllToGammaFitMax        = 800;
Float_t fZllAccFitMin            = 350;
Float_t fZllAccFitMax            = 800;
Float_t fZllRecoEffFitMin        = 300;
Float_t fZllRecoEffFitMax        = 800;
// global vars
Bool_t  fDoRatioStatErr          = false;
//Float_t fZllToPhotonRatio        = 0.0827;
Float_t fZllToPhotonRatio        = 0.0860; // MC value
Float_t fZllToPhotonRatioRelErr  = 0.15; // rel uncertainty
Float_t fZllAccEff               = 0.7838;
Float_t fZllAccEffRelErr         = 0.1;
Float_t fZllRecoEff              = 0.817;
Float_t fZllRecoEffRelErr        = 0.1;
Float_t fBRZnunu                 = 0.2;
Float_t fBRZll                   = 0.06726;
Float_t fLumiCorr                = 4400./4600;
// cuts
Float_t fHTmin                     = 750;
Float_t fHTmax                     = 1E+8;
Float_t fMT2min                    = 0;
Float_t fMT2max                    = 1E+8;


// cutstreams
std::ostringstream* fLogStream     = 0;
std::ostringstream cutStream_PhotonPt;
std::ostringstream cutStreamZll;
std::ostringstream cutStreamGenZll;
std::ostringstream cutStreamGenZllAcc;
std::ostringstream cutStreamGenZllAcc_recoed;
std::ostringstream cutStreamZnunuMT2;
std::ostringstream cutStreamPhotonMT2;
std::ostringstream cutStreamGenPhoton;
std::ostringstream cutStreamGenZnunu;
std::ostringstream fTriggerStream;
std::ostringstream cutStream_PhotonPtforRecoEff;
std::ostringstream cutStream_GenPhotonforRecoEff;

// stuff 
double doubledummy = 0;

TFile * fOutFile=0;
TDirectory *fDir=0;

void DefineCutStreams(float minHT, float maxHT, float minMT2, float maxMT2){
	cutStream_PhotonPt << " " 
	<< "misc.MET >=0"                              << "&&"
	<< "NPhotons==1"                               << "&&"
	<< "photon[0].lv.Pt()>150"                     << "&&"
	<< "photon[0].lv.Pt()<800"                     << "&&"
        << "abs(photon[0].lv.Rapidity())<1.4"          << "&&"
	<< "photon[0].isEGMlooseID==1"                 << "&&"
	<< "photon[0].isEGMlooseRelIso==1"             << "&&"
	<< "misc.CrazyHCAL==0";
	
	cutStream_PhotonPtforRecoEff << " " 
	<< "misc.MET >=0"                              << "&&"
	<< "NPhotons==1"                               << "&&"
	<< "photon[0].lv.Pt()>150"                     << "&&"
        << "abs(photon[0].lv.Rapidity())<2.4"          << "&&"
	<< "photon[0].isEGMlooseID==1"                 << "&&"
	<< "photon[0].isEGMlooseRelIso==1"             << "&&"
	<< "GenPhoton[0].Pt()>150"                     << "&&"
	<< "abs(GenPhoton[0].Rapidity())<2.4"          << "&&"
	<< "misc.CrazyHCAL==0";
	
	cutStream_GenPhotonforRecoEff << " "
	<< "misc.MET>=0"                                          << "&&"
	<< "GenPhoton[0].Pt()>150"                                << "&&"
	<< "abs(GenPhoton[0].Rapidity())<2.4"                     << "&&"
	<< "misc.CrazyHCAL==0";
	
	cutStreamZll << " " 
	<< "misc.MET >=0"                                      << "&&"
	<< "RecoOSDiLeptPt(20,2.4,71,111)>150"                 << "&&"
	<< "RecoOSDiLeptPt(20,2.4,71,111)<800"                 << "&&"
        << "abs(RecoOSDiLeptRapidity(20,2.4,71,111))<1.4"      << "&&"
	<< "misc.CrazyHCAL==0";
	
	cutStreamGenZll << " " 
	<< "misc.MET >=0"                                      << "&&"
	<< "GenDiLeptPt(0,10,71,111,true)>150"                 << "&&"
        << "abs(GenDiLeptRapidity(0,10,71,111,true))<1.4"      << "&&"
	<< "misc.CrazyHCAL==0";
	
	cutStreamGenZnunu << " " 
	<< "misc.MET >=0"                                      << "&&"
	<< "GenZ[0].Pt()>150"                                  << "&&"
        << "abs(GenZ[0].Rapidity())<2.4"                       << "&&"
	<< "misc.CrazyHCAL==0";
	
	cutStreamGenZllAcc << " " 
	<< "misc.MET >=0"                                         << "&&"
	<< "GenDiLeptPt(20,2.4,71,111,true)>150"                  << "&&"
        << "abs(GenDiLeptRapidity(20,2.4,71,1111,true))<1.4"      << "&&"
	<< "misc.CrazyHCAL==0";
	
	cutStreamGenZllAcc_recoed << " " 
	<< "misc.MET >=0"                                         << "&&"
	<< "GenDiLeptPt(20,2.4,71,111,true)>150"                  << "&&"
        << "abs(GenDiLeptRapidity(20,2.4,71,1111,true))<1.4"      << "&&"
	<< "RecoOSDiLeptPt(20,2.4,71,111)>150"                    << "&&"
	<< "misc.CrazyHCAL==0";

	cutStreamGenPhoton << " "
	<< "misc.MET>=0"                                          << "&&"
	<< "GenPhoton[0].Pt()>150"                                << "&&"
	<< "abs(GenPhoton[0].Rapidity())<2.4"                     << "&&"
	<< "misc.CrazyHCAL==0";
	
	cutStreamPhotonMT2 << " " 
	  << "misc.MET>=150"                                               << "&&"
	  << "misc.HT >=" << minHT  << "&&" << "misc.HT<" << maxHT         << "&&"
	  << "misc.MT2>=" << minMT2 << "&&" << "misc.MT2<" << maxMT2       << "&&" 
	  << "(NEles==0  || ele[0].lv.Pt()<10)"                             << "&&"
	  << "(NMuons==0 || muo[0].lv.Pt()<10)"                             << "&&"
	  << "misc.Jet0Pass ==1"                                            << "&&"
	  << "misc.Jet1Pass ==1"                                            << "&&"
	  << "misc.SecondJPt >=100"                                          << "&&"
	  << "misc.PassJetID ==1"                                           << "&&"
	  << "misc.Vectorsumpt < 70"                                        << "&&"
	  << "NPhotons ==1 "                                                << "&&"
//	  << "abs(photon[0].lv.Rapidity())<2.4"                             << "&&"
	  << "photon[0].isEGMlooseID==1"                                    << "&&"
	  << "photon[0].isEGMlooseRelIso==1"                                << "&&"
	  << "rawpfmet[0].Pt()<100"                                         << "&&"
	  << "misc.HBHENoiseFlagIso ==0"                                       << "&&"
	  << "misc.CSCTightHaloID==0"                                       << "&&"
	  << "NJetsIDLoose40 >=3"                                           << "&&"
	  << "misc.MinMetJetDPhi >0.3"                                      << "&&"
	  << "misc.CrazyHCAL==0";
	
	cutStreamZnunuMT2 << " " 
	  << "misc.MET>=150"                                               << "&&"
	  << "misc.HT >=" << minHT  << "&&" << "misc.HT<" << maxHT         << "&&"
	  << "misc.MT2>=" << minMT2 << "&&" << "misc.MT2<" << maxMT2       << "&&" 
	  << "(NEles==0  || ele[0].lv.Pt()<10)"                            << "&&"
	  << "(NMuons==0 || muo[0].lv.Pt()<10)"                            << "&&"
	  << "misc.Jet0Pass ==1"                                           << "&&"
	  << "misc.Jet1Pass ==1"                                           << "&&"
	  << "misc.SecondJPt  >100"                                        << "&&"
	  << "misc.PassJetID ==1"                                          << "&&"
	  << "misc.Vectorsumpt < 70"                                       << "&&"
//	  << "MinGenBosonJetsDR()>0.5"                                     << "&&"
	  << "misc.HBHENoiseFlagIso ==0"                                   << "&&"
	  << "misc.CSCTightHaloID==0"                                      << "&&"
	  << "NJetsIDLoose40 >=3"                                          << "&&"
 	  << "misc.MinMetJetDPhi >0.3"                                     << "&&"
	  << "misc.CrazyHCAL==0";

	fTriggerStream << "";
}

void run_DataGammaJetsZllToZnunu(){

	// defs ---------------------------------------------------------
	gSystem->CompileMacro("../MT2Code/src/MT2Shapes.cc", "k");
	
	// logStream
	fLogStream = new std::ostringstream();
	
	// create dir
	if(!fOutDir.EndsWith("/")) fOutDir += "/";
	char cmd[500];
	sprintf(cmd,"mkdir -p %s", fOutDir.Data());
	system(cmd);

	DefineCutStreams(fHTmin, fHTmax, fMT2min, fMT2max);
	TString filename=fOutDir+"/"+fRootFile;
	fOutFile = new TFile(filename.Data(), "RECREATE");
	fDir     = (TDirectory*) fOutFile;

	// fix output dir
//	if(fHTmax <10000) fOutDir= TString::Format("%s_%d_HT_%d",  fOutDir.Data(), abs(fHTmin),  abs(fHTmax));
//	else              fOutDir= TString::Format("%s_%d_HT_%s",  fOutDir.Data(), abs(fHTmin),  "Inf");
//	if(fMT2max<10000) fOutDir= TString::Format("%s_%d_MT2_%d", fOutDir.Data(), abs(fMT2min), abs(fMT2max));
//	else              fOutDir= TString::Format("%s_%d_MT2_%s", fOutDir.Data(), abs(fMT2min), "Inf");
	
	// log MT2 and HT cuts
	*fLogStream << "------------------------------------------------------------------------------------------------" << endl;
	*fLogStream << "+++ new Znunu with Gamma+jets prediction                                                     +++" << endl;
	*fLogStream << "+++ outputdir: " << fOutDir <<                                                              "+++" << endl; 
	*fLogStream << "------------------------------------------------------------------------------------------------" << endl;
	
	// new prediction ------------------------------------------------------
	Prediction* prediction = new Prediction();
	prediction->fVerbose=fVerbose;
	prediction->fSave   =fSaveResults;


	// Photon Pt
	if(fDoPhotonPtShape){
  	const int gNMT2bins                   = 11;
  	double  gMT2bins[gNMT2bins+1]   = {150, 160, 170, 180, 190, 200, 225, 250, 300, 400, 550, 800}; 	
	prediction->ChPhotonPt = new Channel("PhotonPt", "photon[0].lv.Pt()", cutStream_PhotonPt.str().c_str(), fTriggerStream.str().c_str(),fSamplesPhotonPt);
	prediction->ChPhotonPt->fVerbose =prediction->fVerbose;
//	prediction->ChPhotonPt->GetShapes("PhotonPt", "#gamma Pt", 2, 300, 800);
	prediction->ChPhotonPt->GetShapes("PhotonPt", "#gamma Pt", 100, 150, 800);
//	prediction->ChPhotonPt->GetShapes("PhotonPt", "#gamma Pt", gNMT2bins, gMT2bins);
	}
	
	// Zll Pt
	if(fDoZllPtShape){
  	const int gNMT2bins                   = 11;
  	double  gMT2bins[gNMT2bins+1]   = {150, 160, 170, 180, 190, 200, 225, 250, 300, 400, 550, 800}; 	
	prediction->ChZllPt = new Channel("ZllPt", "RecoOSDiLeptPt(20,2.4,71,111)", cutStreamZll.str().c_str(), fTriggerStream.str().c_str(),fSamplesZllPt);
	prediction->ChZllPt->fVerbose =prediction->fVerbose;
//	prediction->ChZllPt->GetShapes("ZllPt", "Zll Pt", 2, 300, 800);
	prediction->ChZllPt->GetShapes("ZllPt", "Zll Pt", 100, 150, 800);
//	prediction->ChZllPt->GetShapes("ZllPt", "Zll Pt", gNMT2bins, gMT2bins);
	}

	// compute Zll/gamma pt ratio
	if(fDoPhotonPtShape && fDoZllPtShape){
		TH1D* hPhotonToZllPtRatio = prediction->GetRatio(fDoDataZllToPhotonRatio? prediction->ChZllPt->hData    : prediction->ChZllPt->hZJetsToLL, 
				                                 fDoDataZllToPhotonRatio? prediction->ChPhotonPt->hData : prediction->ChPhotonPt->hPhotons, 4);
		DrawHisto(hPhotonToZllPtRatio,hPhotonToZllPtRatio->GetName(),"EX0");
		TString rationame=hPhotonToZllPtRatio->GetName();
		rationame +="_fittedRatio";
		TF1 *f_lin = new TF1(rationame,"pol0(0)", fZllToGammaFitMin , fZllToGammaFitMax);   f_lin->SetLineColor(8);
		if(fDoFits){ 
			hPhotonToZllPtRatio->Fit(rationame,"0L","", fZllToGammaFitMin, fZllToGammaFitMax);    // set al weights to 1
			fZllToPhotonRatio   = f_lin->GetParameter(0);
		} else{
			fZllToPhotonRatio   = prediction->GetLimitedRatio(fDoDataZllToPhotonRatio? prediction->ChZllPt->hData: prediction->ChZllPt->hZJetsToLL,
			                                              fDoDataZllToPhotonRatio? prediction->ChPhotonPt->hData : prediction->ChPhotonPt->hPhotons,
				                                      fZllToGammaFitMin, fZllToGammaFitMax, false, fZllToPhotonRatioRelErr);
			f_lin->SetParameter(0,fZllToPhotonRatio);
		}
  		const int nBins= 3;
		const double Bins[nBins+1] = {150, fZllToGammaFitMin>150?fZllToGammaFitMin:150.0001, fZllToGammaFitMax<800? fZllToGammaFitMax:799.99, 800};
		TH1D* hErrorbar = new TH1D(rationame.Data(), "", nBins, Bins);
		hErrorbar->SetBinContent(2,fZllToPhotonRatio);
		hErrorbar->SetBinError(2,fZllToPhotonRatioRelErr*fZllToPhotonRatio);
		hErrorbar->SetBinContent(1,-10);
		hErrorbar->SetBinError( 1,0);
		hErrorbar->SetFillColor(5);
		hErrorbar->SetFillStyle(3001);
		hErrorbar->Draw("e2same");
		f_lin->Draw("same");
		hPhotonToZllPtRatio->Draw("EX0same");
	}
	
	// GenLevel Zll Pt, no acceptance cuts
	if(fDoGenZllShape){
	prediction->ChGenZllPt = new Channel("GenZllPt", "GenDiLeptPt(0,10,0,1000,true)", cutStreamGenZll.str().c_str(), fTriggerStream.str().c_str(),fSamplesZllPtMConly);
	prediction->ChGenZllPt->fVerbose =prediction->fVerbose;
	prediction->ChGenZllPt->GetShapes("GenZllPt", "GenZll Pt", 8, 0, 800);
	}
	
	// GenLevel Zll Pt, within acceptance
	if(fDoGenAccZllShape){
	prediction->ChGenZllPtAcc = new Channel("GenZllPtAcc", "GenDiLeptPt(20,2.4,71,111,true)", cutStreamGenZllAcc.str().c_str(), fTriggerStream.str().c_str(),fSamplesZllPtMConly);
	prediction->ChGenZllPtAcc->fVerbose =prediction->fVerbose;
	prediction->ChGenZllPtAcc->GetShapes("GenZllPtAcc", "GenZll Pt", 8, 0, 800);
	}
	// Get Acceptance Efficiency
	if(fDoGenZllShape && fDoGenAccZllShape){
		Bool_t binomial =true;
		TH1D* hZllAccEff = prediction->GetRatio(prediction->ChGenZllPtAcc->hZJetsToLL, prediction->ChGenZllPt->hZJetsToLL, 1, binomial);
		DrawHisto(hZllAccEff,hZllAccEff->GetName(),"EX");
//		TString rationame=hZllAccEff->GetName();
//		rationame +="_fittedRatio";
//		TF1 *f_lin = new TF1(rationame,"pol0(0)", fZllAccFitMin ,   fZllAccFitMax);   f_lin->SetLineColor(8);
//		if(fDoFits){ 
//			hZllAccEff->Fit(rationame,"0L","", fZllAccFitMin, fZllAccFitMax);    // set al weights to 1
//			fZllAccEff= f_lin->GetParameter(0);
//		} else{
//			fZllAccEff= prediction->GetLimitedRatio(prediction->ChGenZllPtAcc->hZJetsToLL,prediction->ChGenZllPt->hZJetsToLL,
//				                                fZllAccFitMin, fZllAccFitMax, true, fZllAccEffRelErr);
//			f_lin->SetParameter(0,fZllAccEff);
//		}
//  		const int nBins= 3;
//		const double Bins[nBins+1] = {150, fZllAccFitMin>150?fZllAccFitMin:150.0001, fZllAccFitMax<800? fZllAccFitMax:799.99, 800};
//		TH1D* hErrorbar = new TH1D(rationame.Data(), "", nBins,Bins);
//		hErrorbar->SetBinContent(2,fZllAccEff);
//		hErrorbar->SetBinError(2,fZllAccEffRelErr*fZllAccEff);
//		hErrorbar->SetBinContent(1,-1);
//		hErrorbar->SetBinError(1,0);
//		hErrorbar->SetFillColor(5);
//		hErrorbar->SetFillStyle(3001);
//		hErrorbar->Draw("e2same");
//		f_lin->Draw("same");
//		hZllAccEff->Draw("EX0same");
	}
	
	// GenLevel Zll Pt, within acceptance, OS dilepton recoed
	if(fDoGenAccZllRecoShape){
	prediction->ChGenZllPtAccRecoed = new Channel("GenZllPtAccRecoed", "GenDiLeptPt(20,2.4,71,111,true)", cutStreamGenZllAcc_recoed.str().c_str(), fTriggerStream.str().c_str(),fSamplesZllPtMConly);
	prediction->ChGenZllPtAccRecoed->fVerbose =prediction->fVerbose;
	prediction->ChGenZllPtAccRecoed->GetShapes("GenZllPtAccRecoed", "GenZll Pt", 100, 150, 800);
	}
	if(fDoGenAccZllShape && fDoGenAccZllRecoShape){
		Bool_t binomial =true;
		TH1D* hZllRecoEff = prediction->GetRatio(prediction->ChGenZllPtAccRecoed->hZJetsToLL, prediction->ChGenZllPtAcc->hZJetsToLL, 4, binomial);
		DrawHisto(hZllRecoEff,hZllRecoEff->GetName(),"EX");
		TString rationame=hZllRecoEff->GetName();
		rationame +="_fittedRatio";
		TF1 *f_lin = new TF1(rationame,"pol0(0)", fZllRecoEffFitMin , fZllRecoEffFitMax);   f_lin->SetLineColor(8);
		if(fDoFits){ 
			hZllRecoEff->Fit(rationame,"0L","", fZllRecoEffFitMin, fZllRecoEffFitMax);    // set al weights to 1
			fZllRecoEff= f_lin->GetParameter(0);
		} else{
			fZllRecoEff= prediction->GetLimitedRatio(prediction->ChGenZllPtAccRecoed->hZJetsToLL,prediction->ChGenZllPtAcc->hZJetsToLL,
				                                fZllRecoEffFitMin, fZllRecoEffFitMax, true, fZllRecoEffRelErr);
			f_lin->SetParameter(0,fZllRecoEff);
		}
  		const int nBins= 3;
		const double Bins[nBins+1] = {150, fZllRecoEffFitMin>150?fZllRecoEffFitMin:150.001, fZllRecoEffFitMax<800? fZllRecoEffFitMax:799.99, 800};
		TH1D* hErrorbar = new TH1D(rationame.Data(), "", nBins,Bins);
		hErrorbar->SetBinContent(2,fZllRecoEff);
		hErrorbar->SetBinError(2,fZllRecoEffRelErr*fZllRecoEff);
		hErrorbar->SetBinContent(1,-1);
		hErrorbar->SetBinError(1,0);
		hErrorbar->SetFillColor(5);
		hErrorbar->SetFillStyle(3001);
		hErrorbar->Draw("e2same");
		f_lin->Draw("same");
		hZllRecoEff->Draw("EX0same");
	}
	
	// Photons Hadronic Search MT2
	if(fDoPhotonMT2Shape){
	prediction->ChPhotonsMT2 = new Channel("PhotonsMT2", "photon[0].lv.Pt()", cutStreamPhotonMT2.str().c_str(), fTriggerStream.str().c_str(),fSamplesPhotonPt);
	prediction->ChPhotonsMT2->fVerbose =prediction->fVerbose;
	prediction->ChPhotonsMT2->GetShapes("PhotonsMT2", "MET", 40, 0, 800);
	}

	// Znunu Hadronic Search MT2
	if(fDoZnunuMT2Shape){
	prediction->ChZnunuMT2 = new Channel("ZnunuMT2", "misc.MET", cutStreamZnunuMT2.str().c_str(), fTriggerStream.str().c_str(),fSamplesZnunu);
	prediction->ChZnunuMT2->fVerbose =prediction->fVerbose;
	prediction->ChZnunuMT2->GetShapes("ZnunuMT2", "MET", 40, 0, 800);
	}

	if(fDoZnunuMT2Shape && fDoPhotonMT2Shape){
	TH1D* ZnunuToPhotonMT2Ratio = prediction->GetRatio(prediction->ChZnunuMT2->hZJetsToNuNu, prediction->ChPhotonsMT2->hPhotons, 1);
	DrawHisto(ZnunuToPhotonMT2Ratio,ZnunuToPhotonMT2Ratio->GetName(),"EX0");
	}


	
	// Do Pt spectra comparison plot
	if(fDoPtSpectraComparison && fDoPhotonPtShape && fDoZllPtShape){
		*fLogStream<< "************************* produce pr spectra plot ****************************** " << endl;
		TH1D* hZllMC_cp      = prediction->GetScaledHisto(prediction->ChZllPt->hZJetsToLL, fDoPtSpectraComparisonScaling?(prediction->ChZllPt->hData->Integral())/(prediction->ChZllPt->hZJetsToLL->Integral())    :1, 0, 1);
		TH1D* hZllData_cp    = prediction->GetScaledHisto(prediction->ChZllPt->hData,    1, 0, 1) ;
		TH1D* hPhotonMC_cp   = prediction->GetScaledHisto(prediction->ChPhotonPt->hPhotons,fDoPtSpectraComparisonScaling?(prediction->ChPhotonPt->hData->Integral())/(prediction->ChPhotonPt->hPhotons->Integral()):1, 0, 1);
		TH1D* hPhotonData_cp = prediction->GetScaledHisto(prediction->ChPhotonPt->hData,    1, 0, 1);

		if(fDoPtSpectraComparisonAccCorr){
			TFile *f = new TFile("../RootMacros/ZllAcceptance.root", "OPEN");
			TH1D*  hZllAcc = (TH1D*) f->Get("ZJetsToLL_GenZllPtAcc_ZJetsToLL_GenZllPt_Ratio");
			if (hZllAcc==0) {cout << "WARNING: could not get histo ZJetsToLL_GenZllPtAcc_ZJetsToLL_GenZllPt_Ratio" << endl; exit(1);}
			for(int i=1; i<=hZllMC_cp->GetNbinsX(); ++i){
				if(hZllAcc->GetBinLowEdge(i) != hZllMC_cp->GetBinLowEdge(i)) {cout << "Zll Acc Correction: binnin does not match!!" << endl; exit(1);}
				if(hZllMC_cp  ->GetBinContent(i)<=0) continue;
				double acc_eff   = hZllAcc->GetBinContent(i);
				double orig_mc   = hZllMC_cp     ->GetBinContent(i);
				double orig_data = hZllData_cp   ->GetBinContent(i);
				cout << "bin i " << i << " acc eff " << acc_eff << " orig_mc " << orig_mc << " become " << orig_mc/acc_eff 
				     << " orig_data " << orig_data << " becomes " << orig_data/acc_eff << endl;
				hZllMC_cp   ->SetBinContent(i, orig_mc  /acc_eff);
				hZllData_cp ->SetBinContent(i, orig_data/acc_eff);
			}
			delete f;
		}
		
		hZllMC_cp->SetMarkerStyle(22);
		hZllMC_cp->SetLineColor(kOrange);
		hZllMC_cp->SetMarkerColor(kOrange);
		hZllMC_cp->SetMarkerSize(1.2);
		hZllData_cp->SetMarkerStyle(26);
		hZllData_cp->SetMarkerColor(kBlack);
		hZllData_cp->SetLineColor(kBlack);
		hPhotonMC_cp->SetLineColor(kMagenta+2);
		hPhotonMC_cp->SetMarkerColor(kMagenta+2);
		hPhotonMC_cp->SetMarkerStyle(20);
		hPhotonMC_cp->SetMarkerSize(1.2);
		hPhotonData_cp->SetMarkerStyle(4);
		hPhotonData_cp->SetMarkerColor(kBlack);
		hPhotonData_cp->SetLineColor(kBlack);

		TCanvas *col = new TCanvas("ZllToGammaPtSpectra", "", 0, 0, 500, 500);
		col->SetFillStyle(0);
		col->SetFrameFillStyle(0);
//		col->cd();
		gPad->SetFillStyle(0);
		gStyle->SetPalette(1);
		gPad->SetRightMargin(0.15);

		TLegend* leg2 = new TLegend(0.2, 0.6, 0.5, 0.9);	
		leg2->AddEntry(hPhotonMC_cp,"Gamma+Jets MC","p" );
		leg2->AddEntry(hPhotonData_cp,"Photon Data","p" );
		leg2->AddEntry(hZllMC_cp ,"Zll MC","p" );
		leg2->AddEntry(hZllData_cp,"Zll Data","p" );
		leg2 -> SetFillColor(0);
		leg2 -> SetBorderSize(0);
	//	
		hPhotonMC_cp->SetXTitle("V boson p_{T} (GeV)");
		hPhotonMC_cp->SetYTitle("Events / 7 GeV");
		hPhotonMC_cp->Draw("EX");
		hPhotonData_cp->Draw("EXsame");
		hZllMC_cp->Draw("EXsame");
		hZllData_cp->Draw("EXsame");
		leg2->Draw();
	
		TCanvas *c3 = new TCanvas("ZllToGammaPtSpectraRatio", "", 500, 500);
		TH1D* h_ratio         = (TH1D*) hZllMC_cp      ->Clone("h_ratio");	
		TH1D* h_ratioData     = (TH1D*) hZllData_cp    ->Clone("h_ratioData");	
		TH1D* hPhotonMC_cp2   = (TH1D*) hPhotonMC_cp   ->Clone("hPhotonMC_cp2");
		TH1D* hPhotonData_cp2 = (TH1D*) hPhotonData_cp ->Clone("hPhotonData_cp2");
		h_ratio->Rebin(1); h_ratioData->Rebin(1);
		h_ratio->SetYTitle("#frac{Z(ll)}{#gamma}");
		h_ratio->SetXTitle("boson p_{T} (GeV)");
		h_ratio    ->Divide(hPhotonMC_cp2->Rebin(1));
		h_ratioData->Divide(hPhotonData_cp2->Rebin(1));
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
	}

	if(fDoMT2SpectraCompaison && fDoZnunuMT2Shape && fDoPhotonMT2Shape){
		*fLogStream<<  "************************* MT2 spectra comparison ***************  " << endl;
		TH1D* hZNunuMT2      = prediction->GetScaledHisto(prediction->ChZnunuMT2   ->hZJetsToNuNu, fDoMT2SpectraCompaisonScaling?(prediction->ChPhotonsMT2 ->hData->Integral())/(prediction->ChPhotonsMT2 ->hPhotons->Integral()):1, 0, 1);
		TH1D* hPhotonMT2     = prediction->GetScaledHisto(prediction->ChPhotonsMT2 ->hPhotons,     fDoMT2SpectraCompaisonScaling?(prediction->ChPhotonsMT2 ->hData->Integral())/(prediction->ChPhotonsMT2 ->hPhotons->Integral()):1, 0, 1);
		TH1D* hPhotonDataMT2 = prediction->GetScaledHisto(prediction->ChPhotonsMT2 ->hData, 1, 0, 1);
		
		hZNunuMT2->SetFillStyle(0);
		hZNunuMT2->SetLineWidth(3);

		TCanvas *col = new TCanvas("ZnunuToGammaPtSpectra", "", 0, 0, 500, 500);
		col->SetFillStyle(0);
		col->SetFrameFillStyle(0);
//		col->cd();
		gPad->SetFillStyle(0);
		gStyle->SetPalette(1);
		gPad->SetRightMargin(0.15);
		gPad->SetLogy(1);

		TLegend* leg2 = new TLegend(0.2, 0.6, 0.5, 0.9);	
		leg2->AddEntry(hPhotonMT2,"#gamma+jets MC","f" );
		leg2->AddEntry(hZNunuMT2,"Z(#nu#nu)+jets MC","l" );
		leg2->AddEntry(hPhotonDataMT2,"data","p" );
		leg2 -> SetFillColor(0);
		leg2 -> SetBorderSize(0);
	//	
		hPhotonMT2->SetXTitle("min #Delta #phi(jets,MET)");
		hPhotonMT2->SetYTitle("Events");
		hPhotonMT2->Draw("hist");
		hZNunuMT2->Draw("histsame");
		hPhotonDataMT2->Draw("EXsame");
		leg2->Draw();
//		gPad->RedrawAxis();
		// cout integral above 250 in MT2 and ratio
//		double sumPhotons=0;
//		double sumZnunu  =0;
//		for(int i=1;  i<=hPhotonMT2->GetNbinsX(); ++i){
//			if(hPhotonMT2->GetBinLowEdge(i)>=250){
//				sumPhotons+=hPhotonMT2->GetBinContent(i);
//				sumZnunu  +=hZNunuMT2 ->GetBinContent(i);
//			}
//		}
//		*fLogStream<< ">>> hPhotonMT2: integral above 250: " << sumPhotons << endl;
//		*fLogStream<< ">>> hZNunuMT2 : integral above 250: " << sumZnunu   << endl;
//		*fLogStream<< ">>> -> Ratio  : "                     << sumZnunu/sumPhotons << endl;
		
	}

	if(fDoZnunuGammaGenPtSpectraComparison){
		*fLogStream<< "*************************ZnunuGammaGenPtSpectraComparison**********" << endl;
		// gen photons
		prediction->ChGenPhotonPt = new Channel("GenPhotonPt", "GenPhoton[0].Pt()", cutStreamGenPhoton.str().c_str(), fTriggerStream.str().c_str(),fSamplesPhotonPtMConly);
		prediction->ChGenPhotonPt->fVerbose =prediction->fVerbose;
		prediction->ChGenPhotonPt->GetShapes("GenPhotonPt", "Gen-level #gamma Pt", 50, 150, 800);
		
		// Gen Znunu
		prediction->ChGenZnunuPt = new Channel("GenZnunuPt", "GenZ[0].Pt()", cutStreamGenZnunu.str().c_str(), fTriggerStream.str().c_str(),fSamplesZnunu);
		prediction->ChGenZnunuPt->fVerbose =prediction->fVerbose;
		prediction->ChGenZnunuPt->GetShapes("GenZnunuPt", "Gen-level Z(#nu#nu) Pt", 50, 150, 800);
		
		prediction->DrawHistos(prediction->ChGenPhotonPt->hPhotons, prediction->ChGenZnunuPt->hZJetsToNuNu,
				       "EX0"                              , "EX0same",
				       kViolet+3                          , kBlack,
				       20                                 ,4,
				       "#gamma Pt"                        ,"Z(#nu#nu)");

		TH1D* hGenPhotonToZnunuRatio = prediction->GetRatio(prediction->ChGenZnunuPt->hZJetsToNuNu, prediction->ChGenPhotonPt->hPhotons, 1, false);
		hGenPhotonToZnunuRatio->SetMarkerColor(kBlack);
		hGenPhotonToZnunuRatio->SetLineColor(kBlack);
		hGenPhotonToZnunuRatio->SetMarkerStyle(4);
		DrawHisto(hGenPhotonToZnunuRatio, "GenPhotonToZnunuRatio","EX0", prediction->ChGenPhotonPt);
	}

	if(fDoPhotonRecoEff){
		prediction->ChGenPhotonPtForRecoEff = new Channel("GenPhotonPtForRecoEff", "GenPhoton[0].Pt()", cutStream_GenPhotonforRecoEff.str().c_str(), fTriggerStream.str().c_str(),fSamplesPhotonPtMConly);
		prediction->ChGenPhotonPtForRecoEff->fVerbose =prediction->fVerbose;
		prediction->ChGenPhotonPtForRecoEff->GetShapes("GenPhotonPtForRecoEff", "Gen-level #gamma Pt", 50, 150, 800);
		
		prediction->ChPhotonPtForRecoEff = new Channel("PhotonPtForRecoEff", "GenPhoton[0].Pt()", cutStream_PhotonPtforRecoEff.str().c_str(), fTriggerStream.str().c_str(),fSamplesPhotonPtMConly);
		prediction->ChPhotonPtForRecoEff->fVerbose =prediction->fVerbose;
		prediction->ChPhotonPtForRecoEff->GetShapes("PhotonPtForRecoEff", "Gen-level #gamma Pt", 50, 150, 800);
		
		prediction->DrawHistos(prediction->ChGenPhotonPtForRecoEff->hPhotons, prediction->ChPhotonPtForRecoEff->hPhotons,
				       "EX0"                              , "EX0same",
				       kViolet+3                          , kBlack,
				       20                                 ,4,
				       "#all"                             ,"recoed");
		TH1D* hPhotonRecoEff = prediction->GetRatio(prediction->ChPhotonPtForRecoEff->hPhotons, prediction->ChGenPhotonPtForRecoEff->hPhotons, 1, false);
		hPhotonRecoEff->SetMarkerColor(kBlack);
		hPhotonRecoEff->SetLineColor(kBlack);
		hPhotonRecoEff->SetMarkerStyle(4);
		DrawHisto(hPhotonRecoEff, "PhotonRecoEff","EX0", prediction->ChPhotonPtForRecoEff);
	}

	// Prediction
	if(fDoPrediction){
		*fLogStream<< "************************* Prediction ****************************** " << endl;

		TH1D* MCZnunu = prediction->GetScaledHisto(prediction->ChZnunuMT2->hZJetsToNuNu,fLumiCorr,0);  // scale to lumi 4400
		double MCZnunuEst    = MCZnunu->GetBinContent(1);
		double MCZnunuEstErr = MCZnunu->GetBinError(1);
		delete MCZnunu;
		if(fDoData){ 
		double nGamma    = prediction->ChPhotonsMT2->hData->Integral();
		double nGammaErr = sqrt(nGamma);
		*fLogStream << "********** Data Prediction ***************************************************** " << endl;
		MakePrediction(nGamma, nGammaErr, fZllToPhotonRatio, fZllToPhotonRatioRelErr, fZllAccEff, fZllAccEffRelErr, fZllRecoEff, fZllRecoEffRelErr, MCZnunuEst, MCZnunuEstErr);
		}
		TH1D* hDummy =   prediction->GetScaledHisto(prediction->ChPhotonsMT2->hPhotons, 1, 0);
		float nGammaMC   =hDummy->GetBinContent(1);
		float nGammaMCErr=hDummy->GetBinError(1);
		*fLogStream << "********** MC Prediction ***************************************************** " << endl;
		MakePrediction(nGammaMC, nGammaMCErr, fZllToPhotonRatio, fZllToPhotonRatioRelErr, fZllAccEff, fZllAccEffRelErr, fZllRecoEff, fZllRecoEffRelErr, MCZnunuEst, MCZnunuEstErr);
	}

	// write -----------
	if(fWriteToFile){
		TString logname =fOutDir + ".log"; 
		ofstream f_log (logname.Data(), ios::app);
		f_log << fLogStream->str();
	} else{
		cout << fLogStream->str();
	}

//	fOutFile->Close();

}

// ****************************************** Definition of classes and methods *********************************************************
//
//
void MakePrediction(float nGamma, float nGammaErr, float ZllToPhotonRatio, float ZllToPhotonRatioRelErr, 
		   float ZllAccEff, float ZllAccEffRelErr, float ZllRecoEff, float ZllRecoEffRelErr,
		   float MCZnunuEst, float MCZnunuEstErr){
	*fLogStream << "Prediction with:   " << endl;
	*fLogStream << "  nGamma           " << nGamma            << " pm " << nGammaErr                                 << endl;
	*fLogStream << "  ZllToPhotonRatio " << ZllToPhotonRatio  << " pm " << ZllToPhotonRatioRelErr*ZllToPhotonRatio   << endl;
	*fLogStream << "  ZllAccEff        " << ZllAccEff         << " pm " << ZllAccEffRelErr*ZllAccEff                 << endl;
	*fLogStream << "  ZllRecoEff       " << ZllRecoEff        << " pm " << ZllRecoEff*ZllRecoEffRelErr               << endl;

	double pred         = nGamma*ZllToPhotonRatio/(ZllAccEff*ZllRecoEff)*fBRZnunu/fBRZll*fLumiCorr;
	double prederrstat  = fLumiCorr*sqrt(pow(fBRZnunu/fBRZll,2) *pow(1/(ZllAccEff*ZllRecoEff),2)
		                                                                       *( (pow(nGammaErr,2)       *pow(ZllToPhotonRatio,2))
			                                                                )
			    );
	double prederrsys   = fLumiCorr*sqrt(pow(fBRZnunu/fBRZll,2) *pow(1/(ZllAccEff*ZllRecoEff),2)
		                                                                       *(  (pow(nGamma,2)          *pow(ZllToPhotonRatio*ZllToPhotonRatioRelErr,2))
			                                                                 + (pow(ZllAccEffRelErr,2) *pow(nGamma,2)*pow(ZllToPhotonRatio,2))
			                                                                 + (pow(ZllRecoEffRelErr,2)*pow(nGamma,2)*pow(ZllToPhotonRatio,2))
			                                                                )
			    );
	*fLogStream << "\n" 
		    << "  NZnunu Pred:     " << pred       << " pm (stat)  " << prederrstat << " pm (sys) " << prederrsys << "\n"
		    << "  MC Znunu   :     " << MCZnunuEst << " pm         " << MCZnunuEstErr      << "\n"
	            << "*******************************************************************" << endl;
} 
class Channel {
public:
	Channel(TString name, TString variable, TString cuts, TString trigger, TString samples);
	~Channel();
	TH1D* hQCD; 
	TH1D* hData; 
	TH1D* hPhotons;
	TH1D* hWJets;
	TH1D* hZJetsToLL;
	TH1D* hZJetsToNuNu;
	TH1D* hTop;
	TH1D* hSignal;
	TH1D* hOther;
	TH1D* hTotalSM;
	TString fName;
	TString fCuts;
	TString fTrigger;
	TString fSamples;
	TString fVariable;
	vector<TH1D*> h_shapes;
	void GetShapes(TString SelectionName, TString xtitle, const int nbins, float binmin, float binmax);
	void GetShapes(TString SelectionName, TString xtitle, const int nbins, const double *bins );
	int fVerbose;
	bool fGotShapes;

};

Channel::Channel(TString name, TString variable, TString cuts, TString trigger, TString samples){
	fName=name;
	fCuts=cuts;
	fTrigger=trigger;
	fSamples=samples;
	fVariable=variable;
	hQCD=0;
	hData=0;
	hPhotons=0;
	hWJets=0;
	hZJetsToLL=0;
	hZJetsToNuNu=0;
	hTop=0;
	hSignal=0;
	hOther=0;
	hTotalSM=0;
	fVerbose=4;
	fGotShapes=false;
}
Channel::~Channel(){};

void Channel::GetShapes(TString SelectionName, TString xtitle, const int nbins, float binmin, float binmax){ 
	double bins[nbins];
	bins[0] = binmin;
	for(int i=1; i<=nbins; i++) bins[i] = binmin+i*(binmax-binmin)/nbins;
	GetShapes(SelectionName, xtitle, nbins, bins);
}

void Channel::GetShapes(TString SelectionName, TString xtitle, const int nbins, const double *bins){ 
	MT2Shapes *tA;
	if(fWriteToFile) tA = new MT2Shapes(fOutDir, fRootFile, fLogStream);
	else             tA = new MT2Shapes(fOutDir, fRootFile);
	tA->setVerbose(fVerbose);
	tA->init(fSamples);
	tA->SetPrintSummary(true);
	tA->SetDraw(false);
	tA->SetWrite(false);
  
//                    variable,      cuts,  njet, nlep,  selection_name,      HLT,   xtitle   nbins  bins   
        tA->GetShapes(fVariable,  fCuts,    -1,  -10  , SelectionName,    fTrigger , xtitle , nbins, bins);

	// retrieve shapes
	for(int i=0; i<tA->GetNShapes(); ++i){
		TString name =tA->fh_shapes[i]->GetName();
		if      (name.Contains("QCD_"))        {hQCD         = (TH1D*) tA->fh_shapes[i]->Clone(tA->fh_shapes[i]->GetName()); hQCD->SetDirectory(fDir);}
		else if (name.Contains("PhotonsJets_")){hPhotons     = (TH1D*) tA->fh_shapes[i]->Clone(tA->fh_shapes[i]->GetName()); hPhotons->SetDirectory(fDir);}
		else if (name.Contains("Data_"))       {hData        = (TH1D*) tA->fh_shapes[i]->Clone(tA->fh_shapes[i]->GetName()); hData->SetDirectory(fDir);}
		else if (name.Contains("WJets_"))      {hWJets       = (TH1D*) tA->fh_shapes[i]->Clone(tA->fh_shapes[i]->GetName()); hWJets->SetDirectory(fDir);}
		else if (name.Contains("ZJetsToLL_"))  {hZJetsToLL   = (TH1D*) tA->fh_shapes[i]->Clone(tA->fh_shapes[i]->GetName()); hZJetsToLL->SetDirectory(fDir);}
		else if (name.Contains("ZJetsToNuNu_")){hZJetsToNuNu = (TH1D*) tA->fh_shapes[i]->Clone(tA->fh_shapes[i]->GetName()); hZJetsToNuNu->SetDirectory(fDir);}
		else if (name.Contains("Top_"))        {hTop         = (TH1D*) tA->fh_shapes[i]->Clone(tA->fh_shapes[i]->GetName()); hTop->SetDirectory(fDir);}
		else if (name.Contains("Signal_"))     {hSignal      = (TH1D*) tA->fh_shapes[i]->Clone(tA->fh_shapes[i]->GetName()); hSignal->SetDirectory(fDir);}
		else if (name.Contains("Other_"))      {hOther       = (TH1D*) tA->fh_shapes[i]->Clone(tA->fh_shapes[i]->GetName()); hOther->SetDirectory(fDir);}
	}
	delete tA;
	fGotShapes=true;

	fDir->cd();

	// fix colors
	if(hQCD!=0)        {hQCD    ->SetLineColor(kYellow+1);    hQCD    ->SetFillColor(kYellow+1);      hQCD    ->SetFillStyle(3001);}
	if(hPhotons!=0)    {hPhotons->SetLineColor(kViolet-3);    hPhotons->SetFillColor(kViolet-3);      hPhotons->SetFillStyle(3001);}
	if(hOther!=0)      {hOther  ->SetLineColor(kCyan+2);      hOther  ->SetFillColor(kCyan+2);        hOther  ->SetFillStyle(3001);}
	if(hZJetsToNuNu!=0){hZJetsToNuNu->SetLineColor(kGreen+1); hZJetsToNuNu->SetFillColor(kGreen+1);   hZJetsToNuNu  ->SetFillStyle(3001);}
	if(hSignal!=0)     {hSignal ->SetLineColor(kBlack);       hSignal  ->SetFillColor(kBlack);        hSignal  ->SetFillStyle(3001);}
	if(hZJetsToLL!=0)  {hZJetsToLL->SetLineColor(kOrange);    hZJetsToLL->SetFillColor(kOrange);      hZJetsToLL->SetFillStyle(3001);}
	if(hTop!=0)        {hTop    ->SetLineColor(600);          hTop     ->SetFillColor(600);           hTop->SetFillStyle(3001);}

	if(fSaveResults){
		if(hQCD!=0)         {hQCD->Write();}
		if(hPhotons!=0)     {hPhotons->Write();}
		if(hOther!=0)       {hOther->Write();}
		if(hData!=0)        {hData->Write();}
		if(hTop!=0)         {hTop->Write();}
		if(hZJetsToLL!=0)   {hZJetsToLL->Write();}
		if(hZJetsToNuNu!=0) {hZJetsToNuNu->Write();}
		if(hSignal!=0)      {hSignal->Write();}
	}
	if(fDraw){
		if(hQCD!=0)        DrawHisto(hQCD,           hQCD->GetName(),            "hist", this);
		if(hPhotons!=0)    DrawHisto(hPhotons,       hPhotons->GetName(),        "hist", this);
		if(hOther!=0)      DrawHisto(hOther,         hOther->GetName(),          "hist", this);
		if(hData!=0)       DrawHisto(hData,          hData->GetName(),           "EXO",  this);
		if(hTop!=0)        DrawHisto(hTop,           hTop ->GetName(),           "hist", this);
		if(hZJetsToNuNu!=0)DrawHisto(hZJetsToNuNu,   hZJetsToNuNu->GetName(),    "hist", this);
		if(hSignal!=0)     DrawHisto(hSignal,        hSignal->GetName(),         "hist", this);
		if(hZJetsToLL!=0)  DrawHisto(hZJetsToLL,     hZJetsToLL->GetName(),      "hist", this);
	}
}

class Prediction {
public:
	Prediction();
	~Prediction();
	
	TH1D* GetScaledHisto(TH1D* histo, float scale_fact, float scale_fact_err);
	TH1D* GetScaledHisto(TH1D* histo, float scale_fact, float scale_fact_err, int ngroup);
	TH1D* GetRatio(TH1D* h1, TH1D* h2, int mergenbins=1,  bool binomial=false);
	void  DrawHistos(TH1D* h1orig, TH1D* h2orig, Option_t *drawopt1, Option_t *drawopt2, Int_t color1, Int_t color2, Int_t MarkerStyle1, Int_t MarkerStyle2, TString leg1, TString leg2);
	Channel* ChPhotonPt;
	Channel* ChPhotonPtForRecoEff;
	Channel* ChGenPhotonPtForRecoEff;
	Channel* ChGenPhotonPt;
	Channel* ChZllPt;
	Channel* ChGenZllPt;
	Channel* ChGenZnunuPt;
	Channel* ChGenZllPtAcc;
	Channel* ChGenZllPtAccRecoed;
	Channel* ChPhotonsMT2;
	Channel* ChZnunuMT2;

	bool fSave;
	int fVerbose;
};

Prediction::Prediction(){
	fVerbose =0;
	fSave = false;
}
Prediction::~Prediction(){
}

TH1D* Prediction::GetRatio(TH1D* h1, TH1D* h2, int mergenbins, bool binomial){
	// compute Zll/Photon ratio
	*fLogStream << "------------------------------------------------------------------------" << endl;
	*fLogStream << "Compute ratio between " << h1->GetName() << " and " << h2->GetName()      << endl;           

	
	if(h1==0 || h2==0 ){
		cout << "GetMCZnunuToPhotonRatio: received 0 pointer!" << endl;
		exit(-1);
	}

	// GetScaled histos with only one bin and probagated errors: stat error and error on scale factor
	TH1D *currh1   = GetScaledHisto(h1, 1, 0, mergenbins); 
	TH1D *currh2   = GetScaledHisto(h2, 1, 0, mergenbins); 

	TString name = h1->GetName();
	name         += "_";
	name         += h2->GetName();
	name         += "_Ratio";
	TH1D *ratio = (TH1D*) currh1->Clone(name);
	ratio->Divide(currh1, currh2, 1,1, binomial?"B":""); // binomial errors
	
	ratio->SetLineColor(kBlack);
	ratio->SetMarkerColor(kBlack);
	ratio->SetMarkerStyle(4);

	if(fWriteToFile){
		ratio ->Write();
	}
	delete currh1;
	delete currh2;
	return ratio;
}

TH1D* Prediction::GetScaledHisto(TH1D* histo, float scale_fact, float scale_fact_err){
	// takes histo, scale factor and uncertainty on scale factor
	// and returns scaled and rebinned histo with 1 bin with uncertainty on scale factor propagated.
	TH1D *h        = (TH1D*) histo->Clone(histo->GetName());
	h->Rebin(h->GetNbinsX());
	h->SetBinError(1, sqrt(h->GetBinError(1)*  h->GetBinError(1)   *scale_fact    *scale_fact + 
			  h->GetBinContent(1)*h->GetBinContent(1) *scale_fact_err*scale_fact_err));
	h->SetBinContent(1, h->GetBinContent(1)*scale_fact);
	return h;
}

TH1D* Prediction::GetScaledHisto(TH1D* histo, float scale_fact, float scale_fact_err, int ngroup){
	// takes histo, scale factor and uncertainty on scale factor
	// and returns scaled and rebinned histo with ngroup bins merged into 1 with uncertainty on scale factor propagated.
	if(ngroup>=histo->GetNbinsX()) ngroup=histo->GetNbinsX();
	TH1D *h        = (TH1D*) histo->Clone(histo->GetName());
	h->Rebin(ngroup);
	for(int i=1; i<=h->GetNbinsX(); ++i){
		h->SetBinError(i, sqrt(h->GetBinError(i)*  h->GetBinError(i)   *scale_fact    *scale_fact + 
				  h->GetBinContent(i)*h->GetBinContent(i) *scale_fact_err*scale_fact_err));
		h->SetBinContent(i, h->GetBinContent(i)*scale_fact);
	}
	return h;
}

double Prediction::GetLimitedRatio(TH1D* h1, TH1D* h2, float min, float max, bool binomial, float &relerror){
	if(h1->GetNbinsX() != h2->GetNbinsX()){cout << "GetIntegral ERROR! " << endl; exit(-1);}
	TH1D* h1cp = h1->Clone("h1cp");
	TH1D* h2cp = h2->Clone("h2cp");
	for(int i=1; i<=h1cp->GetNbinsX(); ++i){
		if(h1cp->GetBinLowEdge(i) <min || h1cp->GetBinLowEdge(i)+h1cp->GetBinWidth(i) >max){
			h1cp->SetBinContent(i, 0);
			h1cp->SetBinError(i, 0);
			h2cp->SetBinContent(i, 0);
			h2cp->SetBinError(i, 0);
		}
	}
	TH1D* h1cp_2 = GetScaledHisto(h1cp, 1, 0, h1cp->GetNbinsX());
	TH1D* h2cp_2 = GetScaledHisto(h2cp, 1, 0, h1cp->GetNbinsX());

	TH1D* r      = h1cp_2->Clone("r");
	r->Divide(h1cp_2, h2cp_2, 1, 1, binomial? "B":"");

	*fLogStream << "Getting Limited Ratio for " << h1->GetName()  << " " << h2->GetName() << endl;
	*fLogStream << ">>> ratio = " << r->GetBinContent(1)  << " pm " << r->GetBinError(1) << endl;
	if(fDoRatioStatErr){
		relerror = r->GetBinError(1)/r->GetBinContent(1);
	}
	double  ratio = r->GetBinContent(1);
	delete h1cp;
	delete h2cp;
	delete h1cp_2;
	delete h2cp_2;
	delete r;
	return ratio;
}

void DrawHisto(TH1* h_orig, TString name,  Option_t *drawopt, Channel* channel){
	TH1D* h = (TH1D*)h_orig->Clone(name);
	TString canvname = "canv_"+name;
	TCanvas *col = new TCanvas(canvname, "", 0, 0, 500, 500);
	col->SetFillStyle(0);
	col->SetFrameFillStyle(0);
	gPad->SetFillStyle(0);
	gStyle->SetPalette(1);
	gPad->SetRightMargin(0.15);
	h->DrawCopy(drawopt);
//	gPad->RedrawAxis();
	if(fSaveResults){
		col->Write();
	}
}

void DrawHisto(TH1* h_orig, TString name,  Option_t *drawopt){
	TH1D* h = (TH1D*)h_orig->Clone(name);
	TString canvname = "canv_"+name;
	TCanvas *col = new TCanvas(canvname, "", 0, 0, 500, 500);
	col->SetFillStyle(0);
	col->SetFrameFillStyle(0);
	gPad->SetFillStyle(0);
	gStyle->SetPalette(1);
	gPad->SetRightMargin(0.15);
	h->DrawCopy(drawopt);
}

void Prediction::DrawHistos(TH1D* h1orig, TH1D* h2orig, Option_t *drawopt1, Option_t *drawopt2, Int_t color1, Int_t color2, Int_t MarkerStyle1, Int_t MarkerStyle2, TString leg1, TString leg2){
	TH1D* h1 = (TH1D*)h1orig->Clone("h1");
	TH1D* h2 = (TH1D*)h2orig->Clone("h2");
	TString canvname  = "canv_";
	canvname         +=h1orig->GetName();
	canvname         +="_";
	canvname         +=h2orig->GetName();
	TCanvas *col = new TCanvas(canvname, "", 0, 0, 500, 500);
	col->SetFillStyle(0);
	col->SetFrameFillStyle(0);
	gPad->SetFillStyle(0);
	gStyle->SetPalette(1);
	gPad->SetRightMargin(0.15);
	h1->SetTitle("");
	h2->SetTitle("");
	if(MarkerStyle1!=0) {
		h1->SetMarkerStyle(MarkerStyle1); 
		h1->SetMarkerColor(color1); 
		h1->SetLineColor(color1); 
		h1->SetFillStyle(0);
	}
	if(MarkerStyle2!=0) {
		h2->SetMarkerStyle(MarkerStyle2); 
		h2->SetMarkerColor(color2); 
		h2->SetLineColor(color2); 
		h2->SetFillStyle(0);
	}
	h1->DrawCopy(drawopt1);
	h2->DrawCopy(drawopt2);
	if(leg1!=""){
		TLegend* leg = new TLegend(0.2, 0.6, 0.5, 0.9);	
		leg->AddEntry(h1,leg1,MarkerStyle1==0?"l":"p");
		leg->AddEntry(h2,leg2,MarkerStyle2==0?"l":"p");
		leg->SetFillColor(0);
		leg->SetBorderSize(0);
		leg->Draw();
	}
	if(fSave){
		col->Write();
	}
}
void Close(){
	fFile->Close();
}
