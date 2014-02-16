#include "ZInvEstFromW.hh"

#include "TLorentzVector.h"
#include "TGraphAsymmErrors.h"
#include "helper/Utilities.hh"
//#include "helper/Monitor.hh"

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
#include "TDirectoryFile.h"
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
//#include <vector>
//#include <stdlib.h>
//#include <stdio.h>
//#include <map>
#include <time.h> // access to date/time

#include "helper/Hemisphere.hh"


//utils
#ifndef Utilities_HH
#include "Utilities.hh"
#endif

using namespace std;

//____________________________________________________________________________
ZInvEstFromW::ZInvEstFromW(){
  // Default constructor, no samples are set
}

//____________________________________________________________________________

ZInvEstFromW::ZInvEstFromW(TString outputdir/*, TString outputfile*/){
  // Explicit constructor with output directory and output file
  setOutputDir(outputdir);
  //setOutputFile(outputfile);
}

//____________________________________________________________________________
ZInvEstFromW::~ZInvEstFromW(){
  //fOutputFile->Close();
  //delete fOutputFile;
}

//____________________________________________________________________________

void ZInvEstFromW::init(TString filename){
  if(fVerbose > 0) cout << "------------------------------------" << endl;
  if(fVerbose > 0) cout << "Initializing MassPlotter ... " << endl;
  //Util::SetStyle();
  loadSamples(filename);
 }


float* ZInvEstFromW::Analysis( TString MT2_REGIME,  bool IS_MC, TString outputfile, int MIN_NJETS, int NELE, int NMU, float MT2CUT_LOW, float MT2CUT_HIGH,  TString trigger, TString filter,  bool BENRICH, bool CORRECT_EFF, float * INPUT){
  setOutputFile(outputfile);

  //to be set up here:
  bool CORRECT_MET = true; // correct quantities for Leptons
  bool TP_EFF_MC = false; //use MC TP eff for leptons
  float ELE_PT = 20;
  float MU_PT = 10;

  float JET_PT = 40; //for jet counting
  float BJET_PT = 20; //for jet counting


  //set up by arguments
  //bool STD_SEL = true;  //do standard high MT sel, and no Zinv/lost leptons est
  string ELEorMU = "MU"; //Use MU or ELE

  //Cuts
  map<TString, TString> cutLabels;
  vector <TString> orderedCutLabels;
  //cutLabels["leptSel"] = "Lepton Selection"; orderedCutLabels.push_back("leptSel");
  cutLabels["preSel"] = "All events(jets >= X)"; orderedCutLabels.push_back("preSel");
  cutLabels["dphiSel"] = "Minimum DPhi(MET,jet) > 0.3"; orderedCutLabels.push_back("dphiSel");
  cutLabels["hbheSel"] = "HBHE noise veto"; orderedCutLabels.push_back("hbheSel");
  cutLabels["metSel"] = "MET > 30 GeV"; orderedCutLabels.push_back("metSel");
  cutLabels["vsptSel"] = "VectorSumPt < 70"; orderedCutLabels.push_back("vsptSel");
  cutLabels["jidSel"] = "jets > 50GeV failing PFID event veto"; orderedCutLabels.push_back("jidSel");
  cutLabels["leptVetoSel"] = "Lepton Veto"; orderedCutLabels.push_back("leptVetoSel");
  cutLabels["leptSel"] = "Lepton Selection"; orderedCutLabels.push_back("leptSel");
  cutLabels["mtSel"] = "MT < 100"; orderedCutLabels.push_back("mtSel");
  cutLabels["bSel"] = "BJets"; orderedCutLabels.push_back("bSel");

  if(MT2CUT_LOW == -1){
    cutLabels["mt100Sel"] = "MT2  100 - 150 GeV"; orderedCutLabels.push_back("mt100Sel");
    cutLabels["mt150Sel"] = "MT2  150 - 200 GeV"; orderedCutLabels.push_back("mt150Sel");
    cutLabels["mt200Sel"] = "MT2  200 - 275 GeV";orderedCutLabels.push_back("mt200Sel");
    cutLabels["mt275Sel"] = "MT2  275 - 375 GeV";orderedCutLabels.push_back("mt275Sel");
    cutLabels["mt375Sel"] = "MT2  375 - 500 GeV";orderedCutLabels.push_back("mt375Sel");
    cutLabels["mt500Sel"] = "MT2 > 500 GeV";orderedCutLabels.push_back("mt500Sel");
  }
  
  if(NELE==1 && NMU==0){
    ELEorMU= "ELE";
  }
  else if(NELE==0 && NMU==1){
    ELEorMU = "MU";
  }

  //Efficiencies
  TString EFF_DIR = "/shome/leo/Analysis/Efficiencies/";
  TString MuEff_filename = EFF_DIR+"effMuon.root";
  if(TP_EFF_MC)  MuEff_filename = EFF_DIR+"effMuonMC.root";
  if(ELEorMU == "ELE"){
    MuEff_filename = EFF_DIR+"effElectron.root";
    if(TP_EFF_MC)  MuEff_filename = EFF_DIR+"effElectronMC.root";
  }
  TString BEff_filename = EFF_DIR+"Btag_performance_2011.root";

  if(BENRICH) cout << "---+++ 1b, 1l" << endl;
  else  cout << "---+++ 0b, 1l" << endl;
  cout << "-------------------- " << endl;
  cout << "   * MT2 Regime: " << MT2_REGIME << endl;
  cout << "   * Lept-corrected quantities:  " << CORRECT_MET << endl;
  cout << "   * BEnriched: " << BENRICH << endl;
  cout << "   * Jets: " << MIN_NJETS << endl;
  cout << "   * Nele, NMu: " << NELE <<","<< NMU << endl;
  cout << "-------------------- " << endl;


  //Get efficiencies from TP
  TFile *MuEff_file = TFile::Open(MuEff_filename);
  TH2F * h2_effMu_ptdr=NULL,  * h2_effMu_pteta=NULL;

  if(ELEorMU == "MU"){
    if(TP_EFF_MC){
      h2_effMu_ptdr = (TH2F*) MuEff_file->Get("histo_muon_var1pt_var2dr_mc");
      h2_effMu_pteta = (TH2F*) MuEff_file->Get("histo_muon_var1pt_var2eta_mc");
    }
    else{
      h2_effMu_ptdr = (TH2F*) MuEff_file->Get("histo_muon_var1pt_var2dr_data");
      h2_effMu_pteta = (TH2F*) MuEff_file->Get("histo_muon_var1pt_var2eta_data");
    }
  }
  else if(ELEorMU == "ELE"){
    if(TP_EFF_MC){
      h2_effMu_ptdr = (TH2F*) MuEff_file->Get("histo_electron_var1pt_var2dr_mc");
      h2_effMu_pteta = (TH2F*) MuEff_file->Get("histo_electron_var1pt_var2eta_mc");
    }
    else{
      h2_effMu_ptdr = (TH2F*) MuEff_file->Get("histo_electron_var1pt_var2dr_data");
      h2_effMu_pteta = (TH2F*) MuEff_file->Get("histo_electron_var1pt_var2eta_data");
    }
  }
  
  cout << "   * Using " <<  MuEff_filename << " and " <<  h2_effMu_ptdr->GetName() << " " <<h2_effMu_pteta->GetName() << endl;

  ///// BTag eff
  TFile * BEff_file = TFile::Open(BEff_filename);
  TH2F * h2_mistag_SF,  *h2_mistag_SFerr;
  TH2F * h2_btag_SF,  *h2_btag_SFerr;

  h2_mistag_SF = (TH2F*) BEff_file->Get("demo2/MISTAGSSVHEM_BTAGLEFFCORR");
  h2_mistag_SFerr = (TH2F*)  BEff_file->Get("demo2/MISTAGSSVHEM_BTAGLERRCORR");
  h2_btag_SF = (TH2F*) BEff_file->Get("demo2/BTAGSSVHPT_BTAGLEFFCORR");
  h2_btag_SFerr = (TH2F*) BEff_file->Get("demo2/BTAGSSVHPT_BTAGLERRCORR");

  //The Output File
  TFile *fout = TFile::Open(fOutputFile,"RECREATE");
  
  // Some info
  cout << "   * Cuts for Flow (triggers not included here): <verbatim>" << filter << "</verbatim>" << endl;

  //variables
  map<TString, float> varFloat;

  //loop on samples
  vector<TString> samplesName;
  map<TString, map<TString,TH1F*> > counterH;

  //Needed for W subsampling
  //counter Histos
  createCounterHistos(&counterH, cutLabels, "Other"); samplesName.push_back("Other");
  if(ELEorMU=="ELE"){ createCounterHistos(&counterH, cutLabels, "WtolnuEle"); samplesName.push_back("WtolnuEle");}
  createCounterHistos(&counterH, cutLabels, "WtolnuMu"); samplesName.push_back("WtolnuMu");
  createCounterHistos(&counterH, cutLabels, "WtolnuTau"); samplesName.push_back("WtolnuTau");
  createCounterHistos(&counterH, cutLabels, "MC");
 
  for(size_t i = 0; i < fSamples.size(); ++i){
    bool isThere=false;
    TString tsname = fSamples[i].sname;
    if (tsname=="Wtolnu") continue;
    else if(tsname=="HT-Data") {tsname="data";fSamples[i].sname="data";}
    else if( fSamples[i].name == "DYToNuNu_50_100" ||
	     fSamples[i].name == "DYToNuNu_100_200" ||
	     fSamples[i].name == "DYToNuNu_200_infty" ) {tsname="Zinv";fSamples[i].sname="Zinv";}
    else if( fSamples[i].sname == "Photons" ||  fSamples[i].sname == "DY"  ) { tsname="Other"; fSamples[i].sname="Other"; }

    for(unsigned int k=0; k< samplesName.size(); k++){
      if( samplesName[k] == tsname) {isThere = true; break;}
    }
    if(!isThere)  { 
      samplesName.push_back(tsname);
      varFloat["Lept_eff_TP_err_"+tsname] = 0;
      varFloat["B_eff_err_"+tsname] = 0;
    }
  }
  TString lastS = samplesName.back();
  samplesName.erase(samplesName.begin());//delete Other, to put at the end
  samplesName.pop_back();
  samplesName.push_back("Other");
  samplesName.push_back("MC");
  samplesName.push_back(lastS);
  varFloat["Lept_eff_TP_err_MC"] = 0;
  varFloat["B_eff_err_MC"] = 0;

  //Creating distr. histos
  map<TString, TH1F*> histos;
  createHistos( &histos, samplesName );
  
  
  for(size_t i = 0; i < fSamples.size(); ++i){
    int NSelEvents = 0;

    TString tsname = fSamples[i].sname;
    TString type = fSamples[i].type;

    if(tsname=="HT-Data") tsname="data";

    //Getting NEntries, PUWeight
    float avg_pu_w = 1;
    if(fSamples[i].type!="data"){
      TH1F *h_PUWeights = (TH1F*) fSamples[i].file->Get("h_PUWeights");
      avg_pu_w = h_PUWeights->GetMean();
      fSamples[i].nevents = h_PUWeights->GetEntries();
      delete h_PUWeights;
    }
    Double_t sample_weight = fSamples[i].xsection * fSamples[i].kfact * fSamples[i].lumi / (fSamples[i].nevents* avg_pu_w);
    
    if(fVerbose>1) cout << "Analysis: looping over " << fSamples[i].name << endl;
    if(fVerbose>2) cout << "\t\tOriginal Entries: " << fSamples[i].nevents << endl;
    if(fVerbose>2) cout << "\t\tAverage PU weight: " <<avg_pu_w  << endl;
    if(fVerbose>2) cout << "\t\t sample has weight " << sample_weight << " and " << fSamples[i].tree->GetEntries() << " entries" << endl;
    
    if( fSamples[i].tree->GetEntries()==0 ) continue;

    fMT2tree = new MT2tree();
    fSamples[i].tree->SetBranchAddress("MT2tree", &fMT2tree);
    Long64_t nbytes = 0, nb = 0;
 
    //create counter histos (but for W, done before)
    if(tsname!= "Wtolnu")  createCounterHistos(&counterH, cutLabels, fSamples[i].sname);

    //Filtering TTree
    TString myCuts = filter;
    if( fSamples[i].type=="data") myCuts += " && " + trigger;
    fSamples[i].tree->Draw(">>selList", myCuts);
    TEventList *myEvtList = (TEventList*)gDirectory->Get("selList");
    fSamples[i].tree->SetEventList(myEvtList);
    int counter=0;
    if(fVerbose>2) cout <<  "\t\t Filtering done, size=" <<myEvtList->GetN()  << endl;
    
    bool isW = false;
    if(tsname == "Wtolnu")      isW = true;

    //analysis
    if(myEvtList->GetN() ==0) continue;
    while(myEvtList->GetEntry(counter++) !=-1){
      int jentry = myEvtList->GetEntry(counter-1);
      nb =  fSamples[i].tree->GetEntry(jentry);   nbytes += nb;
      if ( fVerbose>2 && counter % 100000 == 0 )  cout << "+++ Processing event " << counter << endl;

      //Weights
      Double_t weight = 1;
      if (!fMT2tree->misc.isData  ) weight = sample_weight * fMT2tree->pileUp.Weight; // pile-up reweighting for MC

      tsname = fSamples[i].sname;
      if(fMT2tree->misc.isData) tsname="data";
      
      bool isWMU(false), isWELE(false), isWTAU(false);//, isWELSE(false);
      if(isW){
        if( fMT2tree->WDecayMode() == 1 ) {
	  if(ELEorMU=="MU") { tsname="Other"; }
	  else  tsname += "Ele" ;
	  isWELE = true;
        }
        else if( fMT2tree->WDecayMode() == 2 ){
          isWMU =true; 
	  tsname += "Mu" ;
        }
        else if( fMT2tree->WDecayMode() !=0 ){
          tsname += "Tau" ;
	  isWTAU =true;
        }
      }//end W cases
      
      
      //Variables
      float MT2 = fMT2tree->misc.MT2;
      float MET = fMT2tree->misc.MET;
      float HT = fMT2tree->misc.HT;
      float minDPhi = fMT2tree->misc.MinMetJetDPhi;
      int minDPhiIndex =  fMT2tree->misc.MinMetJetDPhiIndex;

      TLorentzVector METlv(0,0,0,0);
      METlv.SetPtEtaPhiM(MET, 0., fMT2tree->misc.METPhi, 0);

      float VSPt = fMT2tree->misc.Vectorsumpt;
      vector<float> px, py, pz, E; //for hemisph.
      int nBTags = 0;
      vector<int> bjets_i;
      vector<TLorentzVector> bjets_lv;
      int nAllBTags = 0;
      vector<int> allbjets_i;
      vector<TLorentzVector> allbjets_lv;
      float BEff = 0, BSqErr = 0;


      //lepton acceptance
      //Acceptance for W
      float V_pt = -1;
      float LeptMC_pt = -1;
      TLorentzVector LeptMC_lv;
      bool inAcceptance(false);
      
      if(isW){
        int PID = 13;
        float minpt = MU_PT;
        if(ELEorMU == "ELE") {
          PID = 11;  minpt = ELE_PT;
        }
        if( fMT2tree->GenLeptFromW(PID,minpt, 2.4, false)) inAcceptance = true;
        V_pt =  fMT2tree->GetGenVPt(24);
        LeptMC_pt = fMT2tree->GetGenLeptPt(0, PID, 24, minpt, 2.4);
	int li = fMT2tree->GetGenLeptIndex(0, PID, 24, minpt, 2.4);
	LeptMC_lv = fMT2tree->genlept[li].lv;
      }
      else if(tsname=="Zinv"){ V_pt = fMT2tree->GetGenVPt(23); } 

      //Jets
      int NJets = 0;
      for(int j=0; j< fMT2tree->NJets; j++){
	if(fMT2tree->jet[j].lv.Pt() > JET_PT) NJets++;
      }
      
      //bjets
      for(int jet=0; jet < fMT2tree->NJets; jet++){
        if(fMT2tree->jet[jet].lv.Pt() < BJET_PT || fabs(fMT2tree->jet[jet].lv.Eta())>2.4 ) continue;

        //BTag
        bool isBTagged = false;
        if(!BENRICH){
          if( fMT2tree->jet[jet].bTagProbSSVHE > 1.74 ) isBTagged=true;
        }
        else{
          if( fMT2tree->jet[jet].bTagProbSSVHP > 2. )  isBTagged=true;
        }

        if(isBTagged){
	  if(fMT2tree->jet[jet].lv.Pt()<500.){
	    nBTags++;
	    bjets_i.push_back(jet);
	    bjets_lv.push_back( fMT2tree->jet[jet].lv);
	  }
	  nAllBTags++;
	  allbjets_i.push_back(jet);
	  allbjets_lv.push_back( fMT2tree->jet[jet].lv);
	}
      }


      //Muons
      int NMuons = 0;
      TLorentzVector theMu;
      int theMu_i = -1;
      for( int mu=0; mu< fMT2tree->NMuons; mu++){
	if( fMT2tree->muo[mu].lv.Pt() < MU_PT) continue;
	theMu = fMT2tree->muo[mu].lv;
	theMu_i = mu;
	NMuons ++;
      }

      //Electrons
      int NEles = 0;
      TLorentzVector theEle;
      int theEle_i=-1;
      for( int el=0; el< fMT2tree->NEles; el++){
	if( fMT2tree->ele[el].lv.Pt() < ELE_PT) continue;

	if(ELEorMU=="ELE"){
	  if(fMT2tree->ele[el].ID90 !=7 ) continue;
	}
	theEle = fMT2tree->ele[el].lv;
	theEle_i=el;
	NEles ++;
      }

      //exactly 1 lepton
      //int NLepts= NMuons+NEles;
      TLorentzVector theLept;
      if(ELEorMU == "ELE") theLept = theEle;
      else  theLept = theMu;

      //MT calculation
      float MT=-999;
      if(NMuons==1) MT = fMT2tree->muo[theMu_i].MT;
      if(NEles==1) MT = fMT2tree->ele[theEle_i].MT;


      //Recomputing quantities 
      //if( CORRECT_MET && (ELEorMU=="MU" && NMuons==1) || (ELEorMU=="ELE" && NEles==1) /*NLepts==1*/){
      if( CORRECT_MET && ((ELEorMU=="MU" && NMuons==1) || (ELEorMU=="ELE" && NEles==1)) ){
        ///// Correcting MET
        double MET_x = MET*cos(fMT2tree->misc.METPhi);
        double MET_y = MET*sin(fMT2tree->misc.METPhi);

        //if(CORRECT_MET){
	MET_x += theLept.Px();
	MET_y += theLept.Py();
	
	MET = sqrt( MET_x*MET_x + MET_y*MET_y);
	METlv.SetPtEtaPhiM(MET, 0.,TMath::ATan2( MET_y,MET_x)  , 0.);
	
	//correction VSPt
	minDPhi = 999;
	TLorentzVector mht = fMT2tree->GetMHTlv(1,20,2.4, false);
	VSPt = (mht-METlv).Pt();
	
	//DPhi change for MT2b
	int NMaxJetsDPhi=9999, NJetsDPhi=0;
	if(MT2_REGIME=="MT2b") NMaxJetsDPhi = 4;
	for(int jet=0; jet < fMT2tree->NJets; jet++){
	  //correctiong dPhi
	  if(fMT2tree->jet[jet].IsGoodPFJet(20,5,0)==false) continue;
	  Double_t dphi = TMath::Abs( fMT2tree->jet[ jet ].lv.DeltaPhi(METlv));
	  if(dphi < minDPhi && NJetsDPhi<NMaxJetsDPhi) minDPhi = dphi;
	  NJetsDPhi++;
	  if(fMT2tree->jet[jet].lv.Pt() < 20. || fabs(fMT2tree->jet[jet].lv.Eta())>2.4 ) continue;
	  
	  //filing LV for Hemis
	  if( fMT2tree->jet[jet].isPATPFIDLoose!=1) continue;
	  px.push_back(fMT2tree->jet[jet].lv.Px());
	  py.push_back(fMT2tree->jet[jet].lv.Py());
	  pz.push_back(fMT2tree->jet[jet].lv.Pz());
	  E .push_back(fMT2tree->jet[jet].lv.E ());
	  
	}//end jets
	
	// get hemispheres (seed 2: max inv mass, association method: default 3 = minimal lund distance)
	if (  px.size()>1){ 
	  Hemisphere* hemi = new Hemisphere(px, py, pz, E, 2, 3);
	  vector<int> grouping = hemi->getGrouping();
	  
	  TLorentzVector pseudojet1(0.,0.,0.,0.);
	  TLorentzVector pseudojet2(0.,0.,0.,0.);
	  
	  for(unsigned int j=0; j<px.size(); ++j){
	    if(grouping[j]==1){
	      pseudojet1.SetPx(pseudojet1.Px() + px[j]);              pseudojet1.SetPy(pseudojet1.Py() + py[j]);
	      pseudojet1.SetPz(pseudojet1.Pz() + pz[j]);              pseudojet1.SetE( pseudojet1.E()  + E[j]);
	    }else if(grouping[j] == 2){
	      pseudojet2.SetPx(pseudojet2.Px() + px[j]);              pseudojet2.SetPy(pseudojet2.Py() + py[j]);
	      pseudojet2.SetPz(pseudojet2.Pz() + pz[j]);              pseudojet2.SetE( pseudojet2.E()  + E[j]);
	    }
	  }
	  delete hemi;
	  
	  MT2 =  fMT2tree->CalcMT2(0, false, pseudojet1, pseudojet2, METlv);
	}
	//}//end CORRECT_MET
      }//end NLepts==1
      

      //Efficiencies correction
      float LeptEff=0,  LeptSqErr=0;
      float TotEff=0, EffSqPercErr = 0;

      if(theLept.Pt()>0 ){
	float MinDR_LeptJet = 9999.;
        //minDR
	
        for(int jet=0; jet < fMT2tree->NJets; jet++){
          if(fMT2tree->jet[jet].lv.Pt() < 20.) continue; // || fabs(fMT2tree->jet[jet].lv.Eta())>5 ) continue;
	  Double_t dr =  fMT2tree->jet[ jet ].lv.DeltaR( theLept);
          if(dr < MinDR_LeptJet)  MinDR_LeptJet = dr;
	}

	
	//here, the efficiencies are in pt,eta bins
        if(theLept.Pt() < 20. ){
          int bin =     h2_effMu_pteta->FindBin(theLept.Pt(), fabs(theLept.Eta()));
          if(bin!=1){
            LeptEff = h2_effMu_pteta->GetBinContent( bin  );
            LeptSqErr = h2_effMu_pteta->GetBinError( bin  );
            LeptSqErr *= LeptSqErr;
          }
          else
            cout << "ERROR: eff bin not found for " << theLept.Pt() << " " << theLept.Eta() << endl;
        }
        else{
          int bin = -1;
          float tmpDr=MinDR_LeptJet, tmpPt=theLept.Pt();

          // for leptons with dR<0.5, I'm taking dR=0.5. Same for pt
          if(MinDR_LeptJet <0.5) tmpDr = 0.5;
          if(theLept.Pt() >200.) tmpPt = 199;
	  
          bin = h2_effMu_ptdr->FindBin(tmpPt, tmpDr );
          if(bin!=-1){
            LeptEff = h2_effMu_ptdr->GetBinContent( bin  );
            LeptSqErr = h2_effMu_ptdr->GetBinError( bin  );
            LeptSqErr *= LeptSqErr;
	  }
	  else
            cout << "ERROR: eff bin not found for " << theLept.Pt() << " " << tmpDr << endl;
	}
	
	TotEff *= LeptEff;
	EffSqPercErr += LeptSqErr/(LeptEff*LeptEff); //sum in quadrature of relative errors 
      }

      


      //selections

      //MT2 sel
      if(MT2 <=  MT2CUT_LOW  || MT2 >   MT2CUT_HIGH) continue;

      TString mysel = "";
      mysel = "preSel";
      if( !(NJets >= MIN_NJETS) ) weight=0;
      counterH[mysel][tsname]->Fill(MT2,weight);
      if(type=="mc")  counterH[mysel]["MC"]->Fill(MT2,weight);
      //Vector boson info
      histos["VectBoson_pt_"+tsname+"_"+mysel]->Fill(V_pt, weight);

      mysel = "dphiSel";
      if( minDPhi <= 0.3 ) weight = 0;
      //if( minDPhi > 0.3 ) weight = 0;
      counterH[mysel][tsname]->Fill(MT2,weight);
      if(type=="mc")  counterH[mysel]["MC"]->Fill(MT2,weight);
      
      mysel = "hbheSel";
      if( fMT2tree->misc.HBHENoiseFlag != 0 )  weight=0;
      if( fMT2tree->misc.CrazyHCAL     != 0 )  weight=0;
      if(fMT2tree->misc.CSCTightHaloID != 0 )  weight=0;
      counterH[mysel][tsname]->Fill(MT2,weight);
      if(type=="mc")  counterH[mysel]["MC"]->Fill(MT2,weight);

      mysel = "metSel";
      if( MET <=30 ) weight=0;
      counterH[mysel][tsname]->Fill(MT2,weight);
      if(type=="mc")  counterH[mysel]["MC"]->Fill(MT2,weight);
	

      mysel = "vsptSel";
      if( VSPt >= 70. )  weight=0;
      counterH[mysel][tsname]->Fill(MT2,weight);
      if(type=="mc")  counterH[mysel]["MC"]->Fill(MT2,weight);


      mysel = "jidSel";
      if( !fMT2tree->misc.PassJetID )  weight=0;
      counterH[mysel][tsname]->Fill(MT2,weight);
      if(type=="mc")  counterH[mysel]["MC"]->Fill(MT2,weight);
      //Acceptance
      if(inAcceptance){
	if(isWELE and ELEorMU=="ELE") histos["MT2_LeptAcc_"+tsname+"_jidSel"]->Fill(MT2,weight);
	else if  (isWMU and ELEorMU=="MU") histos["MT2_LeptAcc_"+tsname+"_jidSel"]->Fill(MT2,weight);
      }
      if( (isW && inAcceptance ) || tsname=="Zinv")histos["VectBoson_pt_"+tsname+"_"+mysel]->Fill(V_pt, weight);


      mysel = "leptVetoSel";
      if(ELEorMU == "MU"){ if( NEles!=0 ) weight=0; }
      else { if( NMuons!=0 ) weight=0; }
      counterH[mysel][tsname]->Fill(MT2,weight);
      if(type=="mc")  counterH[mysel]["MC"]->Fill(MT2,weight);
      //Vector boson info
      if(weight!=0 && ( (isW && inAcceptance ) || tsname=="Zinv" ))histos["VectBoson_pt_"+tsname+"_"+mysel]->Fill(V_pt, weight);
     
      mysel = "leptSel";
      if(ELEorMU == "MU") { if( NMuons!=1 ) weight=0; }
      else { if( NEles!=1 ) weight=0; }
      counterH[mysel][tsname]->Fill(MT2,weight);
      if(type=="mc")  counterH[mysel]["MC"]->Fill(MT2,weight); 
      //filling lepton/btag efficiencies
      if(weight!=0){
        histos["Lept_eff_TP_"+(string)tsname]->Fill( LeptEff, weight);
        varFloat["Lept_eff_TP_err_"+(string)tsname] += LeptSqErr;
        if( type=="mc"){
          histos["Lept_eff_TP_MC"]->Fill( LeptEff, weight);
          varFloat["Lept_eff_TP_err_MC"] += LeptSqErr;
        }

        if(nBTags!=0){
          histos["B_eff_"+(string)tsname]->Fill( BEff, weight);
          varFloat["B_eff_err_"+(string)tsname] += BSqErr;
          if( type=="mc"){
            histos["B_eff_MC"]->Fill( BEff, weight);
            varFloat["B_eff_err_MC"] += BSqErr;
          }
        }
	//also >500GeV
	for(int b=0; b<nAllBTags; b++){
          histos["BJets_pt_"+tsname+"_"+mysel]->Fill( allbjets_lv[b].Pt(), weight);
	}
      }
      
      if(weight!=0){
	histos["Lept_pt_"+mysel+"_"+tsname]->Fill(theLept.Pt(), weight);
	histos["Lept_eta_"+mysel+"_"+tsname]->Fill(theLept.Eta(), weight);
	histos["MT_"+tsname+"_"+mysel]->Fill(MT, weight);
       
	//histos
	histos["MT2_"+tsname+"_"+mysel]->Fill(MT2, weight);
	histos["BJets_N_"+tsname+"_"+mysel]->Fill(nBTags, weight);
	histos["Jets_N_"+tsname+"_"+mysel]->Fill(NJets, weight);
	histos["HT_"+tsname+"_"+mysel]->Fill(HT, weight);
	histos["MET_"+tsname+"_"+mysel]->Fill(MET, weight);
	histos["MinDPhi_"+tsname+"_"+mysel]->Fill(minDPhi, weight);
	
	histos["MinDPhiIndex_"+tsname+"_"+mysel]->Fill(minDPhiIndex, weight);
	histos["MinDPhiPt_"+tsname+"_"+mysel]->Fill( fMT2tree->jet[minDPhiIndex].lv.Pt(), weight);
	histos["MinDPhiEta_"+tsname+"_"+mysel]->Fill( fMT2tree->jet[minDPhiIndex].lv.Eta(), weight);
	histos["MinDPhiPhi_"+tsname+"_"+mysel]->Fill( fMT2tree->jet[minDPhiIndex].lv.Phi(), weight);
	histos["MinDPhiMuDR_"+tsname+"_"+mysel]->Fill( fMT2tree->jet[minDPhiIndex].lv.DeltaR(theMu) , weight);
	
	histos["NVertices_"+tsname+"_"+mysel]->Fill( fMT2tree->pileUp.NVertices , weight);
	histos["PUWeight_"+tsname+"_"+mysel]->Fill(  fMT2tree->pileUp.Weight);
      }
      
      
      histos["caloHT_"+tsname+"_"+mysel]->Fill( fMT2tree->misc.caloHT50_ID, weight);
      histos["VSPt_"+tsname+"_"+mysel]->Fill( VSPt, weight);


      //MT with bsel
      if(BENRICH) { if(nBTags!=0) histos["MT_"+tsname+"_bSel"]->Fill(MT, weight);}
      else { if(nBTags==0) histos["MT_"+tsname+"_bSel"]->Fill(MT, weight);}
      mysel = "mtSel";
      if( MT >= 100 ) weight=0;
      counterH[mysel][tsname]->Fill(MT2,weight);
      if(type=="mc")  counterH[mysel]["MC"]->Fill(MT2,weight);
  

      mysel = "bSel";
      if(BENRICH) { if(nBTags==0) weight=0;}
      else  { if(nBTags!=0) weight=0;}
      counterH[mysel][tsname]->Fill(MT2,weight);
      if(type=="mc")  counterH[mysel]["MC"]->Fill(MT2,weight);

      //Vector boson info
      histos["VectBoson_pt_"+tsname+"_"+mysel]->Fill(V_pt, weight);


      //Histo Filling (before MT2 cuts)

      histos["MT2_"+tsname+"_"+mysel]->Fill(MT2, weight);
      histos["Jets_N_"+tsname+"_"+mysel]->Fill(NJets, weight);
      histos["HT_"+tsname+"_"+mysel]->Fill(HT, weight);
      histos["caloHT_"+tsname+"_"+mysel]->Fill( fMT2tree->misc.caloHT50_ID, weight);
      histos["MET_"+tsname+"_"+mysel]->Fill(MET, weight);
      histos["MinDPhi_"+tsname+"_"+mysel]->Fill(minDPhi, weight);
      histos["VSPt_"+tsname+"_"+mysel]->Fill( VSPt, weight);
      if(weight!=0){
	histos["Lept_pt_"+mysel+"_"+tsname]->Fill(theLept.Pt(), weight);
	histos["Lept_eta_"+mysel+"_"+tsname]->Fill(theLept.Eta(), weight);
	histos["PUWeight_"+tsname+"_"+mysel]->Fill(  fMT2tree->pileUp.Weight);
	
	//PUWeight_avg +=  fMT2tree->pileUp.Weight;
	NSelEvents += 1;
      }
      histos["NVertices_"+tsname+"_"+mysel]->Fill( fMT2tree->pileUp.NVertices , weight);


	
      if(MT2CUT_LOW == -1){
	float oldW=weight;

	mysel = "mt100Sel";
	if( !(MT2 > 100 && MT2<=150) ) weight=0;
	counterH[mysel][tsname]->Fill(MT2,weight);
	if(type=="mc")  counterH[mysel]["MC"]->Fill(MT2,weight);
	
	weight=oldW;
	mysel = "mt150Sel";
	if( !(MT2 > 150 && MT2<=200) ) weight=0;
	counterH[mysel][tsname]->Fill(MT2,weight);
	if(type=="mc")  counterH[mysel]["MC"]->Fill(MT2,weight);
	
	weight=oldW;
	mysel = "mt200Sel";
	if( !(MT2 > 200 && MT2<=275) ) weight=0;
	counterH[mysel][tsname]->Fill(MT2,weight);
	if(type=="mc")  counterH[mysel]["MC"]->Fill(MT2,weight);
	
	weight=oldW;
	mysel = "mt275Sel";
	if( !(MT2 > 275 && MT2<=375) ) weight=0;
	counterH[mysel][tsname]->Fill(MT2,weight);
	if(type=="mc")  counterH[mysel]["MC"]->Fill(MT2,weight);
	
	weight=oldW;
	mysel = "mt375Sel";
	if( !(MT2 > 375 && MT2<=500) ) weight=0;
	counterH[mysel][tsname]->Fill(MT2,weight);
	if(type=="mc")  counterH[mysel]["MC"]->Fill(MT2,weight);
	
	weight=oldW;
	mysel = "mt500Sel";
	if( MT2 <= 500 ) weight=0;
	counterH[mysel][tsname]->Fill(MT2,weight);
	if(type=="mc")  counterH[mysel]["MC"]->Fill(MT2,weight);
      }
    }//end event loop
    
  }//end sample loop


  /////// Print a Latex table
  ////// Print a TWIKI table
  printTwikiTable(  samplesName, cutLabels, orderedCutLabels, counterH );
  cout << endl;
  if(MT2CUT_LOW == -1) printLatexTable(  samplesName, cutLabels, orderedCutLabels, counterH );
  

  //printing effs etc
  float eff, err;
  TString leptLabel = "";
  if(ELEorMU=="ELE")    leptLabel="Ele";
  else if(ELEorMU=="MU")    leptLabel="Mu";
    
  //acceptance
  float lept_acc, lept_acc_err;
  getEfficiency(histos["MT2_LeptAcc_Wtolnu"+leptLabel+"_jidSel"], counterH["jidSel"]["Wtolnu"+leptLabel], &lept_acc, &lept_acc_err);
  histos["MT2_LeptAcc_Wtolnu"+leptLabel+"_jidSel"]->Divide(counterH["jidSel"]["Wtolnu"+leptLabel] );
  //printf("|*W(%snu) acceptance after preselection* | %.2f +- %.2f |\n", leptLabel.Data(), lept_acc, lept_acc_err);

  //VB stat, err
  double z_mc_err;
  float z_mc = histos["VectBoson_pt_Zinv_preSel"]->IntegralAndError(1 ,histos["VectBoson_pt_Zinv_preSel"]->GetNbinsX() , z_mc_err);
  double w_mc_err;
  float w_mc = histos["VectBoson_pt_Wtolnu"+leptLabel+"_preSel"]->IntegralAndError(1 ,histos["VectBoson_pt_Wtolnu"+leptLabel+"_preSel"]->GetNbinsX() , w_mc_err);

  cout << endl;
  cout << "|*T&P efficiencies*|" << endl;
  printBLeptEff("Wtolnu"+leptLabel, histos, varFloat);
  printBLeptEff("Top", histos, varFloat);
  printBLeptEff("MC", histos, varFloat);
  printBLeptEff("data", histos, varFloat);
  cout << endl;
  

  float *result = new float[10];
  /////////////
  // The Estimates
  ////////////
  /// Common inputs
  float zw_ratio_pdf_err_perc = 0.014;
  float lept_acc_pdf_err_perc = 0.01;

  //BTag scale factors
  float ssvhpt_scale = 0.91;
  float ssvhpt_scale_err = sqrt(0.02*0.02 + 0.09*0.09);
  //float ssvhpt_scale_syst = 0.09;
  //  float ssvhpt_scale_err_perc = sqrt( ssvhpt_scale_err* ssvhpt_scale_err +ssvhpt_scale_syst*ssvhpt_scale_syst )/ssvhpt_scale;
  float ssvhpt_scale_err_perc = ssvhpt_scale_err/ssvhpt_scale;

  float ssvhem_scale = 0.95;
  float ssvhem_scale_err = sqrt(0.01*0.01 + 0.10*0.10);
  //  float ssvhem_scale_syst = 0.10;
  float ssvhem_scale_err_perc = ssvhem_scale_err /ssvhem_scale;

  //choosing between DD est and MC closure
  TString dataSample = "data";
  if(IS_MC) dataSample = "MC";


  if(BENRICH){
    //    if(MT2_REGIME=="MT2"){
    //MC prediction for top, after the MT cut
    double mc_pred_err;
    float mc_pred = counterH["mtSel"]["Top"]->IntegralAndError(1, counterH["mtSel"]["Top"]->GetNbinsX()+1, mc_pred_err);
    double data_err;
    float data = counterH["bSel"][dataSample]->IntegralAndError(1,counterH["bSel"][dataSample]->GetNbinsX()+1, data_err );
    float data_err_perc = data_err/data;
    
    //Bkgs from MC
    TH1F* h_mc_bkg;
    if(ELEorMU=="MU") h_mc_bkg = (TH1F*)counterH["bSel"]["WtolnuMu"]->Clone();
    else  h_mc_bkg = (TH1F*)counterH["bSel"]["WtolnuEle"]->Clone(); 
    h_mc_bkg->SetName("h_mc_bkg");
    h_mc_bkg->Add(counterH["bSel"]["WtolnuTau"]);  
    h_mc_bkg->Add(counterH["bSel"]["QCD"]); 
    h_mc_bkg->Add(counterH["bSel"]["Other"]);
    double mc_bkg_err;
    float mc_bkg = h_mc_bkg->IntegralAndError(1 ,h_mc_bkg->GetNbinsX()+1 , mc_bkg_err);
    mc_bkg_err = mc_bkg; //100% error on bkg subtraction

    //b-sel eff from MC for TOP
    float b_eff, b_eff_err;
    getEfficiency(counterH["bSel"]["Top"], counterH["mtSel"]["Top"], &b_eff, &b_eff_err);
    float b_eff_err_perc = b_eff_err/b_eff;

    //TOP est computation
    if(!IS_MC) {
      mc_bkg *= ssvhpt_scale; //SF for MC bkg
      mc_bkg_err = mc_bkg;
      //mc_bkg_err = sqrt(  mc_bkg_err* mc_bkg_err/(mc_bkg*mc_bkg) + ssvhpt_scale_err_perc*ssvhpt_scale_err_perc ); //error is already 100%
      mc_pred *= ssvhpt_scale;
      mc_pred_err = mc_pred*sqrt(  mc_pred_err* mc_pred_err/(mc_pred*mc_pred) + ssvhpt_scale_err_perc*ssvhpt_scale_err_perc );
    }
    float est = data - mc_bkg;
    float est_syst_err=0;
    float est_stat_err = 0;
    //Sanity check on ttbar est
    if(est <= 0) {
      cout << "\n*[WARNING] Negative ttbar estimate, using MC prediction with 100% syst error and no stat err for 1lept,1b*\n" << endl;
      est = -1;//counterH["bSel"]["Top"]->Integral();
      est_syst_err = 0;
      est_stat_err = 0;
    }
    else {
      //BTag SF (data only)
      if(!IS_MC) { b_eff *=ssvhpt_scale; }
      //beff correction (from MC)
      est /= b_eff; 
      
      //errors
      est_syst_err += mc_bkg_err*mc_bkg_err; //MC bkg subtraction, 100% error
      if(!IS_MC) est_syst_err += est*est* ssvhpt_scale_err_perc* ssvhpt_scale_err_perc; // SF errors
      est_syst_err += est*est*b_eff_err_perc*b_eff_err_perc ; //MC errors on b-sel eff
      est_syst_err = sqrt(est_syst_err);
      
      if(data!=0) est_stat_err = est*data_err_perc;
      if(est <= 0) est_stat_err = 0;
      float est_err  =  sqrt( est_stat_err* est_stat_err + est_syst_err*est_syst_err);
      
      //printouts
      printf("|*BKG (1b, 1l, MC)* | %.3f +- %.3f |\n", mc_bkg, mc_bkg_err);
      printf("|*b-tag sel eff (1l, MC)* | %.3f +- %.3f |\n", b_eff, b_eff_err);

      printf("|*Top Estimate (1l) after MT cut* | %.3f +- %.3f ( %.3f stat +- %.3f syst   )|\n", est, est_err, est_stat_err, est_syst_err);
      printf("|*Top MC pred (1l)  after MT cut* | %.3f +- %.3f |\n", mc_pred, mc_pred_err);
      
    }    

    //latex output
    cout << "\n\n" << endl;
    printf("\\begin{table}[h!]\n \\begin{center}\n \\begin{tabular}{l|ccccccccc} \n   \\hline\\hline\n   %\\multirow{2}{*}{$M_{T2}$\\ bin} & \\multicolumn{4}{c}{$H_T$}   \\\\ \n");
    printf("               & $N^{mu}$         & $N^{bkg}$        & data pred.       & MC pred. \\\\ \n ");
    printf("TT $%.0f - %.0f$ & %.2f $\\pm$ %.2f & %.2f $\\pm$ %.2f & %.2f $\\pm$ %.2f &  \\\\ \n", 
	   MT2CUT_LOW, MT2CUT_HIGH,
	   data, data_err,
	   mc_bkg, mc_bkg_err,
	   est, sqrt( est_stat_err* est_stat_err + est_syst_err*est_syst_err ),
	   mc_pred, mc_pred_err
	   );
    printf("   \\hline\\hline \n \\end{tabular} \n \\end{center} \n \\caption{Results for $H_T \\in () \\ \\mathrm{GeV}$.}\\label{table:} \n \\end{table}");
    cout << "\n\n" << endl;
      




    //input for next step
    result[0] = est;
    result[1] = est_stat_err;
    result[2] = est_syst_err;
    result[3] = mc_pred;
    result[4] = mc_pred_err;
    //    }
    //////////
    //MT2b
    if(MT2_REGIME=="MT2b"){
      //here, the only needed input for the estimate is the WJets bsel eff
      //b-sel eff from MC for WJets
      //float b_eff, b_eff_err;
      getEfficiency(counterH["bSel"]["Wtolnu"+leptLabel], counterH["mtSel"]["Wtolnu"+leptLabel], &b_eff, &b_eff_err);
      //float 
      b_eff_err_perc = b_eff_err/b_eff;
      
      //outputs for estimate: btag effs, with SF
      result[5] = b_eff*ssvhpt_scale;
      result[6] = sqrt(b_eff_err_perc*b_eff_err_perc + ssvhpt_scale_err_perc*ssvhpt_scale_err_perc );
    }
  }
  ////////////////////////////////////////
  //ZInv estimate
  /////////////////////////////////////
  else{
    //// Z/W Ratios
    cout << "|*ZW Ratios*|" << endl;
    float zw_ratio, zw_ratio_err;
    getEfficiency(histos["VectBoson_pt_Zinv_preSel"], histos["VectBoson_pt_Wtolnu"+leptLabel+"_preSel"], &eff, &err);
    printf("|*Z/W(%snu) ratio after preselection* | %.2f +- %.2f |\n", leptLabel.Data(), eff, err);
    getEfficiency(histos["VectBoson_pt_Zinv_leptVetoSel"], histos["VectBoson_pt_Wtolnu"+leptLabel+"_leptVetoSel"], &zw_ratio, &zw_ratio_err);
    printf("|*Z/W(%snu) ratio after lepton veto* | %.2f +- %.2f |\n", leptLabel.Data(), zw_ratio, zw_ratio_err);
    getEfficiency(histos["VectBoson_pt_Zinv_jidSel"], histos["VectBoson_pt_Wtolnu"+leptLabel+"_jidSel"], &zw_ratio, &zw_ratio_err);
    printf("|*Z/W(%snu) ratio after JID Veto (used)* | %.2f +- %.2f |\n", leptLabel.Data(), zw_ratio, zw_ratio_err);
    //getEfficiency(counterH["jidSel"]["Zinv"], counterH["jidSel"]["Wtolnu"+leptLabel], &zw_ratio, &zw_ratio_err);
    //printf("|*Z/W(%snu) ratio after JID Veto (used)* | %.2f +- %.2f |\n", leptLabel.Data(), zw_ratio, zw_ratio_err);

    //FIXME TEMP LEO DAMNED BRUNO IS GAY
    //zw_ratio=0.87;
    //zw_ratio_err = 0.02;
 
   
    double data_err;
    float data = counterH["bSel"][dataSample]->IntegralAndError(1,counterH["bSel"][dataSample]->GetNbinsX()+1, data_err );
    float data_err_perc = data_err/data;   
    //    if(MT2_REGIME=="MT2"){
           
    //bkgs from MC
    TH1F* h_mc_bkg;
    //if(leptLabel=="Ele") h_mc_bkg = (TH1F*)counterH["bSel"]["WtolnuMu"]->Clone(); 
    //else 
    h_mc_bkg = (TH1F*)counterH["bSel"]["Other"]->Clone(); h_mc_bkg->SetName("h_mc_bkg");
    h_mc_bkg->SetName("h_mc_bkg");
    h_mc_bkg->Add(counterH["bSel"]["WtolnuTau"]);
    //h_mc_bkg->Add(counterH["bSel"]["Other"]);
    //h_mc_bkg->Add(counterH["bSel"]["VV"]);
    //h_mc_bkg->Add(counterH["bSel"]["QCD"]); //should consider QCD negligible?
    
    double mc_bkg_err;
    float mc_bkg = h_mc_bkg->IntegralAndError(1 ,h_mc_bkg->GetNbinsX() , mc_bkg_err);
    
    float tt_bkg = INPUT[0];
    float tt_bkg_err_stat = INPUT[1];
    float tt_bkg_err_syst = INPUT[2];
    //MC prediction for top 1l
    double mc_tt_bkg_err;
    float mc_tt_bkg = counterH["leptSel"]["Top"]->IntegralAndError(1, counterH["leptSel"]["Top"]->GetNbinsX()+1, mc_tt_bkg_err );
    //MC prediction for top 1l
    double mc_tt0b1l_bkg_err;
    float mc_tt0b1l_bkg = counterH["bSel"]["Top"]->IntegralAndError(1, counterH["bSel"]["Top"]->GetNbinsX()+1, mc_tt0b1l_bkg_err );
    
    //sanity check for ttbar estimate
    //TTBar estimate 
    //you should use tt_bar estimate without bveto rescalin, 1lept 
    float tt_bveto_eff , tt_bveto_eff_err;
    getEfficiency( counterH["bSel"]["Top"], counterH["leptSel"]["Top"], &tt_bveto_eff, &tt_bveto_eff_err);
    float tt_bveto_eff_err_perc = tt_bveto_eff_err/tt_bveto_eff;
    //sanity check for ttbar estimate
    if ( tt_bkg <= 0 ){
      tt_bkg=INPUT[3]*tt_bveto_eff;
      tt_bkg_err_stat = 0;
      tt_bkg_err_syst = INPUT[4];
    } 
    else{
      tt_bkg *= tt_bveto_eff;
      tt_bkg_err_syst = sqrt(tt_bveto_eff*tt_bveto_eff*tt_bkg_err_syst*tt_bkg_err_syst  + tt_bkg*tt_bkg*tt_bveto_eff_err_perc*tt_bveto_eff_err_perc   );// sum in quad of systs
      tt_bkg_err_stat = tt_bveto_eff*tt_bkg_err_stat; // scaling of stat unc
    }
    float tt_bkg_err = sqrt(tt_bkg_err_syst*tt_bkg_err_syst + tt_bkg_err_stat*tt_bkg_err_stat);
    //WJets est
    if(!IS_MC) mc_bkg *= ssvhem_scale;
    float w_est = data - mc_bkg - tt_bkg;
    
    //errors
    //float ssvhem_scale_err_perc = sqrt(ssvhem_scale_err*ssvhem_scale_err + ssvhem_scale_syst*ssvhem_scale_syst) /  ssvhem_scale; //SF
    //float data_err_perc = 1./sqrt(data);
    float w_est_err_stat = sqrt(w_est*w_est*data_err_perc*data_err_perc + tt_bkg_err_stat*tt_bkg_err_stat); // stat error: data and ttbar stat err
    float w_est_err_syst = 0; // syst error
    w_est_err_syst = mc_bkg_err*mc_bkg_err + tt_bkg_err_syst*tt_bkg_err_syst;
    //      w_est_err_syst += w_est*w_est* mt_eff_err_perc* mt_eff_err_perc; //MT correction
    if(!IS_MC) w_est_err_syst += w_est*w_est*ssvhem_scale_err_perc*ssvhem_scale_err_perc;
    w_est_err_syst = sqrt(w_est_err_syst);
    float w_est_err = sqrt(w_est_err_stat*w_est_err_stat + w_est_err_syst*w_est_err_syst );

    cout<<endl;
    printf("|* bVeto eff Top (MC)* | %.2f +- %.2f |\n", tt_bveto_eff, tt_bveto_eff_err);
    printf("|* BKG 0b,1l (MC)* | %.2f +- %.2f |\n", mc_bkg, mc_bkg_err);
    printf("|* BKG 0b,1l (Top)* | %.2f +- %.2f +- %.2f|\n", tt_bkg, tt_bkg_err_stat, tt_bkg_err_syst);
    printf("|* WJets Estimate (0b, 1l)* | %.2f +- %.2f +- %.2f |\n", w_est, w_est_err_stat, w_est_err_syst);
    
    
    //TP eff
    float tp_eff, tp_eff_err_stat, tp_eff_err_syst, tp_eff_err;
    TString TP_SAMPLE = "data";
    
    tp_eff = histos["Lept_eff_TP_"+dataSample]->GetMean();
    tp_eff_err_stat =  histos["Lept_eff_TP_"+TP_SAMPLE]->GetMeanError()*histos["Lept_eff_TP_"+TP_SAMPLE]->GetMeanError(); //stat error
    tp_eff_err_syst = sqrt(varFloat["Lept_eff_TP_err_"+TP_SAMPLE] )/histos["Lept_eff_TP_"+TP_SAMPLE]->GetEntries(); //sys error, sqrt(sum(err**2))/N
    tp_eff_err = sqrt( tp_eff_err_stat* tp_eff_err_stat + tp_eff_err_syst*tp_eff_err_syst );
    
    ///ZInv Est
    //bveto eff
    float wj_bveto_eff , wj_bveto_eff_err;
    getEfficiency( counterH["bSel"]["Wtolnu"+leptLabel], counterH["mtSel"]["Wtolnu"+leptLabel], &wj_bveto_eff, &wj_bveto_eff_err);
    //mt-sel eff from MC for WJets
    float mt_eff, mt_eff_err;
    getEfficiency(counterH["mtSel"]["Wtolnu"+leptLabel], counterH["leptSel"]["Wtolnu"+leptLabel], &mt_eff, &mt_eff_err);
    float mt_eff_err_perc = mt_eff_err/mt_eff;
    
    float z_est = w_est;
    float mc_factor = zw_ratio;
    mc_factor /=  wj_bveto_eff;
    //mc_factor /= lept_acc;
    //mc_factor /= tp_eff;
    mc_factor /= mt_eff;

    z_est *= mc_factor;
    z_est /= lept_acc;
    z_est /= tp_eff;

    if(MT2_REGIME=="MT2b"){
      z_est *= INPUT[5];//multiply for bsel eff, SSVHP
    }
    
    
    /////////////////
    //errors
    ////////////////
    //ZW
    float z_est_err = 0;
    zw_ratio_err *= zw_ratio_err;
    zw_ratio_err +=  zw_ratio*zw_ratio*zw_ratio_pdf_err_perc* zw_ratio_pdf_err_perc;
    zw_ratio_err = sqrt(zw_ratio_err);
    float zw_ratio_err_perc = zw_ratio_err/zw_ratio;
    //TP
    float tp_eff_err_perc = tp_eff_err/tp_eff;
    //ACC
    float lept_acc_err_perc = sqrt( lept_acc_err*lept_acc_err + lept_acc*lept_acc* lept_acc_pdf_err_perc* lept_acc_pdf_err_perc)/lept_acc;
    //DATA
    float z_est_stat_err = z_est*w_est_err_stat/w_est;
    //MC
    float mc_err =  mc_bkg_err;
    
    //SUMMING UP
    float syst_err = zw_ratio_err_perc*zw_ratio_err_perc +  tp_eff_err_perc* tp_eff_err_perc + lept_acc_err_perc*lept_acc_err_perc;
    syst_err *= z_est*z_est;
    syst_err += mc_err*mc_err; // not completely correct
    syst_err += z_est*z_est*(w_est_err_syst/w_est)*(w_est_err_syst/w_est); //syst err on W est
    syst_err += z_est*z_est* mt_eff_err_perc* mt_eff_err_perc; //MT correction
    
    //error on mc_factor
    float mc_factor_err_perc = zw_ratio_err_perc*zw_ratio_err_perc + /*lept_acc_err_perc*lept_acc_err_perc + tp_eff_err_perc*tp_eff_err_perc*/ + mt_eff_err_perc*mt_eff_err_perc;
    mc_factor_err_perc = sqrt(mc_factor_err_perc);
    float mc_factor_err = mc_factor*mc_factor_err_perc;

    if(MT2_REGIME=="MT2b"){
      syst_err += z_est*z_est*INPUT[6]*INPUT[6];
    }
    
    syst_err = sqrt(syst_err);
    
    z_est_err = sqrt( z_est_stat_err*z_est_stat_err + syst_err*syst_err);
    
    cout << endl;
    cout << "|*Eff./Acc.*|" << endl;
    printf("|*Lept. Acceptance* | %.3f +- %.3f |\n", lept_acc, lept_acc_err);
    printf("|*T&P efficiency* | %.3f +- %.3f |\n", tp_eff, tp_eff_err);
    cout << endl;
    cout << "|*MC factor*  | " << endl;
    printf("|*Z/W* | %.3f +- %.3f |\n", zw_ratio, zw_ratio_err);
    printf("|*MT sel eff (1l, MC)* | %.3f +- %.3f |\n", mt_eff, mt_eff_err);
    printf("|*bVeto eff WJets (MC)* | %.2f +- %.2f |\n", wj_bveto_eff, wj_bveto_eff_err);
    if(MT2_REGIME=="MT2b")  printf("|*bSel eff WJets (MC)* | %.2f +- %.2f |\n", INPUT[5], INPUT[5]*INPUT[6]);
    printf("|*MC FACTOR* | %.2f +- %.2f |\n", mc_factor, mc_factor_err);
    
    //MC predictions
    double w_mc_pred_err;
    float w_mc_pred = counterH["bSel"]["Wtolnu"+leptLabel]->IntegralAndError(1, counterH["bSel"]["Wtolnu"+leptLabel]->GetNbinsX()+1, w_mc_pred_err);
    
    cout << "\n\n" << endl;
    printf( "| | *Top est (0b, 1l)* | *Top MC (0b, 1l)* | *WJets est (0b, 1l)* | *WJets MC (0b,1l)* | *ZInv est* | *ZInv MC* |\n" );
    printf( "| *MT2 [%.0f, %.0f)* | %.2f +- %.2f ( %.2f (stat) +- %.2f (syst) ) | %.2f +- %.2f |  %.2f +- %.2f ( %.2f (stat) +- %.2f (syst) ) |  %.2f +- %.2f |  %.2f +- %.2f ( %.2f (stat) +- %.2f (syst) ) | |\n", 
	    MT2CUT_LOW, MT2CUT_HIGH, 
	    //INPUT[0], sqrt(INPUT[1]*INPUT[1] +INPUT[2]*INPUT[2] ), INPUT[1], INPUT[2], 
	    //INPUT[3], INPUT[4],
	    tt_bkg, sqrt(tt_bkg_err_syst* tt_bkg_err_syst + tt_bkg_err_stat*tt_bkg_err_stat),  tt_bkg_err_stat,  tt_bkg_err_syst,
	    //mc_tt_bkg, mc_tt_bkg_err,
	    //INPUT[3], INPUT[4],
	    mc_tt0b1l_bkg, mc_tt0b1l_bkg_err,
	    w_est, w_est_err, w_est_err_stat, w_est_err_syst, 
	    w_mc_pred, w_mc_pred_err,
	    z_est, z_est_err, z_est_stat_err, syst_err  );
    
    
    //    }    

    cout << "\n\n" << endl;
    printf("\\begin{table}[h!]\n \\begin{center}\n \\begin{tabular}{l|ccccccccc} \n   \\hline\\hline\n   %\\multirow{2}{*}{$M_{T2}$\\ bin} & \\multicolumn{6}{c}{$H_T$}   \\\\ \n");
    printf(" & $N^{mu}$  & $N^{W}$    &  $N^{top}$     & $N^{other}$  & $R^{MC}$  & data pred.  & MC pred. \\\\ \n ");
    printf(" $%.0f - %.0f$ & %.2f $\\pm$ %.2f & %.2f $\\pm$ %.2f & %.2f $\\pm$ %.2f & %.2f $\\pm$ %.2f & %.2f $\\pm$ %.2f & %.2f $\\pm$ %.2f &  \\\\", 
	   MT2CUT_LOW, MT2CUT_HIGH,
	   data, data_err,
	   w_est, w_est_err,
	   tt_bkg, tt_bkg_err,
	   mc_bkg, mc_bkg_err,
	   mc_factor, mc_factor_err,
	   z_est, z_est_err
	   );
    printf("   \\hline\\hline \n \\end{tabular} \n \\end{center} \n \\caption{Results for $H_T \\in () \\ \\mathrm{GeV}$.}\\label{table:} \n \\end{table}");
    cout << "\n\n" << endl;

    

  }
  fout->Write();
  fout->Close();

  cout << "Output file in " << outputfile << endl;

  return result;
}













//histograms
void ZInvEstFromW::createHistos(map<TString, TH1F*> *histos, vector<TString> samplesName ){

  TString hname="";
  hname="MT2_genCorr_Wtolnu_jidSel"; (*histos)[hname] = new TH1F(hname,";"+hname,500,0,5000);

  for( vector<TString>::const_iterator s_i=samplesName.begin(); s_i != samplesName.end(); s_i++){
    TString sname = *s_i;

    hname="MT2_"+sname+"_bSel"; (*histos)[hname] = new TH1F(hname,";"+hname,200,0,2000);
    hname="MT2_"+sname+"_leptSel"; (*histos)[hname] = new TH1F(hname,";"+hname,200,0,2000);
    hname="HT_"+sname+"_bSel"; (*histos)[hname] = new TH1F(hname,";"+hname,200,0,2000);
    hname="HT_"+sname+"_leptSel"; (*histos)[hname] = new TH1F(hname,";"+hname,200,0,2000);
    hname="MET_"+sname+"_bSel"; (*histos)[hname] = new TH1F(hname,";"+hname,200,0,2000);
    hname="MET_"+sname+"_leptSel"; (*histos)[hname] = new TH1F(hname,";"+hname,200,0,2000);
    hname="Jets_N_"+sname+"_bSel"; (*histos)[hname] = new TH1F(hname,";"+hname,20,0,20);
    hname="Jets_N_"+sname+"_leptSel"; (*histos)[hname] = new TH1F(hname,";"+hname,20,0,20);
    hname="MinDPhi_"+sname+"_bSel"; (*histos)[hname] = new TH1F(hname,";"+hname,200,0,3.2);
    hname="MinDPhi_"+sname+"_leptSel"; (*histos)[hname] = new TH1F(hname,";"+hname,200,0,3.2);
    hname="MinDPhiIndex_"+sname+"_bSel"; (*histos)[hname] = new TH1F(hname,";"+hname,20,0,20);
    hname="MinDPhiIndex_"+sname+"_leptSel"; (*histos)[hname] = new TH1F(hname,";"+hname,20,0,20);

    hname="MinDPhiPt_"+sname+"_leptSel"; (*histos)[hname] = new TH1F(hname,";"+hname,200,0,2000);
    hname="MinDPhiEta_"+sname+"_leptSel"; (*histos)[hname] = new TH1F(hname,";"+hname,100,-5,5);
    hname="MinDPhiPhi_"+sname+"_leptSel"; (*histos)[hname] = new TH1F(hname,";"+hname,100,-3.2,3.2);
    hname="MinDPhiMuDR_"+sname+"_leptSel"; (*histos)[hname] = new TH1F(hname,";"+hname,100,0,6);


    hname="Jets_pt_"+sname+"_bSel"; (*histos)[hname] = new TH1F(hname,";"+hname,200,0,2000);
    hname="Jets_pt_"+sname+"_leptSel"; (*histos)[hname] = new TH1F(hname,";"+hname,200,0,2000);
    hname="caloHT_"+sname+"_bSel"; (*histos)[hname] = new TH1F(hname,";"+hname,200,0,2000);
    hname="caloHT_"+sname+"_leptSel"; (*histos)[hname] = new TH1F(hname,";"+hname,200,0,2000);
    hname="VSPt_"+sname+"_bSel"; (*histos)[hname] = new TH1F(hname,";"+hname,60,0,300);
    hname="VSPt_"+sname+"_leptSel"; (*histos)[hname] = new TH1F(hname,";"+hname,60,0,300);


    //Vectors Boson info
    hname="VectBoson_pt_"+sname+"_preSel"; (*histos)[hname] = new TH1F(hname,";"+hname,500,0,5000);
    hname="VectBoson_pt_"+sname+"_leptVetoSel"; (*histos)[hname] = new TH1F(hname,";"+hname,500,0,5000);
    hname="VectBoson_pt_"+sname+"_bSel"; (*histos)[hname] = new TH1F(hname,";"+hname,500,0,5000);
    hname="VectBoson_pt_"+sname+"_jidSel"; (*histos)[hname] = new TH1F(hname,";"+hname,500,0,5000);
    hname="MT_"+sname+"_leptSel"; (*histos)[hname] = new TH1F(hname,";"+hname,200,0,1000);
    hname="MT_"+sname+"_bSel"; (*histos)[hname] = new TH1F(hname,";"+hname,200,0,1000);

    
    //Acceptance
    hname="MT2_LeptAcc_"+sname+"_jidSel"; (*histos)[hname] = new TH1F(hname,";"+hname,500,0,5000);
    
    //eff
    hname="Lept_eff_TP_"+sname;           (*histos)[hname] = new TH1F(hname,";"+hname,100,0,1);
    hname="B_eff_"+sname;           (*histos)[hname] = new TH1F(hname,";"+hname,100,0,1);
  
    //Leptons
    hname="Lept_pt_leptSel_"+sname; (*histos)[hname] = new TH1F(hname,";"+hname,200,0,2000);
    hname="Lept_eta_leptSel_"+sname; (*histos)[hname] = new TH1F(hname,";"+hname,100,-3,3);
    hname="Lept_pt_bSel_"+sname; (*histos)[hname] = new TH1F(hname,";"+hname,200,0,2000);
    hname="Lept_eta_bSel_"+sname; (*histos)[hname] = new TH1F(hname,";"+hname,100,-3,3);

    //bjets
    hname="BJets_N_"+sname+"_leptSel";  (*histos)[hname] = new TH1F(hname,";"+hname,10,0,10);
    hname="BJets_pt_"+sname+"_leptSel";  (*histos)[hname] = new TH1F(hname,";"+hname,200,0,2000);

    hname="NVertices_"+sname+"_leptSel"; (*histos)[hname] = new TH1F(hname,";"+hname,50,0,50);
    hname="NVertices_"+sname+"_bSel"; (*histos)[hname] = new TH1F(hname,";"+hname,50,0,50);

    hname="PUWeight_"+sname+"_leptSel"; (*histos)[hname] = new TH1F(hname,";"+hname,100,0,5);
    hname="PUWeight_"+sname+"_bSel"; (*histos)[hname] = new TH1F(hname,";"+hname,100,0,5);


}

  for(map<TString,TH1F*>::iterator h= (*histos).begin(); h!=(*histos).end();h++)
    h->second->Sumw2();

}



//print eff
void ZInvEstFromW::printBLeptEff(TString ssample, map<TString, TH1F*> histos, map<TString, float> varFloat){
  //cout << "--- " << ssample << endl;
  //cout << "T&P efficiency: "<<  histos["Lept_eff_TP_"+ssample]->GetMean() << " +- " <<        histos["Lept_eff_TP_"+ssample]->GetMeanError()
  //     << " +- " << sqrt(       varFloat["Lept_eff_TP_err_"+ssample] /(float) histos["Lept_eff_TP_"+ssample]->GetEntries() ) << endl;
  // cout <<  "BTag efficiency: "<<        histos["B_eff_"+ssample]->GetMean() << " +- " <<      histos["B_eff_"+ssample]->GetMeanError()
  //     << " +- " << sqrt(       varFloat["B_eff_err_"+ssample] /(float)       histos["B_eff_"+ssample]->GetEntries() ) << endl;

  printf("|*T&P eff for %s* | %.3f +- %.3f +- %.3f |\n", ssample.Data(),histos["Lept_eff_TP_"+ssample]->GetMean()  ,
	 histos["Lept_eff_TP_"+ssample]->GetMeanError(),sqrt(       varFloat["Lept_eff_TP_err_"+ssample] /(float) histos["Lept_eff_TP_"+ssample]->GetEntries() ) );


}




//____________________________________________________________________________
void ZInvEstFromW::loadSamples(const char* filename){
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
