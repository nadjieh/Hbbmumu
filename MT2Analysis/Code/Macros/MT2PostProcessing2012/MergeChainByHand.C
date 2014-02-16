#include <iostream>
#include "TChain.h"
#include "TTree.h"
#include "TFile.h"
#include "/shome/haweber/MT2Analysis_8TeV/Code/MT2AnalysisCode/MT2Code/include/MT2tree.hh"

using namespace std;

// MT2PostProcessing does not work for directories beginning with numbers as
// /pnfs/psi.ch/cms/trivcat/store/user//haweber/SUSY/MassTrees/MT2_V02-01-02/20120718_MC_2leptons_JEC_52V9_test/8TeV-WWGJets-FastSim-525-Summer12-v3-StoreResults-InTimePU-START52-V9-v3/
// do this by hand (but no skimming here)

void MergeChainByHand(){

  TChain * tx = new TChain("MassTree");

  string DCAPDIR = "dcap://t3se01.psi.ch:22125/pnfs/psi.ch/cms/trivcat/store/user/";
  string basedir = "/haweber/SUSY/MassTrees/MT2_V02-01-02/20120718_MC_2leptons_JEC_52V9_test/";
  string filedir = "8TeV-WWGJets-FastSim-525-Summer12-v3-StoreResults-InTimePU-START52-V9-v3";
  string  totdir = DCAPDIR + basedir + filedir;
  string  MT2dir = "/shome/haweber/MT2Analysis/MT2trees/MT2_V02-01-02/20120718_MC_2leptons_JEC_52V9_test/NoSkim/";

  Int_t n = 0; // need the slash"/" before output_xx.root
  n+=tx->Add((totdir + (string)"/output_0.root" ).c_str(),0);
  n+=tx->Add((totdir + (string)"/output_1.root" ).c_str(),0);
  n+=tx->Add((totdir + (string)"/output_2.root" ).c_str(),0);
  n+=tx->Add((totdir + (string)"/output_3.root" ).c_str(),0);
  n+=tx->Add((totdir + (string)"/output_4.root" ).c_str(),0);
  n+=tx->Add((totdir + (string)"/output_5.root" ).c_str(),0);
  n+=tx->Add((totdir + (string)"/output_6.root" ).c_str(),0);
  n+=tx->Add((totdir + (string)"/output_7.root" ).c_str(),0);
  n+=tx->Add((totdir + (string)"/output_8.root" ).c_str(),0);
  n+=tx->Add((totdir + (string)"/output_9.root" ).c_str(),0);
  n+=tx->Add((totdir + (string)"/output_10.root").c_str(),0);
  n+=tx->Add((totdir + (string)"/output_11.root").c_str(),0);
  n+=tx->Add((totdir + (string)"/output_12.root").c_str(),0);
  n+=tx->Add((totdir + (string)"/output_13.root").c_str(),0);
  n+=tx->Add((totdir + (string)"/output_14.root").c_str(),0);
  n+=tx->Add((totdir + (string)"/output_15.root").c_str(),0);
  n+=tx->Add((totdir + (string)"/output_16.root").c_str(),0);
  n+=tx->Add((totdir + (string)"/output_17.root").c_str(),0);
  n+=tx->Add((totdir + (string)"/output_18.root").c_str(),0);
  n+=tx->Add((totdir + (string)"/output_19.root").c_str(),0);

  cout << "#files = " << n << endl;

  tx->SetBranchStatus("*",1); //enable all branches
  //tx->SetBranchStatus("*",0); //disable all branches


   TFile *newfile = new TFile((MT2dir+filedir+(string)".root").c_str(),"recreate");
   TTree *newtree = tx->CloneTree(0);

   newtree->CopyEntries(tx);

   newtree->Print();
   newfile->Write();
   delete tx;
   delete newfile;
}
