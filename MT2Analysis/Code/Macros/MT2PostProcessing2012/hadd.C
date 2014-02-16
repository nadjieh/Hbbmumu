//macro to add histogram files
//NOTE: This macro is kept for back compatibility only.
//Use instead the executable $ROOTSYS/bin/hadd
//
//This macro will add histograms from a list of root files and write them
//to a target root file. The target file is newly created and must not be
//identical to one of the source files.
//
//Author: Sven A. Schmidt, sven.schmidt@cern.ch
//Date:   13.2.2001

//This code is based on the hadd.C example by Rene Brun and Dirk Geppert,
//which had a problem with directories more than one level deep.
//(see macro hadd_old.C for this previous implementation).
//
//The macro from Sven has been enhanced by
//   Anne-Sylvie Nicollerat <Anne-Sylvie.Nicollerat@cern.ch>
// to automatically add Trees (via a chain of trees).
//
//To use this macro, modify the file names in function hadd.
//
//NB: This macro is provided as a tutorial.
//    Use $ROOTSYS/bin/hadd to merge many histogram files




#include <string.h>
#include "TChain.h"
#include "TFile.h"
#include "TH1.h"
#include "TTree.h"
#include "TKey.h"
#include "Riostream.h"

TList *FileList;
TFile *Target;

void MergeRootfile( TDirectory *target, TList *sourcelist );


void hadd() {
   // in an interactive ROOT session, edit the file names
   // Target and FileList, then
   // root > .L hadd.C
   // root > hadd()

  string DCAPDIR = "dcap://t3se01.psi.ch:22125/pnfs/psi.ch/cms/trivcat/store/user/";
  string basedir = "/haweber/SUSY/MassTrees/MT2_V02-01-02/20120718_MC_2leptons_JEC_52V9_test/";
  string filedir = "8TeV-WWGJets-FastSim-525-Summer12-v3-StoreResults-InTimePU-START52-V9-v3";
  string  totdir = DCAPDIR + basedir + filedir;
  string  MT2dir = "/shome/haweber/MT2Analysis/MT2trees/MT2_V02-01-02/20120718_MC_2leptons_JEC_52V9_test/NoSkim/";

   Target = TFile::Open( (MT2dir+filedir+(string)".root").c_str(), "RECREATE" );

   FileList = new TList();
   FileList->Add( TFile::Open((totdir + (string)"/output_0.root" ).c_str()) );
   FileList->Add( TFile::Open((totdir + (string)"/output_1.root" ).c_str()) );
   FileList->Add( TFile::Open((totdir + (string)"/output_2.root" ).c_str()) );
   FileList->Add( TFile::Open((totdir + (string)"/output_3.root" ).c_str()) );
   FileList->Add( TFile::Open((totdir + (string)"/output_4.root" ).c_str()) );
   FileList->Add( TFile::Open((totdir + (string)"/output_5.root" ).c_str()) );
   FileList->Add( TFile::Open((totdir + (string)"/output_6.root" ).c_str()) );
   FileList->Add( TFile::Open((totdir + (string)"/output_7.root" ).c_str()) );
   FileList->Add( TFile::Open((totdir + (string)"/output_8.root" ).c_str()) );
   FileList->Add( TFile::Open((totdir + (string)"/output_9.root" ).c_str()) );
   FileList->Add( TFile::Open((totdir + (string)"/output_10.root").c_str()) );
   FileList->Add( TFile::Open((totdir + (string)"/output_11.root").c_str()) );
   FileList->Add( TFile::Open((totdir + (string)"/output_12.root").c_str()) );
   FileList->Add( TFile::Open((totdir + (string)"/output_13.root").c_str()) );
   FileList->Add( TFile::Open((totdir + (string)"/output_14.root").c_str()) );
   FileList->Add( TFile::Open((totdir + (string)"/output_15.root").c_str()) );
   FileList->Add( TFile::Open((totdir + (string)"/output_16.root").c_str()) );
   FileList->Add( TFile::Open((totdir + (string)"/output_17.root").c_str()) );
   FileList->Add( TFile::Open((totdir + (string)"/output_18.root").c_str()) );
   FileList->Add( TFile::Open((totdir + (string)"/output_19.root").c_str()) );

   MergeRootfile( Target, FileList );

}

// Merge all files from sourcelist into the target directory.
// The directory level (depth) is determined by the target directory's
// current level
void MergeRootfile( TDirectory *target, TList *sourcelist ) {

  //  cout << "Target path: " << target->GetPath() << endl;
  TString path( (char*)strstr( target->GetPath(), ":" ) );
  path.Remove( 0, 2 );

  TFile *first_source = (TFile*)sourcelist->First();
  first_source->cd( path );
  TDirectory *current_sourcedir = gDirectory;

  // loop over all keys in this directory
  TIter nextkey( current_sourcedir->GetListOfKeys() );
  TKey *key;
  while (key = (TKey*)nextkey() ) {

    // read object from first source file
    first_source->cd( path );
    TObject *obj = key->ReadObj();

    if ( obj->IsA()->InheritsFrom( "TH1" ) ) {
      // descendant of TH1 -> merge it

      //      cout << "Merging histogram " << obj->GetName() << endl;
      TH1 *h1 = (TH1*)obj;

      // loop over all source files and add the content of the
      // correspondant histogram to the one pointed to by "h1"
      TFile *nextsource = (TFile*)sourcelist->After( first_source );
      while ( nextsource ) {
	
	// make sure we are at the correct directory level by cd'ing to path
	nextsource->cd( path );
	TH1 *h2 = (TH1*)gDirectory->Get( h1->GetName() );
	if ( h2 ) {
	  h1->Add( h2 );
	  delete h2; // don't know if this is necessary, i.e. if 
                     // h2 is created by the call to gDirectory above.
	}

	nextsource = (TFile*)sourcelist->After( nextsource );
      }

    } else if ( obj->IsA()->InheritsFrom( "TDirectory" ) ) {
      // it's a subdirectory

      cout << "Found subdirectory " << obj->GetName() << endl;

      // create a new subdir of same name and title in the target file
      target->cd();
      TDirectory *newdir = target->mkdir( obj->GetName(), obj->GetTitle() );

      // newdir is now the starting point of another round of merging
      // newdir still knows its depth within the target file via
      // GetPath(), so we can still figure out where we are in the recursion
      MergeRootfile( newdir, sourcelist );

    } else {
      // object is of no type that we know or can handle
      cout << "Unknown object type, name: " 
	   << obj->GetName() << " title: " << obj->GetTitle() << endl;
    }

    // now write the merged histogram (which is "in" obj) to the target file
    // note that this will just store obj in the current directory level,
    // which is not persistent until the complete directory itself is stored
    // by "target->Write()" below
    if ( obj ) {
      target->cd();
      obj->Write( key->GetName() );
    }

  } // while ( ( TKey *key = (TKey*)nextkey() ) )

  // save modifications to target file
  target->Write();

}