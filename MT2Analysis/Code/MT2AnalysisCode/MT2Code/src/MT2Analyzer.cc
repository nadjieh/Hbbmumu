#include "MT2Analyzer.hh"

#include "base/TreeAnalyzerBase.hh"
#include "base/TreeReader.hh"
#include "MT2Analysis.hh"

using namespace std;

MT2Analyzer::MT2Analyzer(std::vector<std::string>& fileList) 
	: TreeAnalyzerBase(fileList) {
	fMT2Analysis             = new MT2Analysis(fTR);
	Util::SetStyle();
	removePhoton =false;
	removeZll    =false;
	fID          =-1;  //default process ID
	fbtagFileName="";
	fhadtauFileName="";
	fType1MET    =false;
	fCHSJets     =false;
	fisFastSim   =false;
}

MT2Analyzer::~MT2Analyzer(){
	delete fMT2Analysis;
}

// Method for looping over the tree
void MT2Analyzer::Loop(){
	Long64_t nentries = fTR->GetEntries();
	cout << " total events in ntuples = " << fTR->GetEntries() << endl;
	
	if(fMaxEvents==-1)nentries=fTR->GetEntries();
	if(fMaxEvents > 0){
		nentries=min((Long64_t)fMaxEvents, fTR->GetEntries());
		cout << " only running on first " << nentries << " events" << endl;
	}
	if(nentries ==0) return; // fix in order to prevent assertion error if nentries==0

	// loop over all ntuple entries
    	Long64_t jentry=0;
    	for ( fTR->ToBegin(); !(fTR->AtEnd()) && (jentry<nentries || nentries<0); ++(*fTR) ) {
		PrintProgress(jentry++);
               if ( fCurRun != fTR->Run ) {
        		fCurRun = fTR->Run;
			fMT2Analysis        ->BeginRun(fCurRun);
                  	skipRun = false; // re-initialize
                  	if ( !CheckRun() ) skipRun = true;
			//if(skipRun) cout << "skip run " << fCurRun << endl;
                }
               // Check if new lumi is in JSON file
                if ( !skipRun && fCurLumi != fTR->LumiSection ) {
                  	fCurLumi = fTR->LumiSection;
                  	skipLumi = false; // Re-initialise
                  	if ( !CheckRunLumi() ) skipLumi = true;
			//if(skipLumi) cout << "skip lumi " << fCurLumi << " in run " << fTR->Run << endl;
                }
		if ( !(skipRun || skipLumi) ) {
			fMT2Analysis        ->Analyze();
			double PUWeight = 0;
			//PU mean weight
			if (fPu=="MC2012"){
			  	PUWeight  = fMT2Analysis->GetPUWeight(fTR->PUnumTrueInteractions);
			} else {
			  	PUWeight  = 1;
			}
			
			fMT2Analysis->fH_PUWeights->Fill( PUWeight );
			fMT2Analysis->fH_Events->Fill( 1. );

			if(isScan){
			  fMT2Analysis->fH2_mSugraEvents->Fill(fTR->M0, fTR->M12);
			  fMT2Analysis->fH2_SMSEvents->Fill(fTR->MassGlu, fTR->MassLSP);
			  if(fTR->process>0) fMT2Analysis->fH_mSugraSubProcEvents[fTR->process]->Fill(fTR->M0, fTR->M12);
			}

		}
	}//end loop
}

// Method called before starting the event loop
void MT2Analyzer::BeginJob(TString filename, TString setofcuts, bool isData, string data_PileUp, string mc_PileUp, string JEC){
	fMT2Analysis                    ->ReadCuts(setofcuts);
	fMT2Analysis                    ->SetType(isData);
	fMT2Analysis                    ->SetPUReweighting(fPu, data_PileUp, mc_PileUp);
	fMT2Analysis                    ->SetOutputDir(fOutputDir);
	fMT2Analysis                    ->fVerbose        = fVerbose;
	fMT2Analysis                    ->SetJEC(JEC);
        fMT2Analysis                    ->fRemovePhoton = removePhoton;
        fMT2Analysis                    ->fRemoveZll    = removeZll;
	fMT2Analysis                    ->SetProcessID(fID);
	fMT2Analysis                    ->SetBTagEfficiency(fbtagFileName);
	fMT2Analysis                    ->SetHadTauEfficiency(fhadtauFileName);
	fMT2Analysis             	->doPDF        = doPDF;
        fMT2Analysis             	->isScan        = isScan;
	fMT2Analysis                    ->SetType1MET(fType1MET);
	fMT2Analysis                    ->SetCHSJets(fCHSJets);
	fMT2Analysis                    ->SetFastSim(fisFastSim);


	fMT2Analysis                    ->Begin(filename);


	fMT2Analysis->fH_PUWeights = new TH1F("h_PUWeights",";PU weights",400,0,20);
	fMT2Analysis->fH_Events = new TH1F("h_Events",";Events",10,0,10);
	fMT2Analysis->fH2_mSugraEvents = new TH2F("h_mSugraEvents",";m_{0};m_{1/2}",600,0,3000, 200,0,1000);
	fMT2Analysis->fH2_SMSEvents = new TH2F("h_SMSEvents",";m_{0};m_{1/2}",500,0,2500, 500,0,2500);

	if(isScan){
	  for(int s=1;s<11; s++){
	    ostringstream ss;
	    ss << s; 
	    fMT2Analysis->fH_mSugraSubProcEvents[s] = new TH2F( ("h_mSugraSubProcEvents_"+ss.str()).c_str(),";m_{0};m_{1/2}",150,0,3000, 50,0,1000 );
	  }
	}


}

// Method called after finishing the event loop
void MT2Analyzer::EndJob(){
  fMT2Analysis         ->End();
  cout << " MT2Analyzer::End()                                             " << endl;
  
}
