/*****************************************************************
* root macro to skim MT2trees                                    *
* -> sample and path to shlib given as arguments                 *
*                                                                *
* Pascal Nef                             September 18th, 2011    *  
*****************************************************************/
void MT2treeSkimming(string sample, string shlib, string prefix) {
  	gSystem->Load("libPhysics");
  	gSystem->Load(shlib.c_str());

	string LABEL  = "";
	string file   = sample;
	string outfile = prefix+"/"+sample;

	// log file 
	TString log =sample+".skim.log";
	TString cuts="cuts.skim.log";
	string line="";
	ofstream f_log;
	ofstream f_cuts;
	f_log .open(log.Data(),ios::app);
	f_cuts.open(cuts.Data());
	if(!f_log.is_open()||!f_cuts.is_open()) {cout << "ERROR: cannot open file. " << endl; exit(-1);}

	// cuts --------------------------------------------
  	std::ostringstream cutStream;
	cutStream       << " " 	  
//	  << "misc.MT2 >=50"                                     << "&&" 
//	  << "misc.MET>=30"                                      << "&&"
//	  << "misc.HT > 750 "                                    << "&&"
//	  << "misc.Jet0Pass ==1"                                 << "&&"
//	  << "misc.Jet1Pass ==1"                                 << "&&"
//	  << "misc.SecondJPt  >100"                              << "&&"
//	  << "misc.PassJetID ==1"                                << "&&"
//	  << "misc.Vectorsumpt < 70"                             << "&&"
//	  << "((misc.MinMetJetDPhi >0.3&&NBJets==0)||NBJets>=1)" << "&&"
// Lepton Veto
//	  << "(NEles==0 || ele[0].lv.Pt()<10)"                   << "&&"
//	  << "(NMuons==0 || muo[0].lv.Pt()<10)"                  << "&&"
// Lepton Skim
//	  << "(ele[0].lv.Pt()>10 || muo[0].lv.Pt()>10)"          << "&&"
// LowMT2 ----------------------------
//	  << "misc.LeadingJPt >150"                              << "&&"
//	  << "NBJets >0"                                         << "&&"
//	  << "NJetsIDLoose >=4"                                  << "&&"
// -----------------------------------
//	  << "NJetsIDLoose40 >=3"                                << "&&"
// Photons
	  << "(GenDiLeptPt(0,10,0,1000,false)>=100||GenPhoton[0].Pt()>=100)"  << "&&"
//	  << "(RecoOSDiLeptPt(10,10,0,10000)>=100 ||photon[0].lv.Pt()>=100)"  << "&&"
//	  << "NPhotons >0"  
// Noise
//	  << "misc.HBHENoiseFlag == 0"                           << "&&"
//	  << "misc.CSCTightHaloID==0"                            << "&&"
	  << "misc.CrazyHCAL==0";

	
	TString basecut = cutStream.str();
	string  SEL= "("+basecut+")";
	cout << "Skimming with sel: " << SEL << endl;
	TString cuts_log = basecut.ReplaceAll("&&", "\n");
	f_cuts << "Cuts applied: " <<  cuts_log << endl;
	//  --------------------------------------------


	// files ---------------------------------------
	f_log << "Skimming file: " << file;
	TFile *_file0 = TFile::Open( (file).c_str()); 
	TTree * t = (TTree*) _file0->Get("MassTree");
	TH1F*   h_PUWeights = (TH1F*) _file0->Get("h_PUWeights");
	TH1F*   h_Events    = (TH1F*) _file0->Get("h_Events");
	t->SetMaxTreeSize(19000000000);
	TFile*out = TFile::Open( (outfile).c_str(),"RECREATE");
	TTree *tc = t->CopyTree(SEL.c_str());   
	int nentries = tc->GetEntries();
	f_log << " -> skimmed tree has " << nentries << " entries." <<endl;
	f_log.close();
	f_cuts.close();
	out->Write();
	if(h_PUWeights!=0){
		cout << "writing TH1F h_PUWeights" << endl;
		h_PUWeights->Write();
	}
	if(h_Events!=0){
		cout << "writing TH1F h_Events" << endl;
		h_Events->Write();
	}
	out->Close();
	_file0->Close();
	cout << "Result file: " << outfile << endl;
	// -------------------------------------------------

}
