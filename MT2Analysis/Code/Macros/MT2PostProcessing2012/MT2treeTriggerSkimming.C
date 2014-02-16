/*****************************************************************
* root macro to skim MT2trees                                    *
* -> sample and path to shlib given as arguments                 *
*                                                                *
* Mario Masciovecchio                      February 7th, 2013    *
* from Pascal Nef's MT2treeSkimming.C                            * 
*****************************************************************/
void MT2treeTriggerSkimming(string sample, string trigger, string shlib, string prefix) {
  	gSystem->Load("libPhysics");
	gSystem->Load("libFWCoreFWLite.so");
        gSystem->Load("libDataFormatsFWLite.so");
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
	cutStream       << " && " 	  
/*
//	  << "misc.Event%3==0";
//	  << "misc.MT2>=50"                                      << "&&" 
	  << "misc.MET>=30"                                      << "&&"
	  << "misc.Jet0Pass ==1"                                 << "&&"
	  << "misc.Jet1Pass ==1"                                 << "&&"
//	  << "misc.SecondJPt  >100"                              << "&&"
	  << "misc.PassJetID ==1"                                << "&&"
	  << "misc.Vectorsumpt < 70"                             << "&&"
	  << "(NJetsIDLoose40 +NEles +NMuons)>=2"                << "&&"
//	  << "((misc.MinMetJetDPhi >0.3&&NBJets==0)||NBJets>=1)" << "&&"
//	  << "misc.MinMetJetDPhi4 >0.3"                          << "&&"
// Lepton Veto
//	  << "(NEles+NMuons)==0"                                 << "&&"
//	  << "(NEles==0 || ele[0].lv.Pt()<10)"                   << "&&"
//	  << "(NMuons==0 || muo[0].lv.Pt()<10)"                  << "&&"
// Lepton Skim
//	  << "(ele[0].lv.Pt()>10 || muo[0].lv.Pt()>10)"          << "&&"
// LowMT2 ----------------------------
//	  << "misc.LeadingJPt >150"                              << "&&"
//	  << "NBJets >0"                                         << "&&"
//	  << "NJetsIDLoose >=4"                                  << "&&"
// -----------------------------------
// Signal Veto -----------------------
//	   << "misc.ProcessID!=10"                               << "&&"
// MT2 signal regions ----------------
//	  << "misc.HT>=750&&misc.MET>=30"                        << "&&"//HTskim
//	  << "misc.HT>=450&&misc.HT<=750&&misc.MET>=200"         << "&&"//METskim
	  << "((misc.HT>=450&&misc.HT<=750&&misc.MET>=200)||(misc.HT>=750&&misc.MET>=30))" << "&&"//HT or MET skim
// Stop Skim -------------------------
//	  << "NJetsIDLoose40 >=2"                                << "&&"
//	  << "GetNjets(60,2.4,1) >=4"                            << "&&"
//	  << "(NEles+NMuons) >=2"                                << "&&"
//	  << "NBJetsCSVM     >=1"                                << "&&"
//	  << "NJetsIDLoose40 >=3"                                << "&&"
// Stop Selection -------------------
//	   << "(misc.ProcessID!=10 || Susy.MassLSP>=50)"         << "&&"
//	   << "(misc.ProcessID!=10 || Susy.MassChi<=750)"        << "&&"
//	   << "(misc.ProcessID!=10 || Susy.MassLSP==125)"        << "&&"
//	   << "(misc.ProcessID!=10 || Susy.MassChi==400)"        << "&&"
// Photons
//	  << "(GenDiLeptPt(0,10,0,1000,false)>=100||GenPhoton[0].Pt()>=100)"  << "&&"
//	  << "(RecoOSDiLeptPt(10,10,0,10000)>=100 ||photon[0].lv.Pt()>=100)"  << "&&"
//	  << "NPhotons >0"  
// Noise -- first 6 are official filters
//	  << "(misc.isFastSim || misc.HBHENoiseFlag == 0)"       << "&&"
//	  << "(misc.isFastSim || misc.CSCTightHaloIDFlag == 0)"  << "&&"
//	  << "(misc.isFastSim || misc.hcalLaserEventFlag == 0)"  << "&&"
*/
	<< "(misc.ProcessID==10 || misc.HBHENoiseFlag == 0)"       	<< "&&"
	<< "(misc.ProcessID==10 || misc.CSCTightHaloIDFlag == 0)"  	<< "&&"
	<< "(misc.ProcessID==10 || misc.hcalLaserEventFlag == 0)"  	<< "&&"
	<< "(misc.ProcessID==10 || misc.trackingFailureFlag == 0)"      << "&&"
	<< "(misc.ProcessID==10 || misc.eeBadScFlag == 0)"              << "&&"
	<< "(misc.ProcessID==10 || misc.EcalDeadCellTriggerPrimitiveFlag == 0)"  << "&&" 
// 	<< "(!misc.isData || trigger.HLT_SingleMu==1)"			<< "&&"
// 
	<< "(misc.ProcessID==10 || misc.CrazyHCAL==0)"			<< "&&"
// 	<< "NMuons >=1"							<< "&&"
// 	<< "(Sum$(muo.lv.Pt()>18&&abs(muo.lv.Eta())<1.6)) >=1"		<< "&&"
	<<"(misc.ProcessID==10 || misc.HT>450.)";

// 	//JZB
// 	  << "NJetsIDLoose40 >=3"                                << "&&"
// 	  << "(NEles+NMuons) >=2"                                << "&&"
// 	  << "rawpfmet[0].Pt()>=95"				 << "&&"
// 	  << "(Sum$(ele.lv.Pt()>18&&abs(ele.lv.Eta())<1.6)+Sum$(muo.lv.Pt()>18&&abs(muo.lv.Eta())<1.6)) >=2";

	string basecut = cutStream.str();

	string  SEL= "("+trigger+basecut+")";
// &&"+basecut+")";
	cout << "Skimming with sel: " << SEL << endl;
// 	TString cuts_log = basecut.ReplaceAll("&&", "\n");
// 	f_cuts << "Cuts applied: " <<  cuts_log << endl;
	//  --------------------------------------------


	// files ---------------------------------------
	f_log << "Skimming file: " << file;
	TFile *_file0 = TFile::Open( (file).c_str()); 
	TTree * t = (TTree*) _file0->Get("MassTree");
	t->SetMaxTreeSize(19000000000);
	TFile*out = TFile::Open( (outfile).c_str(),"RECREATE");
	TTree *tc = t->CopyTree(SEL.c_str());   
	int nentries = tc->GetEntries();
	f_log << " -> skimmed tree has " << nentries << " entries." <<endl;
	f_log.close();
	f_cuts.close();
	out->Write();
	out->Close();
	_file0->Close();
	cout << "Result file: " << outfile << endl;
	// -------------------------------------------------

}
