/*****************************************************************
* root macro to skim MT2trees                                    *
* -> sample and path to shlib given as arguments                 *
*                                                                *
* Pascal Nef                             September 18th, 2011    *  
*****************************************************************/
void MT2treeSkimming(string sample, string shlib) {
  	gSystem->Load("libPhysics");
  	gSystem->Load(shlib.c_str());

	string LABEL  = "";
	string file   = sample;
	string outfile = "skimmed/"+sample;

	// log file 
	TString log=sample+".skim.log";
	ofstream f_log;
	f_log.open(log.Data());
	f_log << "Skimming file: " << sample << " with cuts: " << endl;


	// cuts --------------------------------------------
  	std::ostringstream cutStream;
	cutStream       << " " 	  
//	  << "misc.MT2 >= 400"                         << "&&" 
	  << "misc.MET>=30"                            << "&&"
	  << "misc.HT > 400 "                          << "&&"
	  << "misc.caloHT50_ID >600 "                  << "&&"
	  << "misc.Jet0Pass ==1"                       << "&&"
	  << "misc.Jet1Pass ==1"                       << "&&"
	  << "misc.SecondJPt  >100"                    << "&&"
	  << "misc.PassJetID ==1"                      << "&&"
	  << "misc.Vectorsumpt < 70"                   << "&&"
	  << "misc.MinMetJetDPhi >0.3"                 << "&&"
	  << "misc.HBHENoiseFlag == 0"                 << "&&"
	  << "misc.CrazyHCAL==0"                       << "&&"
	  << "(NEles ==0 && NMuons==0)"                << "&&"
// LowMT2 ----------------------------
//	  << "misc.LeadingJPt >150"                    << "&&"
//	  << "NBJets >0"                               << "&&"
//	  << "NJetsIDLoose >=4"                        << "&&"
// -----------------------------------
	  << "NJetsIDLoose >=3";
	
	TString basecut = cutStream.str();
	string  SEL= "("+basecut+")";
	cout << "Skimming with sel: " << SEL << endl;
	TString cuts_log = basecut.ReplaceAll("&&", "\n");
	f_log << cuts_log << endl;
	//  --------------------------------------------


	// files ---------------------------------------
	TFile *_file0 = TFile::Open( (file).c_str()); 
	TTree * t = (TTree*) _file0->Get("MassTree");
	t->SetMaxTreeSize(19000000000);
	TFile*out = TFile::Open( (outfile).c_str(),"RECREATE");
	TTree *tc = t->CopyTree(SEL.c_str());   
	int nentries = tc->GetEntries();
	f_log << "skimmed tree has " << nentries << " entries." <<endl;
	f_log.close();
	out->Write();
	out->Close();
	_file0->Close();
	cout << "Result file: " << outfile << endl;
	// -------------------------------------------------

}
