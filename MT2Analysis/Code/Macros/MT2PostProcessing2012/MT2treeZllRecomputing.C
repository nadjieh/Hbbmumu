/*****************************************************************
* root macro to skim MT2trees                                    *
* -> sample and path to shlib given as arguments                 *
*                                                                *
* Pascal Nef                             September 18th, 2011    *  
*****************************************************************/
void MT2treeZllRecomputing(string sample, string shlib, string prefix) {
  	gSystem->Load("libPhysics");
  	gSystem->Load(shlib.c_str());

	string LABEL  = "";
	string file   = sample;
	TString outname = sample;
	int pos = outname.Last('/');
	outname.Remove(0, pos+1);
	string dummy = prefix.c_str();
	dummy += "/";
	dummy += outname.Data();
	outname = dummy;

	cout << "recomputing file " << sample << endl;
	cout << "outfile:         " << outname << endl;
	// files ---------------------------------------
	TFile *_file0 = TFile::Open( (file).c_str());
	TTree * inTree = (TTree*) _file0->Get("MassTree");
  
	MT2tree *mt2 = new MT2tree();
  	inTree->SetBranchAddress("MT2tree",&mt2);

	TFile *outfile  = new TFile(outname.Data(), "RECREATE");    outfile->cd();
  	TTree *outtree  = new TTree("MassTree","MassTree");
  	outtree->SetDirectory(outfile);
  	MT2tree *outMt2 = new MT2tree();
  	outtree->Branch("MT2tree","MT2tree",&outMt2);
  
	for ( Int_t jentry = 0; jentry < inTree->GetEntries(); jentry++) {
		inTree->GetEntry(jentry);
		if (jentry%10000 == 0) cout << jentry << " events" << endl;

		outMt2 = mt2;
		outMt2->ZllRecalculate();
//		if(outMt2->misc.MT2 <0  )                    continue;
//		if(outMt2->misc.HT  <400)                    continue;
		if(outMt2->RecoOSDiLeptPt(20,2.4,71,111)<150)continue;
		outtree->Fill();
	}
	cout << "--- End of Loop ---   " << endl;
	cout << "--- Writing PU histos " << endl;
	TH1F*   h_PUWeights = (TH1F*) _file0->Get("h_PUWeights");
	TH1F*   h_Events    = (TH1F*) _file0->Get("h_Events");
	outfile->cd();
	outtree->Write();
	if(h_PUWeights!=0){
		cout << "....writing TH1F h_PUWeights" << endl;
		h_PUWeights->Write();
	}
	if(h_Events!=0){
		cout << "....writing TH1F h_Events" << endl;
		h_Events->Write();
	}
	outfile->Close();
	cout << "--- NTuple saved : " << outfile->GetName() << endl;
}
