{

  gROOT->ProcessLine(".x rootlogon.C");
  gSystem->CompileMacro("../MT2Code/src/ScanAnalysis.cc", "k");


  TString outputdir = ".";
  
  TString outname = "histo_mSugra_test";
  
  TString samples   = "../samples/samples_scan.dat";
  int verbose = 0;
  
  ScanAnalysis *tA = new ScanAnalysis(outputdir);
    
  tA->setVerbose(verbose);
  tA->init(samples, "scale_xsection_nlo1.0_m0_m12_10_0_1v1.txt");
  
  std::ostringstream cutStream;
  cutStream << "   misc.Run < 179959 && " 
      //<< "Sum$(jet[].lv.Pt()>40)>2"                            << "&&"
      //<< "misc.MET>=240 && misc.MET<=320"                            << "&&"
      //<< " misc.MT2 > 100"                           << "&&"
      	      << "NJetsIDLoose40 > 2 && "
	    << "misc.HT > 0 "                          << "&&"
      //<< "misc.HT > 950 "                          << "&&"
      //	      << "misc.HT > 750 && misc.HT<=950 "                          << "&&"
      //<< "misc.caloHT50_ID >750 "                  << "&&"
	      << "misc.Jet0Pass ==1"                       << "&&"
	      << "misc.Jet1Pass ==1"                       << "&&"
	      << "misc.PassJetID ==1"                      << "&&"
      //	    << "misc.Vectorsumpt < 70"                   << "&&"
      //<< "misc.MinMetJetDPhi <=0.3"                 << "&&"
    //<< "misc.HBHENoiseFlagIso == 0"                 << "&&"
      
    //      << "misc.CrazyHCAL==0 && "
      //    << "misc.LeadingJPt >150"                    << "&&"
    //      << " misc.CSCTightHaloID==0 && " 
	      << "misc.SecondJPt  >100";//                    << "&&"
    //<< "NBJets >0"                               << "&&"
    //	  << "misc.MT2 > 200 && misc.MT2<400"          << "&&"
    //	  << "misc.MT2 > 100 && misc.MT2<150"          << "&&"
    //	  << "misc.CrazyHCAL==0";
  
  TString cuts = cutStream.str().c_str();
  

  tA->Analysis(outname+".root",cuts);
  
  
  
}
