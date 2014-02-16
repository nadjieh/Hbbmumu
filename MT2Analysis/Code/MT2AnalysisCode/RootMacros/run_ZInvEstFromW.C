{

  TString outputdir = "../Histos/";
  //  TString outname = "ZInv-MT2_0_9999-HT_800-1Mu_pt10";
  bool EFF_CORR=false;
  float MT2_low =  -1;
  float MT2_high = 9999;
  int NEles = 0;
  int NMuons = 1;
  bool IS_MC = false;

  ostringstream MT2_low_s, MT2_high_s;
  MT2_low_s << MT2_low;
  MT2_high_s << MT2_high;
  TString outname = "ZInv-MT2_"+MT2_low_s.str()+"_"+MT2_high_s.str()+"-Jet40-HT_700-1Mu_pt10_4450pb";
  if(IS_MC) outname += "-MC_Closure";

  TString samples   = "../samples/samples_highMT2_lept_4600.dat";
  int verbose = 0;

  gSystem->CompileMacro("../MT2Code/src/ZInvEstFromW.cc", "k");
  gROOT->ProcessLine(".x SetStyle_PRD.C");

  ZInvEstFromW *tA = new ZInvEstFromW(outputdir);
  


  tA->setVerbose(verbose);
  tA->init(samples);
  
  std::ostringstream cutStream;
  cutStream << " " 
    //    	    << "Sum$(jet[].lv.Pt()>40)>2"                            << "&&"
    //<< "misc.MET>=30"                            << "&&"
    //<< " misc.MT2 > 100"                           << "&&"
	    << "jet[2].lv.Pt()>40 && "
	    << "misc.HT > 700 "                          << "&&"
    //<< "misc.caloHT50_ID >800 "                  << "&&"
	    << "misc.Jet0Pass ==1"                       << "&&"
	    << "misc.Jet1Pass ==1"                       << "&&"
    //<< "misc.PassJetID ==1"                      << "&&"
    //	    << "misc.Vectorsumpt < 70"                   << "&&"
    //    << "misc.MinMetJetDPhi >0.3"                 << "&&"
	    << "misc.HBHENoiseFlag == 0"                 << "&&"
	    << "misc.CrazyHCAL==0 && "
 //    << "misc.LeadingJPt >150"                    << "&&"
	    << " misc.CSCTightHaloID==0 && " 
	    << "misc.SecondJPt  >100";//                    << "&&"
  //<<"NJetsIDLoose>3 " << " && "
  //<< "NBJets >0"                               << "&&"
  //<< "misc.PassJetID ==1"                      << "&&"
  //<< "misc.Vectorsumpt < 70"                   << "&&"
  //  << "misc.MinMetJetDPhi >0.3"                 << "&&"
  //  << "misc.HBHENoiseFlag == 0"                 << "&&"
  // "misc.MT2 > 450"                           << "&&"
  //	  << "misc.MT2 > 200 && misc.MT2<400"          << "&&"
  //	  << "misc.MT2 > 100 && misc.MT2<150"          << "&&"
  //	  << "misc.CrazyHCAL==0";
  
  std::ostringstream triggerStream;
  triggerStream << "( "
		<< "(trigger.HLT_HT440_v2 ==1 && misc.Run<161216)" << "||"
		<< "(trigger.HLT_HT450_v2 ==1 && (misc.Run>=161216 && misc.Run< 163269))" << "||"
		<< "(trigger.HLT_HT500_v3 ==1 && (misc.Run>=163269 && misc.Run<=163869))" << "||"
		<< "(trigger.HLT_HT500_v4 ==1 && (misc.Run>=165088 && misc.Run< 165970))" << "||"
		<< "(trigger.HLT_HT550_v5 ==1 && (misc.Run>=165970 && misc.Run<=167043 && misc.Run!=166346))" << "||"
		<< "(trigger.HLT_HT550_v6 ==1 && (misc.Run==166346))" << "||"
		<< "(trigger.HLT_HT550_v7 ==1 && (misc.Run>=167078 && misc.Run< 170249))" << "||"
		<< "(trigger.HLT_HT550_v8 ==1 && (misc.Run>=170249 && misc.Run< 173236))" << "||"
    		<< "(trigger.HLT_HT600_v1 ==1 && (misc.Run>=173236 && misc.Run< 178420))" << ""
    //		<< " || misc.Run< 178420 "
    		<< "|| (trigger.HLT_HT650_v4 ==1 && (misc.Run>=178420 && misc.Run< 179959))" << ""
    		<< "|| (trigger.HLT_PFHT650_v1==1 && misc.Run>=179959)" 
		<< " )";
  TString trigger = triggerStream.str().c_str();
  
  TString cuts = cutStream.str().c_str();
  
  //actual call
  float *ttbarRes = new float[10];
  ttbarRes = tA->Analysis(IS_MC, "histos_"+outname+"_1b.root", 3, NEles, NMuons, MT2_low, MT2_high, trigger, cuts, true, EFF_CORR);
  tA->Analysis(IS_MC, "histos_"+outname+"_0b.root", 3, NEles, NMuons, MT2_low, MT2_high, trigger, cuts, false, EFF_CORR, ttbarRes);

}
