{
  //  TString outputdir = "/17TB/MT2Top/MassPlotssmallb/";
  TString outputdir = "../MassPlots/"; 
  //  TString outputdir = "/";
  TString samples   = "./samples/samplesMineSmearedJetMETNewBTag.dat";
  //TString samples   = "./samples/samplesMine1stTry.dat";

  int verbose =3;

  gSystem->CompileMacro("../MT2Code/src/MassPlotter.cc", "kf");//"k", "f"
  gROOT->ProcessLine(".x SetStyle_PRD.C");

  int gNMT2bins                   = 17;
  double  gMT2bins[gNMT2bins+1]   = {0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 120, 150, 200, 300, 400, 550, 800}; 	
  
//   int gNMT2Bbins                   = 15;
//   double  gMT2Bbins[gNMT2bins+1]   = {0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 120, 150, 200, 300, 500}; 	

  int gNMT2bins_l                   = 14;
  double  gMT2bins_l[gNMT2bins+1]   = {0, 10, 20, 30, 40, 50, 65, 80, 95, 115, 140, 180, 250, 350, 500}; 	

  int gNMT2Bbins = 21;
  double gMT2Bbins[gNMT2Bbins + 1] = {0, 10, 20, 30, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 100, 125, 150, 200, 250, 300, 700};


  MassPlotter *tA = new MassPlotter(outputdir, "MassPlots.root");
  tA->SetSave(false);
  tA->setVerbose(verbose);
  tA->init(samples);
  tA->SetIsPhoton(false);
  
  std::ostringstream cutStream;
  cutStream << " " 
	    << "(misc.ProcessID!=10 || (Susy.MassGlu >= 350 && Susy.MassGlu < 355"<< "&&"
	    << "Susy.MassLSP >= 50 && Susy.MassLSP < 55))"                     << "&&"
    //    << "misc.MT2>125"                                                << "&&"
     	    << "NBJetsCSVT> 0"                                                << "&&"
	    << "misc.MET>30"                                                << "&&"
   // 	    << "misc.HT>1500  "                                  << "&&"
    //          << "NTops > 0"                                                   << "&&"
       	    << "NTaus >0  && tau[0].Isolation3Hits > 1.5"                            << "&&"
    // //	  << "misc.MT2<200&&misc.MT2>150  "                                << "&&"
//	    << "misc.MT2 > 300 "                                               << "&&"
//	    << "misc.MT2<60 "                                               << "&&"
   	//    << "misc.MT2Top>130  "                                              << "&&"
    //	    << "NJetsIDLoose40  >=4"                                         << "&&"
    //	    << "misc.MT2MassiveTop < 1000"                                   << "&&"
 	    << "(NEles==0  || ele[0].lv.Pt()<5)"                             << "&&"
    	    << "(NMuons==0 || muo[0].lv.Pt()<5)"                             << "&&"
    //    << "(NEles!=0  || NMuons!=0)"                                    << "&&"
 	    << "misc.Jet0Pass ==1"                                           << "&&"
  	    << "misc.Jet1Pass ==1"                                           << "&&"
  	    << "misc.PassJetID ==1"                                          << "&&"
    // << "misc.SecondJPt >100"                                         << "&&"
   	    // << "misc.LeadingJPt >150"                                        << "&&"
      	    << "misc.Vectorsumpt < 70"                                       << "&&"
//	    << "(misc.MinMetJetDPhi >0.3||misc.MinMetJetDPhiIndex>3)"        << "&&"
 	    << "misc.MinMetJetDPhi4 >0.3"                                    << "&&"
//	    << "misc.MinMetJetDPhi4 <=0.2"                                   << "&&"
//   	    << "misc.MinMetBJetDPhi >0.7"                                    << "&&"
// 	    << "misc.TopDecayMode  == 16"                                    << "&&"
// 	    << "misc.Event > 7500000   "                                    << "&&"
// 	    << "misc.Event < 10000000   "                                    << "&&"
    //  << "misc.HBHENoiseFlagIso ==0"                               << "&&"
    // // << "misc.CSCTightHaloID==0"                                << "&&"
    //  << "misc.CrazyHCAL==0";
    //	    <<"1.0";

    //HT Trigger
	    // << "misc.HT>750";                             
 	   //  << "misc.HT>725";                             
	    
    //Or between Six and Quad
 	    // <<"((NJets  >=4 && jet[3].lv.Pt() >= 120.0) || (NJets  >=6 && jet[5].lv.Pt() >= 80.0))";
 	    // <<"((NJets  >=4 && jet[3].lv.Pt() >= 100.0) || (NJets  >=6 && jet[5].lv.Pt() >= 70.0))";
 //	    <<"((NJetsIDLoose40  >=5 && jet[4].lv.Pt() >= 95.0) || (NJetsIDLoose40  >=7 && jet[6].lv.Pt() >= 60.0))";
	    // <<"((NJetsIDLoose40  >=5 && jet[4].lv.Pt() >= 90.0) || (NJetsIDLoose40  >=7 && jet[6].lv.Pt() >= 55.0))";
  //New Or
//  	    <<"(((NJetsIDLoose40  >=4 && jet[3].lv.Pt() >= 110.0) || (NJetsIDLoose40  >=6 && jet[5].lv.Pt() >= 75.0))"<<"||"
//  	    <<"((NJetsIDLoose40  >=5 && jet[4].lv.Pt() >= 95.0) || (NJetsIDLoose40  >=7 && jet[6].lv.Pt() >= 60.0)))";
      <<"(((NJetsIDLoose40  >=4 && jet[3].lv.Pt() >= 100.0) || (NJetsIDLoose40  >=6 && jet[5].lv.Pt() >= 65.0))"<<"||"
  	    <<"(NJetsIDLoose40  >=5 && jet[4].lv.Pt() >= 85.0))";
//	    <<"((NJetsIDLoose40  >=1 && jet[0].lv.Pt() >= 20.0)))";


    //    <<"((NJetsIDLoose40  >=4 && jet[3].lv.Pt() >= 120.0))";
 //	    <<"((NJetsIDLoose40  >=6 && jet[5].lv.Pt() >= 80.0) && !(NJetsIDLoose40  >=4 && jet[3].lv.Pt() >= 120.0))";
 //	  << "NJets  >=4 && 0==0";

  std::ostringstream triggerStream;
  triggerStream << "( "
		<<" misc.ProcessID != 0 || ((misc.CrazyHCAL==0 && misc.NegativeJEC==0 " <<"&&"
		<<" misc.CSCTightHaloIDFlag==0 && misc.HBHENoiseFlag==0 " <<"&&"
		<<" misc.hcalLaserEventFlag==0 && misc.trackingFailureFlag==0 " <<"&&"
		<<" misc.eeBadScFlag==0 && misc.EcalDeadCellTriggerPrimitiveFlag==0 )" <<"&&("

 //    << "(trigger.HLT_PFHT650_v5 ==1)" << "||"
//     << "(trigger.HLT_PFHT650_v6 ==1)" << "||"
//     << "(trigger.HLT_PFHT650_v7 ==1)" << "||"
//     << "(trigger.HLT_PFHT650_v8 ==1)" << "||"
//     << "(trigger.HLT_PFHT650_v9 ==1)" << ")";

    // << "(trigger.HLT_QuadJet80_v1 ==1)" << "||"
    // << "(trigger.HLT_QuadJet80_v2 ==1)" << "||"
    // << "(trigger.HLT_QuadJet80_v3 ==1)" << ")";

    // << "(trigger.HLT_SixJet45_v1 ==1)" << "||"
    // << "(trigger.HLT_SixJet45_v2 ==1)" << "||"
    // << "(trigger.HLT_SixJet45_v3 ==1)" << ")";
    
    //Or between Six and Quad
		<< "(((trigger.HLT_QuadJet80_v1 ==1)" << "||"
		<< "(trigger.HLT_QuadJet80_v2 ==1)" << "||"
		<< "(trigger.HLT_QuadJet80_v3 ==1)" << "||"
		<< "(trigger.HLT_QuadJet80_v4 ==1)" << "||"
		<< "(trigger.HLT_QuadJet80_v6 ==1)" << ") &&"
		// << " NJets >=4 "            << "&&"
		// << " jet[3].lv.Pt() >= 120.0 "       << ") || "
		// << " NJets >=4 "            << "&&"
		// << " jet[3].lv.Pt() >= 100.0 "       << ") || "
		// << " NJetsIDLoose40  >=5 "            << "&&"
		// << " jet[4].lv.Pt() >= 95.0 "       << ") || "
		// << " NJetsIDLoose40  >=5 "            << "&&"
		// << " jet[4].lv.Pt() >= 90.0 "       << ") || "
    //New Or
// 		<< " (NJetsIDLoose40  >=4 "            << "&&"
// 		<< " jet[3].lv.Pt() >= 110.0 "       << ") || "
// 		<< " (NJetsIDLoose40  >=5 "            << "&&"
// 		<< " jet[4].lv.Pt() >= 95.0 "       << ")) || "
		<< " ((NJetsIDLoose40  >=4 "            << "&&"
		<< " jet[3].lv.Pt() >= 100.0 "       << ") || "
//              << " ((NJetsIDLoose40  >=5 "            << "&&"
//  		<< " jet[4].lv.Pt() >= 85.0 "       << ") || "
 		<< " (NJetsIDLoose40  >=5 "            << "&&"
 		<< " jet[4].lv.Pt() >= 85.0 "       << "))) || "
//  		<< " (NJetsIDLoose40  >=4 "            << "&&"
//		<< " jet[3].lv.Pt() >= 100.0 "       << "))) || "


		<< "(((trigger.HLT_SixJet45_v1 ==1)" << "||"
		<< "(trigger.HLT_SixJet45_v2 ==1)" << "||"
		<< "(trigger.HLT_SixJet45_v3 ==1)" << "||"
		<< "(trigger.HLT_SixJet45_v4 ==1)" << "||"
		<< "(trigger.HLT_SixJet45_v6 ==1)" << ") &&"
		// << " NJets  >=6 "            << "&&"
		// << " jet[5].lv.Pt() >= 80.0 "       << "))";
		// << " NJets  >=6 "            << "&&"
		// << " jet[5].lv.Pt() >= 0.0 "       << "))";
         	// << " NJetsIDLoose40  >=7 "            << "&&"
		// << " jet[6].lv.Pt() >= 60.0 "       << "))";
         	// << " NJetsIDLoose40  >=7 "            << "&&"
		// << " jet[6].lv.Pt() >= 55.0 "       << "))";
  //New Or
// 		<< " (NJetsIDLoose40  >=6 "            << "&&"
// 		<< " jet[5].lv.Pt() >= 75.0 "       << ") ||"
// 		<< " (NJetsIDLoose40  >=7 "            << "&&"
// 		<< " jet[6].lv.Pt() >= 60.0 "       << ")))))";
 		<< " ((NJetsIDLoose40  >=6 "            << "&&"
 		<< " jet[5].lv.Pt() >= 65.0 "       << ") ||"
//		<< " ((NJetsIDLoose40  >=7 "            << "&&"
//		<< " jet[6].lv.Pt() >= 55.0 "       << ") ||"
//   		<< " (NJetsIDLoose40  >=7 "            << "&&"
//   		<< " jet[6].lv.Pt() >= 55.0 "       << "))))))";
		<< " (NJetsIDLoose40  >=6 "            << "&&"
		<< " jet[5].lv.Pt() >= 65.0 "       << "))))))";



  // 	<< "(trigger.HLT_HT440_v2 ==1 && misc.Run<161216)" << "||"
  // 	<< "(trigger.HLT_HT450_v2 ==1 && (misc.Run>=161216 && misc.Run< 163269))" << "||"
  // 	<< "(trigger.HLT_HT500_v3 ==1 && (misc.Run>=163269 && misc.Run<=163869))" << "||"
  // 	<< "(trigger.HLT_HT500_v4 ==1 && (misc.Run>=165088 && misc.Run< 165970))" << "||"
  // 	<< "(trigger.HLT_HT550_v5 ==1 && (misc.Run>=165970 && misc.Run<=167043 && misc.Run!=166346))" << "||"
  // 	<< "(trigger.HLT_HT550_v6 ==1 && (misc.Run==166346))" << "||"
  // 	<< "(trigger.HLT_HT550_v7 ==1 && (misc.Run>=167078 && misc.Run< 170249))" << "||"
  // 	<< "(trigger.HLT_HT550_v8 ==1 && (misc.Run>=170249 && misc.Run< 173236))" << "||"
  // 	<< "(trigger.HLT_HT600_v1 ==1 && (misc.Run>=173236 && misc.Run< 178420))" << "||"
  // 	<< "(trigger.HLT_HT650_v4 ==1 && (misc.Run>=178420 && misc.Run< 179959))" << ")";
  TString trigger = triggerStream.str().c_str();
  TString cuts = cutStream.str().c_str();
 
  //tA->SetPileUpReweight(false);
  //tA->SetbSFWeights(false);
  //                     variable,              cuts,  njet, nbjet, nlep,     HLT,     xtitle             nbins  bins      flip_order,    log  , comp ,  ratio, stack, overlay
   //tA->makePlot("NBJetsCSVT",              cuts,    -1,  -10 ,-10 ,    trigger , "HT"      ,4,  0, 4,        false,         true ,  true,   true,  true,  true, 1,1);
  //saeid
  //NBJetsCSVT > 0 was removed from the list of the cuts.
  //Be careful, if you have a cut on bjets::
  //NBJetsCSVT >= N ==> nbjets == -N
  //NBJetsCSVT == M ==> nbjets == M
  //saeid
  
  tA->makePlot("misc.MT2",             cuts,    -1,  -1 , -10 ,   trigger , "MT2"     , 10, 0, 350,        false,         true ,  true,   true,  true,  true, 1, 1);
  //tA->TriggerEfficiency(6, 4, 0, 100000000000);
  //tA->SortHighMT2(100.0, 100000000);

  //tA->TopStudy("TTbar",3000);
  //tA->TauContamination(-1, 1000000000, 27, cuts, trigger, 0);
  //tA->vs();
  //tA->Efficiency("SMS");
  //tA->MySmallMakePlot(1000000000);
  //tA->makeSmallCopy(200000000,100);
  //tA->QCD();
  //tA->SpecialMakePlot(10000000000);
  /*
  std::vector< std::string > cuts_v ;
  cuts_v.push_back( " (misc.ProcessID!=10 || (Susy.MassGlu >= 350 && Susy.MassGlu < 355 && Susy.MassLSP >= 50 && Susy.MassLSP < 55))");
  cuts_v.push_back( triggerStream.str() );
  cuts_v.push_back( " ((NJetsIDLoose40  >=4 && jet[3].lv.Pt() >= 100.0) || (NJetsIDLoose40  >=6 && jet[5].lv.Pt() >= 65.0) || (NJetsIDLoose40  >=5 && jet[4].lv.Pt() >= 85.0))" );
  cuts_v.push_back( " (NEles==0  || ele[0].lv.Pt()<5) && (NMuons==0 || muo[0].lv.Pt()<5)" );
  cuts_v.push_back( " misc.Jet0Pass ==1 && misc.Jet1Pass ==1 && misc.PassJetID ==1 " );
  cuts_v.push_back( " misc.MET>30 ");

  cuts_v.push_back( " NBJetsCSVT> 0");
  cuts_v.push_back( " misc.MinMetJetDPhi4 >0.3 " );
  cuts_v.push_back( " misc.MT2>125");
  cuts_v.push_back( " misc.Vectorsumpt < 70 " );
  
  tA->MakeCutFlowTable( cuts_v );*/
}
