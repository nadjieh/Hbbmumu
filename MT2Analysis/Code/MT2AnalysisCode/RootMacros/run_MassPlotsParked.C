{
  //TString outputdir = "../MassPlots/";
  TString outputdir = "/";
  TString samples   = "./samples/samplesMineSmearedJetMETNewBTagParked.dat";
  int verbose =3;

  gSystem->CompileMacro("../MT2Code/src/MassPlotter.cc", "kf");//"k", "f"
  gROOT->ProcessLine(".x SetStyle_PRD.C");

  int gNMT2bins                   = 17;
  double  gMT2bins[gNMT2bins+1]   = {0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 120, 150, 200, 300, 400, 550, 800}; 	
  
  int gNMT2Bbins                   = 15;
  double  gMT2Bbins[gNMT2bins+1]   = {0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 120, 150, 200, 300, 500}; 	

  int gNMT2bins_l                   = 14;
  double  gMT2bins_l[gNMT2bins+1]   = {0, 10, 20, 30, 40, 50, 65, 80, 95, 115, 140, 180, 250, 350, 500}; 	

  MassPlotter *tA = new MassPlotter(outputdir, "MassPlots.root");
  tA->SetSave(false);
  tA->setVerbose(verbose);
  tA->init(samples);
  tA->SetIsPhoton(false);
  
  std::ostringstream cutStream;
  cutStream << " " 
              << "(misc.ProcessID!=10 || (Susy.MassGlu >= 350 && Susy.MassGlu < 355"<< "&&"
              << "Susy.MassLSP >= 50 && Susy.MassLSP < 55))"                     << "&&"
    //	    << "misc.MT2>125"                                                << "&&"
	    << "NBJetsCSVT> 0"                                                << "&&"
	    << "misc.MET>30"                                                << "&&"
   // 	    << "misc.HT>1500  "                                  << "&&"
    //          << "NTops > 0"                                                   << "&&"
    //	    << "NTaus >0  && tau[0].Isolation3Hits > 1.5"                            << "&&"
// 	    << "misc.MT2<60 "                                               << "&&"
	    << "NJetsIDLoose40  >=4"                                         << "&&"
 	    << "(NEles==0  || ele[0].lv.Pt()<5)"                             << "&&"
    	    << "(NMuons==0 || muo[0].lv.Pt()<5)"                             << "&&"
 	    << "misc.Jet0Pass ==1"                                           << "&&"
 	    << "misc.Jet1Pass ==1"                                           << "&&"
 	    << "misc.PassJetID ==1"                                          << "&&"
    	    << "misc.Vectorsumpt < 70"                                       << "&&"
//	    << "(misc.MinMetJetDPhi >0.3||misc.MinMetJetDPhiIndex>3)"        << "&&"
 	    << "misc.MinMetJetDPhi4 >0.3"                                    << "&&"
 // 	    << "misc.MinMetJetDPhi4 <=0.2"                                   << "&&"

	    
    // QuadJet45/50
	    <<"((NJetsIDLoose40  >=4 && jet[3].lv.Pt() >= 80.0) || (NJetsIDLoose40  >=5 && jet[4].lv.Pt() >= 65.0))";



  std::ostringstream triggerStream;
  triggerStream << "( "
		<<" misc.ProcessID != 0 || ((misc.CrazyHCAL==0 && misc.NegativeJEC==0 " <<"&&"
		<<" misc.CSCTightHaloIDFlag==0 && misc.HBHENoiseFlag==0 " <<"&&"
		<<" misc.hcalLaserEventFlag==0 && misc.trackingFailureFlag==0 " <<"&&"
		<<" misc.eeBadScFlag==0 && misc.EcalDeadCellTriggerPrimitiveFlag==0 )" <<"&&("

   //Or between Six and Quad
		<< "(((trigger.HLT_QuadJet50_v1 ==1)" << "||"
		<< "(trigger.HLT_QuadJet50_v2 ==1)" << "||"
		<< "(trigger.HLT_QuadJet50_v3 ==1)" << "||"
		<< "(trigger.HLT_QuadJet50_v5 ==1)" << ") &&"
	
		<< " ((NJetsIDLoose40  >=4 "            << "&&"
		<< " jet[3].lv.Pt() >= 80.0 "       << ") || "
 		<< " (NJetsIDLoose40  >=5 "            << "&&"
 		<< " jet[4].lv.Pt() >= 65.0 "       << "))))))";

  TString trigger = triggerStream.str().c_str();
  TString cuts = cutStream.str().c_str();
 
//   tA->SetPileUpReweight(false);
//   tA->SetbSFWeights(false);
  //                     variable,              cuts,  njet, nbjet, nlep,     HLT,     xtitle             nbins  bins      flip_order,    log  , comp ,  ratio, stack, overlay
  //  tA->makePlot("NBJetsCSVT",              cuts,    -1,  -10 ,-10 ,    trigger , "HT"      ,4,  0, 4,        false,         true ,  true,   true,  true,  true, 1,1);
  //saeid
  //NBJetsCSVT > 0 was removed from the list of the cuts.
  //Be careful, if you have a cut on bjets::
  //NBJetsCSVT >= N ==> nbjets == -N
  //NBJetsCSVT == M ==> nbjets == M
  //saeid

  tA->makePlot("misc.MT2",             cuts,    -1,  -1 , -10 ,   trigger , "MT2"     , 35, 0, 700,        false,         true ,  true,   true,  true,  true, 1, 1);
  //tA->TriggerEfficiency(6, 4, 0, 100000000000);
  //tA->SortHighMT2(100.0, 100000000);

  //tA->TopStudy("TTbar",3000);
  //tA->TauContamination(-1, 1000000000, 27, cuts, trigger, 1);
  //tA->vs();
  //tA->Efficiency("SMS");
  //tA->MySmallMakePlot(1000000000);
  //tA->makeSmallCopy(200000000,100, cuts, trigger);
  //tA->QCD();
  //tA->SpecialMakePlot(10000000000, cuts, trigger);

  /*
  std::vector< std::string > cuts_v ;
  cuts_v.push_back( " (misc.ProcessID!=10 || (Susy.MassGlu >= 350 && Susy.MassGlu < 355 && Susy.MassLSP >= 50 && Susy.MassLSP < 55))");
  cuts_v.push_back( triggerStream.str() );
  cuts_v.push_back( " ((NJetsIDLoose40  >=4 && jet[3].lv.Pt() >= 80.0) || (NJetsIDLoose40  >=5 && jet[4].lv.Pt() >= 65.0))" );
  cuts_v.push_back( " (NEles==0  || ele[0].lv.Pt()<5) && (NMuons==0 || muo[0].lv.Pt()<5)" );
  cuts_v.push_back( " misc.Jet0Pass ==1 && misc.Jet1Pass ==1 && misc.PassJetID ==1 " );
  cuts_v.push_back( " misc.MET>30 ");

  cuts_v.push_back( " NBJetsCSVT> 0");
  cuts_v.push_back( " misc.MinMetJetDPhi4 >0.3 " );
  cuts_v.push_back( " misc.MT2>125");
  cuts_v.push_back( " misc.Vectorsumpt < 70 " );
  
  tA->MakeCutFlowTable( cuts_v );*/

}
