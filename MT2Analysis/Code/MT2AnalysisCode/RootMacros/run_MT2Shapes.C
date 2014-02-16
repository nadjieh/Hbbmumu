{
	TString outputdir = "../MT2Shapes/";
	TString samples   = "./samples/samples_2140_highMT2_oldQCD_fastsim.dat";
	int verbose =3;

	gSystem->CompileMacro("../MT2Code/src/MT2Shapes.cc", "k");
	gROOT->ProcessLine(".x SetStyle_PRD.C");

	int gNMT2bins                   = 19;
	double  gMT2bins[gNMT2bins+1]   = {0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 60, 70, 85, 105, 135, 180, 260, 360, 650}; 	

	int gNMT2bins_l                   = 14;
	double  gMT2bins_l[gNMT2bins+1]   = {0, 10, 20, 30, 40, 50, 65, 80, 95, 115, 140, 180, 250, 350, 500}; 	

	MT2Shapes *tA = new MT2Shapes(outputdir, "MT2Shapes.root");
	tA->setVerbose(verbose);
	tA->init(samples);
  
	std::ostringstream cutStream;
	cutStream << " " 
          << "misc.MET>=30"                            << "&&"
	  << "misc.HT > 400 "                          << "&&"
	  << "misc.caloHT50_ID >600 "                  << "&&"
	  << "misc.Jet0Pass ==1"                       << "&&"
	  << "misc.Jet1Pass ==1"                       << "&&"
//	  << "misc.LeadingJPt >150"                    << "&&"
	  << "misc.SecondJPt  >100"                    << "&&"
//	  << "NBJets >0"                               << "&&"
	  << "misc.PassJetID ==1"                      << "&&"
	  << "misc.Vectorsumpt < 70"                   << "&&"
	  << "misc.MinMetJetDPhi >0.3"                 << "&&"
	  << "misc.HBHENoiseFlag == 0"                 << "&&"
	  << "misc.CrazyHCAL==0";
	
	std::ostringstream triggerStream;
	triggerStream << "( "
	<< "(trigger.HLT_HT440_v2==1 && misc.Run<=161176)" << "||"
	<< "(trigger.HLT_HT450_v2==1 && (misc.Run>=161217 && misc.Run<=163261))" << "||"
	<< "(trigger.HLT_HT500_v3==1 && (misc.Run>=163270 && misc.Run<=163869))" << "||"
	<< "(trigger.HLT_HT550_v4==1 && (misc.Run>=165088 && misc.Run<=165633))" << "||"
	<< "(trigger.HLT_HT550_v5==1 && (misc.Run>=165970 && misc.Run<=167043 && misc.Run!=166346))" << "||"
	<< "(trigger.HLT_HT550_v6==1 && (misc.Run==166346))"                     << "||"
	<< "(trigger.HLT_HT550_v7==1 && (misc.Run>=167078 && misc.Run<=167913))" << "||"
	<< "(trigger.HLT_HT550_v8==1 && (misc.Run>=169561 && misc.Run<=172619))" << " )";
	TString trigger = triggerStream.str().c_str();

  	TString cuts = cutStream.str().c_str();
    
//                    variable,   cuts,  njet, nlep,  selection_name,      HLT,      xtitle      nbins  bins   
        tA->GetShapes("misc.MT2",  cuts,    -3,  0  , "HadronicRegion",    trigger , "M_{T2}"   , 70, 0, 700);
        tA->GetShapes("misc.MT2",  cuts,    -3,  1  , "oneLeptonRegion",   trigger , "M_{T2}"   , 20, 0, 300);
    
}
