
{ TString outputdir = "../MassPlots/";
  
  TString samples = "./samples/samplesMineMuEG.dat";
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
            << "(eleMu.MT2>=0 )" <<"&&"  
   	    << "( eleMu.mu0Ind >=0  )" <<"&&"
 	    << "(eleMu.ele0Ind >=0  )" <<"&&"  //MuEG selections  
            <<"(eleMu.HasNoVetomuoForEleMu)" <<"&&"
            <<"(eleMu.HasNoVetoElecForEleMu)"<<"&&"
            <<"( eleMu.charge ==0)"<<"&&"
            <<"((muo[eleMu.mu0Ind].lv.Pt()>20)||(ele[eleMu.ele0Ind].lv.Pt()>20) )"<<"&&"
            <<"0==0";

  std::ostringstream triggerStream;
  triggerStream << "( "

		<<" misc.ProcessID != 0 || ((misc.CrazyHCAL==0 && misc.NegativeJEC==0 " <<"&&"
		<<" misc.CSCTightHaloIDFlag==0 && misc.HBHENoiseFlag==0 " <<"&&"
		<<" misc.hcalLaserEventFlag==0 && misc.trackingFailureFlag==0 " <<"&&"
		<<" misc.eeBadScFlag==0 && misc.EcalDeadCellTriggerPrimitiveFlag==0 )" <<"&&("

  
     	<< "(trigger.HLT_EMu)"         << ")))";
    
            
  TString trigger = triggerStream.str().c_str();
  TString cuts = cutStream.str().c_str();
 

//    makeplots( TString var, TString maincuts, TString basecuts, int njets, int nbjets, int nleps, TString HLT,
			   //TString xtitle, const int nbins, const double *bins, 
			  // bool flip_order, bool logflag, bool composited, bool ratio, 
			  // bool stacked, bool overlaySUSY, float overlayScale, bool add_underflow, bool saveHistos)
//  tA->makePlot("NBJetsCSVM",          cuts,    -10,  -10 , -10 ,   trigger , "NBJetsCSVM"     ,10,0,10,        false,        true ,  true,   true,  true,  true, 1, 1);
   tA->makePlot("NEles",     cuts,    -10,  -10 , -10 ,   trigger , "n_{e}"            ,10,0,10,        false,          true ,  true,   true,  true,  true, 1, 1);  
  //tA->makePlot("eleMu.mu0Ind",     cuts,    -10,  -10 , -10 ,   trigger , "eleMu.mu0Ind"            ,10,0,10       false,        true ,  true,   true,  true,  true, 1, 1);     
 // tA->makePlot("eleMu.ele0Ind",     cuts,    -10,  -10 , -10 ,   trigger , "eleMu.ele0Ind"            ,10,0,10       false,        true ,  true,   true,  true,  true, 1, 1);     
   //  tA->makePlot("eleMu.MT2",     cuts,    -10,  -10 , -10 ,   trigger , "eleMu.MT2"            ,25,0,1000,        false,        true ,  true,   true,  true,  true, 1, 1);     
 // tA->makePlot("eleMu.METImbalanced",     cuts,    -10,  -10 , -10 ,   trigger , "eleMu.METImbalanced"            ,25,0,1000,        false,        true ,  true,   true,  true,  true, 1, 1);     
  //tA->makePlot("eleMu.MT2Imbalanced",     cuts,    -10,  -10 , -10 ,   trigger , " eleMu.MT2Imbalanced"            ,25,0,1000,        false,        true ,  true,   true,  true,  true, 1, 1);     
   


  



  //tA->TriggerEfficiency(6, 4, 0, 100000000000);jet[3].lv.Pt()
  //tA->SortHighMT2(100.0, 100000000);
  // tA->plotSig("GetMT2elemuon()", cuts, "MT2", 40, 0, 800, true, 1 ,1);

  //tA->TopStudy("TTbar",3000);
  //tA->TauContamination(-1, 1000000000, 27);
  //tA->vs();
  //tA->Efficiency("SMS");
  //tA->MySmallMakePlot(1000000000);
  //tA->makeSmallCopy(200000000,100);
  //tA->QCD();
  //tA->SpecialMakePlot(10000000000);
}
