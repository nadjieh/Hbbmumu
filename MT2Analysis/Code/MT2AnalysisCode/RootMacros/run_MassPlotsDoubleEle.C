
{ TString outputdir = "../MassPlots/";
  
  TString samples = "./samples/samplesMineDoubleEle.dat";
  int verbose =3;
  

   gSystem->CompileMacro("../MT2Code/src/MassPlotter.cc", "kf");//"k", "f"
   gROOT->ProcessLine(".x SetStyle_PRD.C");

  int gNMT2bins                   = 17;
  double  gMT2bins[gNMT2bins+1]   = {0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 120, 150, 200, 300, 400, 550, 800}; 	
  
  int gNMT2Bbins                   = 15;
  double  gMT2Bbins[gNMT2bins+1]   = {0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 120, 150, 200, 300, 500}; 	

  int gNMT2bins_l                   = 14;
  double  gMT2bins_l[gNMT2bins+1]   = {0, 10, 20, 30, 40, 50, 65, 80, 95, 115, 140, 180, 250, 350, 500}; 	

  MassPlotter *tA = new MassPlotter(outputdir, "MassPlots-DoubleEle.root");
  tA->SetSave(false);
  tA->setVerbose(verbose);
  tA->init(samples);
  tA->SetIsPhoton(false);
  
  
 
   std::ostringstream cutStream;
   cutStream << " " 
          
   	    << "(doubleEle.ele0Ind >=0  )" <<"&&"
 	    << "(doubleEle.ele1Ind >=0  )" <<"&&"
           
            << "(doubleEle.MT2>=0 )" <<"&&"            
  
            <<"((ele[doubleEle.ele0Ind].lv.Pt()>20)&&(ele[doubleEle.ele1Ind].lv.Pt()>10) )"<<"&&"
            <<"0==0";

  std::ostringstream triggerStream;
  triggerStream << "( "

		<<" misc.ProcessID != 0 || ((misc.CrazyHCAL==0 && misc.NegativeJEC==0 " <<"&&"
		<<" misc.CSCTightHaloIDFlag==0 && misc.HBHENoiseFlag==0 " <<"&&"
		<<" misc.hcalLaserEventFlag==0 && misc.trackingFailureFlag==0 " <<"&&"
		<<" misc.eeBadScFlag==0 && misc.EcalDeadCellTriggerPrimitiveFlag==0 )" <<"&&("
 

              << "(trigger.HLT_DiElectrons) " << ")))";
            
  TString trigger = triggerStream.str().c_str();
  TString cuts = cutStream.str().c_str();
 

//  tA->makePlot("NBJetsCSVM",          cuts,    -10,  -10 , -10 ,   trigger , "NBJetsCSVM"     ,10,0,10,        false,        true ,  true,   true,  true,  true, 1, 1);
///  tA->makePlot("GetMT2MuEG()",     cuts,    -10,  -10 , -10 ,   trigger , "MT2"            ,20,0,400,        false,        true ,  true,   true,  true,  true, 1, 1);
//  tA->makePlot("doubleTau.signalDoubleTau",     cuts,    -10,  -10 , -10 ,   trigger , "doubleTau"            ,5,-2.5,2.5,        false,        true ,  true,   true,  true,  true, 1, 1);     
//  tA->makePlot("muTau.signalMuTau",     cuts,    -10,  -10 , -10 ,   trigger , "muTau"            ,5,-2.5,2.5,        false,        true ,  true,   true,  true,  true, 1, 1);     
//  tA->makePlot("eleMu.mu0Ind",     cuts,    -10,  -10 , -10 ,   trigger , "eleMu.mu0Ind"            ,5,-2.5,2.5,        false,        true ,  true,   true,  true,  true, 1, 1);     
//  tA->makePlot("eleMu.mu1Ind",     cuts,    -10,  -10 , -10 ,   trigger , "eleMu.mu1Ind"            ,5,-2.5,2.5,        false,        true ,  true,   true,  true,  true, 1, 1);     
 
 tA->makePlot("eleMu.MT2",     cuts,    -10,  -10 , -10 ,   trigger , "eleMu.MT2"            ,20,0,1000,        false,        true ,  true,   true,  true,  true, 1, 1);     

 tA->makePlot("eleMu.METImbalanced",     cuts,    -10,  -10 , -10 ,   trigger , "eleMu.METImbalanced"            ,20,0,1000,        false,        true ,  true,   true,  true,  true, 1, 1);     

 tA->makePlot("eleMu.MT2Imbalanced",     cuts,    -10,  -10 , -10 ,   trigger , " eleMu.MT2Imbalanced"            ,20,0,1000,        false,        true ,  true,   true,  true,  true, 1, 1);     

  // tA->makePlot("doubleEle.NSelEle",     cuts,    -10,  -10 , -10 ,   trigger , "doubleEle.NSelEle"            ,20,0,100,        false,        true ,  true,   true,  true,  true, 1, 1);     
  
   //tA->makePlot("eleMu.MT2",     cuts,    -10,  -10 , -10 ,   trigger , "eleMu.MT2"            ,100,-2,1000,        false,        true ,  true,   true,  true,  true, 1, 1);     


  //tA->makePlot("doubleMu.MT2",     cuts,    -10,  -10 , -10 ,   trigger , "doubleEle.MT2"            ,5,-2.5,2.5,        false,        true ,  true,   true,  true,  true, 1, 1);     

 tA->makePlot("doubleEle.MT2",     cuts,    -10,  -10 , -10 ,   trigger , "doubleEle.MT2"            ,20,0,1000,        false,        true ,  true,   true,  true,  true, 1, 1);     
 
 tA->makePlot("doubleEle.METImbalanced",     cuts,    -10,  -10 , -10 ,   trigger , "doubleEle.METImbalanced"            ,20,0,1000,        false,        true ,  true,   true,  true,  true, 1, 1);     
 
 tA->makePlot("doubleEle.MT2Imbalanced",     cuts,    -10,  -10 , -10 ,   trigger , "doubleEle.MT2Imbalanced"            ,20,0,1000,        false,        true ,  true,   true,  true,  true, 1, 1);     

 tA->makePlot("NEles",     cuts,    -10,  -10 , -10 ,   trigger , "NEles"            ,10,0,10,        false,        true ,  true,   true,  true,  true, 1, 1);     

 tA->makePlot("NMuons",     cuts,    -10,  -10 , -10 ,   trigger , "NMuons"            ,10,0,10,        false,        true ,  true,   true,  true,  true, 1, 1);     

 tA->makePlot("NTaus",     cuts,    -10,  -10 , -10 ,   trigger , "NTaus"            ,10,0,10,        false,        true ,  true,   true,  true,  true, 1, 1);     

 // tA->makePlot("misc.MET",     cuts,    -10,  -10 , -10 ,   trigger , "MET"            ,10,0,10,        false,        true ,  true,   true,  true,  true, 1, 1);     


 // tA->makePlot("doubleEle.MT2",     cuts,    -10,  -10 , -10 ,   trigger , "doubleEle.MT2"            ,100,-2,1000,        false,        true ,  true,   true,  true,  true, 1, 1);     


  // tA->makePlot("GetMT2DoubleMu()",     cuts,    -10,  -10 , -10 ,   trigger , "MT2"            ,20,0,400,        false,        true ,  true,   true,  true,  true, 1, 1); 



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

   


  


