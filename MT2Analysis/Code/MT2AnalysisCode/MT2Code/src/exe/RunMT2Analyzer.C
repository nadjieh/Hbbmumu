// C++ includes
#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>

// ROOT includes
#include <TROOT.h>
#include <TTree.h>
#include <TChain.h>

#include "MT2Analyzer.hh"

#include "FWCore/FWLite/interface/AutoLibraryLoader.h"

using namespace std;

//________________________________________________________________________________________
// Print out usage
void usage( int status = 0 ) {
	cout << "Usage: RunMT2Analyzer [-d dir] [-o filename] [-v verbose] [-j json]                  " << endl;
	cout << "                      [-m set_of_cuts] [-n maxEvents] [-t type]                      " << endl;
	cout << "                      [-p data_PileUp] [-P mc_PileUP]                                " << endl; 
        cout << "                      [-s noPU,MC2012] [-C JEC]                                       " << endl;
	cout << "                      [-w pdf] [-b btag]                                             " << endl;
	cout << "                      [-r photon ] [-i ID ]                                          " << endl;
	cout << "                      [-e type1MET ] [-E CHS ]                                       " << endl;
	cout << "                      [-f FastSim ]                                                  " << endl;
	cout << "                      [-l] file1 [... filen]"                                          << endl;
	cout << "  where:"                                                                              << endl;
	cout << "     dir           is the output directory                                           " << endl;
	cout << "                   default is TempOutput/                                            " << endl;
	cout << "     filename      is the output filename for the MassAnalysis                       " << endl;
	cout << "     verbose       sets the verbose level                                            " << endl;
	cout << "                   default is 0 (quiet mode)                                         " << endl;
	cout << "     json          json file to be read                                              " << endl;
	cout << "     set_of_cuts   optional cuts for MT2Analysis                                     " << endl;
	cout << "     data_PileUp   root file from which the expected # pile-up                       " << endl;
	cout << "                   interactions is read                                              " << endl;
	cout << "     mc_PileUP     root file from which the generated # pile up                      " << endl;
	cout << "                   interactions is read                                              " << endl;
	cout << "     btag          file with btag efficiency histograms                              " << endl;
	cout << "     htau          file with hadronic tau efficiency histograms                      " << endl;
	cout << "     type          data, scan or mc=default                                          " << endl;
	cout << "     photon        add Photon to MET and remove jets/ele matched to photon           " << endl;
	cout << "     ID            ProcessID: 0=data, 1=Znunu, 2=Zll, 3=WJets, 4=Top,                " << endl;
	cout << "                              5=Gamma+jets, 6=QCD ,[empty], 9=Other, 10=Signal       " << endl;
	cout << "     JEC           redo JEC: dir of JEC to be used                                   " << endl;
	cout << "                   /shome/pnef/MT2Analysis/Code/JetEnergyCorrection/[JEC]            " << endl;
	cout << "                   ak5 pf-jets will be corrected with L1FastL2L3 (+RES if type=data) " << endl;
	cout << "     type1MET      use type1 corrected pf-MET (default = false)                      " << endl;
	cout << "     CHS           use CHS - PFjets (default = false)                                " << endl;
	cout << "     FastSim       sample is FastSim (default = false)                               " << endl;//has to be one flag, as scan can be full/fastsim
	cout << "     filen         are the input files (by default: ROOT files)                      " << endl;
	cout << "                   with option -l, these are read as text files                      " << endl;
	cout << "                   with one ROOT file name per line                                  " << endl;
	cout << endl;
	exit(status);
}

//________________________________________________________________________________________
int main(int argc, char* argv[]) {
  	AutoLibraryLoader::enable();

// Default options
	bool isList = false;
	TString outputdir = "TempOutput/";
	TString filename  = "MassTree.root";
	TString setofcuts = "nocuts";
	string  puScenario = "noPU";
  	string  jsonFileName = "";
	string  data_PileUp = "";
	string  mc_PileUp = "";
	string  type = "mc";
	string  JEC  ="";
	string  btagFileName = "";
	string  htauFileName = "";
	bool isData  = false;
	bool isCHSJets = false;
	bool isType1MET = false;
	bool isFastSim = false;
	bool isScan = false;
	int verbose  = 0;
	int maxEvents=-1;
	int ID       =-1;
    	bool removePhoton = false;
    	bool removeZll    = false;
	string photon = "";
	string pdf = "";

// Parse options
	char ch;
	while ((ch = getopt(argc, argv, "s:d:o:v:j:m:n:p:P:t:r:b:u:i:C:w:eEflh?")) != -1 ) {
	  switch (ch) {
	  case 'd': outputdir       = TString(optarg);break;
	  case 'o': filename        = TString(optarg);break;
	  case 'v': verbose         = atoi(optarg);   break;
	  case 'j': jsonFileName    = string(optarg); break;
	  case 'm': setofcuts       = TString(optarg);break;
	  case 'n': maxEvents       = atoi(optarg);   break;
	  case 'p': data_PileUp     = string(optarg); break;
	  case 'P': mc_PileUp       = string(optarg); break;
	  case 't': type            = string(optarg); break;
	  case 'r': photon          = string(optarg); break;
	  case 'b': btagFileName    = string(optarg); break;
	  case 'u': htauFileName    = string(optarg); break;
	  case 'w': pdf             = string(optarg); break;
	  case 'i': ID              = atoi(optarg);   break;
	  case 's': puScenario      = string(optarg); break;
          case 'C': JEC             = string(optarg); break;
	  case 'e': isType1MET      = true; break;
	  case 'E': isCHSJets       = true; break; 
	  case 'f': isFastSim       = true; break; 
	  case 'l': isList          = true; break;
	    //case 'noPU': noPU = true; break;  
	  case '?':
	  case 'h': usage(0); break;
	  default:
	    cerr << "*** Error: unknown option " << optarg << std::endl;
	    usage(-1);
	  }
	}

	argc -= optind;
	argv += optind;

// Check arguments
	if( argc<1 ) {
		usage(-1);
	}
	if      (type=="data") isData =true;
	else if (type=="mc"  ) isData =false;
	else if (type=="scan"  ) isScan =true;
	else    usage(-1);
	if      (photon == "photon") removePhoton=true;
	if      (photon == "Zll"   ) removeZll   =true;
	if      (type=="data" && puScenario!="noPU"){ cout << "ERROR: this is data. don't run PUreweighting" << endl; exit(-1);}
	if      (type=="mc"   && (data_PileUp.length()==0 || mc_PileUp.length()==0) && puScenario!="noPU") {
		cout << "ERROR: need puScenario" << puScenario << " required input files."  << endl; exit(-1);
	}
	setofcuts   ="../../MT2_cuts/"+setofcuts+".dat";
	//	setofcuts   ="/dataLOCAL/MT2Top/CutsMT2/"+setofcuts+".dat";
	if(data_PileUp.length()!=0   ){data_PileUp ="/dataLOCAL/MT2Top/Certification/pileUp_data/"+ data_PileUp;}
	if(mc_PileUp.length()  !=0   ){mc_PileUp   ="/dataLOCAL/MT2Top/Certification/pileUp_mc/"  +   mc_PileUp;}
	if(jsonFileName.length() !=0 ){jsonFileName="/dataLOCAL/MT2Top/Certification/JSONfiles/"  +jsonFileName;}
//	setofcuts   ="/shome/pnef/Projects/CMSAnalysis/MT2Analysis/Code/MT2_cuts/"+setofcuts+".dat";
//	if(data_PileUp.length()!=0   ){data_PileUp ="/shome/pnef/Projects/CMSAnalysis/MT2Analysis/Code/Certification/pileUp_data/"+data_PileUp;}
//	if(mc_PileUp.length()  !=0   ){mc_PileUp   ="/shome/pnef/Projects/CMSAnalysis/MT2Analysis/Code/Certification/pileUp_mc/"  + mc_PileUp;}
//	if(jsonFileName.length() !=0 ){jsonFileName="/shome/pnef/Projects/CMSAnalysis/MT2Analysis/Code/Certification/"            +jsonFileName;}
	if(btagFileName.length() !=0 ){btagFileName="/dataLOCAL/MT2Top/Efficiencies/"            + btagFileName;}
	if(htauFileName.length() !=0 ){htauFileName="/dataLOCAL/MT2Top/Efficiencies/"            + htauFileName;}



	/*

	//      setofcuts   ="/shome/haweber/MT2Analysis_8TeV/Code/CutsMT2/"+setofcuts+".dat";                                                      
        setofcuts = "/shome/paktinat/Top/CMSSW_5_3_7_patch5/src/MT2Analysis_V02-03-02/Code/MT2_Cuts/"+setofcuts+".dat";
        if(data_PileUp.length()!=0   ){data_PileUp ="/shome/haweber/MT2Analysis_8TeV/Code/Certification/pileUp_data/"+ data_PileUp;}
        if(mc_PileUp.length()  !=0   ){mc_PileUp   ="/shome/haweber/MT2Analysis_8TeV/Code/Certification/pileUp_mc/"  +   mc_PileUp;}
        //      if(jsonFileName.length() !=0 ){jsonFileName="/shome/haweber/MT2Analysis_8TeV/Code/Certification/JSONfiles/"  +jsonFileName;}        
        if(jsonFileName.length() !=0 ){jsonFileName="/shome/paktinat/Top/CMSSW_5_3_7_patch5/src/MT2Analysis_V02-03-02/Code/MT2_Cuts/" +jsonFileName\
;}

//      setofcuts   ="/shome/pnef/Projects/CMSAnalysis/MT2Analysis/Code/MT2_cuts/"+setofcuts+".dat";                                                
//      if(data_PileUp.length()!=0   ){data_PileUp ="/shome/pnef/Projects/CMSAnalysis/MT2Analysis/Code/Certification/pileUp_data/"+data_PileUp;}    
//      if(mc_PileUp.length()  !=0   ){mc_PileUp   ="/shome/pnef/Projects/CMSAnalysis/MT2Analysis/Code/Certification/pileUp_mc/"  + mc_PileUp;}     
//      if(jsonFileName.length() !=0 ){jsonFileName="/shome/pnef/Projects/CMSAnalysis/MT2Analysis/Code/Certification/"            +jsonFileName;}   
	if(btagFileName.length() !=0 ){btagFileName="/shome/haweber/MT2Analysis_8TeV/Code/Efficiencies/"            + btagFileName;}
        if(htauFileName.length() !=0 ){htauFileName="/shome/haweber/MT2Analysis_8TeV/Code/Efficiencies/"            + htauFileName;}

	*/








	std::vector<std::string> fileList;
	for(int i = 0; i < argc; i++){
		if( !isList ){
	    		fileList.push_back(argv[i]);
			printf(" Adding file: %s\n",argv[i]);
		} else {
			TString rootFile;
			ifstream is(argv[i]);
			while(rootFile.ReadLine(is) && (!rootFile.IsNull())){
	      			if(rootFile[0] == '#') continue;
				fileList.push_back(rootFile.Data());
				printf(" Adding file: %s\n", rootFile.Data());
	    		}
	  	}
	}


	cout << "--------------" << endl;
	cout << "OutputDir is:                   " << outputdir << endl;
	cout << "Type is:                        " << type << endl;
	cout << "isFastSim                       " << (isFastSim ? "true":"false") << endl;
	cout << "Verbose level is:               " << verbose << endl;
  	cout << "JSON file is:                   " << (jsonFileName.length()>0?jsonFileName:"empty") << endl;
  	cout << "MC_PileUp file:                 " << (mc_PileUp.length()>0?mc_PileUp:"empty") << endl;
  	cout << "Data_PileUp file:               " << (data_PileUp.length()>0?data_PileUp:"empty") << endl;
	cout << "PileUp Scenario:                " << puScenario << endl;
	cout << "ak5-PF have CHS                 " << (isCHSJets? "ENABLED":"DISABLED") << endl;
	cout << "pfMET is:                       " << (isType1MET? "Type1 corrected":"raw") << endl;
  	if(btagFileName.length() !=0){
	cout << "btag file is:                   " << (btagFileName.length()>0?btagFileName:"empty") << endl;
	}
  	if(htauFileName.length() !=0){
	cout << "hadronic tau file is:           " << (htauFileName.length()>0?htauFileName:"empty") << endl;
	}
	cout << "Set of Cuts is:                 " << setofcuts << endl;
	if(JEC.length()!=0){
	cout << "Redo JEC with GlobalTag         " << JEC << endl;
	}
	if(removePhoton){
	cout << "WARNING: Photon is added to MET and jet/ele match to photon is removed!!" << endl;
	}
	if(removeZll   ){
	cout << "WARNING: leptons from Z decay are added to MET (jet overlap already removed)" << endl;
	}
	if(pdf=="pdf"){  cout << "Creating PDF weights" << endl; 
	}
	cout << "--------------" << endl;

	MT2Analyzer *tA = new MT2Analyzer(fileList);
	tA->SetOutputDir(outputdir);
	tA->SetVerbose(verbose);
	tA->SetMaxEvents(maxEvents);
	tA->SetProcessID(ID);
	tA->SetBTagEfficiency(btagFileName);
	tA->SetHadTauEfficiency(htauFileName);
	tA->SetPUReweighting(puScenario);
	tA->SetType1MET(isType1MET);
	tA->SetCHSJets(isCHSJets);
	tA->SetFastSim(isFastSim);
	tA->isScan = isScan;
	tA->removePhoton = removePhoton;
	tA->removeZll    = removeZll;
	if(pdf=="pdf") tA->doPDF=true;
	else tA->doPDF=false;
  	if (jsonFileName!="") tA->ReadJSON(jsonFileName.c_str());
	tA->BeginJob(filename, setofcuts, isData, data_PileUp, mc_PileUp, JEC);
	tA->Loop();
	tA->EndJob();

	delete tA;
	return 0;
}

