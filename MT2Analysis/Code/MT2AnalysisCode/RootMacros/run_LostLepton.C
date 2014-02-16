/*********************************************************************************
*  Pascal, MT2Analysis, May 26th 2011                                            *
*                                                                                *
*********************************************************************************/
#include <iomanip>

// USER INPUT -------------------------------------------
TString outputdir                 = "LostLepton/test/";
TString fTrigger                  = "HT"; // HT or MHT_HT
Bool_t  fPrintTable               = true;
Bool_t  fMakeEfficiencies         = true;
Int_t   fTypeOfEffErr             = 0;
Bool_t  fMakeEfficienciesPresel   = true;  
Bool_t  fMakePrediction           = true;
Bool_t  fIncludeTaus              = true;  
Bool_t  fIncludeTop               = true;  // include top into signal rather than background
Bool_t  fTopEfficencies           = true;  // set to true to take the lepton efficiencies from top
Bool_t  fTopOnly                  = false;  // set to true to treat W as background
Bool_t  fIncludeSingleTop         = false;  // set to true to include single top as signal
Bool_t  fDebug                    = true;
Int_t   fNEventsForEff            = 10000000;
Bool_t  fLLSkim                   = false;
Bool_t  fTagandProbeMC            = false;
Bool_t  fMakeEfficienciesPablo    = false; // set to true in order to use Pablo's efficiencies
Int_t   fVerbose                  = 2;
Bool_t  fMCClosure                = false;
Bool_t  fWeightedProb             = true; // set to true to account for relative contribution of W and Top for prob and prob_err
TString fPablofile_ele_data = "/shome/pnef/SUSY/SUSY_macros/LeptJetMult/efficiencies/effElectron.root"; // pablo's files, currently not used
TString fPablofile_muo_data = "/shome/pnef/SUSY/SUSY_macros/LeptJetMult/efficiencies/effMuon.root";     // pablo's files, currently not used
TString fPablofile_ele_mc   = "/shome/pnef/SUSY/SUSY_macros/LeptJetMult/efficiencies/effElectronMC.root"; // pablo's files, currently not used
TString fPablofile_muo_mc   = "/shome/pnef/SUSY/SUSY_macros/LeptJetMult/efficiencies/effMuonMC.root";     // pablo's files, currently not used

// --------------------------------------------------------
TH2F* WeleEffHistdR;
TH2F* WeleEffHistEta;
TH2F* WeleEffHistdRMC;
TH2F* WeleEffHistEtaMC;
TH2F* WmuoEffHistdR;
TH2F* WmuoEffHistEta;
TH2F* WmuoEffHistdRMC;
TH2F* WmuoEffHistEtaMC;

struct leptonicW{
	double nData;
	double nW;
	double nW_bg;
	double nW_leptveto;
	double TopnW;
	double TopnW_bg;
	double TopnW_leptveto;
	double TopnW_s;
	double TopnW_s_bg;
	double TopnW_s_leptveto;
	double TopnW_t;
	double TopnW_t_bg;
	double TopnW_t_leptveto;
	double TopnW_tW;
	double TopnW_tW_bg;
	double TopnW_tW_leptveto;
	double prob;
	double prob_err;
	double acc;
	double acc_err;
	double rec;
	double rec_err;
	double Topprob;
	double Topprob_err;
	double Topacc;
	double Topacc_err;
	double Toprec;
	double Toprec_err;
	double Top_bg;
	double W_bg;
	double Z_bg;
	double QCD_bg;
	double Other_bg;
	TTree* EffTree;
	vector<double> pt;
	vector<double> eta;
	vector<double> dR;
	vector<double> eff;
	vector<double> err;
	vector<double> ptMC;
	vector<double> etaMC;
	vector<double> dRMC;
	vector<double> effMC;
	vector<double> errMC;
	static const double rel_sys_uncert = 0.05;
	static const double rel_sys_uncert_bg = 1.;
	static const double Wnevents       =  5363746;
	static const double Topnevents     =  3701872.;
	static const double Topnevents_s   =  494967.;
	static const double Topnevents_t   =  474380.;
	static const double Topnevents_tW  =  489417.;
	static const double Wxsection      = 48.49*1.06;
	static const double Topxsection    = 157.5.;
	static const double Topxsection_s  =   1.49;
	static const double Topxsection_t  =  20.54;
	static const double Topxsection_W  =  10.6;
	static const double lumi           = 3161.0;
	double nW_scaled(){ 
		double res=0;
		if(fIncludeTop)        res +=TopnW*lumi*Topxsection/Topnevents;
		if(fIncludeSingleTop)  {
			res +=TopnW_s* lumi*Topxsection_s  /Topnevents_s;
			res +=TopnW_t* lumi*Topxsection_t  /Topnevents_t;
			res +=TopnW_tW*lumi*Topxsection_tW /Topnevents_tW;
		}
		if(!fTopOnly)           res +=nW*lumi*Wxsection/Wnevents;
		return res;
	}
	double nW_bg_scaled(){ 
		double res=0;
		if(fIncludeTop)        res +=TopnW_bg*lumi*Topxsection/Topnevents;
		if(fIncludeSingleTop)  {
			res +=TopnW_bg_s* lumi*Topxsection_s  /Topnevents_s;
			res +=TopnW_bg_t* lumi*Topxsection_t  /Topnevents_t;
			res +=TopnW_bg_tW*lumi*Topxsection_tW /Topnevents_tW;
		}
		if(!fTopOnly)           res +=nW_bg*lumi*Wxsection/Wnevents;
		return res;
	}
	double nW_leptveto_scaled(){ 
		double res=0;
		if(fIncludeTop)        res +=TopnW_leptveto*lumi*Topxsection/Topnevents;
		if(fIncludeSingleTop)  {
			res +=TopnW_leptveto_s* lumi*Topxsection_s  /Topnevents_s;
			res +=TopnW_leptveto_t* lumi*Topxsection_t  /Topnevents_t;
			res +=TopnW_leptveto_tW*lumi*Topxsection_tW /Topnevents_tW;
		}
		if(!fTopOnly)           res +=nW_leptveto*lumi*Wxsection/Wnevents;
		return res;
	}
	double bg(){
		if(fIncludeTop && !fTopOnly) return nW_bg_scaled()            +Z_bg + QCD_bg + Other_bg; 
		else if(fTopOnly)            return nW_bg_scaled() +W_bg      +Z_bg + QCD_bg + Other_bg; 
		else                         return nW_bg_scaled() +Top_bg    +Z_bg + QCD_bg + Other_bg; // add full top as backgroun
	}
	double all_MC(){
		if(fIncludeTop && !fTopOnly) return nW_scaled()               + Z_bg + QCD_bg + Other_bg;
		else if(fTopOnly)            return nW_scaled() + W_bg        + Z_bg + QCD_bg + Other_bg;
		else                         return nW_scaled() + Top_bg      + Z_bg + QCD_bg + Other_bg; 
	}
	double pred(){
		if(fMCClosure){
			return (all_MC()-bg())*(1-prob())/prob();
		}else{
			return (nData - bg())*(1-prob())/prob();
		}
	}
	double pred_error_stat(){
		if(fMCClosure) return 0;
		else return    fabs(sqrt(nData)*(1-prob())/prob());
	}
	double pred_error_sys(){
		if(fMCClosure){
			return sqrt(  pow((all_MC()-bg()) *(1/(prob()*prob())*prob_err_sys()),2)
			             +pow(rel_sys_uncert_bg*bg()*(1-prob())/prob(),2));
		}else{
			return sqrt(  pow((nData-bg()) *(1/(prob()*prob())*prob_err_sys()),2)
			     +pow(rel_sys_uncert_bg*bg()*(1-prob())/prob(),2));
		}
	}
	double prob(){
		if(! fMakeEfficienciesPablo){
			if(fWeightedProb){
				double WEvents   = (nW - nW_bg)      *lumi*Wxsection/Wnevents;
				double TopEvents = (TopnW - TopnW_bg)*lumi*Topxsection/Topnevents;
				return (Topprob*TopEvents + prob*WEvents)/(TopEvents+WEvents);
			}else{
				if(fTopEfficencies)  return Topprob;
				else                 return prob;
			}
		}else if(fMCClosure && fMakeEfficienciesPablo &&!fWeightedProb){
			if(fTopEfficencies)  return TPRecoMC()*Topacc;
			else                 return TPRecoMC()*acc;
		}else if(!fMCClosure && fMakeEfficienciesPablo &&!fWeightedProb){
			if(fTopEfficencies)  return TPReco()*Topacc;
			else                 return TPReco()*acc;
		}else{
			cout << " what do you want? " << endl; exit(1);
		}
	}
	double prob_err_sys(){
		if(fMakeEfficienciesPablo){
			if(fMCClosure) {
				if(fTopEfficencies) {
					return sqrt(pow(Topacc*TPRecoErrMC(),2)+pow(Accerr()*TPRecoMC(),2));
				}else{
					return sqrt(pow(acc*TPRecoErrMC(),2)+pow(Accerr()*TPRecoMC(),2));
				}
			}else   {
				if(fTopEfficencies) { 
					return sqrt(pow(Topacc*TPRecoErr(),2)  +pow(Accerr()*TPReco()  ,2));
				}else {
					return sqrt(pow(acc*TPRecoErr(),2)  +pow(Accerr()*TPReco()  ,2));
				}
			}
		}else{
			if(fWeightedProb){
				double WEvents   = (nW - nW_bg)      *lumi*Wxsection/Wnevents;
				double TopEvents = (TopnW - TopnW_bg)*lumi*Topxsection/Topnevents;
				double err       = (Topprob_err*TopEvents + prob_err*WEvents)/(TopEvents+WEvents);
				return sqrt(pow(rel_sys_uncert*prob(),2)+ pow(err,2));
				
			}else{
				if(fTopEfficencies){
					return sqrt(pow(rel_sys_uncert*Topprob,2)+ pow(Topprob_err,2));
				}else{
					return sqrt(pow(rel_sys_uncert*prob,2)+ pow(prob_err,2));
				}
			}
		}
	}
	double TPRecoMC(){
			double reco_eff=0;
			for(int i=0; i<effMC.size(); ++i) reco_eff+=effMC[i];
			reco_eff=reco_eff/effMC.size();
			return reco_eff;
	}
	double TPRecoErrMC(){
			double reco_err=0;
			for(int i=0; i<errMC.size(); ++i) reco_err+=errMC[i]*errMC[i];
			reco_err=sqrt(reco_err)/errMC.size();
			return sqrt(reco_err*reco_err+ pow(rel_sys_uncert*TPRecoMC(),2));
	}
	double TPReco(){
			double reco_eff=0;
			for(int i=0; i<eff.size(); ++i) reco_eff+=eff[i];
			reco_eff=reco_eff/eff.size();
			return reco_eff;
	}
	double TPRecoErr(){
			double reco_err=0;
			for(int i=0; i<err.size(); ++i) reco_err+=err[i]*err[i];
			reco_err=sqrt(reco_err)/err.size();
			return sqrt(reco_err*reco_err+ pow(rel_sys_uncert*TPReco(),2));
	}
	double Accerr(){
		if(fTopEfficencies) return sqrt(pow(rel_sys_uncert*Topacc,2)+Topacc_err*Topacc_err);
		else                return sqrt(pow(rel_sys_uncert*acc,2)   +acc_err*acc_err);
	}
} Wele, Wmuo;
double nTau_leptveto=0;

//__________________________________
const int gNMT2bins                   = 19;
double    gMT2bins[gNMT2bins+1]   = {0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 60, 70, 85, 105, 135, 180, 260, 360, 500};

const int gNMT2bins_l                   = 14;
double    gMT2bins_l[gNMT2bins+1]   = {0, 10, 20, 30, 40, 50, 65, 80, 95, 115, 140, 180, 250, 350, 500};
//__________________________________



void run_LostLepton(){

	gSystem->Load("libPhysics");
	gSystem->CompileMacro("../MT2Code/src/MassPlotter.cc", "k");

	if( (fTopOnly && !fTopEfficencies) || (!fIncludeTop && fTopOnly) || (fIncludeSingleTop && !fIncludeTop)){
		cout << "make up your mind" << endl; exit(1);
	}

	if(fTagandProbeMC &&!fMakeEfficienciesPablo) exit(1);
        // ------------------------
	if(fMakeEfficiencies)      getEfficiencies( fNEventsForEff);
	if(fMakeEfficienciesPablo) readEfficienciesPablo();
	if(fMakePrediction)   {
		makePlots();
		makePrediction();
	}
	if(fPrintTable) printTable();
        // ------------------------
}

void getEfficiencies(Long64_t nevents){
	
	int verbose = 3;

	MassPlotter *tA = new MassPlotter(outputdir, "MassPlots.root");
	tA->setVerbose(verbose);
	if(fMakeEfficienciesPresel){
		if(fLLSkim) tA->init("samples_LLSkim/samples_WJets_Presel.dat");	
		else        tA->init("samples_LostLepton/samples_WJets_Presel.dat");	
	} else  tA->init("samples_LLSkim/samples_WJets.dat");
		
	tA->PrintWEfficiency(0 ,"W", "ele", nevents, fIncludeTaus);
	tA->PrintWEfficiency(0 ,"W", "muo", nevents, fIncludeTaus);

	Wele.prob          =tA->fWpred.Wenu_prob;
	Wele.prob_err      =tA->fWpred.Wenu_prob_err;
	Wmuo.prob          =tA->fWpred.Wmunu_prob;
	Wmuo.prob_err      =tA->fWpred.Wmunu_prob_err;
	
	Wele.acc           =tA->fWpred.Wenu_acc;
	Wele.acc_err       =tA->fWpred.Wenu_acc_err;
	Wele.rec           =tA->fWpred.Wenu_rec;
	Wele.rec_err       =tA->fWpred.Wenu_rec_err;
	
	Wmuo.acc           =tA->fWpred.Wmunu_acc;
	Wmuo.acc_err       =tA->fWpred.Wmunu_acc_err;
	Wmuo.rec           =tA->fWpred.Wmunu_rec;
	Wmuo.rec_err       =tA->fWpred.Wmunu_rec_err;


	if(fIncludeTop){
		if(fMakeEfficienciesPresel) {
			if(fLLSkim) tA->init("samples_LLSkim/samples_Top_Presel.dat");	
			else        tA->init("samples_LostLepton//samples_Top_Presel.dat");	
		}
		else       tA->init("samples_LostLepton//samples_Top.dat");
		tA->PrintWEfficiency(0 ,"Top", "ele", nevents, fIncludeTaus);
		tA->PrintWEfficiency(0 ,"Top", "muo", nevents, fIncludeTaus);
		Wele.Topacc           =tA->fWpred.TopWenu_acc;
		Wele.Topacc_err       =tA->fWpred.TopWenu_acc_err;
		Wele.Toprec           =tA->fWpred.TopWenu_rec;
		Wele.Toprec_err       =tA->fWpred.TopWenu_rec_err;
		Wmuo.Topacc           =tA->fWpred.TopWmunu_acc;
		Wmuo.Topacc_err       =tA->fWpred.TopWmunu_acc_err;
		Wmuo.Toprec           =tA->fWpred.TopWmunu_rec;
		Wmuo.Toprec_err       =tA->fWpred.TopWmunu_rec_err;
		
		Wele.Topprob          =tA->fWpred.TopWenu_prob;
		Wele.Topprob_err      =tA->fWpred.TopWenu_prob_err;
		Wmuo.Topprob          =tA->fWpred.TopWmunu_prob;
		Wmuo.Topprob_err      =tA->fWpred.TopWmunu_prob_err;

		cout << "Wele prob    " << Wele.prob    << " pm " << Wele.prob_err    << endl;
		cout << "Wele acc     " << Wele.acc     << " pm " << Wele.acc_err     << endl;
		cout << "Wele rec     " << Wele.rec     << " pm " << Wele.rec_err     << endl;
		cout << "TopWele prob " << Wele.Topprob << " pm " << Wele.Topprob_err << endl;
		cout << "TopWele acc  " << Wele.Topacc  << " pm " << Wele.Topacc_err     << endl;
		cout << "TopWele rec  " << Wele.Toprec  << " pm " << Wele.Toprec_err     << endl;
		cout << "Wmuo prob    " << Wmuo.prob    << " pm " << Wmuo.prob_err    << endl;
		cout << "Wmuo acc     " << Wmuo.acc     << " pm " << Wmuo.acc_err     << endl;
		cout << "Wmuo rec     " << Wmuo.rec     << " pm " << Wmuo.rec_err     << endl;
		cout << "TopWmuo prob " << Wmuo.Topprob << " pm " << Wmuo.Topprob_err << endl;
		cout << "TopWmuo acc  " << Wmuo.Topacc  << " pm " << Wmuo.Topacc_err     << endl;
		cout << "TopWmuo rec  " << Wmuo.Toprec  << " pm " << Wmuo.Toprec_err     << endl;
	}else{
		Wele.Topacc           =0;
		Wele.Topacc_err       =0;
		Wele.Toprec           =0;
		Wele.Toprec_err       =0;
		Wmuo.Topacc           =0;
		Wmuo.Topacc_err       =0;
		Wmuo.Toprec           =0;
		Wmuo.Toprec_err       =0;
		
		Wele.Topprob          =0;
		Wele.Topprob_err      =0;
		Wmuo.Topprob          =0;
		Wmuo.Topprob_err      =0;
	}


	delete tA;
}

void makePlots(){
	TString samples;
	if     (fTrigger=="HT")     {
		if(fLLSkim) samples = "samples_LLSkim/samples_MG_V02-01-02_golden_HTonly.dat";
		else        samples = "samples_LostLepton/samples_MG_V02-01-02_golden_HTonly.dat";
	}
	else if(fTrigger=="MHT_HT") samples = "samples_LLSkim/samples_MG_V02-01-02_golden.dat";
	int verbose = 3;

	MassPlotter *tA = new MassPlotter(outputdir, "MassPlots.root");
	tA->setVerbose(verbose);
	tA->init(samples);

	std::ostringstream cutStream;
	cutStream << " " 
          << "misc.MET>=30"                                                  << "&&"
	  << "misc.HT > 800 "                                                << "&&"
	  << "misc.caloHT50_ID >700 "                                        << "&&"
	  << "misc.Jet0Pass ==1"                                             << "&&"
	  << "misc.Jet1Pass ==1"                                             << "&&"
	  << "misc.SecondJPt  >100"                                          << "&&"
	  << "misc.PassJetID ==1"                                            << "&&"
	  << "misc.Vectorsumpt < 70"                                         << "&&"
	  << "NJetsIDLoose >=3"                                              << "&&"
	  << "misc.MinMetJetDPhi >0.3"                                       << "&&"
	  // LowMT2
	  << "misc.LeadingJPt  >150"                                         << "&&"
	  << "NJetsIDLoose >=4"                                              << "&&"
	  << "NBJets >0"                                                     << "&&"
	  // MT2
//	  << "misc.MT2 >200 && misc.MT2 <400"                                << "&&"
//	  << "misc.MT2 > 400"                                                << "&&"
//	  << "misc.MT2 > 80"                                                 << "&&"
	  << "misc.MT2 > 200"                                                 << "&&"
	  // Noise
	  << "misc.HBHENoiseFlag ==0"                                        << "&&"
	  << "misc.CSCTightHaloID==0"                                        << "&&"
	  << "misc.CrazyHCAL==0";

	TString cuts = cutStream.str().c_str();

	//__________________________________________________________________________________________
	// ATTENTION! makePlot adds the underflow and the overflow to the fist and last bin
	//            this is probably not what you want. 
	//            comment out appropriate lines in MassPlotter::makePlot()
	//__________________________________________________________________________________________
	std::ostringstream triggerStream;
	triggerStream << "( "
		<< "(trigger.HLT_HT440_v2==1 && misc.Run<=161176)" << "||"
		<< "(trigger.HLT_HT450_v2==1 && (misc.Run>=161217 && misc.Run<=163261))" << "||"
		<< "(trigger.HLT_HT500_v3==1 && (misc.Run>=163270 && misc.Run<=163869))" << "||"
		<< "(trigger.HLT_HT550_v4==1 && (misc.Run>=165088 && misc.Run<=165633))" << "||"
		<< "(trigger.HLT_HT550_v5==1 && (misc.Run>=165970 && misc.Run<=167043 && misc.Run!=166346))" << "||"
		<< "(trigger.HLT_HT550_v6==1 && (misc.Run==166346))"                     << "||"
		<< "(trigger.HLT_HT550_v7==1 && (misc.Run>=167078 && misc.Run<=167913))" << "||"
		<< "(trigger.HLT_HT550_v8==1 && (misc.Run>=169561 && misc.Run<=173198))" << "||"
		<< "(trigger.HLT_HT600_v1==1 && (misc.Run>=173236 && misc.Run<=173692))" << "||"
		<< "(trigger.HLT_HT600_v1==1 && (misc.Run>=175860 && misc.Run<=177452))" << " )";
  	TString HLT = triggerStream.str().c_str();

	TString ele_cuts = cuts+"&&ele[0].lv.Pt()>10";
	TString muo_cuts = cuts+"&&muo[0].lv.Pt()>10";
	//       samples , variable,                     cuts,      njet, nlep,  HLT, xtitle       , nbins     , bins     , cleaned, log  , comp ,  ratio, stack, overlay
	tA->makePlot("misc.MT2"                        , ele_cuts,  -4, -11 ,  HLT, "M_{T2}"     , 15, 0, 500           ,   false,   true , true,   true,  true,  false , 1);
	tA->makePlot("misc.MT2"                        , muo_cuts,  -4, -13 ,  HLT, "M_{T2}"     , 15, 0, 500           ,   false,   true , true,   true,  true,  false , 1);

	// save number of background events
	Wmuo.Top_bg   =tA->fWpred.Top_bg_mu;
	Wele.Top_bg   =tA->fWpred.Top_bg_e;
	Wmuo.W_bg     =tA->fWpred.W_bg_mu;
	Wele.W_bg     =tA->fWpred.W_bg_e;
	Wmuo.Z_bg     =tA->fWpred.Z_bg_mu;
	Wele.Z_bg     =tA->fWpred.Z_bg_e;
	Wmuo.QCD_bg   =tA->fWpred.QCD_bg_mu;
	Wele.QCD_bg   =tA->fWpred.QCD_bg_e;
	Wmuo.Other_bg =tA->fWpred.Other_bg_mu;
	Wele.Other_bg =tA->fWpred.Other_bg_e;

	// delete pointer
	delete tA;
}

void makePrediction(){
	//                                                                   bg     leptonveto                         
	Wele.nW                        = get_n_events("ele"                , false, false);
	Wele.nW_bg                     = get_n_events("ele"                , true,  false);
	Wele.nW_leptveto               = get_n_events("ele"                , false, true);
	Wele.TopnW                     = get_n_events("eleTop"             , false, false);
	Wele.TopnW_bg                  = get_n_events("eleTop"             , true,  false);
	Wele.TopnW_leptveto            = get_n_events("eleTop"             , false, true);
	if(fIncludeSingleTop){
	Wele.TopnW_s                   = get_n_events("eleTop_s"           , false, false);
	Wele.TopnW_s_bg                = get_n_events("eleTop_s"           , true,  false);
	Wele.TopnW_s_leptveto          = get_n_events("eleTop_s"           , false, true);
	Wele.TopnW_t                   = get_n_events("eleTop_t"           , false, false);
	Wele.TopnW_t_bg                = get_n_events("eleTop_t"           , true,  false);
	Wele.TopnW_t_leptveto          = get_n_events("eleTop_t"           , false, true);
	Wele.TopnW_tW                  = get_n_events("eleTop_tW"          , false, false);
	Wele.TopnW_tW_bg               = get_n_events("eleTop_tW"          , true,  false);
	Wele.TopnW_tW_leptveto         = get_n_events("eleTop_tW"          , false, true);
	}
	Wmuo.nW                        = get_n_events("muo"                , false, false);
	Wmuo.nW_bg                     = get_n_events("muo"                , true,  false);
	Wmuo.nW_leptveto               = get_n_events("muo"                , false, true);
	Wmuo.TopnW                     = get_n_events("muoTop"             , false, false);
	Wmuo.TopnW_bg                  = get_n_events("muoTop"             , true,  false);
	Wmuo.TopnW_leptveto            = get_n_events("muoTop"             , false, true);
	if(fIncludeSingleTop){
	Wmuo.TopnW_s                   = get_n_events("muoTop_s"           , false, false);
	Wmuo.TopnW_s_bg                = get_n_events("muoTop_s"           , true,  false);
	Wmuo.TopnW_s_leptveto          = get_n_events("muoTop_s"           , false, true);
	Wmuo.TopnW_t                   = get_n_events("muoTop_t"           , false, false);
	Wmuo.TopnW_t_bg                = get_n_events("muoTop_t"           , true,  false);
	Wmuo.TopnW_t_leptveto          = get_n_events("muoTop_t"           , false, true);
	Wmuo.TopnW_tW                  = get_n_events("muoTop_tW"          , false, false);
	Wmuo.TopnW_tW_bg               = get_n_events("muoTop_tW"          , true,  false);
	Wmuo.TopnW_tW_leptveto         = get_n_events("muoTop_tW"          , false, true);
	}
	Wele.nData                     = get_n_events("ele_data"           , false, false);
	Wmuo.nData                     = get_n_events("muo_data"           , false, false);
	nTau_leptveto                  = get_n_events("tau_leptveto"       , false, false)*Wele.Wxsection*Wele.lumi/Wele.Wnevents;

	// cout new efficiecies
	if(fMakeEfficienciesPablo ){
		if(fTagandProbeMC){
			double ele_recoeff=0, ele_recoerr=0;
			double muo_recoeff=0, muo_recoerr=0;
			for(int i=0; i<Wele.effMC.size(); ++i){
				ele_recoeff+=Wele.effMC[i];	
				ele_recoerr+=Wele.errMC[i]*Wele.errMC[i];
			}
			for(int i=0; i<Wmuo.effMC.size(); ++i){
				muo_recoeff+=Wmuo.effMC[i];	
				muo_recoerr+=Wmuo.errMC[i]*Wmuo.errMC[i];
			}
			muo_recoeff=muo_recoeff/Wmuo.effMC.size();
			ele_recoeff=ele_recoeff/Wele.effMC.size();
			muo_recoerr=sqrt(muo_recoerr)/Wmuo.effMC.size();
			ele_recoerr=sqrt(ele_recoerr)/Wele.effMC.size();
			cout << "------------------------------" << endl;
			cout << "MC efficiencies:                 " << endl;
			cout << "electrons                     " << endl;
			cout << "W:    reco " << Wele.rec    << " pm " << Wele.rec_err    << endl;
			cout << "WTop: reco " << Wele.Toprec << " pm " << Wele.Toprec_err << endl;
			cout << "T&P:  reco " << ele_recoeff << " pm " << ele_recoerr     << endl;
			cout << "muons                     " << endl;
			cout << "W:    reco " << Wmuo.rec    << " pm " << Wmuo.rec_err    << endl;
			cout << "WTop: reco " << Wmuo.Toprec << " pm " << Wmuo.Toprec_err << endl;
			cout << "T&P:  reco " << muo_recoeff << " pm " << muo_recoerr<< endl;
			cout << "------------------------------" << endl;
		}
		double ele_recoeff=0; double ele_recoerr=0;
		double muo_recoeff=0; double muo_recoerr=0;
		for(int i=0; i<Wele.eff.size(); ++i){
			ele_recoeff+=Wele.eff[i];	
			ele_recoerr+=Wele.err[i]*Wele.err[i];
		}
		for(int i=0; i<Wmuo.eff.size(); ++i){
			muo_recoeff+=Wmuo.eff[i];	
			muo_recoerr+=Wmuo.err[i]*Wmuo.err[i];
		}
		muo_recoeff=muo_recoeff/Wmuo.eff.size();
		muo_recoerr=sqrt(muo_recoerr)/Wmuo.eff.size();
		ele_recoeff=ele_recoeff/Wele.eff.size();
		ele_recoerr=sqrt(ele_recoerr)/Wele.eff.size();
		cout << "Data efficiencies"          << endl;
		cout << "electrons                     " << endl;
		cout << "T&P:  reco " << ele_recoeff << " pm " << ele_recoerr     << endl;
		cout << "muons                     " << endl;
		cout << "T&P:  reco " << muo_recoeff << " pm " << muo_recoerr<< endl;
		cout << "------------------------------" << endl;
	
	}
	makePrint();
}

void makePrint(){

	// __________________________________
	cout << "____________________________________________________"      << endl;
	cout << "W->l nu statistics                                  "      << endl;
	cout << "----------------------------------------------------"      << endl;
	cout << "prob for true W->lv event to be recoed              "      << endl;
	if(!fTopEfficencies){
	cout << " Wele.prob:  " <<Wele.prob 
             << " pm "<<Wele.prob_err_sys() << " (sys)  " 	            << endl;
	}else if(fTopEfficencies){
	cout << " Top Wele.prop: " << Wele.Topprob 
             << " pm "<<Wele.prob_err_sys()   << " (sys)  " 	            << endl;
	}
	if(!fTopEfficencies){
	cout << " Wmunu_prob:  " <<Wmuo.prob 
             << " pm "<<Wmuo.prob_err_sys() << " (sys)  " 	            << endl;
	}else if(fTopEfficencies){
	cout << " Top Wmuo.prop: " << Wmuo.Topprob 
             << " pm "<<Wmuo.prob_err_sys() << " (sys)  " 	            << endl;
	}
	if(fIncludeTop && fTopEfficencies && !fMakeEfficienciesPablo){
	cout << "WARNING: efficiencies for pred taken from Top W->enu"      << endl;
	}else if(! fMakeEfficienciesPablo){
	cout << "WARNING: efficiencies for pred taken from W->enu    "      << endl;
	}else {
	cout << "WARNING: efficiencies from T&P                      "      << endl;
	cout << "Data"                                                      << endl;
	cout << "Wele.TPReco " << Wele.TPReco() << " pm " << Wele.TPRecoErr() << endl;
	cout << "Wmuo.TPReco " << Wmuo.TPReco() << " pm " << Wmuo.TPRecoErr() << endl;
	if(fTagandProbeMC){
	cout << "MC "                                                       << endl;
	cout << "Wele.TPRecoMC " << Wele.TPRecoMC() << " pm " << Wele.TPRecoErrMC() << endl;
	cout << "Wmuo.TPRecoMC " << Wmuo.TPRecoMC() << " pm " << Wmuo.TPRecoErrMC() << endl;
	}}
	cout << "----------------------------------------------------"      << endl;
	cout << "Acc*Reco eff"                                              << endl;
	cout << "Wele.prob() " << Wele.prob() << " pm " << Wele.prob_err_sys() << endl;
	cout << "Wmuo.prob() " << Wmuo.prob() << " pm " << Wmuo.prob_err_sys() << endl;
	cout << "----------------------------------------------------"      << endl;
	cout << "acceptance                                          "      << endl;
	if(fTopEfficencies){
	cout << "Wele.Topacc " << Wele.Topacc << " pm " << Wele.Accerr()     << endl;
	cout << "Wmuo.Topacc " << Wmuo.Topacc << " pm " << Wmuo.Accerr()     << endl;
	}else{
	cout << "Wele.acc " << Wele.acc << " pm " << Wele.Accerr()     << endl;
	cout << "Wmuo.acc " << Wmuo.acc << " pm " << Wmuo.Accerr()     << endl;
	}
	cout << "----------------------------------------------------"      << endl;
	cout << "number of MC W->lnu events                          "      << endl;
	cout << " W->enu:      " << Wele.nW                                 << endl;
	cout << " W->munu:     " << Wmuo.nW                                 << endl;
	cout << "number of true MC W->lnu event (w/ or w/o tau)      "      << endl;
	cout << " nW - nW_bg:  "  << Wele.nW - Wele.nW_bg                   << endl;
	cout << " nW - nW_bg:  "  << Wmuo.nW - Wmuo.nW_bg                   << endl;
	cout << "number of true MC W->lnu events passing l-veto      "      << endl;
	cout << " nW_leptveto: " << Wele.nW_leptveto                        << endl;
	cout << " nW_leptveto: " << Wmuo.nW_leptveto                        << endl;
	if(fIncludeTop){
	cout << "----------------------------------------------------"      << endl;
	cout << "number of MC Top W->lnu events                      "      << endl;
	cout << " Top W->enu:      " << Wele.TopnW                          << endl;
	cout << " Top W->munu:     " << Wmuo.TopnW                          << endl;
	cout << "number of true MC Top W->lnu event (w/ or w/o tau)  "      << endl;
	cout << " Top nW - nW_bg:  "  << Wele.TopnW - Wele.TopnW_bg         << endl;
	cout << " Top nW - nW_bg:  "  << Wmuo.TopnW - Wmuo.TopnW_bg         << endl;
	cout << "number of true MC Top W->lnu events passing l-veto  "      << endl;
	cout << " Top nW_leptveto: " << Wele.TopnW_leptveto                 << endl;
	cout << " Top nW_leptveto: " << Wmuo.TopnW_leptveto                 << endl;
	}
	if(fIncludeSingleTop){
	cout << "----------------------------------------------------"      << endl;
	cout << "number of MC Single Top W->lnu events:s-channel     "      << endl;
	cout << " Top W->enu:      " << Wele.TopnW_s                        << endl;
	cout << " Top W->munu:     " << Wmuo.TopnW_s                        << endl;
	cout << "number of true MC Top W->lnu event (w/ or w/o tau)  "      << endl;
	cout << " Top nW - nW_bg:  "  << Wele.TopnW_s - Wele.TopnW_s_bg     << endl;
	cout << " Top nW - nW_bg:  "  << Wmuo.TopnW_s - Wmuo.TopnW_s_bg     << endl;
	cout << "number of true MC Top W->lnu events passing l-veto  "      << endl;
	cout << " Top nW_leptveto: " << Wele.TopnW_s_leptveto               << endl;
	cout << " Top nW_leptveto: " << Wmuo.TopnW_s_leptveto               << endl;
	cout << "----------------------------------------------------"      << endl;
	cout << "number of MC Single Top W->lnu events:t-channel     "      << endl;
	cout << " Top W->enu:      " << Wele.TopnW_t                        << endl;
	cout << " Top W->munu:     " << Wmuo.TopnW_t                        << endl;
	cout << "number of true MC Top W->lnu event (w/ or w/o tau)  "      << endl;
	cout << " Top nW - nW_bg:  "  << Wele.TopnW_t - Wele.TopnW_t_bg     << endl;
	cout << " Top nW - nW_bg:  "  << Wmuo.TopnW_t - Wmuo.TopnW_t_bg     << endl;
	cout << "number of true MC Top W->lnu events passing l-veto  "      << endl;
	cout << " Top nW_leptveto: " << Wele.TopnW_t_leptveto               << endl;
	cout << " Top nW_leptveto: " << Wmuo.TopnW_t_leptveto               << endl;
	cout << "----------------------------------------------------"      << endl;
	cout << "number of MC Single Top W->lnu events:tW-channel     "     << endl;
	cout << " Top W->enu:      " << Wele.TopnW_tW                       << endl;
	cout << " Top W->munu:     " << Wmuo.TopnW_tW                       << endl;
	cout << "number of true MC Top W->lnu event (w/ or w/o tau)  "      << endl;
	cout << " Top nW - nW_bg:  "  << Wele.TopnW_tW - Wele.TopnW_tW_bg   << endl;
	cout << " Top nW - nW_bg:  "  << Wmuo.TopnW_tW - Wmuo.TopnW_tW_bg   << endl;
	cout << "number of true MC Top W->lnu events passing l-veto  "      << endl;
	cout << " Top nW_leptveto: " << Wele.TopnW_tW_leptveto              << endl;
	cout << " Top nW_leptveto: " << Wmuo.TopnW_tW_leptveto              << endl;
	}
	cout << "----------------------------------------------------"      << endl;
	cout << "MC sanity test                                      "      << endl;
	cout << " electrons:                                         "      << endl;
	cout << "   all events MC     : " << Wele.all_MC()                  << endl;
	cout << "   MC - all bg       : " << Wele.all_MC()-Wele.bg()        << endl;
	cout << "   pred total W->enu : " 
	     <<     (Wele.all_MC()-Wele.bg())/Wele.prob()                   << endl;
	cout << "   pred W lepton veto: " 
	     <<     (Wele.all_MC()-Wele.bg())*(1-Wele.prob())/Wele.prob()   << endl;
	cout << "   true W->enu passing lept veto : "
	     <<     Wele.nW_leptveto_scaled()                               << endl;
	cout << " muons:                                             "      << endl;
	cout << "   all events MC     : " << Wmuo.all_MC()                  << endl;
	cout << "   MC - all bg       : " << Wmuo.all_MC()-Wmuo.bg()        << endl;
	cout << "   pred total W->enu : " 
	     <<     (Wmuo.all_MC()-Wmuo.bg())/Wmuo.prob()                   << endl;
	cout << "   pred W lepton veto: " 
	     <<     (Wmuo.all_MC()-Wmuo.bg())*(1-Wmuo.prob())/Wmuo.prob()   << endl;
	cout << "   true W->enu passing lept veto : "
	     <<     Wmuo.nW_leptveto_scaled()                               << endl;
	cout << "----------------------------------------------------"      << endl;
	cout << "number of events in data                            "      << endl;
	cout << " nWele_data: " << Wele.nData                               << endl;
	cout << " nWmuo_data: " << Wmuo.nData                               << endl;
	cout << "background from other MC                            "      << endl;
	if(!fIncludeTop && !fTopOnly){
	cout << " Top_m   " << Wmuo.Top_bg  << " Top_e   " <<Wele.Top_bg    << endl;
	}
	if(fTopOnly){
	cout << " W_m     " << Wmuo.W_bg    << " W_e     " <<Wele.W_bg      << endl;
	}
	cout << " QCD_m   " << Wmuo.QCD_bg  << " QCD_e   " <<Wele.QCD_bg    << endl;
	cout << " Z_m     " << Wmuo.Z_bg    << " Z_e     " <<Wele.Z_bg      << endl;
	cout << " Other_m " << Wmuo.Other_bg<< " Other_e " <<Wele.Other_bg  << endl;
	if(fIncludeTop && !fTopOnly){
	cout << "bg from Wjets and Top                               "     << endl;
	}
	if(fTopOnly){
	cout << "bg from Top                                         "     << endl;
	}
	if(!fIncludeTop){
	cout << "bg from Wjets                                       "     << endl;
	}
	cout << " muo     " << Wmuo.nW_bg_scaled() 
	     << " ele     " << Wele.nW_bg_scaled()                          << endl;
	cout << "bg subtracted number of event in data               "      << endl;
	cout << " electron: " << Wele.nData-Wele.bg()                       << endl;
	cout << " muon    : " << Wmuo.nData-Wmuo.bg()                       << endl;
	cout << "----------------------------------------------------"      << endl;
	if(fMCClosure){
	cout << "MC CLOSURE TEST "                                          << endl;
	}
	cout << "PREDICTION:                                         "      << endl;
	cout << " number of W->enu passing lepton veto "                           ;
	cout <<    Wele.pred()                                                     ;
	cout << " pm " << Wele.pred_error_stat()       << " (stat)"                ;
	cout << " pm " << Wele.pred_error_sys() << " (sys)"                 << endl;
	cout << " number of W->munu passing lepton veto "                          ;
	cout <<    Wmuo.pred()                                                     ;
	cout << " pm " << Wmuo.pred_error_stat()       << " (stat)"                ;
	cout << " pm " << Wmuo.pred_error_sys()        << " (sys)"          << endl;
	cout << "MC:                                                 "      << endl;
	cout << " MC number of W->enu passing lepton veto  "                       ;
	cout <<    Wele.nW_leptveto_scaled()                                << endl; 
	cout << " MC number of W->munu passing lepton veto "                       ;
	cout <<    Wmuo.nW_leptveto_scaled()                                << endl; 
	cout << "____________________________________________________"      << endl;
	cout << "MC Number of W->tau nu events passing the lept-veto"       << endl;
	if(fIncludeTaus){
	cout << "  leptonic taus *NOT* included!                     "      << endl;
	}else{
	cout << "  leptonic taus included!                           "      << endl;
	}
	cout << "nTau_leptveto " << nTau_leptveto                           << endl;
	cout << "____________________________________________________"      << endl;
//	cout << " ele uncertainty on prob "     << Wele.prob_err/Wele.prob          << endl;
//	cout << "  -> " << (Wele.nData-Wele.bg()) *(1/(Wele.prob()*Wele.prob())*Wele.prob_err)  << endl;
//	cout << " ele uncertainty on prob of  " << Wele.rel_sys_uncert *100 << "%"                                 << endl; 
//	cout << "  -> " << (Wele.nData-Wele.bg()) *(1/(Wele.prob()*Wele.prob())*(Wele.prob*Wele.rel_sys_uncert)) << endl;
//	cout << " ele uncertainty on bg subraction: " << Wele.rel_sys_uncert_bg << " %: " << endl; 
//	cout << "  -> " << Wele.rel_sys_uncert_bg*Wele.bg()*(1-Wele.prob())/Wele.prob()          << endl;
//	cout << " muo uncertainty on prob "     << Wele.prob_err/Wele.prob          << endl;
//	cout << "  -> " << (Wmuo.nData-Wmuo.bg()) *(1/(Wmuo.prob()*Wmuo.prob())*Wmuo.prob_err)  << endl;
//	cout << " muo uncertainty on prob of " << Wmuo.rel_sys_uncert << "%"                                   << endl; 
//	cout << "  -> " << (Wmuo.nData-Wmuo.bg()) *(1/(Wmuo.prob()*Wmuo.prob())*(Wmuo.prob*Wmuo.rel_sys_uncert)) << endl;
//	cout << " muo uncertainty on bg subraction: " << Wmuo.rel_sys_uncert_bg << " %: " << endl; 
//	cout << "  -> " << Wmuo.rel_sys_uncert_bg*Wmuo.bg()*(1-Wmuo.prob())/Wmuo.prob()          << endl;



}
	
double get_n_events(string process, bool background, bool leptveto){
	setStyle();
	TString      basecut   =  " ";
	TChain       *chain    = new TChain("MassTree");
	double       weight    = 1.;
	
	TString path     ="";
	TString path_data="";
	if(fLLSkim) path = "/shome/pnef/SUSY/SUSY_macros/analyzed/RunLeptJetMultAnalyzer/MT2_V00-06-10/20110730_MC_HT300_data_nocuts/SmartSkimmed-LostLepton/";
	else        path = "/shome/pnef/MT2Analysis/MT2trees/MT2_V01-02-01/20111101_data_nocuts_MC_HT400/LostLeptonSkim_MT2gt80/";
	if(fLLSkim) path_data="/shome/pnef/SUSY/SUSY_macros/analyzed/RunLeptJetMultAnalyzer/MT2_V00-06-10/20110730_MC_HT300_data_nocuts/SmartSkimmed-LostLepton2/";
	else        path_data="/shome/pnef/MT2Analysis/MT2trees/MT2_V01-02-01/20111101_data_nocuts_MC_HT400/LostLeptonSkim_MT2gt80/";

	if(process == "ele"){
		chain   ->Add(path+"WJetsToLNu_300_HT_inf_TuneZ2_7TeV-madgraph-tauola_Summer11.root");
		if( background)basecut +="GenLeptFromW(11, 0, 1000, fIncludeTaus)==0  &&";
		if(!leptveto)  basecut +="NEles==1 && NMuons==0 && ele[0].lv.Pt()>10 &&";
		else           basecut +="(NEles==0||ele[0].lv.Pt()<10)&&(NMuons==0||muo[0].lv.Pt()<10) && GenLeptFromW(11, 0, 1000, fIncludeTaus)==1 &&";
	} else if ( process == "eleTop"){
		chain   ->Add(path+"TTJets_TuneZ2_7TeV-madgraph-tauola_Summer11_2.root");
		if( background)basecut +="TopDecayModeResult(11)==0 &&";
		if(!leptveto)  basecut +="NEles==1 && NMuons==0 && ele[0].lv.Pt()>10 &&";
		else           basecut +="(NEles==0||ele[0].lv.Pt()<10)&&(NMuons==0||muo[0].lv.Pt()<10)&&TopDecayModeResult(11)==1 &&";
	} else if ( process == "eleTop_s"){
		chain   ->Add(path+"TToBLNu_TuneZ2_s-channel_7TeV-madgraph_2.root");
		if( background)basecut +="TopDecayModeResult(11)==0 &&";
		if(!leptveto)  basecut +="NEles==1 && NMuons==0 && ele[0].lv.Pt()>10&&";
		else           basecut +="(NEles==0||ele[0].lv.Pt()<10)&&(NMuons==0||muo[0].lv.Pt()<10)&& TopDecayModeResult(11)==1 &&";
	} else if ( process == "eleTop_t"){
		chain   ->Add(path+"TToBLNu_TuneZ2_t-channel_7TeV-madgraph.root");
		if( background)basecut +="TopDecayModeResult(11)==0 &&";
		if(!leptveto)  basecut +="NEles==1 && NMuons==0 && ele[0].lv.Pt()>10&&";
		else           basecut +="(NEles==0||ele[0].lv.Pt()<10)&&(NMuons==0||muo[0].lv.Pt()<10) && TopDecayModeResult(11)==1 &&";
	} else if ( process == "eleTop_tW"){
		chain   ->Add(path+"TToBLNu_TuneZ2_tW-channel_7TeV-madgraph.root");
		if( background)basecut +="TopDecayModeResult(11)==0 &&";
		if(!leptveto)  basecut +="NEles==1 && NMuons==0 && ele[0].lv.Pt()>10&&";
		else           basecut +="(NEles==0||ele[0].lv.Pt()<10)&&(NMuons==0||muo[0].lv.Pt()<10) && TopDecayModeResult(11)==1 &&";
	}else if(process == "muo"){
		chain   ->Add(path+"WJetsToLNu_300_HT_inf_TuneZ2_7TeV-madgraph-tauola_Summer11.root");
		if( background)basecut +="GenLeptFromW(13, 0, 1000, fIncludeTaus)==0  &&";
		if(!leptveto)  basecut +="NEles==0 && NMuons==1 && muo[0].lv.Pt()>10&&";
		else           basecut +="(NEles==0||ele[0].lv.Pt()<10)&&(NMuons==0||muo[0].lv.Pt()<10) && GenLeptFromW(13, 0, 1000, fIncludeTaus)==1 &&";
	} else if ( process == "muoTop"){
		chain   ->Add(path+"TTJets_TuneZ2_7TeV-madgraph-tauola_Summer11_2.root");
		if( background)basecut +="TopDecayModeResult(13)==0 &&";
		if(!leptveto)  basecut +="NEles==0 && NMuons==1 &&muo[0].lv.Pt()>10&&";
		else           basecut +="(NEles==0||ele[0].lv.Pt()<10)&&(NMuons==0||muo[0].lv.Pt()<10) && TopDecayModeResult(13)==1 &&";
	} else if ( process == "muoTop_s"){
		chain   ->Add(path+"TToBLNu_TuneZ2_s-channel_7TeV-madgraph_2.root");
		if( background)basecut +="TopDecayModeResult(13)==0 &&";
		if(!leptveto)  basecut +="NEles==0 && NMuons==1 &&muo[0].lv.Pt()>10&&";
		else           basecut +="(NEles==0||ele[0].lv.Pt()<10)&&(NMuons==0||muo[0].lv.Pt()<10) && TopDecayModeResult(13)==1 &&";
	} else if ( process == "muoTop_t"){
		chain   ->Add(path+"TToBLNu_TuneZ2_t-channel_7TeV-madgraph.root");
		if( background)basecut +="TopDecayModeResult(13)==0 &&";
		if(!leptveto)  basecut +="NEles==0 && NMuons==1 &&muo[0].lv.Pt()>10&&";
		else           basecut +="(NEles==0||ele[0].lv.Pt()<10)&&(NMuons==0||muo[0].lv.Pt()<10) && TopDecayModeResult(13)==1 &&";
	} else if ( process == "muoTop_tW"){
		chain   ->Add(path+"TToBLNu_TuneZ2_tW-channel_7TeV-madgraph.root");
		if( background)basecut +="TopDecayModeResult(13)==0 &&";
		if(!leptveto)  basecut +="NEles==0 && NMuons==1 && &&muo[0].lv.Pt()>10&&";
		else           basecut +="(NEles==0||ele[0].lv.Pt()<10)&&(NMuons==0||muo[0].lv.Pt()<10) && TopDecayModeResult(13)==1 &&";
	} else if ( process == "tau_leptveto"){
		chain   ->Add(path+"WJetsToLNu_300_HT_inf_TuneZ2_7TeV-madgraph-tauola_Summer11.root");
		basecut +="(NEles==0||ele[0].lv.Pt()<10)&&(NMuons==0||muo[0].lv.Pt()<10) && GenLeptFromW(13, 0, 1000, fIncludeTaus)==0 && GenLeptFromW(11, 0, 1000, fIncludeTaus)==0 &&";
	} else if ( process == "ele_data"){
		chain   ->Add(path_data+"/HT-Run2011A-May10ReReco-v1-AOD_3.root");
		chain   ->Add(path_data+"/HT-Run2011A-05Aug2011-v1-AOD_5.root");
		chain   ->Add(path_data+"/HT-Run2011A-PromptReco-v4-AOD_3.root");
		chain   ->Add(path_data+"/HT-Run2011A-PromptReco-v6-AOD_3.root");
		chain   ->Add(path_data+"/HT_Run2011B-PromptReco-v1_AOD.root");
		basecut +="NMuons==0 && NEles==1 && ele[0].lv.Pt()>10 && ";
	} else if ( process == "muo_data"){
		chain   ->Add(path_data+"/HT-Run2011A-May10ReReco-v1-AOD_3.root");
		chain   ->Add(path_data+"/HT-Run2011A-05Aug2011-v1-AOD_5.root");
		chain   ->Add(path_data+"/HT-Run2011A-PromptReco-v4-AOD_3.root");
		chain   ->Add(path_data+"/HT-Run2011A-PromptReco-v6-AOD_3.root");
		chain   ->Add(path_data+"/HT_Run2011B-PromptReco-v1_AOD.root");
		basecut +="NMuons==1 && NEles==0 && muo[0].lv.Pt()>10 &&";
	} else {return -1;}

	basecut  += "misc.MET>=30 && ";
	basecut  += "misc.HT > 800 &&";
	basecut  += "misc.caloHT50_ID >700 &&";
	basecut  += "NJetsIDLoose >=3 && ";
	basecut  += "misc.Jet0Pass==1 && ";
	basecut  += "misc.Jet1Pass==1 && ";
	basecut  += "misc.SecondJPt >100 &&";
	basecut  += "misc.PassJetID ==1 && ";
	basecut  += "misc.Vectorsumpt <70 && ";
	basecut  += "misc.MinMetJetDPhi >0.3 &&";
// Low MT2	
	basecut  += "misc.LeadingJPt >150 &&";
	basecut  += "NBJets >0 &&";
	basecut  += "NJetsIDLoose >=4 && ";
// Mt2 cuts
	basecut  += "misc.MT2 >200 && ";
//	basecut  += "misc.MT2 >80 && ";
//	basecut  += "misc.MT2 >400 && ";
//	basecut  += "misc.MT2 >150 && ";
//	basecut  += "misc.MT2 >100 && misc.MT2<150 && ";
//	basecut  += "misc.MT2 >200 && misc.MT2<400 && ";
//	basecut  += "misc.MT2 >150 && misc.MT2<300 && ";
// Noise
	basecut  += "misc.HBHENoiseFlag == 0 &&";
	basecut  += "misc.CSCTightHaloID == 0 &&";
	basecut  += "misc.CrazyHCAL == 0";



	TString selection;
	if(process == "ele_data" || process == "muo_data"){
		selection = TString::Format("(%f) * (%s)",weight, basecut.Data());
	}else{
		selection = TString::Format("(%.15f*pileUp.Weight) * (%s)",weight, basecut.Data());
	}
 
	double nevents=0;
	if(process == "ele_data" || process == "muo_data" || 
	  (process == "ele"    && !background &&!leptveto&& fTagandProbeMC && !fTopEfficencies)|| 
	  (process == "eleTop" && !background &&!leptveto&& fTagandProbeMC &&  fTopEfficencies)|| 
	  (process == "muo"    && !background &&!leptveto&& fTagandProbeMC && !fTopEfficencies)|| 
	  (process == "muoTop" && !background &&!leptveto&& fTagandProbeMC &&  fTopEfficencies)){
		nevents = GetHistoAndEff(chain, selection, process);
	}else{
		TH1D* h1 = GetHisto(chain, "misc.MET    >>", selection, "h1");	
		cout << "+++ get_n_events: " << process << " background=" << background << " integral " << h1->Integral() << endl; 	
		nevents  = h1->Integral();
		delete h1;
	}


	delete chain;

      	return nevents;	
}

double   GetHistoAndEff(TChain* chainorig, TString basecut, TString process){
	TChain *chain=chainorig->Clone();
	MT2tree *fMT2tree = new MT2tree();
	chain->SetBranchAddress("MT2tree", &fMT2tree);
	Long64_t nentries =  chain->GetEntries();
	Long64_t nbytes = 0, nb = 0;
	int nev =0;
	chain->Draw(">>selList", basecut);

	TEventList *myEvtList = (TEventList*)gDirectory->Get("selList");
	chain->SetEventList(myEvtList);
	int counter=0;
	cout << "+++ Filtering done, size=" <<myEvtList->GetN()  << endl;

	if(myEvtList->GetSize()==0) continue;
	cout << "+++ Passing events for " << process << endl;	    
	while(myEvtList->GetEntry(counter++) !=-1){
		int jentry = myEvtList->GetEntry(counter-1);
		nb =  chain->GetEntry(jentry);   nbytes += nb;
		chain->SetBranchAddress("MT2tree", &fMT2tree);
                  
		if(process=="ele_data"){
			float eff=0; float eff_err=0;
			PabloEffEle(fMT2tree->ele[0].lv.Pt(), fMT2tree->ele[0].lv.Eta(), fMT2tree->LeptJetDR(11,0,0,1), fTypeOfEffErr, eff, eff_err, true);
			Wele.pt .push_back(fMT2tree->ele[0].lv.Pt());
			Wele.eta.push_back(fMT2tree->ele[0].lv.Eta());
			Wele.dR .push_back(fMT2tree->LeptJetDR(11,0,0,1));
			Wele.eff.push_back((double)eff);
			Wele.err.push_back((double)eff_err);
			if(fVerbose>1){
			cout << "++event   " << fMT2tree->misc.Event          << endl;
			cout << "  ele pt  " << fMT2tree->ele[0].lv.Pt()      << endl;
			cout << "  ele eta " << fMT2tree->ele[0].lv.Eta()     << endl;
			cout << "  ele dR  " << fMT2tree->LeptJetDR(11,0,0,1) << endl;
			cout << "  ele eff " << eff << " pm " << eff_err          << endl;
			}
		}else if(process=="ele" || process=="eleTop"){
			float eff=0; float eff_err=0;
			PabloEffEle(fMT2tree->ele[0].lv.Pt(), fMT2tree->ele[0].lv.Eta(), fMT2tree->LeptJetDR(11,0,0,1), fTypeOfEffErr, eff, eff_err, false);
			Wele.ptMC .push_back(fMT2tree->ele[0].lv.Pt());
			Wele.etaMC.push_back(fMT2tree->ele[0].lv.Eta());
			Wele.dRMC .push_back(fMT2tree->LeptJetDR(11,0,0,1));
			Wele.effMC.push_back((double)eff);
			Wele.errMC.push_back((double)eff_err);
			if(fVerbose>1){
			cout << "++event   " << fMT2tree->misc.Event          << endl;
			cout << "  ele pt  " << fMT2tree->ele[0].lv.Pt()      << endl;
			cout << "  ele eta " << fMT2tree->ele[0].lv.Eta()     << endl;
			cout << "  ele dR  " << fMT2tree->LeptJetDR(11,0,0,1) << endl;
			cout << "  ele eff " << eff << " pm " << eff_err          << endl;
			}
		}else if(process=="muo_data"){
			float eff=0; float eff_err=0;
			PabloEffMuo(fMT2tree->muo[0].lv.Pt(), fMT2tree->muo[0].lv.Eta(), fMT2tree->LeptJetDR(13,0,0,1), fTypeOfEffErr, eff, eff_err, true);
			Wmuo.pt .push_back(fMT2tree->muo[0].lv.Pt());
			Wmuo.eta.push_back(fMT2tree->muo[0].lv.Eta());
			Wmuo.dR .push_back(fMT2tree->LeptJetDR(13,0,0,1));
			Wmuo.eff.push_back((double) eff);
			Wmuo.err.push_back((double) eff_err);
			if(fVerbose>1){
			cout << "++event   " << fMT2tree->misc.Event          << endl;
			cout << "  muo pt  " << fMT2tree->muo[0].lv.Pt()      << endl;
			cout << "  muo eta " << fMT2tree->muo[0].lv.Eta()     << endl;
			cout << "  muo dR  " << fMT2tree->LeptJetDR(13,0,0,1) << endl;
			cout << "  muo eff " << eff  << " pm " << eff_err     << endl;
			}
		}else if(process=="muo" || process=="muoTop"){
			float eff=0; float eff_err=0;
			PabloEffMuo(fMT2tree->muo[0].lv.Pt(), fMT2tree->muo[0].lv.Eta(), fMT2tree->LeptJetDR(13,0,0,1), fTypeOfEffErr, eff, eff_err, false);
			Wmuo.ptMC .push_back(fMT2tree->muo[0].lv.Pt());
			Wmuo.etaMC.push_back(fMT2tree->muo[0].lv.Eta());
			Wmuo.dRMC .push_back(fMT2tree->LeptJetDR(13,0,0,1));
			Wmuo.effMC.push_back((double)eff);
			Wmuo.errMC.push_back((double)eff_err);
			if(fVerbose>1){
			cout << "++event   " << fMT2tree->misc.Event          << endl;
			cout << "  ele pt  " << fMT2tree->muo[0].lv.Pt()      << endl;
			cout << "  ele eta " << fMT2tree->muo[0].lv.Eta()     << endl;
			cout << "  ele dR  " << fMT2tree->LeptJetDR(13,0,0,1) << endl;
			cout << "  ele eff " << eff << " pm " << eff_err      << endl;
			}
		}
	}
	return myEvtList->GetN();
}	

TH1D* GetHisto(TChain* chain, TString var, TString basecut, TString name){
	TH1D* h = new TH1D(name,"", 1000, 0, 10000);
	TString varname = var + h->GetName();
	cout << "+++drawing:  " << chain->GetName() << "\n"
	     << "             " << varname << " with cut " << basecut << endl;
	int n = chain->Draw(varname,basecut,"goff");
	return h;
}

void setStyle(){
	gROOT->ProcessLine(".x ~casal/SetStyle_PRD.C");
	gStyle->SetLegendBorderSize(0);
	gStyle->SetFillStyle(0);
	gStyle->SetTextFont(62);
	gStyle->SetTextSize(0.045);
	gStyle->SetPalette(1);
}

void printTable(){

	cout << "*********************************************************************" << endl;
	cout << "\%BEGINLATEX\%"             << endl;
	cout << "\\begin{table}"             << endl
	     << "\\begin{center}"            << endl
             << "\\begin{tabular}{lcc}"      << endl	     
	     << "\\hline\\hline"             << endl
	     << "                               &  $" << "W\\to e \\nu"      << "$ & $" << " W\\to \\mu \\nu$ \\\\"               << endl
	     << "\\hline"                                                                                                         << endl
	     << " $N_{e,\\mu}^{reco}$           &  $" << Wele.nData          << "$ & $" << Wmuo.nData          << "$ \\\\"        << endl
	     << "\\hline"                                                                                                         << endl;
	if(fIncludeTop && !fTopOnly){
	cout << setprecision(4)
	     << " $N_{e,\\mu}^{MC}(W \\& Top)$  &  $" << Wele.nW_scaled()-Wele.nW_bg_scaled()
     	                                                                     << "$ & $" << Wmuo.nW_scaled()-Wmuo.nW_bg_scaled()         
									                                       << "$ \\\\"        << endl
	     << "\\hline"                                                                                                         << endl;
	}
	if(fTopOnly){
	cout << setprecision(4)
	     << " $N_{e,\\mu}^{MC}(Top)$        &  $" << Wele.nW_scaled()-Wele.nW_bg_scaled()
     	                                                                     << "$ & $" << Wmuo.nW_scaled()-Wmuo.nW_bg_scaled()         
									                                       << "$ \\\\"        << endl
	     << "\\hline"                                                                                                         << endl;
		
	}
	if(! fIncludeTop){
	cout << setprecision(4)
	     << " $N_{e,\\mu}^{MC}(W)$          &  $" << Wele.nW_scaled()-Wele.nW_bg_scaled()
     	                                                                     << "$ & $" << Wmuo.nW_scaled()-Wmuo.nW_bg_scaled()         
									                                       << "$ \\\\"        << endl
	     << "\\hline"                                                                                                         << endl;
	}
	if(fIncludeTop && !fTopOnly){
	cout << " $N_{e,\\mu}^{bg}(W \\&Top \\ BG)$ &  $" << Wele.nW_bg_scaled() << "$ & $" << Wmuo.nW_bg_scaled() << "$ \\\\"        << endl;
	}
	if(fTopOnly){
	cout << " $N_{e,\\mu}^{bg}(Top)$        &  $" << Wele.nW_bg_scaled() << "$ & $" << Wmuo.nW_bg_scaled() << "$ \\\\"        << endl
	     << " $N_{e,\\mu}^{bg}(W)  $        &  $" << Wele.W_bg           << "$ & $" << Wmuo.W_bg           << "$ \\\\"        << endl;
	}
	if(!fIncludeTop){
	cout << " $N_{e,\\mu}^{bg}(W)$          &  $" << Wele.nW_bg_scaled() << "$ & $" << Wmuo.nW_bg_scaled() << "$ \\\\"        << endl
	     << " $N_{e,\\mu}^{bg}(Top)$        &  $" << Wele.Top_bg         << "$ & $" << Wmuo.Top_bg         << "$ \\\\"        << endl;
	}
	cout << " $N_{e,\\mu}^{bg}(Z)$          &  $" << Wele.Z_bg           << "$ & $" << Wmuo.Z_bg           << "$ \\\\"        << endl
	     << " $N_{e,\\mu}^{bg}(QCD)$        &  $" << Wele.QCD_bg         << "$ & $" << Wmuo.QCD_bg         << "$ \\\\"        << endl
	     << " $N_{e,\\mu}^{bg}(Other)$      &  $" << Wele.Other_bg       << "$ & $" << Wmuo.Other_bg       << "$ \\\\"        << endl
	     << "\\hline"                                                                                                         << endl;
	cout << " $\\varepsilon      $          &  $" << Wele.prob()         << "$ & $" << Wmuo.prob()       << "$ \\\\"        << endl;
	cout << "\\hline"                                                                                                         << endl
	     << setprecision(3) 
	     << " $W_{e,\\mu}^{pass \\ veto}$ &  $"   << Wele.pred()        
	                                              << " \\pm "             << Wele.pred_error_stat()   
						      << " \\ (stat) \\ \\pm "<< Wele.pred_error_sys()          << " \\ (sys) \\"                 
	                                              << " $ & $"             << Wmuo.pred()        
	                                              << " \\ \\pm "          << Wmuo.pred_error_stat()   
						      << " \\ (stat) \\ \\pm "<< Wmuo.pred_error_sys()          << "\\ (sys) $ \\\\"<< endl   
	     << "\\hline"                                                                                                        << endl
	     << setprecision(4) 
	     << " $W_{e,\\mu}^{pass \\ veto}$ MC & $"<< Wele.nW_leptveto_scaled() << "$ & $" << Wmuo.nW_leptveto_scaled() <<  " $ \\\\ "   << endl 					      
	     << "\\hline\\hline"                                                                                                << endl
	     << "\\end{tabular}"                                                                                                << endl
	     << "\\end{center}"                                                                                                 << endl
	     << "\\end{table}"                                                                                                  << endl
	     << "\%ENDLATEX\%"                                                                                                  << endl
	     << endl;

//	cout << "\%BEGINLATEX\%"                                                                                                << endl
//	     << "\\begin{table}"                                                                                                << endl
//	     << "\\begin{center}"                                                                                               << endl
//             << "\\begin{tabular}{lcc}"                                                                                         << endl	     
//	     << "\\hline\\hline"                                                                                                << endl;
//	if(fIncludeTaus){
//	cout << "                                      &  $" << " W\\to \\tau \\nu, \\ \\tau \\to \\ hadrons$ \\\\"             << endl;
//	}else{
//	cout << "                                      &  $" << " W\\to \\tau \\nu, \\ \\tau \\to \\ anything$  \\\\"           << endl;
//	}	
//	cout << " MC  events  passing  l-veto          &  $" << nTau_leptveto << "$ \\\\"                                       << endl
//	     << "\\hline\\hline"                                                                                                << endl
//	     << "\\end{tabular}"                                                                                                << endl
//	     << "\\end{center}"                                                                                                 << endl
//	     << "\\end{table}"                                                                                                  << endl
//	     << "\%ENDLATEX\%"                                                                                                  << endl
//	     << endl;
//	cout << "*********************************************************************" << endl;
		     

}

void readEfficienciesPablo(){	
	// electros

	cout << "--------------------------------------" << endl;
	cout << "Read Pablo EffTrees " << endl;	
  	TFile* fele   = TFile::Open(fPablofile_ele_data);
  	TFile* feleMC = TFile::Open(fPablofile_ele_mc);
  	if (!fele || !feleMC) {
		cout << "ERROR : Could not open file " << fPablofile_ele_data << " or " fPablofile_ele_mc<< "!"  << endl;
	    	exit(1);
  	}
  
  	fele->GetObject("histo_electron_var1pt_var2dr_data",WeleEffHistdR) ;
  	fele->GetObject("histo_electron_var1pt_var2eta_data",WeleEffHistEta);
  	feleMC->GetObject("histo_electron_var1pt_var2dr_mc",WeleEffHistdRMC) ;
  	feleMC->GetObject("histo_electron_var1pt_var2eta_mc",WeleEffHistEtaMC);
	if(WeleEffHistdR==NULL || WeleEffHistEta==NULL || WeleEffHistdRMC==NULL || WeleEffHistEtaMC==NULL) {cout << " cannot get hist" << endl; exit(-1);}

	// muons
	TFile* fmuo   = TFile::Open(fPablofile_muo_data);
	TFile* fmuoMC = TFile::Open(fPablofile_muo_mc);
  	if (!fmuo || !fmuoMC) {
		cout << "ERROR : Could not open file " << fPablofile_muo_data << " or " << fPablofile_muo_mc << "!"  << endl;
	    	exit(1);
  	}
  
  	fmuo  ->GetObject("histo_muon_var1pt_var2dr_data" ,WmuoEffHistdR) ;
  	fmuo  ->GetObject("histo_muon_var1pt_var2eta_data",WmuoEffHistEta);
  	fmuoMC->GetObject("histo_muon_var1pt_var2dr_mc" ,WmuoEffHistdRMC) ;
  	fmuoMC->GetObject("histo_muon_var1pt_var2eta_mc",WmuoEffHistEtaMC);
	if(WmuoEffHistdR==NULL || WmuoEffHistEta==NULL || WmuoEffHistdRMC==NULL || WmuoEffHistdRMC==NULL) {cout << " cannot get hist" << endl; exit(-1);}

	cout << " done"                                  << endl;
	cout << "--------------------------------------" << endl;
}

void PabloEffEle(float pt_, float eta_, float dr_, int typeOfError, float &eff, float &err, bool data) {
	if(! fMakeEfficienciesPablo) return;

	float efficiency,  e_efficiency1, e_efficiency2;
	float pt, eta, dr, njets;
	float ptwidth, etawidth, drwidth, njetswidth;
	TH2F* H_eta=(data)? WeleEffHistEta->Clone(): WeleEffHistEtaMC->Clone();
	TH2F* H_dR =(data)? WeleEffHistdR ->Clone(): WeleEffHistdRMC ->Clone();
	if(pt_ < 20){
		int bin = 	H_eta->FindBin(pt_, eta_);
		if(bin!=1){ 
			eff = H_eta->GetBinContent( bin  );
			err = H_eta->GetBinError( bin  );
		}
		else exit(-1);
	}
	else{
		int bin = -1;
		float tmpDr=dr_, tmpPt=pt_;
		// for leptons with dR<0.5, I'm taking dR=0.5. Same for pt
		if(dr_ <0.5) tmpDr = 0.5;
		if(pt_ >200) tmpPt = 199;
		bin = H_dR->FindBin(tmpPt, tmpDr );
		if(bin!=-1){ 
			eff = H_dR->GetBinContent( bin  );
			err = H_dR->GetBinError( bin  );
		}
		else exit(-1);
	}

	if(eff==0) {cout << "eff not found " << endl; exit(1);}

}

void PabloEffMuo(float pt_, float eta_, float dr_ , int typeOfError, float &eff, float &err, bool data) {
	if(! fMakeEfficienciesPablo) return;

	float efficiency,  e_efficiency1, e_efficiency2;
	float pt, eta, dr, njets;
	float ptwidth, etawidth, drwidth, njetswidth;
	TH2F* H_eta=(data)? WmuoEffHistEta->Clone(): WmuoEffHistEtaMC->Clone();
	TH2F* H_dR =(data)? WmuoEffHistdR ->Clone(): WmuoEffHistdRMC ->Clone();


	if(pt_ < 20){
		int bin = 	H_eta->FindBin(pt_, eta_);
		if(bin!=1){ 
			eff = H_eta->GetBinContent( bin  );
			err = H_eta->GetBinError( bin  );
		}
		else exit(-1);
	}
	else{
		int bin = -1;
		float tmpDr=dr_;
		float tmpPt=pt_;
		// for leptons with dR<0.5, I'm taking dR=0.5. Same for pt
		if(dr_ <0.5) tmpDr = 0.5;
		if(pt_ >200) tmpPt = 199;
		bin = H_dR->FindBin(tmpPt, tmpDr );
		if(bin!=-1){ 
			eff = H_dR->GetBinContent( bin  );
			err = H_dR->GetBinError( bin  );
		}
		else exit(-1);
	}

	if(eff==0) {cout << "eff not found " << endl; exit(1);}

}
