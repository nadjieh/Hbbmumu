#ifndef MT2Analysis_hh
#define MT2Analysis_hh


#include <cmath>
#include <string>
#include <iostream>
#include <fstream>

#include <numeric>

#include "base/UserAnalysisBase.hh"
#include "base/TreeReader.hh"
#include "helper/Davismt2.h"
#include "helper/TMctLib.h"
#include "helper/Hemisphere.hh"
#include "helper/EventFilterFromListStandAlone.h"
#include "BTagSFWeight.hh"
#include "MT2tree.hh"
#include <TLorentzVector.h>
#include "JetCorrectionUncertainty.h" 
#include "JetCorrectorParameters.h" 
#include "FactorizedJetCorrector.h" 

#include "TH2F.h"
#include "TTree.h"

// #include </shome/leo/Installations/LHAPDF/lhapdf-5.8.4/include/LHAPDF/LHAPDF.h>


static const int gNHemispheres = 4;
static const int gNGenJets     = 20;

class MT2AnalysisJet;


// MT2Analysis Class --------------------------------------------------------------
class MT2Analysis: public UserAnalysisBase{
public:
	MT2Analysis(TreeReader *tr = NULL);
	virtual ~MT2Analysis();
	friend class MT2AnalysisJet;

	void Begin(const char* filename = "Mass_histos.root");
	void Analyze();
	void End();
	void SetType(bool isData=false){fisData=isData;};
	void SetProcessID(int ID){fID=ID;};
	void SetBTagEfficiency(string btagFileName){ fbtagFileName = btagFileName;};
	void SetHadTauEfficiency(string hadtauFileName){ fhadtauFileName = hadtauFileName;};
	void SetType1MET(bool type1MET){fisType1MET = type1MET;};
	void SetCHSJets(bool CHSJets){fisCHSJets = CHSJets;};
	void SetFastSim(bool FastSim){fisFastSim = FastSim;};


  //new 
  void GetElectronIndices();
  void GetMuonIndices();
  void GetTauIndices();
  void FillMT2Elecs();
  void FillMT2Muons();
  void FillMT2Taus();
  //

	// redo JEC
	void SetJEC(string JEC){fJEC=JEC;};
	// parser
	void ReadCuts(const char* SetofCuts);

	// PU
  	enum PUScenario {noPU, MC2012};
	PUScenario fPUScenario;
	void SetPUReweighting(string PU, string data_PileUp, string mc_PileUp){
		if      (PU =="MC2012") {fPUScenario=MC2012;   SetPileUpSrc  (data_PileUp, mc_PileUp);}
		else                     fPUScenario=noPU;
	};

 	 //create pdf weights
        bool doPDF;

  	//is a susy scan?
  	bool isScan;

  	//is a fastsim sample?
  	bool fisFastSim;

  	// remove Photon / Z->ll
  	bool fRemovePhoton;
  	bool fRemoveZll;

	// CHSjets
	bool fisType1MET;

	// Type1 MET
	bool fisCHSJets;
  
  	//Control histos
  	TH1F *fH_PUWeights, *fH_Events ;
  	TH2F *fH2_SMSEvents, *fH2_mSugraEvents;
  	//SUSY subprocess histos
  	TH2F *fH_mSugraSubProcEvents[11];

private:
	// files and trees ----------------------------------------------------------------
	// file for histograms:
	TFile* fHistFile;

        // MT2 mini-tree
        MT2tree* fMT2tree;
        TTree* fATree;

        //Control Histograms
  	//  TH1F *fH_PUWeights ;

	//btagging histograms and files
	TFile *btagfile;
	TH2F  *hbeff;
	TH2F  *hceff;
	TH2F  *hleff;

	//hadronic tau tagging histograms and files
	TFile *hadtaufile;
	TH1F  *htaueff;
	TH1F  *hjeteff;
	TH1F  *heleeff;
	TH1F  *hmuoeff;

        // HCAL Laser standalone filter
         EventFilterFromListStandAlone *myFilter;

	// data members-----------------------------------------
	bool fisData;
	int  fID;
	int  fCounter;
	bool fBasicMT2treeFilled;
	string fJEC;
	string fbtagFileName;
	string fhadtauFileName;


	// vectors for jets and lepton indices
	vector<int> fElecs;
	vector<int> fPhotons;
	vector<int> fTaus;
	vector<int> fMuons;
	vector<bool>fPhotonJetOverlapRemoved;
	// Jets
	vector<MT2AnalysisJet> Jets;
	
	// MT2 and hemisphere
	Davismt2 *fMT2;
	TMctLib  *fMCT;
	Hemisphere *fHemisphere;

  	//PDFs
  	int nPDFs;

	// cut variables
	float fHT;
	float fCaloHT50;
	float fCaloHT50_ID;
	float fCaloMHT30;
	float fCaloMHT30_ID;
	bool  fCrazyHCAL;
	bool  fIsNANObj;
	bool  fNegativeJEC;

	
	//  ---- set of cuts ---
	TString fSetName;
	float fCut_PFMET_min;
	float fCut_HT_min;
	float fCut_caloHT50_min;
	float fCut_caloHT50ID_min;
	float fCut_caloMHT30_min;
	float fCut_caloMHT30ID_min;
	float fCut_JPt_hardest_min;
	float fCut_JPt_second_min;
	float fCut_PtHat_max;
	int   fCut_Run_min;
	int   fCut_Run_max;
	bool  fDoJESUncertainty;
	int   fJESUpDown;
  	int   fCut_NJets40_min;
	int   fCut_NLeptons_min;//only muons and electrons


	// ---- required and vetoed triggers ----
	std::vector<std::string> fRequiredHLT; 
	std::vector<std::string> fVetoedHLT;
	
	typedef std::map <string, bool*> StringBoolMap;
	StringBoolMap fTriggerMap;
	


	// member functions -------------------------------------------------------------
	void BookTree();
	void FillTree();
	void ResetTree(); 
	bool FillMT2TreeBasics();
	void FillMT2treeCalculations();
	void GetLeptonJetIndices();
	bool IsSelectedEvent();
	void InitializeEvent();
	// photons
	bool IsGoodPhotonIsoLoose(int i);
	bool IsGoodPhotonIsoMedium(int i);
	bool IsGoodPhotonIsoTight(int i);
	bool IsGoodPhotonEGMLoose(int i);
	bool IsGoodPhotonEGMMedium(int i);
	bool IsGoodPhotonEGMTight(int i);
	bool IsGoodPhoton(int i);
	float SingleTowerHoverE(int i); // want photon.hadTowOverEm() which we do not store directly
	const float EffAreaChargedHad(float abseta);
	const float EffAreaNeutralHad(float abseta);
	const float EffAreaPhoton(float abseta);

        //Taus
	bool IsGoodTauLooseID(int i);
        bool IsGoodTau(int index);
	//Muons 
	bool IsGoodMT2Muon(const int index);
	float MuPFIso(const int index);
	float MuPFIso04(const int index);
	// Electrons
	bool IsGoodMT2ElectronVetoID(const int index);
  //bool IsGoodMT2ElectronLooseID(const int index);
  //bool IsGoodMT2ElectronMediumID(const int index);

  //bool IsGoodMT2ElectronVetoIDforEleEle(const int index);
	bool IsGoodMT2ElectronSelIDforEleEle(const int index);

	bool IsGoodMT2ElectronVetoIDforEleMu(const int index);
	bool IsGoodMT2ElectronSelIDforEleMu(const int index);

	bool IsGoodMT2ElectronVetoIDforEleTau(const int index);
	bool IsGoodMT2ElectronSelIDforEleTau(const int index);

        bool IsGoodMT2ElectronVetoIDforMuTau(const int index);
        bool IsGoodMT2ElectronVetoIDforTauTau(const int index);

        bool IsGoodMT2ElectronMVANoTrigLoose(const int index);
        bool IsGoodMT2ElectronMVANoTrigTight(const int index);



	const float EffArea(float abseta);
	const float EffArea04(float abseta);
	float ElePFIso(const int index);
	float ElePFIso04(const int index);

	//pfJetID
	bool IsGoodMT2PFJetIDLoose(int index, float ptcut, float absetacut);
	bool IsGoodMT2PFJetIDMedium(int index, float ptcut, float absetacut);
	bool IsGoodMT2PFJetIDTight(int index, float ptcut, float absetacut);
	// jets + MET

        TLorentzVector CorrToMETFromJetSmearing; //saeid 
	TLorentzVector CAJet(int index);
	TLorentzVector MET();
	void Initialize_JetCorrectionUncertainty();
	void Initialize_JetEnergyCorrection();
	JetCorrectionUncertainty *fJecUncPF;   
	JetCorrectionUncertainty *fJecUncCalo; 

	JetCorrectorParameters* fJecCaloL2;
	JetCorrectorParameters* fJecCaloL3; 
	JetCorrectorParameters* fJecCaloRES; 
	JetCorrectorParameters* fJecPFL1;   
	JetCorrectorParameters* fJecPFL2;   
	JetCorrectorParameters* fJecPFL3;   
	JetCorrectorParameters* fJecPFRES;  
	JetCorrectorParameters* fJecrawPFL1;   
	JetCorrectorParameters* fJecrawPFL2;   
	JetCorrectorParameters* fJecrawPFL3;   
	JetCorrectorParameters* fJecrawPFRES;  

	FactorizedJetCorrector* fJetCorrectorPF;
	FactorizedJetCorrector* fMetCorrectorPF;
	FactorizedJetCorrector* fJetCorrectorCalo;

};

// MyJet Helper Class -----------------------------------
class MT2AnalysisJet: public MT2Jet {
public:
	MT2AnalysisJet(int index, TString type, MT2Analysis *ana);
	virtual ~MT2AnalysisJet(){};
	float Pt();
	float Eta();
	float Phi();
	float E();
	float M();

	TString fType;
private:
	TLorentzVector PFJetScaled(TLorentzVector jraw, float area, float rho, int level);
	TLorentzVector CAJetScaled(TLorentzVector jraw, int level);
	float GetPFJEC(float rawpt, float eta, float area, float rho, int level);
	float GetCaloJEC(float rawpt, float eta, int level);
	float GetJECUncertPF(float pt, float eta);
	float GetJECUncertCalo(float pt, float eta);

	MT2Analysis* fAna;

};

#endif
