#ifndef MT2Analyzer_hh
#define MT2Analyzer_hh

#include <TTree.h>
#include <TStyle.h>

#include "base/TreeAnalyzerBase.hh"
#include "base/TreeReader.hh"
#include "MT2Analysis.hh"


class MT2Analyzer : public TreeAnalyzerBase {
public:
	MT2Analyzer(std::vector<std::string>& fileList);
	virtual ~MT2Analyzer();
	void BeginJob(TString filename="MassTree.root" , TString setofcuts="default",
	              bool isData=false, string data_PileUp="", string mc_PileUp="", string JEC="");
	void EndJob();
	void Loop();
	void SetMaxEvents(int a){fMaxEvents=a;};
	void SetProcessID(int ID){fID=ID;};
	void SetBTagEfficiency(string btagFileName){ fbtagFileName = btagFileName;};
	void SetHadTauEfficiency(string hadtauFileName){ fhadtauFileName = hadtauFileName;};
	void SetPUReweighting(string puScenario){fPu = puScenario;};
	void SetType1MET(bool type1MET){fType1MET = type1MET;};
	void SetCHSJets(bool CHSJets){fCHSJets = CHSJets;};
	void SetFastSim(bool FastSim){fisFastSim = FastSim;};
	bool removePhoton;
	bool removeZll;
        bool doPDF;
        bool isScan;
        bool fisFastSim;
	string fPu;
	bool fType1MET;
	bool fCHSJets;
private:
	MT2Analysis             *fMT2Analysis;
  	int fMaxEvents;   
	int fID;
	string fbtagFileName;
	string fhadtauFileName;
};
#endif
