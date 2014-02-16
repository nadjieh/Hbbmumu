#ifndef TopSearch_h
#define TopSearch_h

// Search for top in events

#include <cmath>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
// #include <iostream.h>
#include <TLorentzVector.h>

using namespace std;
	

// ************************************************************
// class TopSearch

// Tries to reconstruct Top candidates by:
//
// first reconstructing Ws from either 2 jets, 1 heavy jet (if recoW1j = 1) or leptonic decays (if recoWlept = 1)
//	if a heavy jet is in a given mass range and is b-tagged, it is considered as a "W1b" = coalescence of a b with a jet from W
//	solutions are kept only if they have a chisquared < a maximum for W in 2 jets (W2j) and in 1 heavy jet (W1j)
//	in addition, if there is a heavy b-tagged jet (W1b) or a single lepton (Wlept) a penalty on the chisq can be given
//		(this is to avoid that these solution overrule systematically the ones with a real chisquared in top reconstruction)
//	in case of ambiguities (overlapping Ws) the best chisquared solution is kept (if ckWoverlap = 1)
//		if ckWoverlap = 0, all solutions will be tried for top reconstruction
//
// then reconstructing Tops from Ws, adding a free jet
//	in case of ambiguities (overlapping Tops) the best top can be chosen as:
//	(can be defined by method SetWithbtag)
//		withbtag = 0: the best chisquared solution
//		withbtag = 1: the one with a correct b-tagged jet (or a W1b) is preferred over the best chisquared (default)
//		withbtag = 2: tops are only reconstructed if the additional jet is b-tagged or there is a W1b
//		withbtag = 3: At least 1 bJet out of the top
// at the end, rejected solutions are removed from the lists ipT and ipW (see below) to provide a coherent interpretation
//
// there are many constants and flags to steer the logic and all of them can be overwritten by a Set method
//
// the event data and the results are stored into several arrays held in the namespace evtData:
//
// evt[i][j] is expected to be input to TopSearch via the namespace evtData
// is defined for each particle i (i < npart) as:
// j = 0: identifier = 1 for jets, 2 for muons, 3 for electrons, 10 for MET
//		10 = MET is created only in pzmissfrW, called from GetWcand, for leptonic W decays
//   = 1: px
//   = 2: py
//   = 3: pz
//   = 4: E
//   = 5: M
//   = 6: pt
//   = 7: eta
//   = 8: phi
//   = 9: btagged (=1) or not (=0)
//		may be incremented in EvtFinal depending on its origin:
//		- correct b from the top, b or W1b (+10)
//		- correct b from the top, b and W1b (+20)
//		- b-tagged decay product of the W in 2j from top (+30)
//		- b-tagged decay product of the W in 1j from top (+40)
//		- b-tagged decay product of lonely W in 2j (+50)
//		- b-tagged decay product of lonely W in 1j (+60)
//		- b-tagged lonely jet (+70)
//
//   = 10: lable in the JetCollection
// Also returned in the namespace evtData are the results:
//
// ipT[i][j] contains the reconstructed Tops (i < npT)
// i = order number of the Top in the list
// j = 0 order number of the W in the list ipW
// j = 1 order number of the 2nd decay particle in the list evt
// j = 2 number of the bjets in a top
//
// ipW[i][j] contains the reconstructed Ws (i < npW)
// i = order number of the W in the list (i < npW)
// j = 0 order number of the 1st decay particle in the list evt
//		may be a completed Missing vector in case of leptonic W decay
//		order number +1000 if W is from Top
// j = 1 order number of the 2nd decay particle in the list evt (W2j)
//		or -999 if W from a single heavy jet compatible with W mass (W1j)
//		or -990 if W from a single heavy jet with b-tag (W1b)
//
// ileft[i] array with tags for all particles (i < npart)
// = 0 if included in a remaining top or W
// = 1 if not used in top nor W (= leftover)
// ************************************************************

class TopSearch {

public:

// Constructor

TopSearch ();

// Destructor

~TopSearch();

	// reconstruction of Ws and tops
	int GoSearch ();
	
	// get the reconstructed Ws
	int GetRecoW (int iW[][2], double minvW[], double chi2W[]) {
		for (int i = 0; i < npW; ++i) {
			iW[i][0] = ipW[i][0];
			iW[i][1] = ipW[i][1];
			minvW[i] = mW[i];
			chi2W[i] = chisqW[i];
		}
			return npW;
	};
	
	// get the reconstructed Tops
	int GetRecoT (int iT[][3], double minvT[], double chi2T[]) {
		for (int i = 0; i < npT; ++i) {
			iT[i][0] = ipT[i][0];
			iT[i][1] = ipT[i][1];
			iT[i][2] = ipT[i][2];
			minvT[i] = mT[i];
			chi2T[i] = chisqT[i];
		}
			return npT;
	};
	
	// set or get the debug level:
	// = 0 if no printout at all
	// = 1 if only final statistics
	// = 2 if also short event per event output
	// = 3 if also longer event per event output
	void SetDebug(int dbg)  {
		debug = dbg;
		if (debug >= 1) cout << " Update: Debug level = " << debug << endl;
	}; 
	int GetDebug(void)  {return debug;}; 
	
	// Define constants for the analysis:
	
	// decide whether to keep only particles within acceptance
	// if = 0: keep all particles independent of acceptance
	// if = 1: keep only particles within eta acceptance (=default)
	void SetInAcc(int in) {
		inAcc = in;
		if (debug >= 1) cout << " Update: Particles only if within eta acceptance" << endl;
	}; 
	
	// to select top, can use chisquared and/or correct b-tag
	// = 0: keep ambiguous solution with best chisquared
	// = 1: btag preferred over chisq
	// = 2: correct b-jet requested for top reconstruction
	void SetWithbtag(int withb) {
		withbtag = withb;
		if (debug >= 1) {
			if (withbtag == 0) cout << " Update: Tops ordered by chisquared only" << endl;
			else if (withbtag == 1) cout << " Update: Correct b-tag for the top preferred over best chisq" << endl;
			else if (withbtag == 2) cout << " Update: Correct b-tag required for the top" << endl;
			else if (withbtag == 3) cout << " Update: At least 1 bJet out of the top" << endl;
		}
	}; 
	
	// to reconstruct Ws with 1 jet or not
	// = 0: do not reconstruct W from 1 jet,
	// = 1: reconstruct W from 1 jet
	void SetRecoW1j(int rW1j) {
		recoW1j = rW1j;
		if (debug >= 1) {
			cout << " Update: W1j reco = " << recoW1j << endl;
		}
	};
	
	// to reconstruct W leptonic decays or not
	// = 0: do not reconstruct W leptonic,
	// = 1: reconstruct W from leptonic
	void SetRecoWlept(int rWlept) {
		recoWlept = rWlept;
		if (debug >= 1) {
			cout << " Update: Wlept reco = " << recoWlept << endl;
		}
	};
	
	// to check W overlaps or not
	// = 0: do not check W overlaps
	// = 1: check W overlaps
	void SetCkWoverlap(int ckWovl) {
		ckWoverlap = ckWovl;
		if (debug >= 1) {
			cout << " Update: Wlept reco = " << ckWoverlap << endl;
		}
	};
	
	// to rescale or not the W 4-vector
	// = 0: no rescaling done,
	// = 1 W 4-vector to be rescaled to the correct W mass
	void SetRescaleW(int scal) {
		rescaleW = scal;
		if (debug >= 1) {
			cout << " Update: Whad 4-vector rescaled = " << rescaleW << endl;
		}
	};

	
	// set maximum value of chisquared to accept a W (default = 7)
	void SetChisqWmax(double mychisqW) {
		chisqWmax = mychisqW;
		if (debug >= 1) {
			cout << " Update: Chisq cut: W = " << chisqWmax << endl;
		}
	};
	
	// set maximum value of chisquared to accept a top (default = 2)
	void SetChisqTmax(double mychisqT) {
		chisqTmax = mychisqT;
		if (debug >= 1) {
			cout << " Update: Chisq cut: top = " << chisqTmax << endl;
		}
	};
	
	
	// set minimum and maximum values of the mass range for a heavy jet to be b+j from a top (default = 40 - 130)
	void SetmassrangeW1b(double mW1bmin, double mW1bmax) {
		massWjbmin = mW1bmin;
		massWjbmax = mW1bmax;
		if (debug >= 1) {
			cout << " Update: W1b mass range: min = " << massWjbmin << ", max = " << massWjbmax << endl;
		}
	};
	
	// Set chisquared penalty for W1b
	void SetChisqWlept(double chi2Wlept) {
		chisqWlept = chi2Wlept;
		if (debug >= 1) {
			cout << " Update: Chisq penalty for W leptonic = " << chisqWlept << endl;
		}
	};
	
	// Set chisquared penalty for W leptonic
	
	// Set chisquared penalty for W1b
	void SetChisqW1b(double chi2W1b) {
		chisqW1b = chi2W1b;
		if (debug >= 1) {
			cout << " Update: Chisq penalty for W1b = " << chisqW1b << endl;
		}
	};
	
	// set maximum value of chisquared to accept a top with leptonic W decay (default = 0.5)
	void SetChisqTWlepmax(double chisqTWlep) {
		chisqTWlepmax = chisqTWlep;
		if (debug >= 1) {
			cout << " Update: Chisq cut: top (W lept) = " << chisqWmax << endl;
		}
	};
	
	// define values for the W and top masses, width, precision
	
	// Delta p / p constant for error propagation (default = 100% = 1.)
	void SetDpbypnom(double dpbyp) {
		dpbypnom = dpbyp;
		if (debug >= 1) {
			cout << " Update: Momentum resolution constant = " << dpbypnom << endl;
		}
	};
	
	// mass of W (default = 80.4 GeV)
	void SetMassW(double mW1) {
		massW = mW1;
		if (debug >= 1) {
			cout << " Update: W mass = " << massW << endl;
		}
	};
	
	// width of W decaying into 2j (default = 2.1 GeV)
	void SetWidthW2j(double gamW2j) {
		widthW2j = gamW2j;
		if (debug >= 1) {
			cout << " Update: W width 2j = " << widthW2j << endl;
		}
	};
	
	// width of W decaying into 1j (default = 10. GeV)
	void SetWidthW1j(double gamW1j) {
		widthW1j = gamW1j;
		if (debug >= 1) {
			cout << " Update: W width 1j = " << widthW1j << endl;
		}
	};
	
	// mass of top (default = 172 GeV)
	void SetMassT(double mT1) {
		massT = mT1;
		if (debug >= 1) {
			cout << " Update: T mass = " << massT << endl;
		}
	};
	
	// width of top (default = 10. GeV)
	void SetWidthT(double gamT) {
		widthT = gamT;
		if (debug >= 1) {
			cout << " Update: T width = " << widthT << endl;
		}
	};

  int Npart(){
    return npart;
  }

  void SetNpart(int i){
    npart = i;
  }

  double Evt(int i, int j){
    return evt[i][j];
  }

  void SetEvt(int i, int j, double myEvt){
    evt[i][j] = myEvt;
  }

  double Mt2(){
    return MT2;
  }
  
  void SetMt2(double i){
    MT2 = i;
  }

  double Met(){
    return MET;
  }
  
  void SetMet(double i){
    MET = i;
  }
  
  double Ht(){
    return HT;
  }

  void SetHt(double i){
    HT = i;
  }

  int Npsjet(int i){
    return npsjet[i];
  }
  
  void SetNpsjet(int i, int j){
    npsjet[i] = j;
  }

  void SetIpsjet(int i, int j, int myIpsjet){
    ipsjet[i][j] = myIpsjet;
  }

  int NpT(){
    return npT;
  }

  void SetNpT(int i){
    npT = i;
  }

  int IpT(int i, int j){
    if(j == 1)
      return evt[ipT[i][j]][10];
    else
      return ipT[i][j];
  }

  /* void SetIpT(int i, int j, int IPT){ */
  /*   ipT[i][j] = IPT; */
  /* } */

  double Mt(int i){
    return mT[i];
  }

  void SetMt(int i, double j){
    mT[i] = j;
  }
  
  double ChisqT(int i){
    return chisqT[i];
  }

  void SetChisqT(int i, double j){
    chisqT[i] = j;
  }

  int NpW(){
    return npW;
  }

  void SetNpW(int i){
    npW = i;
  }

  int IpW(int i, int j){
    if(ipW[i][j] >= 0)
      return evt[ipW[i][j]%1000][10];
    else 
      return ipW[i][j];
  }

  /* void SetIpW(int i, int j, int IPW){ */
  /*   ipW[i][j] = IPW; */
  /* } */

  double Mw(int i){
    return mW[i];
  }

  void SetMw(int i, double j){
    mW[i] = j;
  }
  
  double ChisqW(int i){
    return chisqW[i];
  }

  double ChisqWT(int i){
   //W belong to ith top
    return chisqW[ipT[i][0]];
  }
  
  int nBJets(int i){
   //W belong to ith top
    return ipT[i][2];
  }

  void SetChisqW(int i, double j){
    chisqW[i] = j;
  }

  double PW(int i, int j){
    return pW[i][j];
  }

  void SetPW(int i, int j, double mypW){
    pW[i][j] = mypW;
  }

  int Ileft(int i){
    return ileft[i];
  }

  void SetIleft(int i, int myIleft){
    ileft[i] = myIleft;
  }

  TLorentzVector PT(int i){
    TLorentzVector bJet(evt[ipT[i][1]][1], evt[ipT[i][1]][2], evt[ipT[i][1]][3], evt[ipT[i][1]][4]);
    
    TLorentzVector WBos(pW[ipT[i][0]][1], pW[ipT[i][0]][2], pW[ipT[i][0]][3], pW[ipT[i][0]][4]);

    return (bJet + WBos);
  }
  
  TLorentzVector PW(int i){
    //W belong to ith top
    TLorentzVector WBos(pW[ipT[i][0]][1], pW[ipT[i][0]][2], pW[ipT[i][0]][3], pW[ipT[i][0]][4]);
    return WBos;
  }

  int IsWLeptonic(int i){
    //Is the W in ith top leptonic? if yes, returns the flavor
    int a = 0;
    int WIndex = ipT[i][0];
    int partcleOrder = ipW[WIndex][1];
    if(evt[partcleOrder][0] == 2 || evt[partcleOrder][0] == 3)
      a = evt[partcleOrder][0];
    
    return a;
  }

  int heavyJet(int i){
    //Is there any overlap between jets in the ith top? if yes, return the -990 or -999 (W1b or W1j)
    int h = 0;

 /*    h = -1 * ipW[ipT[i][0]][1]; */
    
/*     if(evt[ipT[i][1]][9] > 0){ */
/*       h *= -1; */
/*     } */
    if(ipW[ipT[i][0]][1] < 0)
      h =  ipW[ipT[i][0]][1];
    else
      h = -1 * (ipW[ipT[i][0]][1] + 100);
    
    if(evt[ipT[i][1]][9] > 0)
      h *= -1;
    
    return h;
  }
  
private:



	static const int ndat = 11;
	static const int npartmax = 100;
	static const int npWmax = 300;
	static const int npTmax = 300;


	int nbJets;

		int npart;							// number of particles in evt
		double evt[npartmax][ndat];			// event data, list of all particles
		double MT2, MET, HT, MHT[6];	
		int npsjet[2], ipsjet[2][npartmax];	// list of jets in pseudojets

		int npT;						// number of reconstructed Tops
		int ipT[npTmax][3];					// 0 is for W, 1 is for additional jet, 2 number of bJets in the top
		double mT[npTmax], chisqT[npTmax];	// invariant mass and chisquared of Top
		int npW;						// number of reconstructed Ws
		int ipW[npWmax][2];					// particles i  reconstructed W
		double mW[npWmax], chisqW[npWmax];	// invariant mass and chisquared of W
		double pW[npWmax][5];				// momentum vector of W (label, px, py, pz, E), used internally
											// the 4-vector may be rescaled to the correct W mass if rescaleW = 1
		int ileft[npartmax];				// list of lonely jets








// private functions:
		virtual int GetnbJets ();
	// reconstructs W candidates
	virtual int GetWcand ();
	// removes overlapping W candidates
	virtual int RemWoverlaps ();
	// reconstructs Top candidates
	virtual int GetTcand ();
	// removes overlapping Top candidates
	virtual int RemToverlaps ();
	// classifies Ws coming from Top and lonely ones
	virtual void CleanWs ();
	// summary print at the end of each event (if debug >= 2)
	virtual void EvtFinal ();
	// writes the final run statistics (if debug >= 1)
	virtual int ttbarStats (int ifl);
	// helper functions
	virtual int IsInHemi (int ind);
	virtual double InvMass2 (double pvect1[], double pvect2[]);
	virtual double InvMass3 (double pvect1[], double pvect2[], double pvect3[]);
	virtual double GetChisqW (int i, int j);
	virtual double GetChisqT (double pWcor[], int i, int j, int k);
	virtual void GetMHT(double MHT[]);
	virtual int pzmissfrW (double lept[], double MET[], double MET1[], double MET2[]);

	// debug flag
	int debug;

	// Definition of constants and flags for the analysis
	int inAcc;				// if =1, keep only particles within eta acceptance
	int withbtag;			// = 1: btag preferred over chisq, = 2: correct b-jet requested for top, = 3: At least 1 bJet out of the top
	int recoW1j;			// = 0: do not reconstruct W from 1 jet, = 1: reconstruct W from 1 jet
	int recoWlept;			// = 0: do not reconstruct W leptonic, = 1: reconstruct W leptonic
	int ckWoverlap;			// = 0: do not check W overlaps, = 1: check W overlaps
	int rescaleW;			// = 0: no rescaling done, = 1 Whad 4-vector rescaled to the correct W mass
	double chisqWmax;		// maximum chisq to accept a W
	double chisqTmax;		// maximum chisq to accept a Top
	double chisqTWlepmax;	// maximum chisq to accept a Top with leptonic W decay
	double massWjbmin;		// minimum mass to consider heavy jet as b+j from top
	double massWjbmax;		// maximum mass to consider heavy jet as b+j from top
	double chisqW1b;		// penalty for chisquared of W = b+j from top
	double chisqWlept;		// penalty for chisquared of W leptonic

	double dpbypnom;		// constant for the jet energy resolution: dpbypnom / sqrt(p)
	double massW;			// W mass
	double widthW2j;		// W width for 2j cases
	double widthW1j;		// W width for 1j cases (from Daniel's hgr)
	double massT;			// Top mass	
	double widthT;			// Top width (default = 10 GeV)
//	double widthT = 1.5;	// Top width (1.5 GeV) used in had top analysis
	
	// arrays for final statistics
	int nread;
	int ntop[3];   // 0j, 1j, 2j
	int naddW[3];
	int nWinT[3];  // 1j, 1b, 2j
	int nWlone[3]; // 1j, 1b, 2j
	int njlone[3]; // 0j, 1j, 2j
	int nbtot[3];  // 0j, 1j, 2j
	int nbinT[3];
	int nbinTbW1b[0];
	int nbinW1jT[3];
	int nbinW2jT[3];
	int nbinW1j[3];
	int nbinW2j[3];
	int nbinJ[3];
	int nhTOK[3];  // 0j, 1j, 2j
	int nhbinT[3];
	int nhWinT[3];
	int nhWOK[3];
	int nhW[3];
	int ntop1[4];  // 1j, 1b, 2j, lep
	int ntop2[10];  // 1j1j, 1j1b, 1j2j, 1jlep, 1b1b, 1b2j, 2blep, 2j2j, 2jlep, leplep
	int nbevt[5];  // 0b, 1b, 2b, 3b, 4b
};
#endif
