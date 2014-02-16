#ifndef MT2DoubleMuon_hh
#define MT2DoubleMuon_hh

// MT2DoubleMu ----------------------------------
class MT2DoubleMuon : public TObject {

public:
  MT2DoubleMuon(){Reset();}
  virtual ~MT2DoubleMuon(){}

  void Reset(){
    mu0Ind            = -1;
    mu1Ind            = -1;
    chargeSum            = -1;
    METImbalanced    = -1;
    METImbalancedPhi = -10;
    MT2            = -1;
    MT2Imbalanced = -1;
    hasNoVetoElec = true; //should be fixed, just a placeholder for now
    hasNoVetoMu   = true; //should be fixed, just a placeholder for now
  }

  Int_t    mu0Ind;
  Int_t    mu1Ind;
  Int_t    chargeSum;
  Float_t  METImbalanced;
  Float_t  METImbalancedPhi;
  Float_t  MT2;
  Float_t  MT2Imbalanced;
  bool     hasNoVetoElec;
  bool     hasNoVetoMu;

  ClassDef(MT2DoubleMuon, 1)
    };
#endif
