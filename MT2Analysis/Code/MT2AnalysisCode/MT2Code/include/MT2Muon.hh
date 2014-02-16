#ifndef MT2Muon_hh
#define MT2Muon_hh

// MT2Muon ----------------------------------
class MT2Muon : public TObject {

public:
  MT2Muon(){Reset();}
  virtual ~MT2Muon(){}

  void Reset(){
    lv.SetPxPyPzE(0, 0, 0, 0);
    MT            = -1;
    Iso           = -1;
    Iso04         = -1;
    Charge        = -10;
    IsTightMuon   = -1;
    PassMu0_MuMu  = -1;
    PassMu1_MuMu  = -1;
    PassMu0_TauMu = -1;
    RejMu1_TauMu = -1;
    PassMu0_EleMu = -1;
    RejMu1_EleMu = -1;
    RejMu_TauTau = -1;
  }
  
  void SetLV(const TLorentzVector v){
    lv = v;
  }

  TLorentzVector lv;

  Float_t  MT;
  Float_t  Iso;
  Float_t  Iso04;
  Int_t    Charge;
  Int_t    IsTightMuon;
  Int_t    PassMu0_MuMu;
  Int_t    PassMu1_MuMu ;
  Int_t    PassMu0_TauMu;
  Int_t    RejMu1_TauMu;
  Int_t    PassMu0_EleMu;
  Int_t    RejMu1_EleMu;
  Int_t    RejMu_TauTau;
  ClassDef(MT2Muon, 9)
};
#endif




