#ifndef MT2Elec_hh
#define MT2Elec_hh

// MT2Elec ----------------------------------
class MT2Elec : public TObject {

public:
  MT2Elec(){Reset();}
  virtual ~MT2Elec(){}

  void Reset(){
  lv.SetPxPyPzE(0, 0, 0, 0);
  MT            = -1;
  Iso           = -1;
  Iso04         = -1;
  Charge        = -10;
  IDMedium      = -1;
  IDLoose       = -1;
  IDVeto        = -1;
  //IDVetoEE      = -1;
  IDSelEE       = -1;

  IDVetoEMU     = -1;
  IDSelEMU      = -1;

  IDVetoETau    = -1; 
  IDSelETau     = -1;


  IDVetoMuTau  = -1;
  IDVetoTauTau = -1;

  PassE0_EE    = -1;
  PassE1_EE    = -1;
  }
  void SetLV(const TLorentzVector v){lv = v;}

  TLorentzVector lv;

  Float_t  MT;
  Float_t  Iso;
  Float_t  Iso04;
  Int_t    Charge;
  Int_t    IDMedium;
  Int_t    IDLoose;
  Int_t    IDVeto;

  //Int_t    IDVetoEE;
  Int_t    IDSelEE;

  Int_t    IDVetoEMU;
  Int_t    IDSelEMU;

  Int_t    IDVetoETau;
  Int_t    IDSelETau;

  Int_t    IDVetoMuTau;
  Int_t    IDVetoTauTau;

  Int_t PassE0_EE;
  Int_t PassE1_EE;

  ClassDef(MT2Elec, 12)
};
#endif
