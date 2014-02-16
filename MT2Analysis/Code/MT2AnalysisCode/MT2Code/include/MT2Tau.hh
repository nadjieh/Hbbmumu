#ifndef MT2Tau_hh
#define MT2Tau_hh

// MT2Tau  ----------------------------------
class MT2Tau : public TObject {

public:
  MT2Tau(){ Reset();}
  virtual ~MT2Tau(){}

  void Reset(){
    lv.SetPxPyPzE(0, 0, 0, 0);
    MT            = -1;
    Charge        = -1;
    JetPt         = -1;
    JetEta        = -9.99;
    JetPhi        = -9.99;
    JetMass       = -1;
    Isolation     = -1;
    ElectronRej   = -1;
    MuonRej       = -1;
    Isolation3Hits= -1;
    IsolationMVA2 = -1;
    ElectronRejMVA3=-1;
    MuonRej2      = -1;
    isLooseID     = 0;
    isLooseID3Hits= 0;
    isLooseIDMVA  = 0;
    PassTau0_TauTau = -1;
    PassTau1_TauTau = -1;
    PassTau_ElTau = -1;
    PassTau_MuTau = -1;
    PassQCDTau0_TauTau = -1;
    PassQCDTau1_TauTau = -1;
  }

  void SetLV(const TLorentzVector v){
    lv = v;
  }
  
  Bool_t IsGoodTau(float minJPt, float maxJEta, int iso, int elerej, int muorej){
    float pt = lv.Pt();
    if( pt < minJPt )               return false; // need to do this first, before getting lv.Eta();
    float eta = lv.Eta();
    if( fabs(eta) > maxJEta)        return false;
    if(iso>=0){//normal: 1 to 4
      if(Isolation<iso)               return false;
    } else if(iso<-10){//mva: -12 to -14
      if(IsolationMVA2<(-iso-10))     return false;
    } else{//3hits: -2 to -4
      if(Isolation3Hits<(-iso))       return false;
    }
    if(elerej>=0){
      if(ElectronRej<elerej)          return false;
    } else{//mva3
      if(ElectronRejMVA3<(-elerej))   return false;
    }
    if(muorej>=0){
      if(MuonRej<muorej)              return false;
    } else{
      if(MuonRej2<(-muorej))          return false;
    }

    return true;
  }

  TLorentzVector lv;

  Float_t  MT;
  Int_t    Charge;
  Float_t  JetPt;
  Float_t  JetEta;
  Float_t  JetPhi;
  Float_t  JetMass;
  Int_t    Isolation;
  Int_t    ElectronRej;
  Int_t    MuonRej;
  Int_t    Isolation3Hits;
  Int_t    IsolationMVA2;
  Int_t    ElectronRejMVA3;
  Int_t    MuonRej2;
  Bool_t   isLooseID;
  Bool_t   isLooseID3Hits;
  Bool_t   isLooseIDMVA;

  Int_t PassTau0_TauTau;
  Int_t PassTau1_TauTau;
  Int_t PassTau_ElTau;
  Int_t PassTau_MuTau;

  Int_t PassQCDTau0_TauTau;
  Int_t PassQCDTau1_TauTau;

  ClassDef(MT2Tau, 5)
    };
#endif
