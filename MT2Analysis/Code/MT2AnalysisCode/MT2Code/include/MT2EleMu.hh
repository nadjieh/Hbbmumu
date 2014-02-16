#ifndef MT2EleMu_hh
#define MT2EleMu_hh

// MT2EleMu ----------------------------------
class MT2EleMu : public TObject {

public:
  MT2EleMu(){Reset();}
  virtual ~MT2EleMu(){} 
 
 

  void Reset(){
    
    mu0Ind            = -1;
    ele0Ind           = -1;
    charge            = -10;
    METImbalanced   = -1.;
    METImbalancedPhi = -10.;
    MT2           = -1.;
    MT2Imbalanced = -1.; 
    HasNoVetoElecForEleMu=true;
    HasNoVetomuoForEleMu=true;
    HasNoVetoElecForEleMu=true;
    isSignalEleMu = false;

   }
  


  Int_t    mu0Ind;
  Int_t    ele0Ind;
  Int_t    charge;
  Float_t  METImbalanced;
  Float_t  METImbalancedPhi;
  Float_t  MT2;
  Float_t  MT2Imbalanced ;
  bool     HasNoVetomuoForEleMu;
  bool     HasNoVetoElecForEleMu;
  bool	   isSignalEleMu;

ClassDef(MT2EleMu, 1)
    }; 
#endif
