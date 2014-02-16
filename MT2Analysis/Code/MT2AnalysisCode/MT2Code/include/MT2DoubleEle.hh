#ifndef MT2DoubleEle_hh
#define MT2DoubleEle_hh

// MT2DoubleEle ----------------------------------
class MT2DoubleEle : public TObject {

public:
  MT2DoubleEle(){Reset();}
  ~MT2DoubleEle(){}

  void Reset(){
    Ele0Ind            = -1;
    Ele1Ind            = -1;
    PairCharge           = -1;
    MT2            = -1;
    MT2Imbalanced = -1;
    METImbalanced = -1;

  }

  Int_t    Ele0Ind;
  Int_t    Ele1Ind;
  Int_t    PairCharge;
  Float_t  MT2;
  Float_t  MT2Imbalanced;
  Float_t  METImbalanced;

  ClassDef(MT2DoubleEle, 1)
    };
#endif

