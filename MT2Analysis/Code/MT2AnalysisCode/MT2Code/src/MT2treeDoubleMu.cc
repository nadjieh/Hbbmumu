#include "MT2tree.hh"

void MT2tree::FillDoubleMu(){
  //{-1,-1} is Default value or no good muon is found or more than two muon is found 

  for(int i=0; i<NMuons; i++){
  
    if  (muo[i].PassMu0_MuMu ==1) //muon i pass first good muon conditions.
      if (doubleMu.mu0Ind==-1)   //Default value
	{ doubleMu.mu0Ind = i;
	  continue;
	} //one goodmuon is found
    
    if  ( muo[i].PassMu1_MuMu ==1) //muon i pass second good muon conditions.
      if (doubleMu.mu1Ind==-1)
	{ doubleMu.mu1Ind = i;
	  continue;
	} //second good muon is found
  }
 
  if  (doubleMu.mu0Ind!=-1 && doubleMu.mu1Ind!=-1){
    
    doubleMu.chargeSum = muo[doubleMu.mu0Ind].Charge + muo[doubleMu.mu1Ind].Charge;//check two good muon charge
  
    //cout<<" Indices for Double Muon "<<doubleMu.mu0Ind<<" "<<doubleMu.mu1Ind<<endl;
  
    doubleMu.MT2     = CalcMT2(0, 0,muo[doubleMu.mu0Ind].lv,muo[doubleMu.mu1Ind].lv, pfmet[0]);
 
    TLorentzVector met = -(muo[doubleMu.mu0Ind].lv + muo[doubleMu.mu1Ind].lv);
    //cout<<"met.Pt() "<<met.Pt()<<endl;
    doubleMu.METImbalanced = met.Pt();
    doubleMu.METImbalancedPhi = met.Phi();
    doubleMu.MT2Imbalanced = CalcMT2(0, 0,muo[doubleMu.mu0Ind].lv,muo[doubleMu.mu1Ind].lv, met);
  }
}
