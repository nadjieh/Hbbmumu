#include "MT2tree.hh"

void MT2tree::FillEleMu(){
  //{-1,-1} is Default value or no good muon is found or more than two muon is found 

  for(int i=0; i<NMuons; i++)
    {
      if  (muo[i].PassMu0_EleMu ==1) //muon i pass first good muon conditions.
	if (eleMu.mu0Ind==-1){  //Default value
	  eleMu.mu0Ind = i;
	  continue;
	} //one good muon is

      if  (muo[i].RejMu1_EleMu ==1){
	eleMu.HasNoVetomuoForEleMu=false;
      }
    }
    

  for(int i=0; i<NEles; i++){    
    if  (ele[i].IDSelEMU) //muon i pass first good muon conditions.
      if (eleMu.ele0Ind==-1) {  //Default value
	eleMu.ele0Ind = i; 
	continue;
      }//one good ele is 
    
    if (ele[i].IDVetoEMU){
      eleMu.HasNoVetoElecForEleMu=false;
    }
  }

for(int i=0; i<NEles; i++){    
    if  (ele[i].IDSelEMU) //muon i pass first good muon conditions.
      if (eleMu.ele0Ind==-1) {  //Default value
	eleMu.ele0Ind = i; 
	continue;
      }//one good ele is 
    
    if (ele[i].IDVetoEMU){
      eleMu.HasNoVetoElecForEleMu=false;
    }
  }

for(int i=0; i<NEles; i++){    
    if  (ele[i].IDSelEMU) //muon i pass first good muon conditions.
      if (eleMu.ele0Ind==-1) {  //Default value
	eleMu.ele0Ind = i; 
	continue;
      }//one good ele is 
    
    if (ele[i].IDVetoEMU){
      eleMu.HasNoVetoElecForEleMu=false;
    }
  }
for(int i=0; i<NEles; i++){    
    if  (ele[i].IDSelEMU) //muon i pass first good muon conditions.
      if (eleMu.ele0Ind==-1) {  //Default value
	eleMu.ele0Ind = i; 
	continue;
      }//one good ele is 
    
    if (ele[i].IDVetoEMU){
      eleMu.HasNoVetoElecForEleMu=false;
    }
  }

   
 

  if(eleMu.ele0Ind != -1 && eleMu.mu0Ind !=-1){
 

    
    eleMu.charge = muo[eleMu.mu0Ind].Charge + ele[eleMu.ele0Ind].Charge; //check two good muon charge


  
    if(fVerbose > 3){ 
      cout<<" Indices for Ele Muon "<<eleMu.mu0Ind<<" "<<eleMu.ele0Ind<<endl;
      
      cout<<" muo ";
      muo[eleMu.mu0Ind].lv.Print();
      
      cout<<" ele ";
      ele[eleMu.ele0Ind].lv.Print();
    }  
  
    eleMu.MT2    = CalcMT2(0, 0,muo[eleMu.mu0Ind].lv,ele[eleMu.ele0Ind].lv, pfmet[0]);
 
    TLorentzVector met = -(muo[eleMu.mu0Ind].lv +ele[eleMu.ele0Ind].lv);
    if(fVerbose > 3) cout<<"met.Pt() "<<met.Pt()<<endl;
    eleMu.METImbalanced = met.Pt();
    eleMu.METImbalancedPhi = met.Phi();
    eleMu.MT2Imbalanced = CalcMT2(0, 0,muo[eleMu.mu0Ind].lv,ele[eleMu.ele0Ind].lv, met);
    if(fVerbose > 3) cout<<"  eleMu.MT2Imbalanced  "<<eleMu.MT2Imbalanced <<endl;
  }
}

