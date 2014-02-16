#ifndef BTagWeight_HH
#define BTagWeight_HH

#include <vector>
#include <stdlib.h>
#include <stdio.h>
#include <map>
#include <string>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <math.h>

#include "TString.h"
#include "TF1.h"
#include "TH1.h"
#include "TH1D.h"
#include "TFile.h"

using namespace std;

inline Float_t getBTagEventWeightError(float &werr, vector<float> jetEff, vector<float> jetEffErr, vector<float> jetSF, vector<float> jetSFErr, int nBTags){
  // jetEff contains mean efficiencies per jet (average sample eff per flavour)
  // jetSF contains SF per jet (pt, eta, flavour dependent)

  int taggers = 1; //we use only one tagger

  int njets=jetEff.size();
  std::vector<int> comb(jetEff.size());
  for(int i=0;i < njets; i++) comb[i]=0;
  int idx=0;
  int max=taggers+1; //force sorted tagging //1 << taggers;
  float pMC=0;
  float pData=0;
  float pMCErr=0;
  float pDataErr=0;
  if(njets==0) return 0.;
  while(comb[jetEff.size()-1] < max)    {
    
    std::vector<int> tags;
    for(int j=0;j<taggers;j++) tags.push_back(0);
    
    float mc=1.;
    float data=1.;
    float mcerr = 0.;
    float dataerr = 0.;
    for(size_t j=0;j<njets;j++) // loop on jets
      {
	// if none tagged, take the 1-eff SF for the loosest:
	float tagMc = 1.-jetEff[j];
	float tagData = 1.-jetEff[j]*jetSF[j];
	if(comb[j]> 0) //if at least one tagged take the SF for the tightest tagger -> we do have the same
	  {
	     tagMc = jetEff[j];
	     tagData = jetEff[j]*jetSF[j];
	  }
	for(size_t k=0;k< taggers; k++ ) // loop on taggers
	  {
	    bool tagged = (comb[j] > k) ; ///((comb[j] >> k) & 0x1) == 1;
	    if(tagged) tags[k]++;
	  }
	mc*=tagMc;       
	data*=tagData;       
      }
	//need 2 for loops since I need mc, data
    for(size_t j=0;j<njets;j++) // loop on jets
      { //get factor error
	// if none tagged, take the 1-eff SF for the loosest:
	float tagMc = 1.-jetEff[j];
	float tagData = 1.-jetEff[j]*jetSF[j];
	float tagMcErr = jetEffErr[j];
	float tagDataErr = sqrt(pow(jetEffErr[j]*jetSF[j],2) + pow(jetEff[j]*jetSFErr[j],2));
	if(comb[j]> 0) //if at least one tagged take the SF for the tightest tagger -> we do have the same
	  {
	     tagMcErr = jetEffErr[j];
	     tagDataErr = sqrt(pow(jetEffErr[j]*jetSF[j],2) + pow(jetEff[j]*jetSFErr[j],2));
	     tagMc = jetEff[j];
	     tagData = jetEff[j]*jetSF[j];
	  }
	if(tagMc!=0 && tagData!=0) mcerr+= pow(mc*tagMcErr/tagMc,2);       
	if(tagMc!=0 && tagData!=0) dataerr+= pow(data*tagDataErr/tagData,2);       
	//if(tagMc==0 || tagData==0) cout << "data or mc has no weight, ignore jet " << j << " (tagmc " << tagMc << ", tagData " << tagData << ")" << endl;
      }
      
    if( tags[0] >= abs(nBTags) && nBTags<0)
      {
	//  std::cout << mc << " " << data << " " << data/mc << " " << pData/pMC << endl;
	pMC+=mc;
	pData+=data;
	pMCErr += mcerr;
	pDataErr += dataerr;
	//n    std::cout << std::endl<< "mc, data,ratioThis,ratioTot " <<  mc << " " << data << " " << data/mc << " " << pData/pMC << endl;
      }
    else if( tags[0] == abs(nBTags) && nBTags>0)
      {
	pMC+=mc;
	pData+=data;
	pMCErr += mcerr;//are squared errors
	pDataErr += dataerr;//are squared errors
      }
    while (comb[idx] == max -1  && idx+1 < jetEff.size()) idx++; // find first jets for which we did not already test all configs 
    comb[idx]++;  // test a new config for that jet
    for(int i=0;i<idx;i++) { comb[i]=0; } // reset the tested configs for all previous jets
    idx=0;
  }
  if(pMC==0) {
    werr = 0.;
    return 0;
  } 
  werr =  sqrt(pDataErr/(pMC*pMC) + pMCErr*pow(pData/(pMC*pMC),2));
  return pData/pMC;

  
}


inline Float_t getBTagEventWeight(vector<float> jetEff, vector<float> jetSF, int nBTags=1){
  // jetEff contains mean efficiencies per jet (average sample eff per flavour)
  // jetSF contains SF per jet (pt, eta, flavour dependent)

  int taggers = 1; //we use only one tagger

  int njets=jetEff.size();
  std::vector<int> comb(jetEff.size());
  for(int i=0;i < njets; i++) comb[i]=0;
  int idx=0;
  int max=taggers+1; //force sorted tagging //1 << taggers;
  float pMC=0;
  float pData=0;
  if(njets==0) return 0.;
  while(comb[jetEff.size()-1] < max)    {
    
    std::vector<int> tags;
    for(int j=0;j<taggers;j++) tags.push_back(0);
    
    float mc=1.;
    float data=1.;
    for(size_t j=0;j<njets;j++) // loop on jets
      {
	// if none tagged, take the 1-eff SF for the loosest:
	float tagMc = 1.-jetEff[j];
	float tagData = 1.-jetEff[j]*jetSF[j];
	if(comb[j]> 0) //if at least one tagged take the SF for the tightest tagger -> we do have the same
	  {
	    //int k=comb[j]-1;
	    //tagMc=jets[j][k].eff;
	    //tagData=jets[j][k].eff*jets[j][k].sf;
	     tagMc = jetEff[j];
	     tagData = jetEff[j]*jetSF[j];
	     /*
	     if(comb[j]< taggers) //if at least one tagged take the SF for the tightest tagger
	       {
		 int k1=comb[j];
		 tagMc*=1-jets[j][k1].eff/jets[j][k].eff;
		 tagData*=1-jets[j][k1].eff/jets[j][k].eff*jets[j][k1].sf/jets[j][k].sf;
		 
	       }
	     */
	  }
	
	
	for(size_t k=0;k< taggers; k++ ) // loop on taggers
	  {
	    bool tagged = (comb[j] > k) ; ///((comb[j] >> k) & 0x1) == 1;
	    if(tagged) tags[k]++;
	  }
	
	//	cout << "- " << tagMc << " " << tagData << endl;

	mc*=tagMc;       
	data*=tagData;       
      }
      
    if( tags[0] >= nBTags)
      {
	//  std::cout << mc << " " << data << " " << data/mc << " " << pData/pMC << endl;
	pMC+=mc;
	pData+=data;
	//n    std::cout << std::endl<< "mc, data,ratioThis,ratioTot " <<  mc << " " << data << " " << data/mc << " " << pData/pMC << endl;
      }
    while (comb[idx] == max -1  && idx+1 < jetEff.size()) idx++; // find first jets for which we did not already test all configs 
    // next combination
    comb[idx]++;  // test a new config for that jet
    for(int i=0;i<idx;i++) { comb[i]=0; } // reset the tested configs for all previous jets
    idx=0;
  }
  if(pMC==0) return 0; 
  return pData/pMC;
  
  
}


//NBtags <0 => request for >=nBTags
inline Float_t getBTagEventWeight2( vector< vector<float> >  jetEff, vector< vector<float> > jetSF, int nBTags){
  //template <class Filter> float BTagWeight::weight(vector<vector<JetInfo> >jets)
  //{
  // jetEff contains []
  
  int taggers = 1; //we use only one tagger
  unsigned int njets=jetEff.size();
  std::vector<unsigned int> comb(jetEff.size());
 for(size_t i=0;i < njets; i++) comb[i]=0;
 unsigned int idx=0;
 unsigned int max=taggers+1; //force sorted tagging //1 << taggers;
 float pMC=0;
 float pData=0;
 if(njets==0) return 0.;
 while(comb[njets-1] < max)
 {
   std::vector<int> tags;
   for(size_t j=0;j<taggers;j++) tags.push_back(0);

   float mc=1.;
   float data=1.;
   for(size_t j=0;j<njets;j++) // loop on jets
    {
  // if none tagged, take the 1-eff SF for the loosest:
     float tagMc = 1.-jetEff[j][0];
     float tagData = 1.-jetEff[j][0]*jetSF[j][0];
     if(comb[j]> 0) //if at least one tagged take the SF for the tightest tagger
     {
   	   int k=comb[j]-1;
           tagMc=jetEff[j][k];
           tagData=jetEff[j][k]*jetSF[j][k];

     if(comb[j]< taggers) //if at least one tagged take the SF for the tightest tagger
     {
           int k1=comb[j];
           tagMc*=1-jetEff[j][k1]/jetEff[j][k];
           tagData*=1-jetEff[j][k1]/jetEff[j][k]*jetSF[j][k1]/jetSF[j][k];

     }
     }

     for(size_t k=0;k< taggers; k++ ) // loop on taggers
      {
       bool tagged = (comb[j] > k) ; ///((comb[j] >> k) & 0x1) == 1;
       if(tagged) tags[k]++;
      }

    mc*=tagMc;       
    data*=tagData;       
   
  }
   if( (tags[0] == abs(nBTags) && nBTags>0) || (tags[0] >= abs(nBTags) && nBTags<0 ) /*Filter::filter(tags)*/)
   {
    pMC+=mc;
    pData+=data;
   }
   while (comb[idx] == max -1  && idx+1 < njets) idx++; // find first jets for which we did not already test all configs 
   // next combination
  comb[idx]++;  // test a new config for that jet
  for(size_t i=0;i<idx;i++) { comb[i]=0;  } // reset the tested configs for all previous jets
  idx=0;
 }
  if(pMC==0) return 0; 
  return pData/pMC;
}

















//inline float* getBTagSF(TString tagger,float pt, float eta ){//OLD
inline float getBTagSF(float &err, TString tagger,float pt, float eta){//NEW
  float ptmin[14] = {30, 40, 50, 60, 70, 80, 100, 120, 160, 210, 260, 320, 400, 500};
  float ptmax[14] = {40, 50, 60, 70, 80,100, 120, 160, 210, 260, 320, 400, 500, 670};
  float SFb=0;
  float SFb_error[14];
  if(tagger=="SSVHEM"){
    //Tagger: SSVHEM within 30 < pt < 670 GeV, abs(eta) < 2.4, x = pt
    SFb = 0.896462*((1.+(0.00957275*pt))/(1.+(0.00837582*pt)));
    float tSFb_error[14] = { 0.0316234,      0.0310149,      0.02381,      0.0223228,      0.023461,      0.0202517,      0.0156249,    0.0214799,      0.0399369,      0.0416666,      0.0431031,      0.0663209,      0.0687731,      0.0793305 };
    memcpy(SFb_error, tSFb_error,140);  
  }
  else if(tagger=="SSVHPT"){
    //Tagger: SSVHPT within 30 < pt < 670 GeV, abs(eta) < 2.4, x = pt
    SFb = 0.422556*((1.+(0.437396*pt))/(1.+(0.193806*pt)));
    float tSFb_error[] = {
      0.0403485,
      0.0396907,
      0.0291837,
      0.0325778,
      0.0335716,
      0.0255023,
      0.0300639,
      0.0253228,
      0.0409739,
      0.043561,
      0.0458427,
      0.0763302,
      0.0781752,
      0.108927 };
    memcpy(SFb_error, tSFb_error,140);
    
  }

  float SFErr  = 0;
  for(int i=0;i<14;i++){
    if( pt >= ptmin[i] && pt <= ptmax[i])
      SFErr = SFb_error[i];
  }
  if(pt>ptmax[13]){
    float ptt = 670.;
    if(tagger=="SSVHEM"){
	SFb = 0.896462*((1.+(0.00957275*ptt))/(1.+(0.00837582*ptt)));
	SFErr = 2.* SFb_error[13];
   }
   else if(tagger=="SSVHPT"){
	SFb = 0.422556*((1.+(0.437396*ptt))/(1.+(0.193806*ptt)));
	SFErr = 2.* SFb_error[13];
	
   }
  }
  if(pt<ptmin[0]){
    float ptt = 30.;
    if(tagger=="SSVHEM"){
	SFb = 0.896462*((1.+(0.00957275*ptt))/(1.+(0.00837582*ptt)));
	SFErr = 0.12;
   }
   else if(tagger=="SSVHPT"){
	SFb = 0.422556*((1.+(0.437396*ptt))/(1.+(0.193806*ptt)));
	SFErr = 0.12;
   }
  }

  //float *res = new float[2];//OLD
  //res[0] = SFb; res[1] = SFErr;//OLD
  //return res;//OLD
  err = SFErr;
  return SFb;
}


//inline float* getMistagSF(TString tagger, float pt, float eta){//OLD
inline float getMistagSF(float &err, TString tagger, float pt, float eta, int up=0){//NEW

  eta = fabs(eta);
  TF1* SFl, *SFlmin, *SFlmax;
  if( tagger == "SSVHEM" && eta>= 0.0 && eta< 0.8) {
      SFl = new TF1("SFlight","((0.86318+(0.000801639*x))+(-1.64119e-06*(x*x)))+(2.59121e-10*(x*(x*x)))", 20.,670.);
      SFlmin = new TF1("SFlightMin","((0.790364+(0.000463086*x))+(-4.35934e-07*(x*x)))+(-9.08296e-10*(x*(x*x)))", 20.,670.);
      SFlmax = new TF1("SFlightMax","((0.935969+(0.0011402*x))+(-2.84645e-06*(x*x)))+(1.42654e-09*(x*(x*x)))", 20.,670.);
    }
  
  //if( tagger == "SSVHEM" && sEtamin == "0.0" && sEtamax == "2.4") {
  //    SFl = new TF1("SFlight","((0.890254+(0.000553319*x))+(-1.29993e-06*(x*x)))+(4.19294e-10*(x*(x*x)))", 20.,670.);
  //    SFlmin = new TF1("SFlightMin","((0.817099+(0.000421567*x))+(-9.46432e-07*(x*x)))+(1.62339e-10*(x*(x*x)))", 20.,670.);
  //    SFlmax = new TF1("SFlightMax","((0.963387+(0.000685092*x))+(-1.65343e-06*(x*x)))+(6.76249e-10*(x*(x*x)))", 20.,670.);
  //  }
  
  if( tagger == "SSVHEM" && eta>= 0.8 && eta< 1.6) {
      SFl = new TF1("SFlight","((0.958973+(-0.000269555*x))+(1.381e-06*(x*x)))+(-1.87744e-09*(x*(x*x)))", 20.,670.);
      SFlmin = new TF1("SFlightMin","((0.865771+(-0.000279908*x))+(1.34144e-06*(x*x)))+(-1.75588e-09*(x*(x*x)))", 20.,670.);
      SFlmax = new TF1("SFlightMax","((1.0522+(-0.000259296*x))+(1.42056e-06*(x*x)))+(-1.999e-09*(x*(x*x)))", 20.,670.);
    }
  if( tagger == "SSVHEM" && eta>=1.6 && eta < 2.41)    {
      SFl = new TF1("SFlight","((0.923033+(-0.000898227*x))+(4.74565e-06*(x*x)))+(-6.11053e-09*(x*(x*x)))", 20.,670.);
      SFlmin = new TF1("SFlightMin","((0.828021+(-0.000731926*x))+(4.19613e-06*(x*x)))+(-5.81379e-09*(x*(x*x)))", 20.,670.);
      SFlmax = new TF1("SFlightMax","((1.01812+(-0.00106483*x))+(5.29518e-06*(x*x)))+(-6.40728e-09*(x*(x*x)))", 20.,670.);
    }
  if( tagger == "SSVHPT" && eta>= 0.0 && eta < 2.41)    {
      SFl = new TF1("SFlight","((0.97409+(0.000646241*x))+(-2.86294e-06*(x*x)))+(2.79484e-09*(x*(x*x)))", 20.,670.);
      SFlmin = new TF1("SFlightMin","((0.807222+(0.00103676*x))+(-3.6243e-06*(x*x)))+(3.17368e-09*(x*(x*x)))", 20.,670.);
      SFlmax = new TF1("SFlightMax","((1.14091+(0.00025586*x))+(-2.10157e-06*(x*x)))+(2.41599e-09*(x*(x*x)))", 20.,670.);
    }
  
  float SF = SFl->Eval(pt);
  float SFmin = fabs(SFlmin->Eval(pt) - SF);
  float SFmax = fabs(SFlmax->Eval(pt) - SF);
  float SFErr = SFmin;
  if(SFmin < SFmax) SFErr = SFmax;
  if(pt>670.){
     SF = SFl->Eval(670.);
     SFmin = fabs(SFlmin->Eval(670.) - SF);
     SFmax = fabs(SFlmax->Eval(670.) - SF);
     SFErr = 2.*SFmin;
     if(SFmin < SFmax) SFErr = 2.*SFmax;
  }
 // float * SFArr = new float[2];//OLD
  delete SFl; delete SFlmin; delete SFlmax;
 // SFArr[0] = SF; SFArr[1] = SFErr;//OLD
 // return SFArr;//OLD
   if(up==1){
      SF = SFmax + SF;
   } else if(up==(-1)){
      SF = SFmin + SF;
   }
   err = SFErr;//NEW
   return SF;//NEW
}


#endif