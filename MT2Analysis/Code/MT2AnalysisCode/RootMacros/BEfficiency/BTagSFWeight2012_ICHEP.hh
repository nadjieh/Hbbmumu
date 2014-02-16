#ifndef BTagSFWeight_HH
#define BTagSFWeight_HH

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

// Calculate the MC weight due to b-tag SF and error due to statistical uncertainty of the b efficiency
// NBtags <0  => request for >=nBTags
// NBtags >=0 => request for ==nBTags
// The error here is only due to eff-MC-statistics, not due to SF uncertainty
// if you don't want to calculate the error use getBTagEventWeight
inline Float_t getBTagEventWeightError(float &werr, vector<float> jetEff, vector<float> jetEffErr, vector<float> jetSF, int nBTags){
	if(jetEff.size()!=jetSF.size()){
		cout << "efficiency vector and SF vector have different sizes, return -1" << endl;
		return -1;
	}
	bool calcerror = true;
	if(jetEffErr.size() != jetEff.size()){
		calcerror = false;
		jetEffErr = jetEff;//make dummy copy
	}
  unsigned int taggers = 1; //we use only one tagger
  unsigned int njets=jetEff.size();
  std::vector<int> comb(jetEff.size());
  for(unsigned int i=0;i < njets; i++) comb[i]=0;
  int idx=0;
  int max=taggers+1; //force sorted tagging //1 << taggers;
  float pMC=0;
  float pData=0;
  float pMCErr=0;
  float pDataErr=0;
  if(njets==0) return 0.;
  while(comb[jetEff.size()-1] < max)    {
    std::vector<int> tags;
    for(unsigned int j=0;j<taggers;j++) tags.push_back(0);
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
	float tagDataErr = jetEffErr[j]*jetSF[j];
	if(comb[j]> 0) //if at least one tagged take the SF for the tightest tagger -> we do have the same
	  {
	     tagMcErr = jetEffErr[j];
	     tagDataErr = jetEffErr[j]*jetSF[j];
	     tagMc = jetEff[j];
	     tagData = jetEff[j]*jetSF[j];
	  }
	if(tagMc!=0 && tagData!=0) mcerr+= pow(mc*tagMcErr/tagMc,2);       
	if(tagMc!=0 && tagData!=0) dataerr+= pow(data*tagDataErr/tagData,2);       
      }
    if( (tags[0] >= abs(nBTags) && nBTags<0)  || (tags[0]>=0 && nBTags==-99) )
      {
	//  std::cout << mc << " " << data << " " << data/mc << " " << pData/pMC << endl;
	pMC+=mc;
	pData+=data;
	pMCErr += mcerr;
	pDataErr += dataerr;
	//    std::cout << std::endl<< "mc, data,ratioThis,ratioTot " <<  mc << " " << data << " " << data/mc << " " << pData/pMC << endl;
      }
    else if( tags[0] == abs(nBTags) && nBTags>=0)
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
    if(calcerror) werr = 0.;
    return 0;
  } 
  if(calcerror) werr =  sqrt(pDataErr/(pMC*pMC) + pMCErr*pow(pData/(pMC*pMC),2));
  return pData/pMC;

  
}

// Calculate the MC weight due to b-tag SF
// NBtags <0  => request for >=nBTags
// NBtags >=0 => request for ==nBTags
// no errors are calculated here
inline Float_t getBTagEventWeight(vector<float> jetEff, vector<float> jetSF, int nBTags=-1){
	if(jetEff.size()!=jetSF.size()){
		cout << "efficiency vector and SF vector have different sizes, return -1" << endl;
		return -1;
	}
  unsigned int taggers = 1; //we use only one tagger -- if using two taggers look up https://twiki.cern.ch/twiki/pub/CMS/BTagWeight/BTagWeight4.h
  unsigned int njets=jetEff.size();
  std::vector<int> comb(jetEff.size());
  for(unsigned int i=0;i < njets; i++) comb[i]=0;
  int idx=0;
  int max=taggers+1; //force sorted tagging //1 << taggers;
  float pMC=0;
  float pData=0;
  if(njets==0) return 0.;
  while(comb[jetEff.size()-1] < max)    {
    std::vector<int> tags;
    for(unsigned int j=0;j<taggers;j++) tags.push_back(0);
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
	mc*=tagMc;       
	data*=tagData;       
      }
   if( (tags[0] == abs(nBTags) && nBTags>=0) || (tags[0] >= abs(nBTags) && nBTags<0 ) || (tags[0]>=0 && nBTags==-99) )
      {
	//  std::cout << mc << " " << data << " " << data/mc << " " << pData/pMC << endl;
	pMC+=mc;
	pData+=data;
	//  std::cout << std::endl<< "mc, data,ratioThis,ratioTot " <<  mc << " " << data << " " << data/mc << " " << pData/pMC << endl;
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

// gets the b and c SF for jets
// err contains the error on SF for b jets
// note: for c jets need to inflate the error by 2
inline float getBtagSF(float &err, TString tagger,float pt, float eta){
  float ptmin[14] = {30, 40, 50, 60, 70, 80, 100, 120, 160, 210, 260, 320, 400, 500};
  float ptmax[14] = {40, 50, 60, 70, 80,100, 120, 160, 210, 260, 320, 400, 500, 670};
  float SFb=0;
  float SFb_error[14];
  bool overflow  = false;
  bool underflow = false;
  eta = fabs(eta);
  if(pt> ptmax[13]){ overflow = true; pt = 670;}
  if(pt< ptmin[0] ){underflow = true; pt = 30;}
  if(tagger=="CSVM"){
	SFb = 0.6981*((1.+(0.414063*pt))/(1.+(0.300155*pt)));
	float tSFb_error[14] = {	0.0295675,	0.0295095,	0.0210867,	0.0219349,	0.0227033,	0.0204062,	0.0185857,	0.0256242,	0.0383341,	0.0409675,	0.0420284,	0.0541299,	0.0578761,	0.0655432 };
    memcpy(SFb_error, tSFb_error,140);
  } else if(tagger=="CSVT"){
	SFb = 0.901615*((1.+(0.552628*pt))/(1.+(0.547195*pt)));
	float tSFb_error[14] = {	0.0364717,	0.0362281,	0.0232876,	0.0249618,	0.0261482,	0.0290466,	0.0300033,	0.0453252,	0.0685143,	0.0653621,	0.0712586,	0.094589,	0.0777011,	0.0866563 };
    memcpy(SFb_error, tSFb_error,140);
  } else if(tagger=="JPM"){
	SFb = 0.90806*((1.+(0.000236997*pt))/(1.+(5.49455e-05*pt)));
	float tSFb_error[14] = {	0.0352594,	0.0353008,	0.0299008,	0.0276606,	0.0292312,	0.0336607,	0.0284701,	0.029544,	0.0358872,	0.0367869,	0.0375048,0.0597367,	0.0653152,	0.074242 };
    memcpy(SFb_error, tSFb_error,140);
  } else if(tagger=="JPT"){
	SFb = 0.835882*((1.+(0.00167826*pt))/(1.+(0.00120221*pt)));
	float tSFb_error[14] = {	0.0475813,	0.0472359,	0.0378328,	0.0334787,	0.034681,	0.0398312,	0.0481646,	0.0392262,	0.0463086,	0.0534565,	0.0545823,	0.102505,	0.113198,	0.138116 };
    memcpy(SFb_error, tSFb_error,140);
  }
  float SFErr  = 0;
  for(int i=0;i<14;i++){
    //prescription for 2012: inflate error by 1.5
    if( pt >= ptmin[i] && pt <= ptmax[i] && !overflow && !underflow)
      SFErr = SFb_error[i] * 1.5;

  }
  if(overflow){
    //prescription for 2012: inflate error by 1.5, inflate by 2 for overflow
	SFErr = 2.* SFb_error[13] * 1.5;
  }
  if(underflow){
    //prescription for 2012: inflate error by 1.5, 0.12 constant from 2011
	SFErr = 0.12 * 1.5;
  }

  err = SFErr;
  return SFb;
}

// get light SF due to b-tagging per jet
// if up=-1, get lower bound function
// if up=+1, get upper bound function
inline float getMistagSF(TString tagger, float pt, float eta, int up=0){

  bool overflow  = false;
  bool underflow = false;//is always false
  if(pt> 670){ overflow = true; pt = 670;}
  eta = fabs(eta);
  TF1* SFl, *SFlmin, *SFlmax, *SFcorrection;
	if( tagger == "CSVM" && eta>= 0.0 && eta < 0.8)
	{
	SFl = new TF1("SFlight","((1.06182+(0.000617034*x))+(-1.5732e-06*(x*x)))+(3.02909e-10*(x*(x*x)))", 20.,670.);
	SFlmin = new TF1("SFlightMin","((0.972455+(7.51396e-06*x))+(4.91857e-07*(x*x)))+(-1.47661e-09*(x*(x*x)))", 20.,670.);
	SFlmax = new TF1("SFlightMax","((1.15116+(0.00122657*x))+(-3.63826e-06*(x*x)))+(2.08242e-09*(x*(x*x)))", 20.,670.);
	}
	if( tagger == "CSVM" && eta>= 0.8 && eta < 1.6)
	{
	SFl = new TF1("SFlight","((1.111+(-9.64191e-06*x))+(1.80811e-07*(x*x)))+(-5.44868e-10*(x*(x*x)))", 20.,670.);
	SFlmin = new TF1("SFlightMin","((1.02055+(-0.000378856*x))+(1.49029e-06*(x*x)))+(-1.74966e-09*(x*(x*x)))", 20.,670.);
	SFlmax = new TF1("SFlightMax","((1.20146+(0.000359543*x))+(-1.12866e-06*(x*x)))+(6.59918e-10*(x*(x*x)))", 20.,670.);
	}
	if( tagger == "CSVM" && eta>= 1.6 && eta < 2.41)
	{
	SFl = new TF1("SFlight","((1.08498+(-0.000701422*x))+(3.43612e-06*(x*x)))+(-4.11794e-09*(x*(x*x)))", 20.,670.);
	SFlmin = new TF1("SFlightMin","((0.983476+(-0.000607242*x))+(3.17997e-06*(x*x)))+(-4.01242e-09*(x*(x*x)))", 20.,670.);
	SFlmax = new TF1("SFlightMax","((1.18654+(-0.000795808*x))+(3.69226e-06*(x*x)))+(-4.22347e-09*(x*(x*x)))", 20.,670.);
	}
	if( tagger == "CSVT" && eta>= 0.0 && eta < 2.41)
	{
	SFl = new TF1("SFlight","((0.948463+(0.00288102*x))+(-7.98091e-06*(x*x)))+(5.50157e-09*(x*(x*x)))", 20.,670.);
	SFlmin = new TF1("SFlightMin","((0.899715+(0.00102278*x))+(-2.46335e-06*(x*x)))+(9.71143e-10*(x*(x*x)))", 20.,670.);
	SFlmax = new TF1("SFlightMax","((0.997077+(0.00473953*x))+(-1.34985e-05*(x*x)))+(1.0032e-08*(x*(x*x)))", 20.,670.);
	}
	if( tagger == "JPM"  && eta>= 0.0 && eta < 0.8)
	{
	SFl = new TF1("SFlight","((0.970028+(0.00118179*x))+(-4.23119e-06*(x*x)))+(3.61065e-09*(x*(x*x)))", 20.,670.);
	SFlmin = new TF1("SFlightMin","((0.840326+(0.000626372*x))+(-2.08293e-06*(x*x)))+(1.57604e-09*(x*(x*x)))", 20.,670.);
	SFlmax = new TF1("SFlightMax","((1.09966+(0.00173739*x))+(-6.37946e-06*(x*x)))+(5.64527e-09*(x*(x*x)))", 20.,670.);
	}
	if( tagger == "JPM" && eta>= 1.6 && eta < 2.41)
	{
	SFl = new TF1("SFlight","((0.918387+(0.000898595*x))+(-2.00643e-06*(x*x)))+(1.26486e-09*(x*(x*x)))", 20.,670.);
	SFlmin = new TF1("SFlightMin","((0.790843+(0.000548016*x))+(-6.70941e-07*(x*x)))+(1.90355e-11*(x*(x*x)))", 20.,670.);
	SFlmax = new TF1("SFlightMax","((1.0459+(0.00124924*x))+(-3.34192e-06*(x*x)))+(2.51068e-09*(x*(x*x)))", 20.,670.);
	}
	if( tagger == "JPM" && eta>= 1.6 && eta < 2.41)
	{
	SFl = new TF1("SFlight","((0.790103+(0.00117865*x))+(-2.07334e-06*(x*x)))+(1.42608e-09*(x*(x*x)))", 20.,670.);
	SFlmin = new TF1("SFlightMin","((0.667144+(0.00105593*x))+(-1.43608e-06*(x*x)))+(5.24039e-10*(x*(x*x)))", 20.,670.);
	SFlmax = new TF1("SFlightMax","((0.913027+(0.00130143*x))+(-2.71061e-06*(x*x)))+(2.32812e-09*(x*(x*x)))", 20.,670.);
	}
	if( tagger == "JPT" && eta>= 0.0 && eta < 2.41)
	{
	SFl = new TF1("SFlight","((0.831392+(0.00269525*x))+(-7.33391e-06*(x*x)))+(5.73942e-09*(x*(x*x)))", 20.,670.);
	SFlmin = new TF1("SFlightMin","((0.671888+(0.0020106*x))+(-5.03177e-06*(x*x)))+(3.74225e-09*(x*(x*x)))", 20.,670.);
	SFlmax = new TF1("SFlightMax","((0.990774+(0.00338018*x))+(-9.63606e-06*(x*x)))+(7.73659e-09*(x*(x*x)))", 20.,670.);
	}
	//for overflow use inclusive bins
	if( tagger == "CSVM" && overflow)
	{
	SFl = new TF1("SFlight","((1.04318+(0.000848162*x))+(-2.5795e-06*(x*x)))+(1.64156e-09*(x*(x*x)))", 20.,670.);
	SFlmin = new TF1("SFlightMin","((0.962627+(0.000448344*x))+(-1.25579e-06*(x*x)))+(4.82283e-10*(x*(x*x)))", 20.,670.);
	SFlmax = new TF1("SFlightMax","((1.12368+(0.00124806*x))+(-3.9032e-06*(x*x)))+(2.80083e-09*(x*(x*x)))", 20.,670.);
	}

	if( tagger == "JPM" && overflow)
	{
	SFl = new TF1("SFlight","((0.871294+(0.00215201*x))+(-6.77675e-06*(x*x)))+(5.79197e-09*(x*(x*x)))", 20.,670.);
	SFlmin = new TF1("SFlightMin","((0.7654+(0.00149792*x))+(-4.47192e-06*(x*x)))+(3.67664e-09*(x*(x*x)))", 20.,670.);
	SFlmax = new TF1("SFlightMax","((0.977076+(0.00280638*x))+(-9.08158e-06*(x*x)))+(7.9073e-09*(x*(x*x)))", 20.,670.);
	}

	//correction function for 2012;
	if(tagger=="CSVM") SFcorrection = new TF1("SFcorrection","1.10422 + -0.000523856*x + 1.14251e-06*x*x", 20.,670.);
	if(tagger=="CSVT") SFcorrection = new TF1("SFcorrection","1.19275 + -0.00191042*x + 2.92205e-06*x*x" , 20.,670.);
	if(tagger=="JPM" ) SFcorrection = new TF1("SFcorrection","0.964487 + 0.00134038*x + -1.43995e-06*x*x", 20.,670.);
	if(tagger=="JPT" ) SFcorrection = new TF1("SFcorrection","1.06452 + 0.00080354*x + -1.34109e-06*x*x" , 20.,670.);

  //in 2012: correct with SFcorrection;
  float corr = SFcorrection->Eval(pt);
  float SF = corr* SFl->Eval(pt);
  float SFmin = corr* SFlmin->Eval(pt);
  float SFmax = corr* SFlmax->Eval(pt);
  if(overflow){
   //error is inflated by 2 --> difference between SF and SFmin/max is doubled
   // new SFmin = SF - 2*(SF-SFmin); new SFmax = SF + 2*(SFmax-SF);
   SFmin = 2.*SFmin - SF;
   SFmax = 2.*SFmax - SF;
  }
  //error prescription of 2012: relative inflate uncertainty by 1.5 for tight WP:
  // relative w.r.t. SF i.e. inflate by (1.5 * SF)
  // old rel error = (SF-SFmin)/SF or (SFmax-SF)/SF
  // new rel error = [(SF-SFmin)/SF]_new = 1.5 [(SF-SFmin)/SF]_old
  // --> SF-SFmin = 1.5 (SF-SFmin_old) --> SFmin = 1.5 SFmin_old - 0.5 SF
  // new SFmin = SF - (1.5*SF)*(SF-SFmin); new SFmax = SF + (1.5*SF)*(SFmax-SF);
   if(tagger=="CSVT" || tagger=="JPT"){
     //SFmin = SF - (1.5*SF)*(SF-SFmin);
     //SFmax = SF + (1.5*SF)*(SFmax-SF);
     SFmin = 1.5 * SFmin - 0.5 * SF;
     SFmax = 1.5 * SFmax - 0.5 * SF;
   }


  delete SFl; delete SFlmin; delete SFlmax, delete SFcorrection;
   if(up==1) return SFmax;
   else if(up==(-1)) return SFmin;
   return SF;
}

// get fastsim correction to jet SF
inline float FastSimCorrectionFactor(float &err, TString tagger, int flavour, float pt, float eta, string scan=""){
  float ptmin[14] = {30, 40, 50, 60, 70, 80, 100, 120, 160, 210, 260, 320, 400, 500};
  float ptmax[14] = {40, 50, 60, 70, 80,100, 120, 160, 210, 260, 320, 400, 500, 670};
  eta = fabs(eta);
  flavour = abs(flavour);
  float CFb[14];
  float CFb_error[14];
  float CFb_syst[14];
  bool overflow  = false;
  bool underflow = false;
  eta = fabs(eta);
  if(pt> ptmax[13]){ overflow = true; pt = 670;}
  if(pt< ptmin[0] ){underflow = true; pt = 30;}
  if(tagger=="CSVM"){
    if(flavour==5){
	float tCFb[14] = {0.982194,0.980998,0.992014,0.994472,0.996825,0.999822,1.00105,1.00023,0.991994,0.979123,0.947207,0.928006,0.874260,0.839610};
	float tCFb_error[14] = {0.00253112,0.00296453,0.00113963,0.00128363,0.00232566,0.00232353,0.00219086,0.00156856,0.00322279,0.00400414,0.00737465,0.0105033,0.0171706,0.0344172};
	float tCFb_T1_syst[14] = {0.0305974,0.0251517,0.0205015,0.0187029,0.0138344,0.0155380,0.0153906,0.0210581,0.0175900,-0.00234255,0.000241935,-0.0287645,-0.0472476,-0.0841584};
	float tCFb_T1bbbb_syst[14] = {0.00897405,0.00984249,0.00694051,0.00454724,0.00505632,0.00173861,0.00184828,0.00124377,-0.00265479,-0.0100402,-0.0112412,-0.0261436,-0.0221387,-0.0377308};
	float tCFb_T1tttt_syst[14] = {0.0112096,0.0127103,0.0107696,0.0105987,0.0102283,0.00953639,0.0107003,0.0118546,0.00837368,0.000790179,-0.00111371,-0.0146178,-0.00818416,-0.0197257};
	float tCFb_T2_syst[14] = {0.0197094,0.0218538,0.00671038,0.00481349,0.00234514,-0.00960910,-0.00872135,-0.0109075,-0.0185559,-0.0352550,-0.0374648,-0.0604555,-0.0752100,-0.0999645};
	float tCFb_T2bb_syst[14] = {0.0125569,0.0119411,0.0100657,0.0106521,0.00982046,0.00745928,0.00802320,0.00942034,0.00741357,0.00160137,0.00219074,-0.00892913,0.00172952,-0.000213087};
	float tCFb_T2bw_syst[14] = {0.0111744,0.0112791,0.00760594,0.00597137,0.00484192,0.00301468,0.00359970,0.00540084,0.00215334,-0.00427964,-0.00468144,-0.0184798,-0.0110016,-0.0187086};
	float tCFb_T2tt_syst[14] = {0.00574604,0.00677246,0.00509557,0.00374240,0.00314873,0.000637591,-0.000242591,-4.16636e-05,-0.00292352,-0.00581479,-0.000461876,-0.00676391,0.00488830,3.05474e-05};
	float tCFb_T3w_syst[14] = {0.0197131,0.0171196,0.0159192,0.0127636,0.0132435,0.00963777,0.00937313,0.00896174,0.00418186,-0.00353286,-0.00389037,-0.0171415,-0.0120094,-0.0215860};
	if(scan=="T1")     memcpy(CFb_syst, tCFb_T1_syst,    140);
	if(scan=="T1bbbb") memcpy(CFb_syst, tCFb_T1bbbb_syst,140);
	if(scan=="T1tttt") memcpy(CFb_syst, tCFb_T1tttt_syst,140);
	if(scan=="T2")     memcpy(CFb_syst, tCFb_T2_syst,    140);
	if(scan=="T2bb")   memcpy(CFb_syst, tCFb_T2bb_syst,  140);
	if(scan=="T2bW")   memcpy(CFb_syst, tCFb_T2bw_syst,  140);
	if(scan=="T2tt")   memcpy(CFb_syst, tCFb_T2tt_syst,  140);
	if(scan=="T3W")    memcpy(CFb_syst, tCFb_T3w_syst,   140);
	memcpy(CFb, tCFb,140);
	memcpy(CFb_error, tCFb_error,140);
    }//flavour==5
    else if(flavour==4){
	float tCFc[14] = {0.988545,0.981714,1.00946,1.01591,1.02810,1.02195,1.02590,1.01936,0.991228,0.955343,0.944433,0.917282,0.935018,1.06375};
	float tCFc_error[14] = {0.00746259,0.00661831,0.00968682,0.00751322,0.00675507,0.00562821,0.00862890,0.00768003,0.0188981,0.0261163,0.0450601,0.0448453,0.148805,0.177157}; // stat + PU
	float tCFc_T1_syst[14] = {-0.00104785,-0.00216204,0.00839430,0.00884008,0.00204188,0.00832790,0.0169768,0.0229713,0.0157189,-0.00730190,-0.0692086,-0.108517,-0.137035,-0.181932};
	float tCFc_T1bbbb_syst[14] = {-0.0169490,-0.0109324,-0.0173578,-0.0226300,-0.0354031,-0.0380664,-0.0406916,-0.0448566,-0.0634652,-0.0916214,-0.142743,-0.168372,-0.179460,-0.223442};
	float tCFc_T1tttt_syst[14] = {-0.00769350,0.00246567,0.00672805,0.00625175,0.0121922,0.0183616,0.0224260,0.0350031,0.0361672,0.0372230,0.0116431,0.0207569,0.0382855,0.0252644};
	float tCFc_T2_syst[14] = {-0.0147888,-0.00520468,-0.00901467,-0.0194454,-0.00635600,-0.00759417,-0.00953454,-0.0174082,-0.0184701,-0.0257653,-0.0740010,-0.0899951,-0.0860117,-0.0738075};
	float tCFc_T2bb_syst[14] = {-0.0126963,0.00240847,-0.0237588,-0.0202803,-0.0362858,-0.0296324,-0.0417801,-0.0566426,-0.0675621,-0.0768022,-0.141505,-0.160204,-0.199828,-0.237504};
	float tCFc_T2bw_syst[14] = {-0.0177944,-0.0168491,-0.0145971,-0.0171311,-0.0170042,-0.0143744,-0.0160470,-0.0149559,-0.0172561,-0.0137762,-0.0203696,0.00322482,0.0229054,0.0400957};
	float tCFc_T2tt_syst[14] = {-0.0183669,-0.0125071,-0.0174156,-0.0164738,-0.0167200,-0.0149260,-0.0180894,-0.0154648,-0.0141536,-0.0119079,-0.0206974,0.000753522,0.0221000,0.0209901};
	float tCFc_T3w_syst[14] = {-0.0126896,-0.00609615,-0.00314039,-0.00273418,-0.00209022,0.000352204,0.000533044,0.00463945,-0.000409096,-0.00550145,-0.0442329,-0.0519994,-0.0384817,-0.0126860};
	if(scan=="T1")     memcpy(CFb_syst, tCFc_T1_syst,    140);
	if(scan=="T1bbbb") memcpy(CFb_syst, tCFc_T1bbbb_syst,140);
	if(scan=="T1tttt") memcpy(CFb_syst, tCFc_T1tttt_syst,140);
	if(scan=="T2")     memcpy(CFb_syst, tCFc_T2_syst,    140);
	if(scan=="T2bb")   memcpy(CFb_syst, tCFc_T2bb_syst,  140);
	if(scan=="T2bW")   memcpy(CFb_syst, tCFc_T2bw_syst,  140);
	if(scan=="T2tt")   memcpy(CFb_syst, tCFc_T2tt_syst,  140);
	if(scan=="T3W")    memcpy(CFb_syst, tCFc_T3w_syst,   140);
	memcpy(CFb, tCFc,140);
	memcpy(CFb_error, tCFc_error,140);
    }//flavour==4
    else{//light jets
      if(eta<=1.2){
	float tCFudsg[14] = {1.21878,1.28615,1.37535,1.38966,1.40320,1.49835,1.44308,1.58198,1.55687,1.65790,1.90233,1.92259,2.66174,3.08688};
	float tCFudsg_error[14] = {0.0182686,0.0373732,0.0461870,0.0288973,0.0333528,0.0513836,0.0420353,0.106627,0.0658359,0.117285,0.185533,0.214071,0.487274,0.871502}; 
	float tCFudsg_T1_syst[14] = {0.386488,0.524838,0.679631,0.682134,0.731629,0.757618,0.695844,0.724127,0.623855,0.660598,0.829350,0.905624,1.40528,1.87998};
	float tCFudsg_T1bbbb_syst[14] = {0.269423,0.377897,0.460441,0.456512,0.470195,0.481113,0.438144,0.464775,0.347294,0.411220,0.550301,0.623299,1.14485,1.53694};
	float tCFudsg_T1tttt_syst[14] = {0.118477,0.162964,0.223318,0.220063,0.222306,0.267305,0.222287,0.283804,0.252221,0.324747,0.527015,0.659528,1.19317,1.50547};
	float tCFudsg_T2_syst[14] = {0.411044,0.558016,0.670788,0.776938,0.802923,0.895418,0.806768,0.812508,0.671626,0.715303,0.865253,0.889535,1.36529,1.83958};
	float tCFudsg_T2bb_syst[14] = {0.419325,0.559732,0.664588,0.701051,0.760269,0.876007,0.749539,0.827054,0.657627,0.702294,0.858618,0.837998,1.36137,1.75727};
	float tCFudsg_T2bw_syst[14] = {0.249750,0.331114,0.370544,0.380683,0.375024,0.433907,0.370687,0.430421,0.397210,0.479439,0.675053,0.815746,1.36142,1.86164};
	float tCFudsg_T2tt_syst[14] = {0.241447,0.297617,0.365921,0.372697,0.378869,0.434225,0.385061,0.452832,0.412124,0.498940,0.675028,0.813003,1.31961,1.57929};
	float tCFudsg_T3w_syst[14] = {0.287246,0.388381,0.480550,0.504640,0.531340,0.572774,0.532622,0.586227,0.529575,0.600174,0.792410,0.882505,1.42788,1.91256};
	if(scan=="T1")     memcpy(CFb_syst, tCFudsg_T1_syst,    140);
	if(scan=="T1bbbb") memcpy(CFb_syst, tCFudsg_T1bbbb_syst,140);
	if(scan=="T1tttt") memcpy(CFb_syst, tCFudsg_T1tttt_syst,140);
	if(scan=="T2")     memcpy(CFb_syst, tCFudsg_T2_syst,    140);
	if(scan=="T2bb")   memcpy(CFb_syst, tCFudsg_T2bb_syst,  140);
	if(scan=="T2bW")   memcpy(CFb_syst, tCFudsg_T2bw_syst,  140);
	if(scan=="T2tt")   memcpy(CFb_syst, tCFudsg_T2tt_syst,  140);
	if(scan=="T3W")    memcpy(CFb_syst, tCFudsg_T3w_syst,   140);
	memcpy(CFb, tCFudsg,140);
	memcpy(CFb_error, tCFudsg_error,140);
      }
      else{
	float tCFudsg[14] = {1.46970,1.48732,1.69024,1.64494,1.79297,1.90760,1.99867,2.21659,2.20103,2.42645,2.67594,4.24735,3.98979,15.0457};
	float tCFudsg_error[14] = {0.104716,0.0392025,0.106315,0.115751,0.106807,0.0642086,0.138742,0.182345,0.169922,0.297889,0.320088,0.927736,1.24666,15.1860};
	float tCFudsg_T1_syst[14] = {0.686015,0.720420,0.991786,1.02970,1.22030,1.34286,1.64405,2.09951,2.92008,4.41435,6.25081,9.86965,11.9982,25.3907};
	float tCFudsg_T1bbbb_syst[14] = {0.618506,0.635309,0.794568,0.803646,0.886742,0.988190,1.10888,1.31924,1.49423,1.96107,2.36883,3.67770,4.80525,11.2853};
	float tCFudsg_T1tttt_syst[14] = {0.772327,0.874528,1.19814,1.24806,1.49608,1.73841,2.00430,2.54257,3.27898,4.35726,5.31846,7.44186,9.19039,15.6896};
	float tCFudsg_T2_syst[14] = {0.549720,0.580865,0.765356,0.788296,0.913463,1.03193,1.19510,1.40819,1.77895,2.69320,3.44912,5.91765,9.20944,18.0392};
	float tCFudsg_T2bb_syst[14] = {0.569517,0.548840,0.763820,0.726966,0.927079,0.959964,1.09951,1.29263,1.39602,1.97896,2.41141,3.67147,4.17557,11.7192};
	float tCFudsg_T2bw_syst[14] = {0.697431,0.759470,1.03429,1.05697,1.23279,1.38067,1.52550,1.89634,2.29738,2.87713,3.64427,5.54452,6.93274,13.9094};
	float tCFudsg_T2tt_syst[14] = {0.694594,0.753930,1.01105,1.02488,1.18455,1.36025,1.57676,1.87545,2.24691,2.81635,3.46050,5.75946,6.89900,15.8855};
	float tCFudsg_T3w_syst[14] = {0.773925,0.839802,1.08844,1.16056,1.35051,1.53349,1.76781,2.23936,2.95149,4.32639,5.98244,9.13821,12.3203,24.2016};
	if(scan=="T1")     memcpy(CFb_syst, tCFudsg_T1_syst,    140);
	if(scan=="T1bbbb") memcpy(CFb_syst, tCFudsg_T1bbbb_syst,140);
	if(scan=="T1tttt") memcpy(CFb_syst, tCFudsg_T1tttt_syst,140);
	if(scan=="T2")     memcpy(CFb_syst, tCFudsg_T2_syst,    140);
	if(scan=="T2bb")   memcpy(CFb_syst, tCFudsg_T2bb_syst,  140);
	if(scan=="T2bW")   memcpy(CFb_syst, tCFudsg_T2bw_syst,  140);
	if(scan=="T2tt")   memcpy(CFb_syst, tCFudsg_T2tt_syst,  140);
	if(scan=="T3W")    memcpy(CFb_syst, tCFudsg_T3w_syst,   140);
	memcpy(CFb, tCFudsg,140);
	memcpy(CFb_error, tCFudsg_error,140);
      }
    }//light jets
  }//CSVM
  if(tagger=="CSVT"){
    if(flavour==5){
	float tCFb[14] = {0.968117,0.967822,0.978278,0.981281,0.987679,0.986590,0.990246,0.984504,0.967024,0.940042,0.873019,0.850847,0.769561,0.650192};
	float tCFb_error[14] = {0.00223422,0.00367427,0.00145554,0.00337572,0.00344106,0.00591257,0.00218050,0.00472939,0.00353119,0.00739502,0.0193330,0.0158257,0.0306048,0.0603701}; // stat + PU
	float tCFb_T1_syst[14] = {0.0404417,0.0338855,0.0268602,0.0154959,0.0112958,0.00860152,0.00871485,0.00643355,-0.0189310,-0.0548930,-0.0645469,-0.110725,-0.135682,-0.123968};
	float tCFb_T1bbbb_syst[14] = {0.0129219,0.0137310,0.0113178,0.00933427,0.00921154,0.00691691,0.00647020,0.00674507,-0.00834304,-0.0302220,-0.0412711,-0.0777479,-0.0823319,-0.0629126};
	float tCFb_T1tttt_syst[14] = {0.0151952,0.0154848,0.0141975,0.0152467,0.0151888,0.0164594,0.0190932,0.0238398,0.0151870,-0.00265440,-0.0101969,-0.0403102,-0.0398199,-0.0128917};
	float tCFb_T2_syst[14] = {0.0335889,0.0293072,0.00701576,0.00240481,0.00541189,-0.0178958,-0.0182847,-0.0249884,-0.0600149,-0.0852849,-0.108664,-0.139017,-0.160629,-0.136853};
	float tCFb_T2bb_syst[14] = {0.0189191,0.0207380,0.0180343,0.0183147,0.0184177,0.0161918,0.0175651,0.0203825,0.0112045,-0.0102938,-0.0205988,-0.0525127,-0.0512335,-0.0153682};
	float tCFb_T2bw_syst[14] = {0.0163777,0.0155786,0.0114890,0.0101303,0.00891703,0.00743636,0.00856078,0.0129112,0.00174296,-0.0183686,-0.0265031,-0.0612646,-0.0610615,-0.0313782};
	float tCFb_T2tt_syst[14] = {0.00601965,0.00787908,0.00701177,0.00568053,0.00516674,0.00290364,0.00112420,0.00258505,-0.00883520,-0.0173537,-0.00851473,-0.0216886,-0.00751114,0.0275140};
	float tCFb_T3w_syst[14] = {0.0284607,0.0245074,0.0227164,0.0191328,0.0193200,0.0166523,0.0169420,0.0180767,0.00327435,-0.0184413,-0.0281253,-0.0619053,-0.0645325,-0.0373348};
	if(scan=="T1")     memcpy(CFb_syst, tCFb_T1_syst,    140);
	if(scan=="T1bbbb") memcpy(CFb_syst, tCFb_T1bbbb_syst,140);
	if(scan=="T1tttt") memcpy(CFb_syst, tCFb_T1tttt_syst,140);
	if(scan=="T2")     memcpy(CFb_syst, tCFb_T2_syst,    140);
	if(scan=="T2bb")   memcpy(CFb_syst, tCFb_T2bb_syst,  140);
	if(scan=="T2bW")   memcpy(CFb_syst, tCFb_T2bw_syst,  140);
	if(scan=="T2tt")   memcpy(CFb_syst, tCFb_T2tt_syst,  140);
	if(scan=="T3W")    memcpy(CFb_syst, tCFb_T3w_syst,   140);
	memcpy(CFb, tCFb,140);
	memcpy(CFb_error, tCFb_error,140);
    }//flavour==5
    else if(flavour==4){
	float tCFc[14] = {0.960959,0.973876,0.984323,0.996344,1.02418,0.985580,0.994745,0.970489,0.914155,0.872072,0.945289,0.783816,0.942773,0.527354};
	float tCFc_error[14] = {0.0155733,0.0121900,0.0131678,0.0113739,0.0213937,0.0123294,0.0153230,0.0156350,0.0409568,0.0654966,0.112785,0.187795,0.331301,0.162462}; // stat + PU
	float tCFc_T1_syst[14] = {-0.0150188,-0.0463591,-0.0315633,-0.0258283,-0.0327155,-0.0267018,-0.0148436,-0.0318809,-0.0359985,-0.0996782,-0.199801,-0.181507,-0.214658,-0.0237959};
	float tCFc_T1bbbb_syst[14] = {-0.0401338,-0.0623124,-0.0551870,-0.0723321,-0.0831102,-0.0953026,-0.0870973,-0.100656,-0.138856,-0.176190,-0.267756,-0.235408,-0.300116,-0.0866552};
	float tCFc_T1tttt_syst[14] = {-0.00762913,-0.0119695,0.00505692,0.00110212,0.00860715,0.0268639,0.0194188,0.0441512,0.0267859,0.0139230,-0.0531917,0.0474467,0.0777428,0.324803};
	float tCFc_T2_syst[14] = {-0.0210332,-0.0259317,-0.0103462,-0.0366923,-0.0190886,0.00934875,-0.00387061,-0.000786334,-0.0159257,-0.0524397,-0.116956,-0.0820041,-0.0804845,0.241786};
	float tCFc_T2bb_syst[14] = {-0.0376442,-0.0346557,-0.0518969,-0.0535199,-0.0908062,-0.0842518,-0.0974765,-0.105002,-0.117947,-0.151280,-0.254247,-0.216998,-0.284635,-0.0142322};
	float tCFc_T2bw_syst[14] = {-0.0175963,-0.0338167,-0.0177681,-0.0222040,-0.0251838,-0.0129322,-0.0205147,-0.0147394,-0.0377501,-0.0349037,-0.0639912,0.0514183,0.113051,0.453520};
	float tCFc_T2tt_syst[14] = {-0.0179060,-0.0272217,-0.0173259,-0.0195665,-0.0235815,-0.0123463,-0.0248161,-0.0109724,-0.0342969,-0.0314718,-0.0592955,0.0535336,0.124457,0.367309};
	float tCFc_T3w_syst[14] = {-0.0114200,-0.0202786,-0.00171709,-0.00232449,-0.00353348,0.00844371,-0.000867724,0.0146573,-0.0153649,-0.0254736,-0.0965763,-0.0415070,-0.00707200,0.310379};
	if(scan=="T1")     memcpy(CFb_syst, tCFc_T1_syst,    140);
	if(scan=="T1bbbb") memcpy(CFb_syst, tCFc_T1bbbb_syst,140);
	if(scan=="T1tttt") memcpy(CFb_syst, tCFc_T1tttt_syst,140);
	if(scan=="T2")     memcpy(CFb_syst, tCFc_T2_syst,    140);
	if(scan=="T2bb")   memcpy(CFb_syst, tCFc_T2bb_syst,  140);
	if(scan=="T2bW")   memcpy(CFb_syst, tCFc_T2bw_syst,  140);
	if(scan=="T2tt")   memcpy(CFb_syst, tCFc_T2tt_syst,  140);
	if(scan=="T3W")    memcpy(CFb_syst, tCFc_T3w_syst,   140);
	memcpy(CFb, tCFc,140);
	memcpy(CFb_error, tCFc_error,140);
    }//flavour==4
    else{//light jets
      if(eta<=1.2){
	float tCFudsg[14] = {1.24890,1.35145,1.37205,1.32472,1.39976,1.45884,1.59912,1.58971,1.30841,1.55936,1.28346,2.21265,2.06927,2.88109};
	float tCFudsg_error[14] = {0.0751438,0.0651619,0.0604241,0.0726285,0.0968158,0.0931768,0.163039,0.187749,0.198200,0.465354,0.339473,1.07079,1.07723,2.53188}; // stat + PU
	float tCFudsg_T1_syst[14] = {0.541349,0.749937,0.899876,0.913299,1.04832,1.00971,0.948618,0.828252,0.461639,0.515000,0.339278,0.834556,0.922787,1.72162};
	float tCFudsg_T1bbbb_syst[14] = {0.262307,0.382883,0.463978,0.397656,0.450687,0.470582,0.514196,0.478655,0.229171,0.438639,0.200500,0.754263,0.887022,1.60800};
	float tCFudsg_T1tttt_syst[14] = {0.107974,0.146424,0.175264,0.125604,0.148329,0.175837,0.216767,0.212357,0.0737581,0.268330,0.220413,0.810560,0.954838,1.67344};
	float tCFudsg_T2_syst[14] = {0.514257,0.722493,0.824915,1.30112,1.38608,1.41843,1.04881,0.975846,0.539212,0.606663,0.415293,0.763953,0.812614,1.59538};
	float tCFudsg_T2bb_syst[14] = {0.516994,0.673320,0.844278,0.855887,1.08265,1.10285,1.20266,0.989551,0.686524,0.782166,0.602174,0.973691,0.905448,1.66941};
	float tCFudsg_T2bw_syst[14] = {0.284863,0.382788,0.346198,0.313557,0.329289,0.356919,0.347319,0.354319,0.164351,0.357468,0.301092,0.901424,1.13151,1.92501};
	float tCFudsg_T2tt_syst[14] = {0.258683,0.321056,0.326149,0.318511,0.344513,0.360544,0.380937,0.390651,0.208644,0.402611,0.327629,0.955151,1.09279,1.85737};
	float tCFudsg_T3w_syst[14] = {0.402833,0.489644,0.564344,0.587424,0.657689,0.628729,0.663659,0.617975,0.328334,0.460770,0.340796,0.833781,0.972206,1.83632};
	if(scan=="T1")     memcpy(CFb_syst, tCFudsg_T1_syst,    140);
	if(scan=="T1bbbb") memcpy(CFb_syst, tCFudsg_T1bbbb_syst,140);
	if(scan=="T1tttt") memcpy(CFb_syst, tCFudsg_T1tttt_syst,140);
	if(scan=="T2")     memcpy(CFb_syst, tCFudsg_T2_syst,    140);
	if(scan=="T2bb")   memcpy(CFb_syst, tCFudsg_T2bb_syst,  140);
	if(scan=="T2bW")   memcpy(CFb_syst, tCFudsg_T2bw_syst,  140);
	if(scan=="T2tt")   memcpy(CFb_syst, tCFudsg_T2tt_syst,  140);
	if(scan=="T3W")    memcpy(CFb_syst, tCFudsg_T3w_syst,   140);
	memcpy(CFb, tCFudsg,140);
	memcpy(CFb_error, tCFudsg_error,140);
      }
      else{
	float tCFudsg[14] = {1.67634,1.70105,1.75999,1.78459,2.19343,2.73199,3.49277,2.58863,2.48824,4.01723,3.86956,0.000456049,2.30988,0.000855693};
	float tCFudsg_error[14] = {0.222165,0.161403,0.112342,0.275101,0.364229,0.330588,1.00953,0.404417,1.07731,2.65686,3.18286,5.25051e-05,2.38652,0.000438728}; // stat + PU
	float tCFudsg_T1_syst[14] = {1.12564,1.13767,1.41534,1.51698,2.21567,2.46153,3.52618,3.66658,6.82062,9.53548,11.6888,11.6385,37.8436,55.0746};
	float tCFudsg_T1bbbb_syst[14] = {0.966686,0.911206,1.01981,1.05643,1.36584,1.83112,2.31922,1.90143,3.29097,3.21403,4.94706,3.62533,10.3981,5.51893};
	float tCFudsg_T1tttt_syst[14] = {1.07420,1.13792,1.37766,1.43242,1.81959,2.25043,2.85233,3.01381,3.90803,7.90416,7.84521,8.13804,21.4383,90.5637};
	float tCFudsg_T2_syst[14] = {0.914312,1.00641,1.19010,1.21705,1.54337,1.79515,2.79849,2.58661,3.30393,5.87701,6.11172,9.81007,66.2097,110.301};
	float tCFudsg_T2bb_syst[14] = {0.952269,0.947869,1.07177,1.07862,1.61536,1.89066,2.38571,2.09688,2.12828,3.69276,4.06967,3.12065,4.10128,-0.649336};
	float tCFudsg_T2bw_syst[14] = {1.06254,1.07746,1.35255,1.45264,1.89413,2.11516,2.95536,2.60233,3.42582,4.82154,7.41456,6.31043,6.52546,13.4704};
	float tCFudsg_T2tt_syst[14] = {1.01487,1.07556,1.33002,1.38201,1.77121,2.12163,2.80742,2.57413,3.43778,4.60502,5.42910,4.53808,6.23334,23.3936};
	float tCFudsg_T3w_syst[14] = {1.17544,1.23456,1.48954,1.57960,2.18616,2.55298,3.25695,3.47147,5.16497,7.76216,11.6997,10.8320,19.7364,64.8392};
	if(scan=="T1")     memcpy(CFb_syst, tCFudsg_T1_syst,    140);
	if(scan=="T1bbbb") memcpy(CFb_syst, tCFudsg_T1bbbb_syst,140);
	if(scan=="T1tttt") memcpy(CFb_syst, tCFudsg_T1tttt_syst,140);
	if(scan=="T2")     memcpy(CFb_syst, tCFudsg_T2_syst,    140);
	if(scan=="T2bb")   memcpy(CFb_syst, tCFudsg_T2bb_syst,  140);
	if(scan=="T2bW")   memcpy(CFb_syst, tCFudsg_T2bw_syst,  140);
	if(scan=="T2tt")   memcpy(CFb_syst, tCFudsg_T2tt_syst,  140);
	if(scan=="T3W")    memcpy(CFb_syst, tCFudsg_T3w_syst,   140);
	memcpy(CFb, tCFudsg,140);
	memcpy(CFb_error, tCFudsg_error,140);
      }
    }//light jets
  }//CSVT

  float CFvalue  = 0;
  float CFvalue_err  = 0;
  float CFvalue_syst  = 0;
  for(int i=0;i<14;i++){
    if( pt >= ptmin[i] && pt <= ptmax[i])
      CFvalue = CFb[i];
      CFvalue_err = CFb_error[i];
      CFvalue_syst = CFb_syst[i];
  }
  if(overflow){
	CFvalue_err = 2.*CFvalue_err;
	CFvalue_syst= 2.*CFvalue_syst;
  }
  if(scan=="T1"||scan=="T1bbbb"||scan=="T1tttt"||scan=="T2"||scan=="T2bb"||scan=="T2bW"||scan=="T2tt"||scan=="T3W") 
	CFvalue_err = sqrt(CFvalue_err*CFvalue_err + CFvalue_syst*CFvalue_syst);
  err = CFvalue_err;
  return CFvalue;
}


//the monster function with gets you SF + total error due to all possible systematics
//uses getBTagEventWeightError and getBTagEventWeight from above
inline Float_t getBTagEventWeightErrorTotal(float &werr, vector<float> jetEff, vector<float> jetEffErr, vector<float> jetSF, vector<float> jetSFBCdown, vector<float> jetSFBCup, vector<float> jetSFlightdown, vector<float> jetSFlightup, vector<float> jetSFFSdown, vector<float> jetSFFSup, int nBTags){
	float statefferr = 0;
//	float weight = getBTagEventWeightError(statefferr, jetEff, jetEffErr, jetSF, nBTags);
	vector<float> jetEffup; jetEffup.clear(); vector<float> jetEffdown; jetEffdown.clear();
	for(unsigned int i = 0; i<jetEff.size(); ++i){
		jetEffup.push_back(jetEff[i]+jetEffErr[i]);
		jetEffdown.push_back(jetEff[i]-jetEffErr[i]);
	}
	float weight = getBTagEventWeightError(jetEff, jetSF, nBTags);
	float staterrup = getBTagEventWeightError(jetEffup, jetSF, nBTags);
	float staterrdown = getBTagEventWeightError(jetEffdown, jetSF, nBTags);

	float bcerrorup = fabs(weight-getBTagEventWeight(jetEff, jetSFBCup, nBTags));// b error
	float bcerrordown = fabs(weight-getBTagEventWeight(jetEff, jetSFBCdown, nBTags));// b error
	float lighterrorup = fabs(weight-getBTagEventWeight(jetEff, jetSFlightup, nBTags));// mistag error
	float lighterrordown = fabs(weight-getBTagEventWeight(jetEff, jetSFlightdown, nBTags));// mistag error
	float FSerrorup = fabs(weight-getBTagEventWeight(jetEff, jetSFFSup, nBTags));//fastsim error
	float FSerrordown = fabs(weight-getBTagEventWeight(jetEff, jetSFFSdown, nBTags));//fastsim error

	float bcerr, lerr, fserr, sterr;
//	sterr = statefferr;
	if(staterrup>staterrdown) sterr = staterrup;
	else sterr = staterrdown;
	if(bcerrorup>bcerrordown) bcerr = bcerrorup;
	else bcerr = bcerrordown;
	if(lighterrorup>lighterrordown) lerr = lighterrorup;
	else lerr = lighterrordown;
	if(FSerrorup>FSerrordown) fserr = FSerrorup;
	else fserr = FSerrordown;
	werr = sqrt(sterr*sterr + bcerr*bcerr + lerr*lerr + fserr*fserr);
	return weight;
}

//the monster function with gets you SF + total error due to all possible systematics
// is double sided error up and down
//uses getBTagEventWeightError and getBTagEventWeight from above
inline Float_t getBTagEventWeightTwoSidedErrorTotal(float &werrup, float &werrdown, vector<float> jetEff, vector<float> jetEffErr, vector<float> jetSF, vector<float> jetSFBCdown, vector<float> jetSFBCup, vector<float> jetSFlightdown, vector<float> jetSFlightup, vector<float> jetSFFSdown, vector<float> jetSFFSup, int nBTags){
	float statefferr = 0;
//	float weight = getBTagEventWeightError(statefferr, jetEff, jetEffErr, jetSF, nBTags);
//	float staterrdown = statefferr; float staterrup = statefferr;
	vector<float> jetEffup; jetEffup.clear(); vector<float> jetEffdown; jetEffdown.clear();
	for(unsigned int i = 0; i<jetEff.size(); ++i){
		jetEffup.push_back(jetEff[i]+jetEffErr[i]);
		jetEffdown.push_back(jetEff[i]-jetEffErr[i]);
	}
	float weight = getBTagEventWeightError(jetEff, jetSF, nBTags);
	float staterrup = getBTagEventWeightError(jetEffup, jetSF, nBTags);
	float staterrdown = getBTagEventWeightError(jetEffdown, jetSF, nBTags);

	float bcerrorup = fabs(weight-getBTagEventWeight(jetEff, jetSFBCup, nBTags));// b error
	float bcerrordown = fabs(weight-getBTagEventWeight(jetEff, jetSFBCdown, nBTags));// b error
	float lighterrorup = fabs(weight-getBTagEventWeight(jetEff, jetSFlightup, nBTags));// mistag error
	float lighterrordown = fabs(weight-getBTagEventWeight(jetEff, jetSFlightdown, nBTags));// mistag error
	float FSerrorup = fabs(weight-getBTagEventWeight(jetEff, jetSFFSup, nBTags));//fastsim error
	float FSerrordown = fabs(weight-getBTagEventWeight(jetEff, jetSFFSdown, nBTags));//fastsim error

	werrup   = sqrt(staterrup  *staterrup   + bcerrorup  *bcerrorup   + lighterrorup  *lighterrorup   + FSerrorup  *FSerrorup  );
	werrdown = sqrt(staterrdown*stateffdown + bcerrordown*bcerrordown + lighterrordown*lighterrordown + FSerrordown*FSerrordown);
	return weight;
}



#endif
