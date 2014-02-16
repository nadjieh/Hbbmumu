//#include <stdio.h>
{
#include <stdio.h>

    gROOT->ProcessLine(".x SetStyle_PRD.C");

  TString DIR = "../Histos/";

  TString B="1b";
  TString outF_name = "histos_ZInv-MT2_-1_9999-Jet40-HT_800-1Mu_pt10_4450pb_"+B;

  TFile *infile = TFile::Open(DIR+outF_name+".root");
  
  const int totS = 100;
  TString names[totS], legs[totS];
  TCanvas *c[100];
  int nc=0;
  double Ws[totS];
  int colors[totS];
  int nf = 0;
  double LUMI = 1;
  double QCD_K = 1.1216;

  
  
  colors[nf] = 0; legs[nf] = "data"; names[nf++] = "data"; 
  colors[nf] = 401;  legs[nf] = "QCD"; names[nf++] = "QCD";
  colors[nf] = 416;  legs[nf] = "Drell-Yan"; names[nf++] = "DY";
  colors[nf] = 418; legs[nf] = "W(#tau#nu)+jets"; names[nf++] = "WtolnuTau";
  colors[nf] = 602; legs[nf] = "Top";  names[nf++] = "Top";
  colors[nf] = 419; legs[nf] = "W(#mu#nu)+jets"; names[nf++] = "WtolnuMu";
  colors[nf] = 417; legs[nf] = "W(e#nu)+jets"; names[nf++] = "WtolnuEle";
  
  
  TH1F * MT2_bSel[totS], * MT2_leptSel[totS];
  TH1F *Lept_pt_leptSel[totS], *Lept_eta_leptSel[totS];
  TH1F * BJets_N_leptSel[totS];
  TH1F* Jets_N_leptSel[totS];
  TH1F* MET_leptSel[totS];
  TH1F* HT_leptSel[totS];
  TH1F* MinDPhi_leptSel[totS];

  TH1F *Lept_pt_bSel[totS], *Lept_eta_bSel[totS];
  TH1F* Jets_N_bSel[totS];
  TH1F* MET_bSel[totS];
  TH1F* HT_bSel[totS];
  TH1F* MinDPhi_bSel[totS];


  for(int i=0; i< nf; i++){
      
    infile->GetObject("MT2_"+names[i]+"_bSel", MT2_bSel[i] );
    if(i==0) MT2_bSel[i]->SetMarkerStyle(20);  MT2_bSel[i]->SetFillColor( colors[i]); MT2_bSel[i]->Rebin(2); MT2_bSel[i]->Scale(LUMI);
    MT2_bSel[i]->GetXaxis()->SetRangeUser(0,1000);

    //infile->GetObject("MT2_"+names[i]+"_bSel", MT2_bSel[i] );
    //if(i==0) MT2_bSel[i]->SetMarkerStyle(20);  MT2_bSel[i]->SetFillColor( colors[i]); MT2_bSel[i]->Rebin(2); MT2_bSel[i]->Scale(LUMI);

    /*        infile->GetObject("Lept_pt_bSel_"+names[i], Lept_pt_bSel[i] );
    if(i==0) Lept_pt_bSel[i]->SetMarkerStyle(20);  Lept_pt_bSel[i]->SetFillColor( colors[i]); Lept_pt_bSel[i]->Rebin(5); Lept_pt_bSel[i]->Scale(LUMI);
    
    infile->GetObject("Lept_eta_bSel_"+names[i], Lept_eta_bSel[i] );
    if(i==0) Lept_eta_bSel[i]->SetMarkerStyle(20);  Lept_eta_bSel[i]->SetFillColor( colors[i]); Lept_eta_bSel[i]->Rebin(2); Lept_eta_bSel[i]->Scale(LUMI);
    */
    infile->GetObject("Jets_N_"+names[i]+"_bSel", Jets_N_bSel[i] );
    if(i==0) Jets_N_bSel[i]->SetMarkerStyle(20);  Jets_N_bSel[i]->SetFillColor( colors[i]); Jets_N_bSel[i]->Rebin(1); Jets_N_bSel[i]->Scale(LUMI);
    Jets_N_bSel[i]->GetXaxis()->SetRangeUser(2,10);

    infile->GetObject("MET_"+names[i]+"_bSel", MET_bSel[i] );
    if(i==0) MET_bSel[i]->SetMarkerStyle(20);  MET_bSel[i]->SetFillColor( colors[i]); MET_bSel[i]->Rebin(4); MET_bSel[i]->Scale(LUMI);
    MET_bSel[i]->GetXaxis()->SetRangeUser(0,1000);

    infile->GetObject("HT_"+names[i]+"_bSel", HT_bSel[i] );
    if(i==0) HT_bSel[i]->SetMarkerStyle(20);  HT_bSel[i]->SetFillColor( colors[i]); HT_bSel[i]->Rebin(4); HT_bSel[i]->Scale(LUMI);
    HT_bSel[i]->GetXaxis()->SetRangeUser(700,2000);

    infile->GetObject("MinDPhi_"+names[i]+"_bSel", MinDPhi_bSel[i] );
    if(i==0) MinDPhi_bSel[i]->SetMarkerStyle(20);  MinDPhi_bSel[i]->SetFillColor( colors[i]); MinDPhi_bSel[i]->Rebin(4); MinDPhi_bSel[i]->Scale(LUMI);
    



    infile->GetObject("MT2_"+names[i]+"_leptSel", MT2_leptSel[i] );
    if(i==0) MT2_leptSel[i]->SetMarkerStyle(20);  MT2_leptSel[i]->SetFillColor( colors[i]); MT2_leptSel[i]->Rebin(2); MT2_leptSel[i]->Scale(LUMI);
    MT2_leptSel[i]->GetXaxis()->SetRangeUser(0,1000);

    infile->GetObject("Lept_pt_leptSel_"+names[i], Lept_pt_leptSel[i] );
    if(i==0) Lept_pt_leptSel[i]->SetMarkerStyle(20);  Lept_pt_leptSel[i]->SetFillColor( colors[i]); Lept_pt_leptSel[i]->Rebin(5); Lept_pt_leptSel[i]->Scale(LUMI);
    Lept_pt_leptSel[i]->GetXaxis()->SetRangeUser(0,1000);

    infile->GetObject("Lept_eta_leptSel_"+names[i], Lept_eta_leptSel[i] );
    if(i==0) Lept_eta_leptSel[i]->SetMarkerStyle(20);  Lept_eta_leptSel[i]->SetFillColor( colors[i]); Lept_eta_leptSel[i]->Rebin(2); Lept_eta_leptSel[i]->Scale(LUMI);

    infile->GetObject("BJets_N_"+names[i]+"_leptSel", BJets_N_leptSel[i] );
    if(i==0) BJets_N_leptSel[i]->SetMarkerStyle(20);  BJets_N_leptSel[i]->SetFillColor( colors[i]); BJets_N_leptSel[i]->Rebin(1); BJets_N_leptSel[i]->Scale(LUMI);
    BJets_N_leptSel[i]->GetXaxis()->SetRangeUser(0,6);

    infile->GetObject("Jets_N_"+names[i]+"_leptSel", Jets_N_leptSel[i] );
    if(i==0) Jets_N_leptSel[i]->SetMarkerStyle(20);  Jets_N_leptSel[i]->SetFillColor( colors[i]); Jets_N_leptSel[i]->Rebin(1); Jets_N_leptSel[i]->Scale(LUMI);
    Jets_N_leptSel[i]->GetXaxis()->SetRangeUser(2,10);

    infile->GetObject("MET_"+names[i]+"_leptSel", MET_leptSel[i] );
    if(i==0) MET_leptSel[i]->SetMarkerStyle(20);  MET_leptSel[i]->SetFillColor( colors[i]); MET_leptSel[i]->Rebin(4); MET_leptSel[i]->Scale(LUMI);
    MET_leptSel[i]->GetXaxis()->SetRangeUser(0,1000);

    infile->GetObject("HT_"+names[i]+"_leptSel", HT_leptSel[i] );
    if(i==0) HT_leptSel[i]->SetMarkerStyle(20);  HT_leptSel[i]->SetFillColor( colors[i]); HT_leptSel[i]->Rebin(4); HT_leptSel[i]->Scale(LUMI);
    HT_leptSel[i]->GetXaxis()->SetRangeUser(700,2000);

    infile->GetObject("MinDPhi_"+names[i]+"_leptSel", MinDPhi_leptSel[i] );
    if(i==0) MinDPhi_leptSel[i]->SetMarkerStyle(20);  MinDPhi_leptSel[i]->SetFillColor( colors[i]); MinDPhi_leptSel[i]->Rebin(4); MinDPhi_leptSel[i]->Scale(LUMI);

  }



  TLegend *l = new TLegend(0.7,0.7,0.9,0.9);
  l->AddEntry(  MT2_bSel[0], legs[0], "p");
  for(int i=1; i<nf; i++) l->AddEntry(  MT2_bSel[i], legs[i], "fp");
  l->SetFillColor(000000);
  l->SetBorderSize(0);
  

  //bsel
  if(B=="1b"){
    c[nc++] = new TCanvas(outF_name+"-MT2_bSel",outF_name+"-MT2_bSel" );  gPad->SetLogy();
    THStack s_MT2_bSel("s_MT2_bSel","s_MT2_bSel");
    for(int i=1; i<nf; i++){     s_MT2_bSel.Add(MT2_bSel[i]);  }
    MT2_bSel[0].GetXaxis()->SetTitle("M_{T2} (GeV)"); MT2_bSel[0]->GetYaxis()->SetTitle("events"); 
    MT2_bSel[0]->Draw("ep");   s_MT2_bSel.Draw("hsames");  MT2_bSel[0]->Draw("samese");   l->Draw();
    /*
      c[nc++] = new TCanvas(outF_name+"-Lept_pt_bSel",outF_name+"-Lept_pt_bSel" );  gPad->SetLogy();
      THStack s_Lept_pt_bSel("s_Lept_pt_bSel","s_Lept_pt_bSel");
      for(int i=1; i<nf; i++){    s_Lept_pt_bSel.Add(Lept_pt_bSel[i]); }
      Lept_pt_bSel[0].GetXaxis()->SetTitle("M_{T2} (GeV)"); Lept_pt_bSel[0]->GetYaxis()->SetTitle("events");
      Lept_pt_bSel[0]->Draw("ep");  s_Lept_pt_bSel.Draw("hsames");  Lept_pt_bSel[0]->Draw("samese");  l->Draw();
    */
    c[nc++] = new TCanvas(outF_name+"-Jets_N_bSel",outF_name+"-Jets_N_bSel" );  gPad->SetLogy();
    THStack s_Jets_N_bSel("s_Jets_N_bSel","s_Jets_N_bSel");
    for(int i=1; i<nf; i++){    s_Jets_N_bSel.Add(Jets_N_bSel[i]); }
    Jets_N_bSel[0].GetXaxis()->SetTitle("NJets"); Jets_N_bSel[0]->GetYaxis()->SetTitle("events");
    Jets_N_bSel[0]->Draw("ep");  s_Jets_N_bSel.Draw("hsames");  Jets_N_bSel[0]->Draw("samese");  l->Draw();
    
    c[nc++] = new TCanvas(outF_name+"-MET_bSel",outF_name+"-MET_bSel" );  gPad->SetLogy();
    THStack s_MET_bSel("s_MET_bSel","s_MET_bSel");
    for(int i=1; i<nf; i++){    s_MET_bSel.Add(MET_bSel[i]); }
    MET_bSel[0].GetXaxis()->SetTitle("MET"); MET_bSel[0]->GetYaxis()->SetTitle("events");
    MET_bSel[0]->Draw("ep");  s_MET_bSel.Draw("hsames");  MET_bSel[0]->Draw("samese");  l->Draw();
    
    c[nc++] = new TCanvas(outF_name+"-HT_bSel",outF_name+"-HT_bSel" );  gPad->SetLogy();
    THStack s_HT_bSel("s_HT_bSel","s_HT_bSel");
    for(int i=1; i<nf; i++){    s_HT_bSel.Add(HT_bSel[i]); }
    HT_bSel[0].GetXaxis()->SetTitle("HT"); HT_bSel[0]->GetYaxis()->SetTitle("events");
    HT_bSel[0]->Draw("ep");  s_HT_bSel.Draw("hsames");  HT_bSel[0]->Draw("samese");  l->Draw();
    
    c[nc++] = new TCanvas(outF_name+"-MinDPhi_bSel",outF_name+"-MinDPhi_bSel" );  gPad->SetLogy();
    THStack s_MinDPhi_bSel("s_MinDPhi_bSel","s_MinDPhi_bSel");
    for(int i=1; i<nf; i++){    s_MinDPhi_bSel.Add(MinDPhi_bSel[i]); }
    MinDPhi_bSel[0].GetXaxis()->SetTitle("MinDPhi"); MinDPhi_bSel[0]->GetYaxis()->SetTitle("events");
    MinDPhi_bSel[0]->Draw("ep");  s_MinDPhi_bSel.Draw("hsames");  MinDPhi_bSel[0]->Draw("samese");  l->Draw();
  }
  //leptSel
  else{
    c[nc++] = new TCanvas(outF_name+"-MT2_leptSel",outF_name+"-MT2_leptSel" );   gPad->SetLogy();
    THStack s_MT2_leptSel("s_MT2_leptSel","s_MT2_leptSel");
    for(int i=1; i<nf; i++){    s_MT2_leptSel.Add(MT2_leptSel[i]);  }
    MT2_leptSel[0].GetXaxis()->SetTitle("M_{T2} (GeV)"); MT2_leptSel[0]->GetYaxis()->SetTitle("events");
    MT2_leptSel[0]->Draw("ep");   s_MT2_leptSel.Draw("hsames");  MT2_leptSel[0]->Draw("samese");  l->Draw();
    
    c[nc++] = new TCanvas(outF_name+"-Lept_pt_leptSel",outF_name+"-Lept_pt_leptSel" );  gPad->SetLogy();
    THStack s_Lept_pt_leptSel("s_Lept_pt_leptSel","s_Lept_pt_leptSel");
    for(int i=1; i<nf; i++){    s_Lept_pt_leptSel.Add(Lept_pt_leptSel[i]);  }
    Lept_pt_leptSel[0].GetXaxis()->SetTitle("Lepton p_{T} (GeV)"); Lept_pt_leptSel[0]->GetYaxis()->SetTitle("events");
    Lept_pt_leptSel[0]->Draw("ep");  s_Lept_pt_leptSel.Draw("hsames");  Lept_pt_leptSel[0]->Draw("samese");  l->Draw();
    
    c[nc++] = new TCanvas(outF_name+"-BJets_N_leptSel",outF_name+"-BJets_N_leptSel" );  gPad->SetLogy();
    THStack s_BJets_N_leptSel("s_BJets_N_leptSel","s_BJets_N_leptSel");
    for(int i=1; i<nf; i++){    s_BJets_N_leptSel.Add(BJets_N_leptSel[i]);  }
    BJets_N_leptSel[0].GetXaxis()->SetTitle("NBJets"); BJets_N_leptSel[0]->GetYaxis()->SetTitle("events");
    BJets_N_leptSel[0]->Draw("ep");  s_BJets_N_leptSel.Draw("hsames");  BJets_N_leptSel[0]->Draw("samese");  l->Draw();
    
    c[nc++] = new TCanvas(outF_name+"-Jets_N_leptSel",outF_name+"-Jets_N_leptSel" );  gPad->SetLogy();
    THStack s_Jets_N_leptSel("s_Jets_N_leptSel","s_Jets_N_leptSel");
    for(int i=1; i<nf; i++){    s_Jets_N_leptSel.Add(Jets_N_leptSel[i]);  }
    Jets_N_leptSel[0].GetXaxis()->SetTitle("NJets"); Jets_N_leptSel[0]->GetYaxis()->SetTitle("events");
    Jets_N_leptSel[0]->Draw("ep");  s_Jets_N_leptSel.Draw("hsames");  Jets_N_leptSel[0]->Draw("samese");  l->Draw();
    
    c[nc++] = new TCanvas(outF_name+"-MET_leptSel",outF_name+"-MET_leptSel" );  gPad->SetLogy();
    THStack s_MET_leptSel("s_MET_leptSel","s_MET_leptSel");
    for(int i=1; i<nf; i++){    s_MET_leptSel.Add(MET_leptSel[i]);  }
    MET_leptSel[0].GetXaxis()->SetTitle("MET"); MET_leptSel[0]->GetYaxis()->SetTitle("events");
    MET_leptSel[0]->Draw("ep");  s_MET_leptSel.Draw("hsames");  MET_leptSel[0]->Draw("samese");  l->Draw();
    
    c[nc++] = new TCanvas(outF_name+"-HT_leptSel",outF_name+"-HT_leptSel" );  gPad->SetLogy();
    THStack s_HT_leptSel("s_HT_leptSel","s_HT_leptSel");
    for(int i=1; i<nf; i++){    s_HT_leptSel.Add(HT_leptSel[i]);  }
    HT_leptSel[0].GetXaxis()->SetTitle("HT"); HT_leptSel[0]->GetYaxis()->SetTitle("events");
    HT_leptSel[0]->Draw("ep");  s_HT_leptSel.Draw("hsames");  HT_leptSel[0]->Draw("samese");  l->Draw();
    
    c[nc++] = new TCanvas(outF_name+"-MinDPhi_leptSel",outF_name+"-MinDPhi_leptSel" );  gPad->SetLogy();
    THStack s_MinDPhi_leptSel("s_MinDPhi_leptSel","s_MinDPhi_leptSel");
    for(int i=1; i<nf; i++){    s_MinDPhi_leptSel.Add(MinDPhi_leptSel[i]);  }
    MinDPhi_leptSel[0].GetXaxis()->SetTitle("MinDPhi"); MinDPhi_leptSel[0]->GetYaxis()->SetTitle("events");
    MinDPhi_leptSel[0]->Draw("ep");  s_MinDPhi_leptSel.Draw("hsames");  MinDPhi_leptSel[0]->Draw("samese");  l->Draw();
  }
  


  TString outPlotDir = DIR+"/plots/"+outF_name;
  system( ( "mkdir -p "+outPlotDir).Data() );
  TFile *fout = TFile::Open(outPlotDir+"/plots_"+outF_name+".root","RECREATE");
  TString htmlBody = "";
  for(int i=0; i<nc; i++){
    c[i]->Write();
    TString iName = outPlotDir+"/"+c[i]->GetName()+".png";
    c[i]->SaveAs(iName);
    TString image = c[i]->GetName();
    htmlBody += "<img src=\""+image+".png\" />";
  }

  TString htmlFile = outPlotDir+"/index.html";
  FILE * pFile = fopen ( htmlFile.Data() ,"w");
  TString  html = "<html>\n<body>\n<h1>File:"+outF_name+"</h1><hr>"+htmlBody+"\n</body>\n</html>\n";
  fprintf(pFile, "%s", html.Data() );
  fclose(pFile);

  fout->Write();
  fout->Close();
  
}
