void SetStyle_PRD(){

  TStyle *RootStyle = new TStyle("Root-Style","Single Top Style for PRD");
  
#ifdef __CINT__
  TStyle *GloStyle;
  GloStyle = gStyle;                          // save the global style reference
  gStyle = RootStyle;
#endif


  // Canvas

  RootStyle->SetCanvasColor     (0);
  RootStyle->SetCanvasBorderMode(0);

  // Pads

  RootStyle->SetPadColor       (0);
  RootStyle->SetPadBorderMode  (0);
  RootStyle->SetPadBottomMargin(0.17);
  RootStyle->SetPadTopMargin   (0.07);
  RootStyle->SetPadLeftMargin  (0.18);
  RootStyle->SetPadRightMargin (0.08);
  RootStyle->SetPadTickX       (0);
  RootStyle->SetPadTickY       (0);

  // Frames

  RootStyle->SetFrameLineWidth (3);
  RootStyle->SetFrameFillColor (0);
  RootStyle->SetFrameBorderMode(0);
  RootStyle->SetFrameBorderSize(0);


  // Histograms

   RootStyle->SetHistFillColor(0);
   RootStyle->SetHistLineWidth(2);

  // Functions

  //RootStyle->SetFuncColor(1);
  //RootStyle->SetFuncStyle(0);
  //RootStyle->SetFuncWidth(2);

  //Legends 

  //RootStyle->SetLegendBorderSize(0);
  //RootStyle->SetFillStyle(0);
  //RootStyle->SetTextFont(62); 
  //RootStyle->SetTextSize(0.045);

  // Labels, Ticks, and Titles

  //RootStyle->SetTickLength ( 0.015,"X");
  RootStyle->SetTitleSize  ( 0.070,"X");
  RootStyle->SetTitleOffset( 1.100,"X");
  RootStyle->SetLabelSize  ( 0.070,"X");
  RootStyle->SetNdivisions ( 505  ,"X");

//RootStyle->SetTickLength ( 0.015,"Y");
  RootStyle->SetTitleSize  ( 0.070,"Y");
  RootStyle->SetTitleOffset( 1.300,"Y");
  RootStyle->SetLabelSize  ( 0.070,"Y");
  RootStyle->SetNdivisions ( 505  ,"Y");

//   RootStyle->SetTickLength ( 0.015,"Z");
//   RootStyle->SetTitleSize  ( 0.060,"Z");
//   RootStyle->SetTitleOffset( 1.100,"Z");
//   RootStyle->SetLabelOffset( 0.015,"Z");
//   RootStyle->SetLabelSize  ( 0.050,"Z");
//   RootStyle->SetLabelFont  ( 42   ,"Z");
//   RootStyle->SetTitleFont  ( 42   ,"Z");
//   RootStyle->SetNdivisions ( 707   ,"Z");


  RootStyle->SetTitleBorderSize  (0);
  RootStyle->SetTitleFillColor  (0);  
//RootStyle->SetTitleFont  (42);
//RootStyle->SetTitleColor  (1);

  RootStyle->SetLineWidth  (2);

  // Options

  RootStyle->SetOptFit        (0111);
  RootStyle->SetOptStat       (1);
  RootStyle->SetStatBorderSize(0);
  RootStyle->SetStatColor     (0);
  RootStyle->SetStatW(0.5);
//  RootStyle->SetMarkerStyle(20);
//RootStyle->SetMarkerSize(1.25);

//RootStyle->SetPalette(1);

}
