void rootlogon() {
  gSystem->SetIncludePath(" -I../MT2Code/include/ -I../TESCO/include/");
  gSystem->Load("libPhysics");
  gSystem->Load("../MT2Code/shlib/libDiLeptonAnalysis.so");
  gROOT->ProcessLine(".x SetStyle_PRD.C");
}
