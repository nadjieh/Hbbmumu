void rootlogon() {
  gSystem->SetIncludePath(" -I./include/ -I../TESCO/include/ -I${VO_CMS_SW_DIR}/slc5_amd64_gcc462/cms/cmssw/CMSSW_5_3_7/src/");
  //  gSystem->Load("/public/V_CMSSW64/slc5_amd64_gcc462/cms/cmssw/CMSSW_5_3_7/lib/slc5_amd64_gcc462/libFWCoreFWLite.so");
  gSystem->Load("$(VO_CMS_SW_DIR)/slc5_amd64_gcc462/cms/cmssw/CMSSW_5_3_7/lib/$(SCRAM_ARCH)/libFWCoreFWLite.so");
  gSystem->Load("$(VO_CMS_SW_DIR)/slc5_amd64_gcc462/cms/cmssw/CMSSW_5_3_7/lib/$(SCRAM_ARCH)/libFWCoreUtilities.so");
  gSystem->Load("$(VO_CMS_SW_DIR)/slc5_amd64_gcc462/cms/cmssw/CMSSW_5_3_7/lib/$(SCRAM_ARCH)/libDataFormatsCommon.so");
  gSystem->Load("$(VO_CMS_SW_DIR)/slc5_amd64_gcc462/cms/cmssw/CMSSW_5_3_7/lib/$(SCRAM_ARCH)/libDataFormatsFWLite.so");

  gSystem->Load("libPhysics");
  gSystem->Load("shlib/libDiLeptonAnalysis.so");
  gROOT->ProcessLine(".x SetStyle_PRD.C");
}
