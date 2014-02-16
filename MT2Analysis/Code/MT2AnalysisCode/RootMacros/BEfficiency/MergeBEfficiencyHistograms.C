#include "TString.h"
#include <iostream>



using namespace std;

void MergeBEfficiencyHistograms(){


	TString directory = "/shome/haweber/MT2Analysis_8TeV/Code/Efficiencies/InputHistos/";
	TString oldfile1  = "BEfficiencies_histos_ge2jets_0l_CSVM_HT750orMET200_MinDPhiPt40_onlyQCD.root";
	TString oldfile2  = "BEfficiencies_histos_ge2jets_0l_CSVM_HT750orMET200_MinDPhiPt40_noQCD.root";
	TString newfile   = "BEfficiencies_histos_ge2jets_0l_CSVM_HT750orMET200_MinDPhiPt40.root";

	cout << "type following three lines as command:" << endl;
	cout << "cd " << directory << endl;
	cout << "hadd " << newfile << " " << oldfile1 << " " << oldfile2 << endl;
	cout << "cd -" << endl;

}