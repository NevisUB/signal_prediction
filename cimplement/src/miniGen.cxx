#include <iostream>
#include <cstring>
#include <vector>
#include <iterator>
#include <algorithm>
#include <getopt.h>

#include "TClass.h"
#include "TTree.h"
#include "TFile.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "THStack.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TGaxis.h"
#include "TChain.h"
#include "TROOT.h"
#include "TColor.h"

#include "boone_fnc.h"
#include "miniClass.h"

#define no_argument 0
#define required_argument 1
#define optional_argument 2


int main(int argc, char* argv[])
{


	std::string rootlocation = "/rootfiles";
	int iarg = 0;
	opterr=1;
	int index; 
	/*************************************************************
	 *************************************************************
	 *		Command Line Argument Reading
	 ************************************************************
	 ************************************************************/

	const struct option longopts[] = 
	{
		{"file", 		required_argument, 	0, 'f'},
		{"help",		no_argument,	0, 'h'},
		{0,			no_argument, 		0,  0},
	};


	while(iarg != -1)
	{
		iarg = getopt_long(argc,argv, "f:h", longopts, &index);

		switch(iarg)
		{
			case 'f':
				rootlocation = optarg;
				break;
			case '?':
			case 'h':
				std::cout<<"Allowed arguments:"<<std::endl;
				std::cout<<"\t-f\t--file\t\tInput directy location of rootfiles"<<std::endl;
				return 0;
		}

	}

	std::string prefix = "output";
	std::string postfix = ".root";	

	gROOT->ProcessLine("#include <vector>");

	TH1::SetDefaultSumw2(true);
	TGaxis::SetMaxDigits(3);

	float fEnergy = 0;
	float fCosTheta = 0;
	float fWeight = 1;

	bool fPassOsc = 0;

	int fEventNumber = 0;
	int fNuType = 0;
	int fNDecay = 0;
	int fNFSP = 0;
	int fNHit = 0;
	int fNUANCEChan = 0;
	int fNuParentID = 0;

	float fNuMomX = 0;
	float fNuMomY = 0;
	float fNuMomZ = 0;
	float fNuMomT = 0;

	std::vector<int> * fNCER = 0;
	std::vector<int> * fNSCI = 0;
	std::vector<int> * fFSPType = 0;
	std::vector<float> * fTime = 0;
	std::vector<float> * fVertexX = 0;
	std::vector<float> * fVertexY = 0;
	std::vector<float> * fVertexZ = 0;
	std::vector<float> * fMomX = 0;
	std::vector<float> * fMomY = 0;
	std::vector<float> * fMomZ = 0;
	std::vector<float> * fMomT = 0;
	std::vector<float> * fTimePMT = 0;

	float fNuanceNuEneIn = 0;
	float fNuanceNuEneOut = 0;
	float fNuanceTargFermiMom = 0;
	float fNuanceFourMomTransfer = 0;
	float fNuanceEnergyTransfer = 0;
	float fNuanceX = 0;
	float fNuanceY = 0;
	float fNuancePhiCM = 0;
	float fNuanceFermiMomDir = 0;

	TChain * chain = new TChain("MiniBooNE_CCQE","chain");
	std::string fileName =rootlocation+"/"+prefix+"*"+postfix;
	chain->Add(fileName.c_str());

	chain->SetMakeClass(1);

	chain->SetBranchAddress("Energy",&fEnergy);
	chain->SetBranchAddress("CosTheta",&fCosTheta);
	chain->SetBranchAddress("Weight",&fWeight);
	chain->SetBranchAddress("PassOsc",&fPassOsc);
	chain->SetBranchAddress("EventNumber",&fEventNumber);
	chain->SetBranchAddress("NuType",&fNuType);
	chain->SetBranchAddress("NDecay",&fNDecay);
	chain->SetBranchAddress("NFSP",&fNFSP);
	chain->SetBranchAddress("NHit",&fNHit);
	chain->SetBranchAddress("NuMomX",&fNuMomX);
	chain->SetBranchAddress("NuMomY",&fNuMomY);
	chain->SetBranchAddress("NuMomZ",&fNuMomZ);
	chain->SetBranchAddress("NuMomT",&fNuMomT);
	chain->SetBranchAddress("NCER",&fNCER);
	chain->SetBranchAddress("NSCI",&fNSCI);
	chain->SetBranchAddress("FSPType",&fFSPType);
	chain->SetBranchAddress("VertexX",&fVertexX);
	chain->SetBranchAddress("VertexY",&fVertexY);
	chain->SetBranchAddress("VertexZ",&fVertexZ);
	chain->SetBranchAddress("Time",&fTime);
	chain->SetBranchAddress("MomX",&fMomX);
	chain->SetBranchAddress("MomY",&fMomY);
	chain->SetBranchAddress("MomZ",&fMomZ);
	chain->SetBranchAddress("MomT",&fMomT);
	chain->SetBranchAddress("NuParentID",&fNuParentID);
	chain->SetBranchAddress("TimePMT",&fTimePMT);
	chain->SetBranchAddress("NUANCEChan",&fNUANCEChan);
	chain->SetBranchAddress("NuanceNuEneIn",&fNuanceNuEneIn);
	chain->SetBranchAddress("NuanceNuEneOut",&fNuanceNuEneOut);
	chain->SetBranchAddress("NuanceTargFermiMom",&fNuanceTargFermiMom);
	chain->SetBranchAddress("NuanceFourMomTransfer",&fNuanceFourMomTransfer);
	chain->SetBranchAddress("NuanceEnergyTransfer",&fNuanceEnergyTransfer);
	chain->SetBranchAddress("NuanceX",&fNuanceX);
	chain->SetBranchAddress("NuanceY",&fNuanceY);
	chain->SetBranchAddress("NuancePhiCM",&fNuancePhiCM);
	chain->SetBranchAddress("NuanceFermiMomDir",&fNuanceFermiMomDir);

	/**************************************************************
	 *		Enough nonsense
	 *
	 * ************************************************************/




	//Which of the selections are you defining as "Signal"
	std::string signal = "CCQE_nue";	
	//The rest..
	std::vector<std::string> v_selection_names = {signal,"CCQE_numu","CCPIP","NCPI0","CHPI0","DELTA","NCQE_PI","Other"};
	//What variables do you want to plot
	std::vector<std::string> v_variable_names = {"Enu","Evis","CosTheta","Eqe","tHit","-Q^2"};
	//The upper and lower bounds for said variable histgrams
	std::vector<double> lower = {0.3,0.3,-1,0.3,0,0};
	std::vector<double> higher = {1.5,1.5,1,1.5,10000,2};
	//The colors (can have different pass and fail, but meh)
	std::vector<int> cols = {kRed-6, kBlue-7, kMagenta+3,kOrange+1,kMagenta-3,kCyan+1, kOrange+2  ,kGray};


	//Initilise the class with all that
	miniClass myclass(v_variable_names, v_selection_names, signal, lower, higher);


	//Some simple counting variables
	int Ne=0;
	double NeW = 0.0;
	int Nm=0;
	double NmW =0.0;
	int Nw0=0;

	double gw = 1.0;// 1.0/3.0;

	long entries = chain->GetEntries();
	for (long entry = 0; entry < entries; ++entry) {
		chain->GetEntry(entry);

		if (entry % 5000 == 0) {
			std::cout << "At entry: " << entry << " out of " << entries << std::endl;
		}
	
		//Dont lke negative or 0 weights..
		if(fWeight<=0){Nw0++;continue;};


		if(fWeight<=0 && fPassOsc){
			std::cerr<<"ERROR; pass osc and negative wight"<<std::endl;
			exit(EXIT_FAILURE);
		}

		if(fNuType ==NUE || fNuType ==NUEBAR){Ne++;NeW+=fWeight ;};
		if(fNuType ==NUMU || fNuType ==NUMUBAR){Nm++;NmW += fWeight;};

		//Whats the background Type? Needs to correspond to the above selection types.
		int IBKGD = calcIBKGD(fNuType, fNUANCEChan);
		std::string sIBKGD = stringCalcIBKGD( fNuType, fNUANCEChan);


		//Fill all the respective histograms
		double Mn = 1;
		myclass.Fill("Evis", sIBKGD, fPassOsc, fEnergy/1e3   ,fWeight*gw);
		myclass.Fill("Enu", sIBKGD, fPassOsc, fNuMomT  ,fWeight*gw);
		myclass.Fill("Eqe", sIBKGD, fPassOsc, Mn*fEnergy/1e3/(Mn-fEnergy/1e3 +fEnergy/1e3*fCosTheta )  ,fWeight*gw);
		myclass.Fill("CosTheta", sIBKGD, fPassOsc, fCosTheta ,fWeight*gw);
		myclass.Fill("-Q^2", sIBKGD, fPassOsc, -fNuanceFourMomTransfer/1e6  ,fWeight*gw);
		myclass.Fill("tHit", sIBKGD, fPassOsc, fNHit  ,fWeight*gw);
	//	myclass.Fill("weight", sIBKGD, fPassOsc, fWeight*gw  ,1.0);
		
		myclass.Fill2D("Evis",sIBKGD,fPassOsc, fNuMomT, fEnergy/1e3, fWeight*gw);
		myclass.Fill2D("CosTheta",sIBKGD,fPassOsc, fNuMomT, fCosTheta, fWeight*gw);
		myclass.Fill2D("Eqe",sIBKGD,fPassOsc, fNuMomT,  Mn*fEnergy/1e3/(Mn-fEnergy/1e3 +fEnergy/1e3*fCosTheta ), fWeight*gw); 
		myclass.Fill2D("-Q^2",sIBKGD,fPassOsc, fNuMomT, -fNuanceFourMomTransfer/1e6 ,fWeight*gw);
		myclass.Fill2D("tHit",sIBKGD,fPassOsc, fNuMomT, fNHit, fWeight*gw);
		

	}
	//Set the colors.
	myclass.setColors(cols,cols);
	myclass.writeOut("out.root");

	std::cout<<"COUNTING: "<<entries<<" Ne: "<<Ne<<" NeW: "<<NeW<<" Nm: "<<Nm<<" NmW: "<<NmW<<std::endl;
	std::cout<<"ZeroWeights: "<<Nw0<<" out of "<<entries<<std::endl;

	delete chain;
	return 0;
}


