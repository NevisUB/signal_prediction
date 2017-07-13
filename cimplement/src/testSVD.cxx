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
#include "TMatrixT.h"
#include "TVectorT.h"

#include "miniSVD.h"

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

	std::cout<<"Starting "<<std::endl;	
	TFile *f2 = TFile::Open("out.root");
	TH1D * truth_MC = (TH1D*)f2->Get("Truth_Enu_Reco_Eqe_t;1");
	TH1D * Hbini = (TH1D*)f2->Get("Truth_Enu_Reco_Eqe_r;1");
	TH2D * response_MC = (TH2D*)f2->Get("Truth_Enu_Reco_Eqe_pass;1");



	std::cout<<"Starting (obs) "<<std::endl;	
	TFile *f = TFile::Open("fout.root");
	TMatrixT<double> * A = (TMatrixT<double>*)f->Get("adet;1");	
	TVectorT<double>* b = (TVectorT<double>*)f->Get("bini;1");	
	TVectorT<double>* xini = (TVectorT<double>*)f->Get("xini;1");	
	std::cout<<"Done (obs) "<<std::endl;	
	f->Close();
	std::cout<<"Closing is the problem?"<<std::endl;



	TH1D * LEE = (TH1D*)Hbini->Clone("LEE");
	for(int i=0; i<=LEE->GetNbinsX(); i++){
		LEE->SetBinContent(i,0.0);
		if(i<3 && i != 0){
			LEE->SetBinContent(i,33);
		}
	}

		miniSVD * mylee =  new miniSVD(response_MC, LEE, truth_MC,10);
		TH1D * hmm = mylee->unfold();
		hmm->Print();
		hmm->SetLineColor(kBlack);
		hmm->SetLineWidth(3);


	std::vector<TH1D*> vec_reco;
	std::vector<TH1D*> v_ans;


	TRandom3  * rang = new TRandom3(201);

	int NN = 6;
	for(int k=0; k<NN; k++){	
		std::cout<<"Generating possion #: "<<k<<std::endl;	
		TH1D * bdat =(TH1D *)Hbini->Clone("bdat");

		for(int i=0; i<=bdat->GetNbinsX(); i++){
			bdat->SetBinContent(i, rang->Poisson(bdat->GetBinContent(i)));
		}
		vec_reco.push_back(bdat);
	}

	TCanvas *c2 = new TCanvas();
	c2->cd();
	for(int k=0;k<NN;k++){
		
			std::cout<<"On kth universe "<<k<<std::endl;	
			miniSVD * mysvd =  new miniSVD(response_MC, vec_reco.at(k) ,truth_MC, 2);
			mysvd->unfold();	

		
	//	TSVDUnfold *ts = new TSVDUnfold(vec_reco.at(k), bini, xini, Adet);
		//TH1D *ans =ts->Unfold(3);
		TH1D * ans = mysvd->unfold();
		//ans->Print();
	//	TH1D * d = ts->GetD();
	//	d->Draw("same");
		v_ans.push_back(ans);
		for(int i=0; i<=ans->GetNbinsX()+1; i++){
			//std::cout<<ans->GetBinContent(i)<<std::endl;
		}
	}

	std::vector<int> tmp_colors ={kRed,kBlue,kGreen,kOrange,kCyan,kViolet};
	TCanvas * c = new TCanvas();
	c->Divide(2,2);
	c->cd(1);
	truth_MC->SetMarkerStyle(21);
	truth_MC->SetLineColor(kBlack);
	truth_MC->Draw();
	for(int i= 0; i< NN; i++){
			v_ans[i]->SetLineColor(tmp_colors[i]);
			v_ans[i]->SetLineWidth(2);

		v_ans[i]->Draw("hist same");
	}
	truth_MC->Draw("same ap");


	c->cd(2);
	Hbini->SetMarkerStyle(21);
	Hbini->SetLineColor(kBlack);
	Hbini->Draw();
	for(int i= 0; i< NN; i++){
			vec_reco[i]->SetLineColor(tmp_colors[i]);
			vec_reco[i]->SetLineWidth(3);

		vec_reco[i]->Draw("hist same");
	}
	Hbini->Draw("same ap");

	c->cd(3);
		std::cout<<"LEE unfolded: "<<std::endl;
		for(int i=0; i<=hmm->GetNbinsX(); i++){
			std::cout<<i<<" "<<hmm->GetBinContent(i)<<std::endl;
		}

		hmm->Draw("hist");	
	c->cd(4);
		LEE->Draw("hist");


	c->SaveAs("foo.pdf","pdf");	
	f2->Close();
	
	return 0;
}


