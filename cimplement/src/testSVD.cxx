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
	TH1D * Hxini = (TH1D*)f2->Get("Truth_Enu_Reco_Eqe_t;1");
	TH1D * Hbini = (TH1D*)f2->Get("Truth_Enu_Reco_Eqe_r;1");
	TH2D * HAdet = (TH2D*)f2->Get("Truth_Enu_Reco_Eqe_pass;1");



	std::cout<<"Starting (obs) "<<std::endl;	
	TFile *f = TFile::Open("fout.root");
	TMatrixT<double> * A = (TMatrixT<double>*)f->Get("adet;1");	
	TVectorT<double>* b = (TVectorT<double>*)f->Get("bini;1");	
	TVectorT<double>* xini = (TVectorT<double>*)f->Get("xini;1");	
	std::cout<<"Done (obs) "<<std::endl;	
	f->Close();
	std::cout<<"Closing is the problem?"<<std::endl;

	
	
	
	std::vector<TH1D*> v_bdat;
	std::vector<TH1D*> v_ans;


	TRandom3  * rang = new TRandom3();

	int NN = 6;
	for(int k=0; k<NN; k++){	
		std::cout<<"Generating possion #: "<<k<<std::endl;	
		TH1D * bdat =(TH1D *)Hbini->Clone("bdat");

		for(int i=0; i<=bdat->GetNbinsX(); i++){
			bdat->SetBinContent(i, rang->Poisson(bdat->GetBinContent(i)));
		}
		v_bdat.push_back(bdat);
	}

	TCanvas *c2 = new TCanvas();
	c2->cd();
	for(int k=0;k<NN;k++){
		
			std::cout<<"On kth universe "<<k<<std::endl;	
			miniSVD * mysvd =  new miniSVD(HAdet,v_bdat.at(k),Hxini);
			mysvd->init(0.0001);
			mysvd->unfold();	

		
	//	TSVDUnfold *ts = new TSVDUnfold(v_bdat.at(k), bini, xini, Adet);
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
	c->Divide(2,1);
	c->cd(1);
	Hxini->SetMarkerStyle(21);
	Hxini->SetLineColor(kBlack);
	Hxini->Draw("ap");
	for(int i= 0; i< NN; i++){
			v_ans[i]->SetLineColor(tmp_colors[i]);
			v_ans[i]->SetLineWidth(2);

		v_ans[i]->Draw("hist same");
	}
	Hxini->Draw("same ap");


	c->cd(2);
	Hbini->SetMarkerStyle(21);
	Hbini->SetLineColor(kBlack);
	Hbini->Draw("ap");
	for(int i= 0; i< NN; i++){
			v_bdat[i]->SetLineColor(tmp_colors[i]);
			v_bdat[i]->SetLineWidth(3);

		v_bdat[i]->Draw("hist same");
	}
	Hbini->Draw("same ap");

	c->SaveAs("foo.pdf","pdf");	
	f2->Close();
	
	return 0;
}


