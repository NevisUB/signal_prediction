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
#include "miniBBB.h"

#define no_argument 0
#define required_argument 1
#define optional_argument 2


int main(int argc, char* argv[])
{

	gStyle->SetOptStat(0);
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
	TH1D * truth_MC_pass = (TH1D*)f2->Get("Enu_CCQE_nue_pass;1");


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
		if(i<3 && i !=0 ){
			LEE->SetBinContent(i,33);
		}
	}

	std::vector<TH1D*> vec_reco;
	std::vector<TH1D*> v_ans_svd;
	std::vector<TH1D*> v_ans_bbb;
	std::vector<TH1D*> v_ans_inv;

	int USEK = 2;

	TRandom3  * rang = new TRandom3(0);
	std::cout<<"TRandom3 Seed: "<<rang->GetSeed()<<std::endl;
	int NN = 5;
	for(int k=0; k<NN; k++){	
		std::cout<<"Generating possion #: "<<k<<std::endl;	
		TH1D * bdat =(TH1D *)Hbini->Clone("bdat");

		for(int i=0; i<=bdat->GetNbinsX(); i++){
			bdat->SetBinContent(i, rang->Poisson(bdat->GetBinContent(i)));
		}
		vec_reco.push_back(bdat);
	}

	TCanvas *c2 = new TCanvas("a","a",1200,600);
	c2->Divide(2,1);
	c2->cd(1);
	TH1D * tt = (TH1D*)truth_MC->Clone("tt");
	TH1D * rr = (TH1D*)Hbini->Clone("rr");
	tt->SetMarkerStyle(21);
	rr->SetMarkerStyle(21);
	tt->SetMarkerColor(kBlue-6);
	rr->SetMarkerColor(kRed-7);
	tt->SetLineWidth(2);
	rr->SetLineWidth(2);
	tt->SetLineColor(kBlue-6);
	rr->SetLineColor(kRed-7);
	tt->GetXaxis()->SetTitle("Energy [GeV]");
	tt->GetYaxis()->SetTitle("Events");
	tt->GetXaxis()->SetTitleSize(0.04);
	tt->GetXaxis()->SetLabelSize(0.04);
	tt->GetYaxis()->SetLabelSize(0.04);
	tt->GetYaxis()->SetTitleOffset(1.45);
	tt->SetMinimum(0);
	tt->SetTitle("Intrinsic #nu_{e} CCQE only");
	tt->SetMaximum(170);

	TLegend * leg = new TLegend(0.58,0.6,0.89,0.89);
	leg->AddEntry(tt,"True E_{#nu} MC","lep");
	leg->AddEntry(rr,"Reco E_{QE} MC","lep");
	leg->SetFillStyle(0);
//	leg->SetBorderSize(0.1);

	tt->Draw();
	rr->Draw("same ap");
	leg->Draw();

	c2->cd(2)->SetRightMargin(0.175);
	TH2D * prob_MC = (TH2D*)response_MC->Clone("prob");
	for(int i=0; i<=truth_MC->GetNbinsX();i++){
	for(int j=0; j<=truth_MC->GetNbinsX();j++){
		prob_MC->SetBinContent(i,j ,response_MC->GetBinContent(i,j)/truth_MC->GetBinContent(j));
	}
	}
	prob_MC->GetYaxis()->SetTitle("Truth E_{#nu} [GeV]");
	prob_MC->GetXaxis()->SetTitle("Reco E_{QE} [GeV]");
	prob_MC->SetTitle("Response Matrix");
	prob_MC->GetXaxis()->SetTitleOffset(1.1);
	prob_MC->GetYaxis()->SetTitleOffset(1.1);

	prob_MC->Draw("colz");
	

	c2->SaveAs("mod1.pdf","pdf");

	for(int k=0;k<NN;k++){

		std::cout<<"On kth universe "<<k<<std::endl;	
	
		miniBBB * myinv =  new miniBBB(response_MC, vec_reco.at(k), truth_MC, truth_MC_pass, Hbini);
		miniBBB * mybbb =  new miniBBB(response_MC, vec_reco.at(k), truth_MC, truth_MC_pass, Hbini);
		miniSVD * mysvd =  new miniSVD(response_MC, vec_reco.at(k) ,truth_MC, 3);


		TH1D * ans = mybbb->unfold(USEK);
		TH1D * ans2 = mysvd->unfold();
		TH1D *ans3 = myinv->unfold(2000);

		v_ans_bbb.push_back(ans);
		v_ans_svd.push_back(ans2);
		v_ans_inv.push_back(ans3);
	}

	
	std::vector<int> tmp_colors ={kOrange-3, kMagenta-3, kGreen-6,kBlue-7,kRed-7};
	TCanvas * c = new TCanvas("","",1200,800);

	c->Divide(2,2);

	c->cd(1);
	Hbini->SetMarkerStyle(21);
	Hbini->SetLineColor(kBlack);
	Hbini->Draw();
	Hbini->SetMaximum(100);
	Hbini->SetTitle("Reconstructed E_{QE}");
	Hbini->GetXaxis()->SetTitle("E_{QE} [GeV]");
	Hbini->GetYaxis()->SetTitle("Events");
	Hbini->GetXaxis()->SetTitleSize(0.04);
	Hbini->GetXaxis()->SetLabelSize(0.04);
	Hbini->GetYaxis()->SetLabelSize(0.04);

	TLegend * leg_reco = new TLegend(0.6,0.7,0.89,0.89);
	leg_reco->SetBorderSize(0);
	leg_reco->AddEntry(Hbini,"Reco E_{QE} MC","lep");
	leg_reco->Draw();


	for(int i= 0; i< NN; i++){
		vec_reco[i]->SetLineColor(tmp_colors[i]);
		vec_reco[i]->SetLineWidth(1.5);
		vec_reco[i]->SetTitle("Reconstucted E_{QE}");
		vec_reco[i]->Draw("hist same");
	}
	Hbini->Draw("same ap");

	c->cd(2)->SetBottomMargin(20);
	v_ans_inv[0]->Draw("hist");
	v_ans_inv[0]->SetMaximum(200);
	v_ans_inv[0]->GetXaxis()->SetTitle("E_{#nu} [GeV]");
	v_ans_inv[0]->GetYaxis()->SetTitle("Events");
	v_ans_inv[0]->SetLineColor(tmp_colors[0]);
	v_ans_inv[0]->SetLineWidth(1.5);
	v_ans_inv[0]->SetTitle("Unfolded Truth: Matrix Inversion");
	v_ans_inv[0]->GetXaxis()->SetTitleSize(0.04);
	v_ans_inv[0]->GetXaxis()->SetLabelSize(0.04);
	v_ans_inv[0]->GetYaxis()->SetLabelSize(0.04);


	truth_MC->SetMarkerStyle(21);
	truth_MC->SetLineColor(kBlack);
	for(int i= 1; i< NN; i++){
		v_ans_inv[i]->SetLineColor(tmp_colors[i]);
		v_ans_inv[i]->SetLineWidth(1.5);
		v_ans_inv[i]->Draw("hist same");
	}
	TLegend * leg_inv = new TLegend(0.6,0.7,0.89,0.89);
	leg_inv->AddEntry(truth_MC,"True E_{#nu} MC","lep");
	leg_inv->SetFillStyle(0);
	leg_inv->SetBorderSize(0);
	leg_inv->Draw();

	truth_MC->SetMarkerStyle(21);
	truth_MC->SetLineColor(kBlack);
	truth_MC->Draw("same ap");
	
	c->cd(3);
	v_ans_svd[0]->Draw("hist");
	v_ans_svd[0]->SetMaximum(200);
	v_ans_svd[0]->GetXaxis()->SetTitle("E_{#nu} [GeV]");
	v_ans_svd[0]->GetYaxis()->SetTitle("Events");
	v_ans_svd[0]->SetLineColor(tmp_colors[0]);
	v_ans_svd[0]->SetLineWidth(1.5);
	v_ans_svd[0]->GetXaxis()->SetTitleSize(0.04);
	v_ans_svd[0]->GetXaxis()->SetLabelSize(0.04);
	v_ans_svd[0]->GetYaxis()->SetLabelSize(0.04);


	v_ans_svd[0]->SetTitle("Unfolded Truth: SVD Unfolding, k=3");
	truth_MC->SetMarkerStyle(21);
	truth_MC->SetLineColor(kBlack);
	for(int i= 1; i< NN; i++){
		v_ans_svd[i]->SetLineColor(tmp_colors[i]);
		v_ans_svd[i]->SetLineWidth(1.5);
		v_ans_svd[i]->Draw("hist same");
	}

	TLegend * leg_svd = new TLegend(0.6,0.7,0.89,0.89);
	leg_svd->AddEntry(truth_MC,"True E_{#nu} MC","lep");
	leg_svd->SetBorderSize(0);
	leg_svd->Draw();
	truth_MC->SetMarkerStyle(21);
	truth_MC->SetLineColor(kBlack);
	truth_MC->Draw("same ap");


	c->cd(4);
	v_ans_bbb[0]->Draw("hist");
	v_ans_bbb[0]->SetMaximum(200);
	v_ans_bbb[0]->GetXaxis()->SetTitle("E_{#nu} [GeV]");
	v_ans_bbb[0]->GetYaxis()->SetTitle("Events");
	v_ans_bbb[0]->SetLineColor(tmp_colors[0]);
	v_ans_bbb[0]->SetLineWidth(1.5);
	v_ans_bbb[0]->SetTitle("Unfolded Truth: D'Agostini Iteration, k=3");
	v_ans_bbb[0]->GetXaxis()->SetTitleSize(0.04);
	v_ans_bbb[0]->GetXaxis()->SetLabelSize(0.04);
	v_ans_bbb[0]->GetYaxis()->SetLabelSize(0.04);


	truth_MC->SetMarkerStyle(21);
	truth_MC->SetLineColor(kBlack);
	for(int i= 1; i< NN; i++){
		v_ans_bbb[i]->SetLineColor(tmp_colors[i]);
		v_ans_bbb[i]->SetLineWidth(1.5);
		v_ans_bbb[i]->Draw("hist same");
	}
	TLegend * leg_bbb = new TLegend(0.6,0.7,0.89,0.89);
	leg_bbb->AddEntry(truth_MC,"True E_{#nu} MC","lep");
	leg_bbb->SetBorderSize(0);
	leg_bbb->Draw();

	truth_MC->SetMarkerStyle(21);
	truth_MC->SetLineColor(kBlack);
	truth_MC->Draw("same ap");



	c->SaveAs("foo.pdf","pdf");	
	f2->Close();

	return 0;
}


