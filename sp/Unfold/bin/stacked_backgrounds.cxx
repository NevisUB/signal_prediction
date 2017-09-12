#include "Unfold/Core/SPIO.h"

#include "Combined/CombinedTypes.h"
#include "Combined/CombinedUtil.h"

#include "TApplication.h"
#include "TCanvas.h"
#include "TLine.h"
#include "THStack.h"
#include "TLegend.h"

#include <iostream>

int main(int argc, char** argv) {
	std::string fname(argv[1]);
	std::string fname_dirt(argv[2]);
	std::string fname_data(argv[3]);

	std::string param = "RecoEnuQE";
	//std::vector<double> bins_reco = {200,300.,375.,475.,550.,675.,800.,950.,1100.,1300.,1500.,3000.};
	//std::vector<double> bins_reco = {150,200,300.,375.,475.,550.,675.,800.,950.,1100.,1300.,1500.,3000.};
	//std::vector<double> bins_reco = {150,200,250,300.,350,400,450,500,550,600,650,700,750,800.,850,900,950.,1000,1100.,1300.,1500.,3000.};
	std::vector<double> bins_reco = {150,175,200,225,250,275,300.,325,350,375,400,425,450,475,500,525,550,575,600,650,700,750,800.,850,900,950.,1000,1100.,1300.,1500.,3000.};
	
	std::vector<double> summed_mc(bins_reco.size()-1,0.0);

	double pot_scale = 6.46/41.10;
	//Published total backgrounds
	std::vector<double> minibkg = {180.80171,108.22448,120.03353,63.887782,89.806966,67.249431,69.855878,57.014477,51.846417,38.586738,69.381391};

	//for(auto& v : bins_reco) v /= 1000.0;


	std::cout<<"Loading Data"<<std::endl;
	TFile * f_data = new TFile( (fname_data).c_str());
	TTree * data = (TTree*)f_data->Get("MiniBooNE_CCQE");
	float d_Weight,d_Energy,d_CosTheta,d_RecoEnuQE;
	data->SetBranchAddress("Weight",&d_Weight);
	data->SetBranchAddress("Energy",&d_Energy);
	data->SetBranchAddress("CosTheta",&d_CosTheta);
	data->SetBranchAddress("RecoEnuQE", &d_RecoEnuQE);

	TH1D * h_data = new TH1D("h_data","h_data", bins_reco.size()-1,  &bins_reco[0] );

	for(int i=0; i< data->GetEntries();i++){
		data->GetEntry(i);
		h_data->Fill(1000*d_RecoEnuQE, d_Weight);
	}

	std::vector<double>summed_data;
	std::vector<double>summed_excess;

	for(int i=0; i<bins_reco.size()-1;i++){
		summed_data.push_back(h_data->GetBinContent(i+1));		
	}




	sp::SPIO spio;
	spio.set_verbosity((sp::msg::Level_t)0);

	//
	// Get the vector of histograms for each background type
	//
	std::vector<int> nicer_colors = {kRed,kOrange-7, kRed-7, kOrange-2, kGreen-2, kGreen+1, kGreen-8, kGray,kGray};


	auto th1d_v = spio.gen_background(fname,fname_dirt,param,bins_reco);

	std::vector<double> summed_mc_err(bins_reco.size()-1,0.0);
	TH1D* sum_bkg = (TH1D*)th1d_v.at(0).Clone("sum");
	sum_bkg->Reset();
	sum_bkg->Sumw2();
	  TList *list = new TList;
  	  for(auto &v: th1d_v){
		list->Add(&v);
	}
  	sum_bkg->Merge(list);
	sum_bkg->Scale(pot_scale);
	for(int i=0; i<bins_reco.size()-1; i++){
		summed_mc_err.at(i) = sum_bkg->GetBinError(i+1);
	}






	TApplication app("app", 0, 0);


	THStack ths("ths","");
	TLegend tl(0.5,0.6,0.89,0.89);
	tl.SetLineColor(kWhite);

	for(size_t bkgd_id = (size_t) sp::kBKGD_MAX-1 ; bkgd_id >0 ; --bkgd_id) {
		if ((sp::StackedBkgdType_t)bkgd_id == sp::kBKGD_INVALID) continue;
		//if ((sp::StackedBkgdType_t)bkgd_id == sp::kBKGD_DIRT)    continue; //Mark, added dirt

		auto& th1d = th1d_v[bkgd_id];

		for(size_t bin_id=1; bin_id < bins_reco.size(); ++bin_id) {

			summed_mc.at(bin_id-1) +=th1d.GetBinContent(bin_id)*pot_scale;


			auto dx = bins_reco.at(bin_id) - bins_reco.at(bin_id-1);

			auto modifier = 1.0 / (dx );
			modifier *= pot_scale; // POT normalize

			auto bin_content = th1d.GetBinContent(bin_id) * modifier;
			auto bin_error   = th1d.GetBinError(bin_id) * modifier;

			th1d.SetBinContent(bin_id, bin_content);
			th1d.SetBinError(bin_id, bin_error);
		}
		th1d.SetLineWidth(2);
		th1d.SetLineColor(kBlack);
		th1d.SetFillColor( nicer_colors.at(bkgd_id)  );
		ths.Add(&th1d);
		tl.AddEntry(&th1d, sp::StackedBkgd2String((sp::StackedBkgdType_t)bkgd_id).c_str(), "f");
	}

	tl.AddEntry(h_data,"MiniBooNE Data 6.46e20 POT","lp");

	for(int i=0; i< summed_mc.size();i++){
		//std::cout<<"Sum: "<<bins_reco.at(i)<<" Ours: "<<summed_mc.at(i)<<" Ratio: "<<summed_mc.at(i)/minibkg.at(i)<<std::endl;
	}



	std::vector<double> excess_err;
	for(int i=0; i<bins_reco.size()-1;i++){
		summed_excess.push_back( summed_data.at(i)-summed_mc.at(i));
		excess_err.push_back(sqrt(summed_data.at(i) +pow(summed_mc_err.at(i),2) )/h_data->GetBinWidth(i+1) );
	}
	

	TCanvas *c = new TCanvas("bkg","bkg",800,800);
	TPad *p1 =  new TPad("p1","p1",0,0.3,1,1.0);
	p1->SetBottomMargin(0);
	p1->Draw();
	p1->cd();

	ths.Draw("hist");
	ths.GetXaxis()->SetRangeUser(bins_reco.front()  ,1000);
	ths.SetMaximum(3);
	ths.GetYaxis()->SetTitle("Events/MeV");

	h_data->SetMarkerColor(kBlack);
	h_data->SetMarkerStyle(20);
	h_data->SetLineColor(kBlack);
	h_data->SetLineWidth(2);
	h_data->SetMarkerSize(1);
	h_data->Scale(1,"width");	
	h_data->Draw("E1 same");

	tl.Draw("sames");
	//	c.Update();
	//	c.Modified();


	c->cd();
	TPad * p2 = new TPad("p2","p2",0,0.05,1,0.3);
	p2->SetTopMargin(0);
	p2->SetBottomMargin(0.2);
//	p2->SetGrid();
	p2->Draw();
	p2->cd();

	TH1D *excess = (TH1D*)h_data->Clone("excess");
	TH1D *bkg  = (TH1D*)(ths.GetStack()->Last());
	excess->SetTitle("");
	excess->Sumw2();
	excess->SetStats(0);
	excess->Add(bkg,-1);
	excess->SetError(&excess_err[0]);
	excess->Draw("ep");
	excess->GetXaxis()->SetRangeUser(bins_reco.front(), 1000);
	excess->GetXaxis()->SetTitle("E^{CCQE} [MeV]");
	excess->GetYaxis()->SetTitle("Excess/MeV");

	excess->GetYaxis()->SetTitleSize(20);
	excess->GetYaxis()->SetTitleFont(65);
	excess->GetYaxis()->SetTitleOffset(1.55);
	excess->GetYaxis()->SetLabelFont(65); // Absolute font size in pixel (precision 3)
	excess->GetYaxis()->SetLabelSize(15);

	// X axis ratio plot settings
	excess->GetXaxis()->SetTitleSize(20);
	excess->GetXaxis()->SetTitleFont(65);
	excess->GetXaxis()->SetTitleOffset(4.);
	excess->GetXaxis()->SetLabelFont(65); // Absolute font size in pixel (precision 3)
	excess->GetXaxis()->SetLabelSize(15);
	
	TLine * l = new TLine(bins_reco.front(),0, 1100,0);
	l->Draw();

	p2->Update();
	c->Update();
	c->Modified();

	c->SaveAs("binning_scheme.pdf","pdf");




	app.Run();




	return 0;
}
