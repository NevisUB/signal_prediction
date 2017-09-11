#include "Unfold/Core/SPIO.h"
#include "Unfold/Algo/UnfoldAlgoDAgnostini.h"
#include "Unfold/Algo/UnfoldAlgoInverse.h"
#include "Unfold/Algo/UnfoldAlgoSVD.h"
#include "Unfold/Algo/ModelNueCCQE.h"

#include "Combined/CombinedTypes.h"
#include "Combined/CombinedUtil.h"

#include "TCanvas.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TLine.h"
#include "TLatex.h"

int main(int argc, char** argv) {

	TRandom3 * rangen = new TRandom3(0);

	gStyle->SetOptStat(0);
	sp::SPIO a;
	a.set_verbosity((sp::msg::Level_t)0);

	std::string core = "/home/mark/work/uBooNE/lee_unfolding/";

	a.add_mc_in_file(core+"rootfiles/filterd_ccqe_nue_nuebar.root");
	//a.add_mc_in_file("/home/vgenty/signal/simplifyTreeOsc/filtered_ccqe_nue_nuebar/filterd_ccqe_nue_nuebar.root");
	a.set_mc_tree_name("MiniBooNE_CCQE");

	sp::ModelNueCCQE model;
	a.set_model(&model);

	a.initialize();

	std::vector<std::string> var_truth = {"NuMomT"};
	std::vector<std::string> var_reco = {"RecoEnuQE"};
	std::vector<double> bins_reco = {200,  300.,  375. , 475.,  550.,  675.,  800.,  950.,  1100.  ,1300. , 1500. , 2000};
	//	std::vector<double> bins_reco = {140,200.,250,300.,350 ,400,450,500., 550., 600, 700., 800., 950.,   1100., 1300., 1500., 3000.};
	std::vector<double> bins_truth = {200,250,300,350,400,450, 500, 550, 650,  800.,1000,2000};
	//std::vector<double> bins_truth = {140,200.0,250,300,350,400,300,320,340,360,380,400,440,480,500, 550, 600., 650,   800.,1000,3000};
	int N_bins_reco = bins_reco.size()-1;
	int N_bins_truth = bins_truth.size()-1;
	std::cout<<"# True Bins: "<<N_bins_truth<<" # Reco Bins: "<<N_bins_reco<<std::endl;

	a.add_true_parameter(var_truth,bins_truth,sp::kOP_GEV2MEV);
	a.add_reco_parameter(var_reco,bins_reco,sp::kOP_GEV2MEV);

	std::cout << "int_response matrix" << std::endl;
	a.init_response_matrix();

	std::cout << "fill resp" << std::endl;
	a.fill_responses();

	std::cout << "write" << std::endl;
	a.write_unfold_file();

	std::cout << "Beginning" << std::endl;

	sp::UnfoldAlgoDAgnostini alg; 
	//sp::UnfoldAlgoSVD alg; 

	alg.set_verbosity((sp::msg::Level_t)0);
	alg.Initialize(&(a.Responses().front()));


	double pot_scale = 6.46/41.10;


	//Got to tihnk more bout flows
	std::vector<double> miniobs = {232,  156,  156,   79,   81,   70,   63,   65,   62,   34,   70};
	std::vector<double> minibkg = {180.80171,108.22448,120.03353,63.887782,89.806966,67.249431,69.855878,57.014477,51.846417,38.586738,69.381391};
	double Nsignal = 0;

	std::cout<<"Loading in MC and Dirt."<<std::endl;
	std::string fname = core + "rootfiles/filtered_passosc.root";
	std::string fname_dirt = core + "/rootfiles/merged_filtered_out_osc_mc_dirt.root";
	std::string param = var_reco.at(0);
	std::vector<double> summed_mc(N_bins_reco,0.0);

	std::vector<double> sigerr(N_bins_reco,0);
	//
	// Get the vector of histograms for each background type
	//
	sp::SPIO spio;
	spio.set_verbosity((sp::msg::Level_t)0);

	std::cout<<"Generating Background"<<std::endl;
	auto th1d_v = spio.gen_background(fname,fname_dirt,param,bins_reco);

	for(auto & v :th1d_v){
		for(int i=0; i< N_bins_reco;i++){
			summed_mc.at(i) += v.GetBinContent(i+1)*pot_scale;
		}

	}
	/*
	   for(size_t bkgd_id = 0; bkgd_id < (size_t) sp::kBKGD_MAX; ++bkgd_id) {
	   if ((sp::StackedBkgdType_t)bkgd_id == sp::kBKGD_INVALID) continue;
	   std::cout<<"Loading up: "<<sp::StackedBkgd2String((sp::StackedBkgdType_t)bkgd_id)<<std::endl;
	   auto& th1d = th1d_v[bkgd_id];

	   for(size_t bin_id=1; bin_id < bins_reco.size(); ++bin_id) {

	   summed_mc.at(bin_id-1) += th1d.GetBinContent(bin_id)*pot_scale;

	   }
	   }
	 */


	std::cout<<"Loading Data"<<std::endl;
	TFile * f_data = new TFile( (core+"rootfiles/output_osc_data_detail_1.root").c_str());
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



	/*std::cout<<"Have loaded Data and MC"<<std::endl;
	  for(int i=0; i< summed_mc.size();i++){
	  std::cout<<"Data || Ours: "<<h_data->GetBinContent(i+1)<<" Published: "<<miniobs.at(i)<<std::endl;
	  std::cout<<"MC || Ours: "<<summed_mc.at(i)<<" Published: "<<minibkg.at(i)<<std::endl;
	  }*/


	//changed briefly... 
	TVectorD mini_signal(N_bins_reco);
	for(int i=0; i<N_bins_reco; i++){
		//This one is the old working one
		//mini_signal(i) = (miniobs.at(i)-minibkg.at(i))+alg.r(i)*pot_scale    ;
		mini_signal(i) = (h_data->GetBinContent(i+1)-summed_mc.at(i))+alg.r(i)*pot_scale    ;


		Nsignal += mini_signal(i);
		sigerr.at(i) =sqrt(h_data->GetBinContent(i+1)+summed_mc.at(i)) ;

		if(sigerr.at(i)!=sigerr.at(i)) {
			std::cout<<"Failure, should fabs"<<std::endl;
			exit(EXIT_FAILURE);

		}
	}

	TMatrixD sigcorr(N_bins_reco, N_bins_reco);
	sigcorr.Zero();
	for(int i=0; i<N_bins_reco;i++){
		sigcorr(i,i)=h_data->GetBinContent(i+1)+summed_mc.at(i);
		std::cout<<"MYD(i,i) "<<sigcorr(i,i)<<std::endl;
	}
	alg.Setd(&mini_signal);
	alg.SetD(&sigcorr);
	std::cout<<"juuust setting D"<<std::endl;
	alg.Unfold();


	//	alg.SetRegularization(3);
	//	alg.TestUnfolding("CCQE_poison_unfold_test");













	TH1D ans =  alg.GetHistU();
	TH1D sig = alg.GetHistD();

	TH1D sig_events = alg.GetHistD(); //pure for chi^2 calc

	TH1D reco = alg.GetHistR();
	TH1D truth = alg.GetHistT();
	sig.SetLineColor(kRed); 
	ans.SetLineColor(kBlue); 
	truth.SetLineColor(kBlack);
	reco.SetLineColor(kBlack); 

	sig.SetMarkerStyle(21);
	sig.SetError(&sigerr[0]);
	sig.Scale(1,"width");
	ans.Scale(1,"width");
	truth.Scale(pot_scale,"width");
	reco.Scale(pot_scale,"width");

	TFile *f = new TFile("CCQE_ans.root","RECREATE");







	TCanvas *cb = new TCanvas("a","a",1200,1200);
	cb->Divide(2,2);
	cb->cd(1);//->SetLogy();
	TH1D tt = truth;
	TH1D rr = reco;
	tt.SetMarkerStyle(21);
	rr.SetMarkerStyle(21);
	tt.SetMarkerColor(kBlue-6);
	rr.SetMarkerColor(kRed-7);
	tt.SetLineWidth(2);
	rr.SetLineWidth(2);
	tt.SetLineColor(kBlue-6);
	rr.SetLineColor(kRed-7);
	tt.GetXaxis()->SetTitle("Energy (Truth or Reco) [MeV]");
	tt.GetYaxis()->SetTitle("Events/MeV");
	tt.GetXaxis()->SetTitleSize(0.04);
	tt.GetXaxis()->SetLabelSize(0.04);
	tt.GetYaxis()->SetLabelSize(0.04);
	tt.GetYaxis()->SetTitleOffset(1.45);
	tt.SetMinimum(0);
	tt.SetTitle("Intrinsic #nu_{e} CCQE only");
	tt.SetMaximum(5);
	tt.GetXaxis()->SetRange(1,10);

	TLegend * legr = new TLegend(0.58,0.6,0.89,0.89);
	legr->AddEntry(&tt,"True E_{#nu} MC","lep");
	legr->AddEntry(&rr,"Reco E_{QE} MC","lep");
	legr->SetFillStyle(0);
	//	legr->SetBorderSize(0.1);

	tt.Draw();
	rr.Draw("same ap");
	legr->Draw();

	cb->cd(2)->SetRightMargin(0.175);
	TH2D prob_MC = alg.GetHistA();

	prob_MC.GetYaxis()->SetTitle("Truth E_{#nu} Bin");
	prob_MC.GetXaxis()->SetTitle("Reco E_{QE} Bin");
	prob_MC.SetTitle("Response Matrix");
	prob_MC.GetXaxis()->SetTitleOffset(1.1);
	prob_MC.GetYaxis()->SetTitleOffset(1.1);

	prob_MC.Draw("colz");
	prob_MC.GetXaxis()->SetRange(1,10);
	prob_MC.GetYaxis()->SetRange(1,10);

	cb->cd(3);
	TH1D eff = alg.GetHistEff();

	std::vector<double> errEff(alg.n_t);
	for(int b=0; b<errEff.size();b++){
		errEff.at(b)=sqrt(alg.Ep(b,b));  ;
	}
	eff.SetError(&errEff[0]);

	eff.GetYaxis()->SetTitle("Efficiency");
	eff.GetXaxis()->SetTitle("Truth Bin");
	eff.Draw("E1");
	eff.SetMinimum(0);
	eff.SetMaximum(0.4);


	cb->SaveAs("CCQE_signals_response.pdf","pdf");



	//return 0;












	///////////////////////////////////////////////////////////////////////

	TCanvas *c = new TCanvas("c","c",1200,800);



	sig.SetTitle("Reconstructed E^{QE}" );
	sig.GetYaxis()->SetTitle("Events/MeV");
	sig.GetXaxis()->SetTitle("Reco E^{QE}_{#nu} [MeV]");
	sig.SetLineColor(kRed-7);
	sig.SetMarkerColor(kRed-7);
	sig.SetMarkerStyle(21);
	sig.SetLineWidth(2);
	sig.Draw("e1");
	sig.GetXaxis()->SetRange(1,10);

	reco.SetFillColor(kGray);
	reco.SetMarkerColor(kBlack);
	reco.SetMarkerStyle(21);
	//	reco.Draw("same e2");
	reco.Draw("same hist");
	sig.Draw("same");

	TLegend * leg1 = new TLegend(0.58,0.6,0.89,0.89);
	leg1->AddEntry(&sig,"6.46e20 POT #nu-mode data","lep");
	//leg1->AddEntry(&reco,"MC intrinsic #nu_{e} CCQE","lef");
	leg1->AddEntry(&reco,"MC intrinsic #nu_{e} CCQE","l");
	leg1->SetFillStyle(0);
	leg1->SetBorderSize(0.0);
	leg1->Draw();







	c->Write();
	c->SaveAs("CCQE_signal.pdf","pdf");


	TCanvas *c2 =  new TCanvas("c1","c1",1200,1200);
	TCanvas *cu2 =  new TCanvas("cu2","cu2",1200,1200);
	TCanvas *cU =  new TCanvas("c2","c2",3000,3000);
	TCanvas *cR =  new TCanvas("cR","cR",3000,3000);
	TCanvas *cBias = new TCanvas("cBias","cBias",3000,3000);

	c2->Divide(3,3);
	cU->Divide(3,3);
	cR->Divide(3,3);
	cBias->Divide(3,3);


	std::vector<int> kreg = {1,2,3,4,5,6,7,8,9};
	std::vector<int> cols = {kBlue-7,kGreen-6,kRed-7, kOrange-3, kMagenta-3, kGreen+3, kBlue-7,kGreen-6,kRed-7};
	std::vector<TH1D> us(kreg.size());
	std::vector<TH2D> US(kreg.size());
	std::vector<TH1D> uR(kreg.size());
	std::vector<TH1D> uBias(kreg.size());
	std::vector<TLegend*> leg(kreg.size());
	std::vector<TH1D> us_stat(kreg.size());
	//std::vector<int> kreg = {1,2,5,8,100};

	//Some basic Lcurve and Chi^2 analysis on refolded spectra
	TCanvas *c_chicurve = new TCanvas("lc","lc,",1200,1200);	
	std::vector<double> refold_avg_chi;
	std::vector<double> refold_k;

	std::vector<double> bias_avg;


	double N_random_refold = 10.0;
	int N_mc_bias = 10;


	for(int k=0; k<kreg.size(); k++){
		std::string nam = "Reg " +std::to_string(kreg.at(k)) ;

		c2->cd(k+1)->SetLeftMargin(0.1);
		alg.SetRegularization(kreg.at(k));
		alg.Unfold();

		//Section to re-fold and calc numbers of events.
		TVectorD refolded = alg.A*alg.u;
		double Nrefold = 0;
		for(int i=0; i<alg.n_r; i++){
			Nrefold+=refolded(i);
		}
		std::cout<<std::setprecision(10)<<"REFOLD || On k="<<kreg.at(k)<<" and refolded-number-of-events is "<<Nrefold<<" . Data is "<<Nsignal<<std::endl;





		std::vector<double> errA(alg.n_t,0.0);
		std::vector<double> errD(alg.n_t,0.0);
		std::vector<double> err(alg.n_t,0.0);
		std::vector<double> errS(alg.n_t,0.0);

		for(int i=0;i<err.size(); i++){
			errD.at(i) = sqrt(fabs(alg.U(i,i)));
			errA.at(i) = sqrt(fabs(alg.UA(i,i)));
			errS.at(i) = sqrt(fabs(alg.u(i)));
			err.at(i) = sqrt(fabs(alg.UA(i,i)+alg.U(i,i)));
			std::cout<<"u= "<<alg.u(i)<<" +/-(stat) "<<errS.at(i)<<" +/-(D) "<<errD.at(i)<<" +/-(A) "<<errA.at(i)<<" +/-(AD) "<<err.at(i)<<" bias "<<alg.b(i)<<" +/- "<<sqrt(alg.B(i,i))<<std::endl;
		}


		us.at(k) = alg.GetHistU();
		us.at(k).SetTitle(  nam.c_str() );
		us.at(k).GetYaxis()->SetTitle("Events/MeV");
		us.at(k).GetXaxis()->SetTitle("True E_{#nu} [MeV]");
		us.at(k).SetError(&err[0]);
		us.at(k).SetLineColor(kBlack);
		us.at(k).SetLineWidth(2);
		us.at(k).SetFillColor(cols.at(k));
		us.at(k).Scale(1,"width");
		us.at(k).Draw("E2");

		us.at(k).SetMinimum(0);
		us.at(k).GetXaxis()->SetRangeUser(bins_truth.front(),1000);
		us.at(k).SetMaximum(10);//has to be after!

		us_stat.at(k) = alg.GetHistU();
		us_stat.at(k).SetError(&errS[0]);
		us_stat.at(k).SetLineColor(kBlack);
		us_stat.at(k).SetLineWidth(1);
		us_stat.at(k).Scale(1,"width");
		us_stat.at(k).Draw("same E1");

		truth.SetLineColor(kBlack);
		truth.SetLineWidth(2);
		truth.SetMarkerStyle(21);
		truth.SetMarkerColor(kBlack);
		truth.Draw("same hist");

		leg.at(k) = new TLegend(0.5,0.65,0.89,0.89);
		leg.at(k)->AddEntry(&us.at(k),"Raw excess","lef");
		leg.at(k)->AddEntry(&truth,"MC True #nu_{e} CCQE","l");
		leg.at(k)->SetFillStyle(0);
		leg.at(k)->SetBorderSize(0);
		leg.at(k)->Draw();

		if(k==1){
			cu2->cd();
			us_stat.at(k).SetFillColor(kGreen+2);
			us_stat.at(k).SetLineColor(kGreen+2);
			us_stat.at(k).SetMarkerColor(kBlack);
			us_stat.at(k).SetMarkerStyle(29);
			us.at(k).Draw("E2");
			us_stat.at(k).Draw("same E2");
			us_stat.at(k).Draw("same P");
			truth.Draw("same hist");
			leg.at(k)->Draw();
			us.at(k).GetXaxis()->SetRangeUser(bins_truth.front(),1000);
		}


		cU->cd(k+1);//->SetLogz();
		US.at(k) = alg.GetCovU();
		US.at(k).SetTitle(  nam.c_str() );
		US.at(k).GetYaxis()->SetTitle("True E_{#nu} Bin a");
		US.at(k).GetXaxis()->SetTitle("True E_{#nu} Bin b");
		US.at(k).Draw("colz");

		cR->cd(k+1);
		uR.at(k) = alg.GetHistRefold();
		uR.at(k).SetTitle(nam.c_str());
		double ch2 = uR.at(k).Chi2Test(&sig_events,"CHI2,WW");

		double Nrefold_2 = uR.at(k).GetSumOfWeights();
		uR.at(k).Scale(1,"width");
		uR.at(k).SetMarkerStyle(29);
		uR.at(k).SetMarkerColor(kBlack);
		uR.at(k).SetLineColor(kBlack);
		uR.at(k).SetMarkerSize(2);


		//RandomRefolding && MC bias calculation
		std::vector<double> bias(alg.n_t,0);
		std::vector<double> avgu(alg.n_t,0);

		double avg_refold_chi = 0;
		for(int k=0;k<(int)N_random_refold;k++){
			TH1D tmp = alg.GetHistRandRefold(rangen);

			double temp_chi = tmp.Chi2Test(&sig_events,"CHI2,WW");

			avg_refold_chi +=temp_chi;
			std::cout<<"RandRefold: "<<k<<" chi^2: "<<temp_chi<<std::endl;


			//and MC bias calc
			TVectorD tmpV(alg.n_r);
			for(int i=0; i<alg.n_r; i++){
				tmpV(i)=tmp.GetBinContent(i+1);
			}

			auto tmp_alg = alg;
			tmp_alg.Setd(&tmpV);
			tmp_alg.Unfold();
			for(int a=0; a<alg.n_t; a++){
				avgu.at(a) += tmp_alg.u(a)/N_random_refold;
			}

		}

		for(int a=0;a<alg.n_t; a++){
			bias.at(a) = avgu.at(a)-alg.u(a);
			std::cout<<"ThisBias bin:"<<a<<" "<<bias.at(a)<<" "<<alg.b(a)<<std::endl;
		}

		double avg_bias = 0;//Just average the bins		
		for(auto b: bias){
			avg_bias += b/((double)N_bins_truth);
		}

		bias_avg.push_back(pow(fabs(avg_bias),2.0));
		std::cout<<"ThisBias Avg: "<<avg_bias<<std::endl;





		//Step 1. Get random (but slightly different "true" spectra)


		/* This is multi-truth biasing
		   TVectorD pois_u = alg.u;
		   std::vector<double> bias(alg.n_t,0.0);
		   for(int i=0; i< N_mc_bias; i++){

		   for(int a=0; a<alg.n_t; a++){
		   pois_u(a) = rangen->Poisson( alg.u(a) ); 
		   }


		   TVectorD re = alg.A*pois_u;
		   auto tmp_alg = alg;
		   tmp_alg.Setd(&re);
		   tmp_alg.Unfold();

		   for(int a=0; a<alg.n_t; a++){
		   double thisbias = tmp_alg.u(a)-pois_u(a);
		   bias.at(a) += (thisbias)/((double)N_mc_bias);
		   }			
		   }
		   double avg_bias = 0;			
		   for(auto b: bias){
		   avg_bias += b/((double)N_bins_truth);
		   }

		   bias_avg.push_back(pow(fabs(avg_bias),2.0));
		   std::cout<<"ThisBias Avg: "<<avg_bias<<std::endl;
		 */

		//This is multi 


		std::cout<<"RandRefoldAverage: "<<avg_refold_chi/N_random_refold<<std::endl;
		refold_avg_chi.push_back(avg_refold_chi/(N_random_refold*(double)N_bins_reco));
		refold_k.push_back((double)kreg.at(k));

		std::string s_ch2 = "#chi^{2}/ndof : " + std::to_string(ch2)+"/"+std::to_string(N_bins_reco);
		TLegend * legR = new TLegend(0.3,0.7,0.89,0.89);
		legR->SetHeader(s_ch2.c_str() );
		legR->AddEntry(&sig, ("MiniBooNE Excess: " + std::to_string(Nsignal)).c_str() ,"lp");
		legR->AddEntry(&uR.at(k), ("Refolded Excess: "+std::to_string(Nrefold_2)).c_str(),"l");
		uR.at(k).Draw("hist");
		sig.Draw("same");
		legR->Draw();
		uR.at(k).GetXaxis()->SetRangeUser(bins_truth.front(),1000);
		uR.at(k).SetMaximum(1.2);




		//Bias plot filling
		cBias->cd(k+1);
		uBias.at(k) = alg.GetHistBias();	
		uBias.at(k).SetTitle(  nam.c_str() );
		uBias.at(k).GetYaxis()->SetTitle("Bias/MeV");
		uBias.at(k).GetXaxis()->SetTitle("True E_{#nu} [MeV]");
		uBias.at(k).SetLineColor(kBlack);
		uBias.at(k).SetLineWidth(2);
		uBias.at(k).SetFillColor(cols.at(k));
		uBias.at(k).Scale(1,"width"); // should it be bias /MeV?
		uBias.at(k).Draw("E2");
		uBias.at(k).GetXaxis()->SetRangeUser(bins_truth.front(),1000);
/* At this point the standard
deviations of the biases are approximately equal to the biases themselves, and
therefore any further bias reduction would introduce as much error as it removes.
*/
	}

	alg.TestRegularization("CCQE_lcurves", 1,11,11);
	//alg.SetDirectRegularization(1000);
	//alg.TestRegularization("CCQE_lcurves_direct", -8,8,200);  //-20 13

	c_chicurve->cd();
	c_chicurve->SetLogy();

	TGraph *g_chi = new TGraph(kreg.size(), &refold_k[0], &refold_avg_chi[0] );
	g_chi->SetTitle("Average #chi^{2}/ndof ");
	g_chi->GetXaxis()->SetTitle("Regularisation Parameter k");
	g_chi->GetYaxis()->SetTitle("#chi^{2}/ndof ");
	g_chi->SetMarkerSize(2);
	g_chi->SetMarkerStyle(2);
	g_chi->Draw("ACP");
	TLine l_chicurve(refold_k.front(),1,refold_k.back(),1);
	l_chicurve.Draw();
	c_chicurve->Write();
	c_chicurve->SaveAs("CCQE_chicurve.pdf","pdf");

	TCanvas * c_bias = new TCanvas("cbias","cbias",1200,1200);
	c_bias->cd();
	c_bias->SetLogy();

	TGraph *g_bias = new TGraph(kreg.size(), &refold_k[0], &bias_avg[0] );
	g_bias->SetTitle("Average Bias MC ");
	g_bias->GetXaxis()->SetTitle("Regularisation Parameter k");
	g_bias->GetYaxis()->SetTitle("<bias> ");
	g_bias->SetMarkerSize(2);
	g_bias->SetMarkerStyle(2);
	g_bias->Draw("ACP");
	//TLine l_bias(refold_k.front(),1,refold_k.back(),1);
	//l_bias.Draw();
	c_bias->Write();
	c_bias->SaveAs("CCQE_bias.pdf","pdf");





	c2->Write();
	cU->Write();
	cBias->Write();
	cu2->SaveAs("CCQE_unfolded.pdf","pdf");
	c2->SaveAs("CCQE_reg_vary.pdf","pdf");
	cU->SaveAs("CCQE_reg_corr.pdf","pdf");
	cR->SaveAs("CCQE_reg_refold.pdf","pdf");
	cBias->SaveAs("CCQE_reg_bias.pdf","pdf");


	TCanvas *cr = new TCanvas();
	TH1D ratio = us.at(1);	
	ratio.Divide(&truth);
	ratio.SetFillColor(kBlue-7);
	ratio.GetYaxis()->SetTitle("Ratio to MiniBooNE MC Central Value");
	ratio.SetMarkerStyle(5);
	ratio.SetMarkerSize(1);
	ratio.Draw("e2");
	ratio.SetMaximum(22);
	ratio.GetXaxis()->SetRangeUser(bins_truth.front(),1000);

	TLegend * leg2 = new TLegend(0.58,0.6,0.89,0.89);
	leg2->AddEntry(&ratio,"Intrinsic #nu_{e} CCQE Model","lf");
	leg2->SetFillStyle(0);
	leg2->SetBorderSize(0.0);
	leg2->Draw();

	TLine *line = new TLine(bins_truth.front(),1,1000,1);
	line->SetLineColor(kBlack);
	line->SetLineStyle(2);
	line->Draw();


	TCanvas *cr_scaled = new TCanvas();
	TFile * fr = new TFile((core+"cross_section_generator_difference/ratio_out.root").c_str());
	TGraph * gr = (TGraph*)fr->Get("Graph");
	TH1D * ratio_scaled = (TH1D*)ratio.Clone("scaled");

	std::vector<double> bincenter;
	std::vector<double> binval;


	for(int i=0; i<ratio_scaled->GetNbinsX()+2; i++){
		ratio_scaled->SetBinContent(i+1, ratio.GetBinContent(i+1)/gr->Eval(0.001*ratio.GetBinCenter(i+1) )  );
		bincenter.push_back( ratio.GetBinCenter(i+1) );
		binval.push_back(  ratio.GetBinContent(i+1)/gr->Eval(0.001*ratio.GetBinCenter(i+1) ));

		std::cout<<"Ratio: "<<ratio.GetBinContent(i+1)<<" "<<ratio.GetBinCenter(i+1)<<" "<<gr->Eval(0.001*ratio.GetBinCenter(i+1) )<<std::endl; 
		std::cout<<"Res: "<<ratio_scaled->GetBinContent(i+1)<<std::endl;
	}



	ratio_scaled->Draw("e2");
	ratio_scaled->SetMaximum(10);
	ratio_scaled->GetXaxis()->SetRangeUser(bins_truth.front(),1000);

	leg2->Draw();
	line->Draw();
	cr_scaled->SaveAs("CCQE_model_ratio_scaleo.pdf","pdf");





	for(int i=0; i<=ratio.GetNbinsX()+2; i++){
		std::cout<<"Ratio "<<ratio.GetBinContent(i)<<" "<<ratio_scaled->GetBinContent(i)<<std::endl;
	}
	cr->Write();
	cr->SaveAs("CCQE_model_ratio.pdf","pdf");

	f->Close();

	TFile* f_graph_out =  new TFile("CCQE_final_tgraph.root","RECREATE");
	f_graph_out->cd();
	TGraph * graph_out = new TGraph(bincenter.size(), &bincenter[0], &binval[0]  );
	graph_out->Write();
	f_graph_out->Close();





	return 0;
}
