#include "Unfold/Core/SPIO.h"
#include "Unfold/Algo/UnfoldAlgoDAgnostini.h"
#include "Unfold/Algo/UnfoldAlgoInverse.h"
#include "Unfold/Algo/UnfoldAlgoSVD.h"
#include "Unfold/Algo/ModelNueCCQE.h"
#include "Unfold/Algo/ModelNCDelta.h"

#include "Combined/CombinedTypes.h"
#include "Combined/CombinedUtil.h"

#include "TCanvas.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TLine.h"
#include "TLatex.h"

#include <getopt.h>
#define no_argument 0
#define required_argument 1
#define optional_argument 2


#define ALGO_SVD 0
#define ALGO_DAG 1

#define MODEL_CCQE 0
#define MODEL_DELTARES 1

int main(int argc, char** argv) {
	int c;
	int index; 
	int iarg = 0;
	opterr=1;


	int ALGO = ALGO_SVD;
	int MODEL = MODEL_CCQE;

	std::string base_name = "CCQE_SVD";
	std::string algo_name = "SVD";
	std::string model_name = "CCQE";


	const struct option longopts[] = 
	{
		{"algo",	 	required_argument, 	0, 'a'},
		{"model",		required_argument,	0, 'm'},
		{0,			no_argument, 		0,  0}
	};
	const char * tin;

	while(iarg != -1)
	{
		iarg = getopt_long(argc,argv, "a:m:", longopts, &index);

		switch(iarg)
		{
			case 'a':
				tin = optarg;//way to read in null terminated c strings and compare to known list
				if(strcmp(tin, "svd")==0){  ALGO = ALGO_SVD; algo_name = "SVD";}
				if(strcmp(tin, "dag")==0) { ALGO = ALGO_DAG; algo_name = "DAG";}
				break;
			case 'm':
				tin = optarg;//way to read in null terminated c strings and compare to known list
				if(strcmp(tin, "ccqe")==0){ MODEL = MODEL_CCQE; model_name = "CCQE";}
				if(strcmp(tin, "deltares")==0) {  MODEL = MODEL_DELTARES; model_name = "DELTARES";}
				break;

			case '?':
				std::cout<<"Abandon hope all ye who enter this value. "<<std::endl;
			case 'h':
				std::cout<<"******************************************"<<std::endl;
				std::cout<<"Allowed arguments:"<<std::endl;
				std::cout<<"\t-m, --model\t\t\tSets the model, ncdelta or ccqe [default ccqe]"<<std::endl;
				std::cout<<"\t-a, --algo\t\t\tSets the unfolding algo, svd or dag [default svd]"<<std::endl;
				std::cout<<"******************************************"<<std::endl;

				return 0;
		}

	}

	//Define naming scheme
	base_name = model_name +"_"+algo_name;





	TRandom3 * rangen = new TRandom3(0);

	gStyle->SetOptStat(0);
	sp::SPIO a;
	a.set_verbosity((sp::msg::Level_t)0);


	std::string core = "/home/mark/work/uBooNE/lee_unfolding/";
	std::string bkg_file = "rootfiles/filtered_passosc.root";
	std::string data_file = "rootfiles/output_osc_data_detail_1.root";	
	std::string dirt_file =  "rootfiles/merged_filtered_out_osc_mc_dirt.root";


	std::string use_file;
	switch(MODEL)
	{
		case MODEL_CCQE:
			use_file = "rootfiles/filterd_ccqe_nue_nuebar.root";
			break;
		case MODEL_DELTARES:
			use_file= "rootfiles/filtered_nc_delta.root";
			break;
	}


	a.add_mc_in_file(core +use_file);
	a.set_mc_tree_name("MiniBooNE_CCQE");


	sp::ModelNueCCQE CCmodel;
	sp::ModelNCDelta NCmodel;

	switch(MODEL)
	{
		case MODEL_CCQE:
			a.set_model(&CCmodel); 
			break;
		case MODEL_DELTARES:
			a.set_model(&NCmodel);
			break;
	}


	a.initialize();

	std::vector<std::string> var_truth = {"NuMomT"};
	std::vector<std::string> var_reco = {"RecoEnuQE"};



	std::vector<double> bins_truth ;//=  {300,   475.,  550.,  675.,  800.,  950.,  1100.  ,1300. , 1500.,2500 };
	std::vector<double> bins_reco ;//=  {200 ,300,  375. , 475.,  550.,  675.,  800.,  950.,  1100.  ,1300. , 1500. , 2000};


	switch(MODEL)
	{
		case MODEL_CCQE:
			bins_reco = {200 ,300,  375. , 475.,  550.,  675.,  800.,  950.,  1100.  ,1300. , 1500.,1750 , 2000,2500};
			bins_truth = {200,250,300,350,400,450, 500,600,  800.,1000,1500,2000,2500};
			break;
		case MODEL_DELTARES:
			bins_truth =  {300,  475.,  550.,  675.,  800.,  950.,  1100.  ,1300. , 1500. ,2500};
			bins_reco =  {200 ,300,  375. , 475.,  550.,  675.,  800.,  950.,  1100.  ,1300. , 1500. , 2000};
			break;
	}

	//double max_plot_bin_truth  = 1000; //CCQE
	double max_plot_bin_truth = bins_truth.back() ;//2000;
	double max_plot_bin_reco = bins_reco.back() ;//1500;


	int N_bins_reco = bins_reco.size()-1;
	int N_bins_truth = bins_truth.size()-1;

	int N_bins_lee = 5; 


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


	sp::UnfoldAlgoBase * alg2;

	switch(ALGO)
	{
		case ALGO_SVD:
			alg2 = new sp::UnfoldAlgoSVD();
			break;
		case ALGO_DAG:
			alg2 = new sp::UnfoldAlgoDAgnostini();
			break;
	}


	//	sp::UnfoldAlgoSVD * alg2 = new sp::UnfoldAlgoSVD();


	alg2->set_verbosity((sp::msg::Level_t)0);
	alg2->Initialize(&(a.Responses().front()));


	double pot_scale = 6.46/41.10;


	/******************************* all TH1's that we need ******************************/
	TH1D * h_excess = new TH1D("h_excess","h_excess",bins_reco.size()-1,&bins_reco[0]);
	TH1D * h_data = new TH1D("h_data","h_data", bins_reco.size()-1,  &bins_reco[0] );



	//Got to tihnk more bout flows
	std::vector<double> miniobs = {232,  156,  156,   79,   81,   70,   63,   65,   62,   34,   70};
	std::vector<double> minibkg = {180.80171,108.22448,120.03353,63.887782,89.806966,67.249431,69.855878,57.014477,51.846417,38.586738,69.381391};
	double Nsignal = 0;

	std::cout<<"Loading in MC and Dirt."<<std::endl;
	std::string fname = core + bkg_file;
	std::string fname_dirt = core + dirt_file;
	std::string param = var_reco.at(0);
	std::vector<double> summed_mc(N_bins_reco,0.0);
	std::vector<double> summed_mc_err(N_bins_reco,0.0);

	std::vector<double> sigerr(N_bins_reco,0);
	//
	// Get the vector of histograms for each background type
	//

	sp::SPIO spio;
	spio.set_verbosity((sp::msg::Level_t)0);
	//spio.set_model(&model);

	std::cout<<"Generating Background"<<std::endl;
	auto th1d_v = spio.gen_background(fname,fname_dirt,param,bins_reco);


	TH1D* sum_bkg = (TH1D*)th1d_v.at(0).Clone("sum");
	sum_bkg->Reset();
	sum_bkg->Sumw2();
	TList *list = new TList;
	for(auto &v: th1d_v){
		list->Add(&v);
	}
	sum_bkg->Merge(list);
	sum_bkg->Scale(pot_scale);

	for(auto & v :th1d_v){
		for(int i=0; i< N_bins_reco;i++){
			summed_mc.at(i) += v.GetBinContent(i+1)*pot_scale;
		}

	}
	for(int i=0; i<N_bins_reco; i++){
		std::cout<<"ERRORCHECK: sumw2: "<<sum_bkg->GetBinContent(i+1)<<" +/- "<<sum_bkg->GetBinError(i+1)<<" obv: "<<summed_mc.at(i)<<" +/- "<<sqrt(summed_mc.at(i))<<std::endl;
		summed_mc_err.at(i) = sum_bkg->GetBinError(i+1);
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
	TFile * f_data = new TFile( (core+data_file).c_str());
	TTree * data = (TTree*)f_data->Get("MiniBooNE_CCQE");
	float d_Weight,d_Energy,d_CosTheta,d_RecoEnuQE;
	data->SetBranchAddress("Weight",&d_Weight);
	data->SetBranchAddress("Energy",&d_Energy);
	data->SetBranchAddress("CosTheta",&d_CosTheta);
	data->SetBranchAddress("RecoEnuQE", &d_RecoEnuQE);


	for(int i=0; i< data->GetEntries();i++){
		data->GetEntry(i);
		h_data->Fill(1000*d_RecoEnuQE, d_Weight);
	}






	TVectorD mini_signal(N_bins_reco);
	for(int i=0; i<N_bins_reco; i++){
		//This one is the old working one
		//mini_signal(i) = (miniobs.at(i)-minibkg.at(i))+alg2->r(i)*pot_scale    ;

		h_excess->SetBinContent(i+1, h_data->GetBinContent(i+1)-summed_mc.at(i));

		mini_signal(i) = (h_data->GetBinContent(i+1)-summed_mc.at(i))+alg2->r(i)*pot_scale    ;


		Nsignal += mini_signal(i);
		sigerr.at(i) =sqrt(h_data->GetBinContent(i+1)+ pow(summed_mc_err.at(i),2.0) ) ;

		if(sigerr.at(i)!=sigerr.at(i)) {
			std::cout<<"Failure, should fabs"<<std::endl;
			exit(EXIT_FAILURE);

		}



	}

	TMatrixD sigcorr(N_bins_reco, N_bins_reco);
	sigcorr.Zero();
	for(int i=0; i<N_bins_reco;i++){
		sigcorr(i,i)=h_data->GetBinContent(i+1)+pow(summed_mc_err.at(i),2.0);
		std::cout<<"MYD(i,i) "<<sigcorr(i,i)<<std::endl;
	}

	alg2->Setd(&mini_signal);
	alg2->SetD(&sigcorr);
	std::cout<<"juuust setting D"<<std::endl;
	alg2->Unfold();


	//	alg2->SetRegularization(3);
	//	alg2->TestUnfolding("CCQE_poison_unfold_test");













	TH1D ans =  alg2->GetHistU();
	TH1D sig =  alg2->GetHistD();

	TH1D sig_events = alg2->GetHistD(); //pure for chi^2 calc

	TH1D reco = alg2->GetHistR();
	TH1D truth = alg2->GetHistT();
	sig.SetLineColor(kRed); 
	ans.SetLineColor(kBlue); 
	truth.SetLineColor(kBlack);
	reco.SetLineColor(kBlack); 

	//silly SetError expects underflow
	sigerr.insert ( sigerr.begin() , 0 );	

	sig.SetMarkerStyle(21);
	sig.SetError(&sigerr[0]);
	sig.Scale(1,"width");
	ans.Scale(1,"width");
	truth.Scale(pot_scale,"width");
	reco.Scale(pot_scale,"width");

	TFile *f = new TFile("CCQE_ans.root","RECREATE");








	/************************* Plot 3: Grid of input, responce and efficiency **********************/



	TCanvas *c_responce = new TCanvas("c_response","c_response",1200,1200);
	c_responce->Divide(2,2);
	c_responce->cd(1);

	TH1D tt = truth;
	tt.SetMarkerStyle(21);
	tt.SetMarkerColor(kBlue-6);
	tt.SetLineWidth(2);
	tt.SetLineColor(kBlue-6);
	tt.GetXaxis()->SetTitle("Truth [MeV]");
	tt.GetYaxis()->SetTitle("Events/MeV");
	tt.GetXaxis()->SetTitleSize(0.04);
	tt.GetXaxis()->SetLabelSize(0.04);
	tt.GetYaxis()->SetLabelSize(0.04);
	tt.GetYaxis()->SetTitleOffset(1.45);
	tt.SetMinimum(0);
	tt.SetTitle("True Variable");
	//tt.SetMaximum(5);
	tt.GetXaxis()->SetRange(1,10);

	TLegend * legr = new TLegend(0.58,0.6,0.89,0.89);
	legr->AddEntry(&tt,"True E_{#nu} MC","lep");
	legr->SetFillStyle(0);
	//	legr->SetBorderSize(0.1);

	tt.Draw();
	//rr.Draw("same ap");
	legr->Draw();


	c_responce->cd(2);
	TH1D rr = reco;
	rr.SetMarkerStyle(21);
	rr.SetMarkerColor(kRed-7);
	rr.SetLineWidth(2);
	rr.SetLineColor(kRed-7);
	rr.GetXaxis()->SetTitle("Reco [MeV]");
	rr.GetYaxis()->SetTitle("Events/MeV");
	rr.GetXaxis()->SetTitleSize(0.04);
	rr.GetXaxis()->SetLabelSize(0.04);
	rr.GetYaxis()->SetLabelSize(0.04);
	rr.GetYaxis()->SetTitleOffset(1.45);
	rr.SetMinimum(0);
	rr.SetTitle("Reconstructed Variable");
	//rr.SetMaximum(5);
	rr.GetXaxis()->SetRange(1,10);




	TLegend * legr2 = new TLegend(0.58,0.6,0.89,0.89);
	legr2->AddEntry(&rr,"Reco E_{QE} MC","lep");
	legr2->SetFillStyle(0);
	rr.Draw();
	legr2->Draw();





	c_responce->cd(4)->SetRightMargin(0.175);
	TH2D prob_MC = alg2->GetHistA();
	gStyle->SetPalette(kInvertedDarkBodyRadiator);

	prob_MC.GetYaxis()->SetTitle("Truth E_{#nu} Bin");
	prob_MC.GetXaxis()->SetTitle("Reco E_{QE} Bin");
	prob_MC.SetTitle("Response Matrix");
	prob_MC.GetXaxis()->SetTitleOffset(1.1);
	prob_MC.GetYaxis()->SetTitleOffset(1.1);

	prob_MC.Draw("colz");
	prob_MC.GetXaxis()->SetRange(1,10);
	prob_MC.GetYaxis()->SetRange(1,10);

	c_responce->cd(3);
	TH1D eff = alg2->GetHistEff();

	std::vector<double> errEff(alg2->n_t+1,0.0);
	for(int b=1; b<errEff.size();b++){
		errEff.at(b)=sqrt(alg2->Ep(b-1,b-1));  ;
	}
	eff.SetError(&errEff[0]);

	eff.GetYaxis()->SetTitle("Efficiency");
	eff.GetXaxis()->SetTitle("Truth Bin");
	eff.Draw("E1");
	eff.SetMinimum(0);
	eff.SetMaximum(0.4);


	c_responce->SaveAs( (base_name+"_eff_response.pdf").c_str(),"pdf");



	//return 0;


	/************************ Plot 0: The data and MC in chosen binning **********************/
	TCanvas *c_obsv = new TCanvas("c_obsv","c_obsv",1200,800);	
	c_obsv->cd();

	h_data->Scale(1,"width");
	h_data->SetLineColor(kRed);
	h_data->SetMarkerStyle(3);
	h_data->SetMarkerColor(kRed);	

	sum_bkg->GetYaxis()->SetTitle("Events/MeV");
	sum_bkg->GetXaxis()->SetTitle("Reco E^{QE}_{#nu} [MeV]");

	sum_bkg->SetTitle("Observed data and MC bkg");
	sum_bkg->Scale(1,"width");	
	sum_bkg->SetFillColor(kGray);
	sum_bkg->SetLineColor(kBlack);


	sum_bkg->Draw("hist");
	h_data->Draw("same E1");
	sum_bkg->SetMaximum(2.5);

	c_obsv->SaveAs( (base_name+ "_obsv.pdf").c_str(),"pdf");



	/************************* Plot 1: Just the excess, more of a check **********************/
	TCanvas *c_excess = new TCanvas("c_excess","c_excess",1200,800);	
	c_excess->cd();
	h_excess->SetLineColor(kGreen+3);
	h_excess->Scale(1,"width");

	h_excess->SetTitle("MiniBooNE LEE");
	h_excess->Draw("E1");
	h_excess->GetYaxis()->SetTitle("Events/MeV");
	h_excess->GetXaxis()->SetTitle("Reco E^{QE}_{#nu} [MeV]");

	c_excess->SaveAs((base_name +"_excess.pdf").c_str(),"pdf");






	/************************* Plot 2: Signal over chosen background **********************/

	TCanvas *c_signal = new TCanvas("c_signal","c_signal",1200,800);

	sig.SetTitle("Reconstructed E^{QE}" );
	sig.GetYaxis()->SetTitle("Events/MeV");
	sig.GetXaxis()->SetTitle("Reco E^{QE}_{#nu} [MeV]");
	sig.SetLineColor(kRed-7);
	sig.SetMarkerColor(kRed-7);
	sig.SetMarkerStyle(21);
	sig.SetLineWidth(2);
	sig.Draw("e1");
	sig.GetXaxis()->SetRangeUser(bins_reco.front(),max_plot_bin_reco);

	sig.SetMinimum(0);

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


	c_signal->Write();
	c_signal->SaveAs((base_name+"_signal.pdf").c_str(),"pdf");


	/************************* Plot 3s**********************/

	TCanvas *c_vary =  new TCanvas("c1","c1",1200,1200);
	TCanvas *c_unfold =  new TCanvas("c_unfold","c_unfold",1200,1200);
	TCanvas *c_refold =  new TCanvas("c_refold","c_refold",3000,3000);
	TCanvas *c_bias = new TCanvas("c_bias","c_bias",3000,3000);
	TCanvas *c_frac_cov =  new TCanvas("c_frac_cov","c_frac_cov",3000,3000);
	TCanvas *c_full_cov =  new TCanvas("c_full_cov","c_full_cov",3000,3000);
	TCanvas *c_corr =  new TCanvas("c_vary","c_vary",3000,3000);

	c_vary->Divide(3,3);
	c_refold->Divide(3,3);
	c_bias->Divide(3,3);
	c_corr->Divide(3,3);
	c_frac_cov->Divide(3,3);
	c_full_cov->Divide(3,3);

	double N_random_refold = 10.0;
	int N_mc_bias = 10;
	int N_chols = 10;



	std::vector<int> kreg = {1,2,3,4,5,6,7,8};
	std::vector<int> cols = {kBlue-7,kGreen-6,kRed-7, kOrange-3, kMagenta-3, kGreen+3, kBlue-7,kGreen-6,kRed-7};
	std::vector<TH1D> us(kreg.size());
	std::vector<TMatrixT<double>> UU(kreg.size());
	std::vector<TH2D> U_full_cov(kreg.size());
	std::vector<TH2D> U_frac_cov(kreg.size());
	std::vector<TH2D> U_corr(kreg.size());
	

	std::vector<TH1D> uR(kreg.size());
	std::vector<TH1D> uBias(kreg.size());
	std::vector<TLegend*> leg(kreg.size());
	std::vector<TH1D> us_stat(kreg.size());

	std::vector<std::vector<TH1D>> u_chols( kreg.size(), std::vector<TH1D>(N_chols)   );
	std::vector<std::vector<TH1D>> u_chols_ratio( kreg.size(), std::vector<TH1D>(N_chols)   );

	//std::vector<int> kreg = {1,2,5,8,100};

	//Some basic Lcurve and Chi^2 analysis on refolded spectra
	TCanvas *c_chicurve = new TCanvas("lc","lc,",1200,1200);	
	std::vector<double> refold_avg_chi;
	std::vector<double> refold_k;

	std::vector<double> bias_avg;




	/**********************************************************
	 *
	 *		Some rudimentary Reg Choice
	 *
	 **********************************************************/

	int best_reg=-1;
	std::vector<int> v_best_reg;







	

	std::vector<std::vector<int>> base_cols;
	std::vector<int> tri_cols = {kRed, kMagenta,kBlue,kCyan,kGreen,kYellow};
	std::vector<int> rec_cols = {kOrange, kPink,kViolet,kAzure,kTeal,kSpring};	
	

	for(int k=0; k<kreg.size(); k++){
		std::string nam = "Reg " +std::to_string(kreg.at(k)) ;

		c_vary->cd(k+1)->SetLeftMargin(0.1);
		alg2->SetRegularization(kreg.at(k));
		alg2->Unfold();

		//Section to re-fold and calc numbers of events.
		TVectorD refolded = alg2->A*alg2->u;
		double Nrefold = 0;
		for(int i=0; i<alg2->n_r; i++){
			Nrefold+=refolded(i);
		}
		std::cout<<std::setprecision(10)<<"REFOLD || On k="<<kreg.at(k)<<" and refolded-number-of-events is "<<Nrefold<<" . Data is "<<Nsignal<<std::endl;



		std::vector<double> errA(alg2->n_t+1,0.0);
		std::vector<double> errD(alg2->n_t+1,0.0);
		std::vector<double> err(alg2->n_t+1,0.0);
		std::vector<double> errS(alg2->n_t+1,0.0);

		for(int i=1;i<err.size(); i++){
			errD.at(i) = sqrt(fabs(alg2->U(i-1,i-1)));
			errA.at(i) = sqrt(fabs(alg2->UA(i-1,i-1)));
			errS.at(i) = sqrt(fabs(alg2->u(i-1)));
			err.at(i) = sqrt(fabs(alg2->UA(i-1,i-1)+alg2->U(i-1,i-1)));
			std::cout<<"u= "<<alg2->u(i-1)<<" +/-(stat) "<<errS.at(i-1)<<" +/-(D) "<<errD.at(i-1)<<" +/-(A) "<<errA.at(i-1)<<" +/-(AD) "<<err.at(i-1)<<" bias "<<alg2->b(i-1)<<" +/- "<<sqrt(fabs(alg2->B(i-1,i-1) ) )<<std::endl;
		}


		us.at(k) = alg2->GetHistU();
		UU.at(k) = alg2->U;
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
		us.at(k).GetXaxis()->SetRangeUser(bins_truth.front(),max_plot_bin_truth);
		us.at(k).SetMaximum(10);//has to be after!

		us_stat.at(k) = alg2->GetHistU();
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


		// Cholosky Decomp
		for(int i=0; i<N_chols; i++){
				std::cout<<"On Chol Deco #: "<<i<<std::endl;
				u_chols.at(k).at(i) = alg2->SampleCovarianceU();
				u_chols.at(k).at(i).SetLineWidth(1);
				int rnd_col = rec_cols.at(rangen->Integer(rec_cols.size()));
				rnd_col = rnd_col+ rangen->Integer(20)-10;
				u_chols.at(k).at(i).SetLineColor(rnd_col);
			
				u_chols_ratio.at(k).at(i) = u_chols.at(k).at(i);
				u_chols_ratio.at(k).at(i).Scale(1,"width");
				u_chols_ratio.at(k).at(i).Divide(&truth); 
		}



		c_full_cov->cd(k+1);//->SetLogz();
		U_full_cov.at(k) = alg2->GetCovU();
		U_full_cov.at(k).SetTitle(  nam.c_str() );
		U_full_cov.at(k).GetYaxis()->SetTitle("True E_{#nu} Bin a");
		U_full_cov.at(k).GetXaxis()->SetTitle("True E_{#nu} Bin b");
		U_full_cov.at(k).Draw("colz");

		//Fractional 
		c_frac_cov->cd(k+1);//->SetLogz();
		U_frac_cov.at(k) = U_full_cov.at(k);
			for(int a=0; a< alg2->n_t; a++){
				for(int b=0; b< alg2->n_t; b++){
					U_frac_cov.at(k).SetBinContent(a+1,b+1,  U_full_cov.at(k).GetBinContent(a+1,b+1)/(alg2->u(a)* alg2->u(b)    ));	
				}
			}

		U_frac_cov.at(k).SetTitle(  nam.c_str() );
		U_frac_cov.at(k).GetYaxis()->SetTitle("True E_{#nu} Bin a");
		U_frac_cov.at(k).GetXaxis()->SetTitle("True E_{#nu} Bin b");
		U_frac_cov.at(k).Draw("colz");

		//Correlation 
		c_corr->cd(k+1);//->SetLogz();
		U_corr.at(k) = U_full_cov.at(k);
			for(int a=0; a< alg2->n_t; a++){
				for(int b=0; b< alg2->n_t; b++){
					U_corr.at(k).SetBinContent(a+1,b+1,  U_full_cov.at(k).GetBinContent(a+1,b+1)/(sqrt(U_full_cov.at(k).GetBinContent(a+1,a+1)*U_full_cov.at(k).GetBinContent(b+1,b+1)   ) ));	
				}
			}
		U_corr.at(k).SetTitle(  nam.c_str() );
		U_corr.at(k).GetYaxis()->SetTitle("True E_{#nu} Bin a");
		U_corr.at(k).GetXaxis()->SetTitle("True E_{#nu} Bin b");
		U_corr.at(k).Draw("colz");



		c_refold->cd(k+1);
		uR.at(k) = alg2->GetHistRefold();
		uR.at(k).SetTitle(nam.c_str());
		double ch2 =0;// uR.at(k).Chi2Test(&sig_events,"CHI2,WW");
		double ch2_lee = 0;
		double ch2_max = 0;

		for(int i=1; i < N_bins_reco; i++){
			double tmp = pow(uR.at(k).GetBinContent(i)-sig_events.GetBinContent(i),2.0)/sigcorr(i-1,i-1);
			ch2 += tmp;
			if(i<N_bins_lee) ch2_lee +=tmp;
			if(tmp > ch2_max && i <N_bins_reco-1) ch2_max = tmp;
			std::cout<<"CHI k: "<<k<<" bin: "<<i<<" "<<uR.at(k).GetBinContent(i)<<" "<<sig_events.GetBinContent(i)<<" tmp: "<<tmp<<" chi: "<<ch2<<std::endl;		
		}



		double Nrefold_2 = uR.at(k).GetSumOfWeights();
		uR.at(k).Scale(1,"width");
		uR.at(k).SetMarkerStyle(29);
		uR.at(k).SetMarkerColor(kBlack);
		uR.at(k).SetLineColor(kBlack);
		uR.at(k).SetMarkerSize(2);


		//RandomRefolding && MC bias calculation
		std::vector<double> bias(alg2->n_t,0);
		std::vector<double> avgu(alg2->n_t,0);

		double avg_refold_chi = 0;
		for(int j=0;j<(int)N_random_refold;j++){
			TH1D tmp = alg2->GetHistRandRefold(rangen);

			double temp_chi =0;// tmp.Chi2Test(&sig_events,"CHI2,WW");

			for(int i=1; i < N_bins_reco; i++){
				double thischi = pow(tmp.GetBinContent(i)-sig_events.GetBinContent(i),2.0)/sigcorr(i,i);
				temp_chi += thischi;
				//std::cout<<"CHI k: "<<k<<" bin: "<<i<<" "<<uR.at(k).GetBinContent(i)<<" "<<sig_events.GetBinContent(i)<<" tmp: "<<tmp<<" chi: "<<ch2<<std::endl;		
			}




			avg_refold_chi +=temp_chi;
			std::cout<<"RandRefold: "<<k<<" chi^2: "<<temp_chi<<std::endl;


			//and MC bias calc
			TVectorD tmpV(alg2->n_r);
			for(int i=0; i<alg2->n_r; i++){
				tmpV(i)=tmp.GetBinContent(i+1);
			}

			/*
			   auto tmp_alg = *alg2;
			   tmp_alg.Setd(&tmpV);
			   tmp_alg.Unfold();
			   for(int a=0; a<alg2->n_t; a++){
			   avgu.at(a) += tmp_alg.u(a)/N_random_refold;
			   }
			 */ //currently broke

		}

		for(int a=0;a<alg2->n_t; a++){
			bias.at(a) = avgu.at(a)-alg2->u(a);
			std::cout<<"ThisBias bin:"<<a<<" "<<bias.at(a)<<" "<<alg2->b(a)<<std::endl;
		}

		double avg_bias = 0;//Just average the bins		
		for(auto b: bias){
			avg_bias += b/((double)N_bins_truth);
		}

		bias_avg.push_back(pow(fabs(avg_bias),2.0));
		std::cout<<"ThisBias Avg: "<<avg_bias<<std::endl;


		bool bias_zero = true;
		double bias_insig = 100;
		//bias consistent with zero
		for(int a=0; a<alg2->n_t; a++){
			std::cout<<"BINOUT: "<<a<<" B "<<alg2->B(a,a)<<" sqrt(B) "<<sqrt(alg2->B(a,a))<<std::endl;

			double this_bias_insig = alg2->b(a)/sqrt(alg2->B(a,a));	
			double this_bias_zero = fabs(alg2->b(a)) - sqrt(alg2->B(a,a)) ;
			if(this_bias_insig < bias_insig) bias_insig = this_bias_insig;
			if(fabs(alg2->b(a)) >= sqrt(alg2->B(a,a))    ) bias_zero = false;
		}


		//Step 1. Get random (but slightly different "true" spectra)


		/* This is multi-truth biasing
		   TVectorD pois_u = alg2->u;
		   std::vector<double> bias(alg2->n_t,0.0);
		   for(int i=0; i< N_mc_bias; i++){

		   for(int a=0; a<alg2->n_t; a++){
		   pois_u(a) = rangen->Poisson( alg2->u(a) ); 
		   }


		   TVectorD re = alg2->A*pois_u;
		   auto tmp_alg2->= alg2->
		   tmp_alg2->Setd(&re);
		   tmp_alg2->Unfold();

		   for(int a=0; a<alg2->n_t; a++){
		   double thisbias = tmp_alg2->u(a)-pois_u(a);
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

		std::string s_ch2 = "#chi^{2}/ndof : " + std::to_string(ch2)+"/"+std::to_string(N_bins_reco) + " || " + std::to_string(ch2_lee)+"/"+std::to_string(N_bins_lee) +" || " + std::to_string(ch2_max)+"/"+std::to_string(1);
		;
		TLegend * legR = new TLegend(0.3,0.7,0.89,0.89);
		legR->SetHeader(s_ch2.c_str() );
		legR->AddEntry(&sig, ("MiniBooNE Excess: " + std::to_string(Nsignal)).c_str() ,"lp");
		legR->AddEntry(&uR.at(k), ("Refolded Excess: "+std::to_string(Nrefold_2)).c_str(),"l");
		uR.at(k).Draw("hist");
		sig.Draw("same");
		legR->Draw();
		uR.at(k).GetXaxis()->SetRangeUser(bins_reco.front(),max_plot_bin_reco);
		uR.at(k).SetMaximum(1.2);
		uR.at(k).SetMinimum(0.0);




		//Bias plot filling
		c_bias->cd(k+1);
		uBias.at(k) = alg2->GetHistBias();	
		uBias.at(k).SetMarkerStyle(2);
		uBias.at(k).SetMarkerSize(2);
		uBias.at(k).SetTitle(  nam.c_str() );
		uBias.at(k).GetYaxis()->SetTitle("Bias");
		uBias.at(k).GetXaxis()->SetTitle("True E_{#nu} [MeV]");
		uBias.at(k).SetLineColor(kBlack);
		uBias.at(k).SetLineWidth(2);
		uBias.at(k).SetFillColor(cols.at(k));
		//uBias.at(k).Scale(1,"width"); // should it be bias /MeV? doubt it really
		uBias.at(k).Draw("E2");
		uBias.at(k).GetXaxis()->SetRangeUser(bins_truth.front(),max_plot_bin_truth);
		/* At this point the standard
		   deviations of the biases are approximately equal to the biases themselves, and
		   therefore any further bias reduction would introduce as much error as it removes.
		 */


		std::cout<<"REGCHECK: "<<kreg.at(k)<<" X^2/ndof < 1: "<<ch2/N_bins_reco<<"  bias_insig < 1: "<<fabs(bias_insig)<<" is bias zero: "<<bias_zero<<std::endl;
		if(ch2/N_bins_reco < 1 && fabs(bias_insig) < 1 && bias_zero){
			v_best_reg.push_back(k);
		}


	}




	if(v_best_reg.size()==0){
		best_reg = 1;
	}else{
		best_reg = v_best_reg.front();
	}

	std::cout<<"REGCHECK: chosen :"<<kreg.at(best_reg)<<std::endl;



	c_unfold->cd();
	us_stat.at(best_reg).SetFillColor(kGreen+2);
	us_stat.at(best_reg).SetLineColor(kGreen+2);
	us_stat.at(best_reg).SetMarkerColor(kBlack);
	us_stat.at(best_reg).SetMarkerStyle(29);
	us_stat.at(best_reg).SetMarkerSize(2);
	us.at(best_reg).Draw("E2");
	us_stat.at(best_reg).Draw("same E2");
	us_stat.at(best_reg).Draw("same P");
	truth.Draw("same hist");
	leg.at(best_reg)->Draw();
	us.at(best_reg).GetXaxis()->SetRangeUser(bins_truth.front(), max_plot_bin_truth);




	/***********************************************************************
				Draw Colerated pairs!
	 ***********************************************************************/
	TCanvas *c_chol = new TCanvas("chol","chol",1200,800);
	c_chol->cd();



	int use_chol = 2;//best_reg;
	us.at(use_chol).Draw("E1");
	truth.Draw("same hist");
	us.at(use_chol).GetXaxis()->SetRangeUser(bins_truth.front(), max_plot_bin_truth);

	for(int i=0; i< N_chols; i++){
		u_chols.at(use_chol).at(i).Scale(1,"width");
		u_chols.at(use_chol).at(i).Draw("same hist");
	}
	

	c_chol->SaveAs((base_name+"_chol.pdf").c_str(),"pdf");



	/***********************************************************************
				Draw Colerated pairs! Ratio
	 ***********************************************************************/
	TCanvas *c_chols_ratio = new TCanvas("chols_ratio","chols_ratio",1200,800);
	c_chols_ratio->cd();



	int use_chols_ratio = 2;//best_reg;
		u_chols_ratio.at(use_chols_ratio).at(0).Draw("hist");
		u_chols_ratio.at(use_chols_ratio).at(0).GetXaxis()->SetRangeUser(bins_truth.front(), max_plot_bin_truth);
		u_chols_ratio.at(use_chols_ratio).at(0).SetMaximum(8);
		u_chols_ratio.at(use_chols_ratio).at(0).SetMinimum(0);
	
	for(int i=1; i< N_chols; i++){
		u_chols_ratio.at(use_chols_ratio).at(i).Draw("same hist");
	}
	
	TLine *line_chol_ratio = new TLine(bins_truth.front(),1,max_plot_bin_truth,1);
	line_chol_ratio->SetLineColor(kBlack);
	line_chol_ratio->SetLineStyle(2);
	line_chol_ratio->Draw();



	c_chols_ratio->SaveAs((base_name+"_chols_ratio.pdf").c_str(),"pdf");











	alg2->TestRegularization("CCQE_lcurves", 1,6,5);
	//alg2->SetDirectRegularization(1000);
	//alg2->TestRegularization("CCQE_lcurves_direct", -8,8,200);  //-20 13

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
	c_chicurve->SaveAs((base_name+"_chicurve.pdf").c_str(),"pdf");

	TCanvas * c_mcbias = new TCanvas("cmcbias","cmcbias",1200,1200);
	c_mcbias->cd();
	c_mcbias->SetLogy();

	TGraph *g_bias = new TGraph(kreg.size(), &refold_k[0], &bias_avg[0] );
	g_bias->SetTitle("Average Bias MC ");
	g_bias->GetXaxis()->SetTitle("Regularisation Parameter k");
	g_bias->GetYaxis()->SetTitle("<bias> ");
	g_bias->SetMarkerSize(2);
	g_bias->SetMarkerStyle(2);
	g_bias->Draw("ACP");
	//TLine l_bias(refold_k.front(),1,refold_k.back(),1);
	//l_bias.Draw();
	c_mcbias->Write();
	c_mcbias->SaveAs((base_name+"_bias.pdf").c_str(),"pdf");



	c_vary->Write();
	c_corr->Write();
	c_full_cov->Write();
	c_frac_cov->Write();
	c_bias->Write();

	c_unfold->SaveAs((base_name+"_unfolded.pdf").c_str(),"pdf");
	c_vary->SaveAs((base_name+"_reg_vary.pdf").c_str(),"pdf");
	c_corr->SaveAs((base_name+"_reg_corr.pdf").c_str(),"pdf");
	c_frac_cov->SaveAs((base_name+"_reg_frac_cov.pdf").c_str(),"pdf");
	c_full_cov->SaveAs((base_name+"_reg_full_cov.pdf").c_str(),"pdf");
	c_refold->SaveAs((base_name+"_reg_refold.pdf").c_str(),"pdf");
	c_bias->SaveAs((base_name+"_reg_bias.pdf").c_str(),"pdf");


	/*	Not finished yet
		TCanvas *c_sample =  new TCanvas();
		c_sample->cd();
		TVectorD result(n_t); 
		std::vector<TH1D*> sample_list;
		us.at(i)
		for(int i=0; i< 4; i++){
		sample_list.at(i) = alg2->SampleCovarianceU(&result);
		sample_list.at(i)->Draw();		
		}
	 */

	TCanvas *cr = new TCanvas();
	TH1D ratio = us.at(best_reg);	
	ratio.Divide(&truth);
	ratio.SetFillColor(kBlue-7);
	ratio.GetYaxis()->SetTitle("Ratio to MiniBooNE MC Central Value");
	ratio.SetMarkerStyle(5);
	ratio.SetMarkerSize(1);
	ratio.Draw("e2");
	ratio.SetMaximum(8);
	ratio.GetXaxis()->SetRangeUser(bins_truth.front(),max_plot_bin_truth);
	ratio.GetYaxis()->SetTitleOffset(1.1);

	TLegend * leg2 = new TLegend(0.58,0.6,0.89,0.89);
	leg2->AddEntry(&ratio,"Intrinsic #nu_{e} CCQE Model","lf");
	leg2->SetFillStyle(0);
	leg2->SetBorderSize(0.0);
	leg2->Draw();

	TLine *line = new TLine(bins_truth.front(),1,max_plot_bin_truth,1);
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
	ratio_scaled->GetXaxis()->SetRangeUser(bins_truth.front(),max_plot_bin_truth);

	leg2->Draw();
	line->Draw();
	cr_scaled->SaveAs((base_name+"_model_ratio_scaleo.pdf").c_str(),"pdf");




	for(int i=0; i<=ratio.GetNbinsX()+2; i++){
		std::cout<<"Ratio "<<ratio.GetBinContent(i)<<" "<<ratio_scaled->GetBinContent(i)<<std::endl;
	}
	cr->Write();
	cr->SaveAs((base_name+"_model_ratio.pdf").c_str(),"pdf");

	f->Close();

	TFile* f_graph_out =  new TFile((base_name+"_final_tgraph.root").c_str(),"RECREATE");
	f_graph_out->cd();

	ratio.Write();
	TGraph * graph_out = new TGraph(bincenter.size(), &bincenter[0], &binval[0]  );
	graph_out->Write();
	f_graph_out->Close();





	return 0;
}
