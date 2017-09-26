#ifndef __UNFOLDALGOBASE_CXX_
#define __UNFOLDALGOBASE_CXX__

#include "UnfoldAlgoBase.h"



namespace sp {

	/**********************************************************************
	 *
	 *		Initializer!	
	 *
	 * *******************************************************************/


	void UnfoldAlgoBase::Initialize(const Response * response_in){


		SP_DEBUG()<<"Starting on response: "<<response_in->_name<<std::endl;

		n_t = response_in->_true_param->_hist.GetNbinsX();
		n_r = response_in->_reco_param->_hist.GetNbinsX();

		SP_DEBUG() << "n_t : " << n_t << " & n_r: " << n_r << std::endl;
		hist_u = (TH1D*)response_in->_true_param->_hist.Clone("unfold");
		hist_r = (TH1D*)response_in->_reco_param->_hist.Clone("reco");

		u.ResizeTo(n_t);
		b.ResizeTo(n_t);
		ep.ResizeTo(n_t);

		d.ResizeTo(n_r);

		D.ResizeTo(n_r,n_r); 
		D.Zero();

		U.ResizeTo(n_t,n_t); 
		U.Zero();

		UD.ResizeTo(n_t,n_t); 
		UD.Zero();

		UA.ResizeTo(n_t,n_t); 
		UA.Zero();

		SP_DEBUG()<<"Initilizing truth MC."<<std::endl;
		t.ResizeTo(n_t);
		t = hist2TVec(&response_in->_true_param->_hist);
		SP_DEBUG()<<"Initilized truth MC, dimension: "<<t.GetNrows()<<std::endl;

		for(int i =0; i<n_t; i++){
			if(t(i)==0){
				SP_WARNING()<<"Truth has a 0 @ position: "<<i<<" of "<<t.GetNrows()<<std::endl;
			}else{
				SP_DEBUG()<<"Truth has "<<t(i)<<" @ position: "<<i<<" of "<<t.GetNrows()<<std::endl;
			}
		}


		SP_DEBUG()<<"Initilizing reco MC."<<std::endl;
		r.ResizeTo(n_r);
		r = hist2TVec(&response_in->_reco_param->_hist);
		SP_DEBUG()<<"Initilized reco MC, dimension: "<<r.GetNrows()<<std::endl;

		if(response_in->_true_param->_hist.GetBinContent(0)!=0){
			SP_WARNING()<<"Truth MC has filled underflow bin of value: "<< response_in->_true_param->_hist.GetBinContent(0)<<std::endl;
		}
		if(response_in->_true_param->_hist.GetBinContent(n_t+1)!=0){
			SP_WARNING()<<"Truth MC has filled overflow bin of value: "<< response_in->_true_param->_hist.GetBinContent(n_t+1)<<std::endl;
		}
		if(response_in->_reco_param->_hist.GetBinContent(0)!=0){
			SP_WARNING()<<"Reco MC has filled underflow bin of value: "<< response_in->_reco_param->_hist.GetBinContent(0)<<std::endl;
		}
		if(response_in->_reco_param->_hist.GetBinContent(n_r+1)!=0){
			SP_WARNING()<<"Reco MC has filled overflow bin of value: "<< response_in->_reco_param->_hist.GetBinContent(n_r+1)<<std::endl;
		}







		SP_DEBUG()<<"Initilizing Number of Events."<<std::endl;
		N.ResizeTo(n_r,n_t);
		SP_DEBUG()<<"input reponse has dimensions: "<<response_in->_response.GetNrows()<<" "<<response_in->_response.GetNcols()<<std::endl;
		N = response_in->_response;
		SP_DEBUG()<<"|| Initilized Number of events N, dimension: "<<N.GetNrows()<<" "<<N.GetNcols()<<std::endl;

		SP_DEBUG()<<"Initilizing Response."<<std::endl;
		A.ResizeTo(n_r,n_t);

		for(int i=0;i<n_r; i++){
			for(int a=0; a<n_t; a++){
				A(i,a) = N(i,a) / t(a);
				if(A(i,a)>1){
					SP_CRITICAL()<<"A has entries > 1."<<std::endl;
					SP_CRITICAL()<<i<<" "<<a<<" A:"<<A(i,a)
						<<" N:"<<N(i,a)<<" t:"<<t(a)
						<<" Nin:"<<response_in->_response(i,a)
						<<" tin: "<<response_in->_true_param->_hist.GetBinContent(a)
						<<" r: "<<r(i)<<" rin: "<<response_in->_reco_param->_hist.GetBinContent(i)<<std::endl;

					throw sperr();
				}
			}
		}
		SP_DEBUG()<<"Initilized Response A, dimension: "<<A.GetNrows()<<" "<<A.GetNcols()<<std::endl;

		SP_DEBUG()<<"Calculating a global efficiency vector."<<std::endl;
		ep.Zero();
		for(int b=0; b<n_t; b++){
			for(int l=0; l<n_r; l++){
				ep(b) += A(l,b) ;
			}
			SP_DEBUG()<<"Efficiency is: "<<ep(b)<<" at bin "<<b<<" "<<std::endl;
			if(ep(b)==0){
				SP_CRITICAL()<<"Efficiency is: "<<ep(b)<<" at bin #"<<b<<" with t: "<<t(b)<<std::endl;
				for(int l=0; l<n_r; l++){
					SP_CRITICAL()<<"bin: "<<l<<" N: "<<N(l,b)<<" r: "<<r(l)<<" A: "<<A(l,b)<<" t: "<<t(b)<<std::endl;
				}
				SP_CRITICAL()<<"Efficiency cannot be zero! Suggestion to merge last (first) true bins to solve."<<std::endl;
				throw sperr();
			}
		}



		covA.resize(n_r, std::vector<std::vector<double>>(n_t, std::vector<double>(n_r,0)));//thre indicies, as we have ignored inter truth correlations i think`

		SP_DEBUG()<<"Calculating covariance on Response. MUST Double CHECK formula."<<std::endl;
		//This is necessary to propagate the errors on A through..
		for(int r=0;r<n_r; r++){
			for(int a=0;a<n_t; a++){
				for(int s=0;s<n_r; s++){
					if(r==s){
						covA.at(r).at(a).at(s) = 1.0/t(a) *A(r,a)*(1-A(r,a)); 
					}else{
						covA.at(r).at(a).at(s) = - 1.0/t(a)*A(r,a)*A(s,a); 
					}
				}
			}
		}


		Ep.ResizeTo(n_t,n_t);
		Ep.Zero();
		SP_DEBUG()<<"Calculating covariance on global efficieincy."<<std::endl;
		for(int a=0;a<n_t; a++){
			for(int i=0; i< n_r; i++){
				for(int j=0; j< n_r; j++){
					Ep(a,a) += covA.at(i).at(a).at(j);
				}
			}
		}






		this->_Initialize_();
	};

	/**********************************************************************
	 *
	 *			Unfolding functions, not sure these work
	 *
	 * *******************************************************************/

	void UnfoldAlgoBase::Unfold(const TVectorD * data, const TMatrixD * Din){
		D.ResizeTo(n_r,n_r);
		D = *Din;

		d.ResizeTo(n_r);
		d = *data;

		this->Unfold();
	}

	void UnfoldAlgoBase::Unfold(const TVectorD * data){
		d.Zero();
		d = *data;

		D.ResizeTo(n_r,n_r);
		D.Zero();

		for(int i=0; i<n_r; i++) D(i,i)=d(i);

		this->Unfold();
	}


	/**********************************************************************
	 *
	 *			Setting functions
	 *
	 * *******************************************************************/


	void UnfoldAlgoBase::Setd(TVectorD * din){
		if( din->GetNrows() != n_r){
			SP_CRITICAL()<<"Passed a D_in that is size: "<<din->GetNrows()<<" but n_r is "<<n_r<<std::endl;
			throw sperr();
		}

		d.ResizeTo(n_r);
		for(int i=0; i<n_r;i++) 
			d(i) = (*din)(i);

		if (this->logger().level() == msg::kDEBUG) {
			SP_DEBUG() << "Set d ==>" << std::endl;
			d.Print();
		}

	}

	void UnfoldAlgoBase::SetD(TMatrixD * din){
		if( din->GetNrows() != n_r){
			SP_CRITICAL()<<"Passed a D_in that is size: "<<din->GetNrows()<<" but n_r is actually"<<n_r<<std::endl;
			throw sperr();
		}

		D.ResizeTo(n_r,n_r);
		for(int i=0; i<n_r;i++){
			for(int j=0; j<n_r;j++){
				D(i,j) = (*din)(i,j);
			}
		}

		if (this->logger().level() == msg::kDEBUG) {
			SP_DEBUG() << "Set D ==>" << std::endl;
			D.Print();
		}

	}


	void UnfoldAlgoBase::SetRegularization(double reg){
		SP_DEBUG()<<"Setting regularization parameter to "<<reg<<std::endl;
		regularization = reg;
	}


	void UnfoldAlgoBase::SetSeed(int inseed){
		SP_DEBUG() << "Setting seed to " << inseed << std::endl;
		seed = inseed;
	}


	/**********************************************************************
	 *
	 *			Getting functions
	 *
	 * *******************************************************************/


	double UnfoldAlgoBase::GetRegularization(){
		return regularization;
	}


	TH2D UnfoldAlgoBase::GetHistA(){
		TH2D tmp("tmp","tmp",n_r,0,n_r-1,n_t, 0, n_t-1 );
		for(int i=0; i<n_r; i++){
			for(int a=0; a<n_t; a++){
				tmp.SetBinContent(i+1 ,a+1, A(i,a));
			}
		}

		return tmp;
	}



	TH1D UnfoldAlgoBase::GetHistU(){

		for(int i=0; i<n_t; i++)
			hist_u->SetBinContent(i+1, u(i));

		return *hist_u;
	}

	TH1D UnfoldAlgoBase::GetHistBias(){
		std::vector<double> err;
		//need to push one back as overflow
		err.push_back(0);
		for(int a=0; a<n_t; a++){
			err.push_back( sqrt(B(a,a)));
		}

		for(int a=0; a<n_t; a++){
			hist_u->SetBinContent(a+1, b(a));
		}
		hist_u->SetError(&err[0]);

		return *hist_u;
	}



	TH2D UnfoldAlgoBase::GetCovU(){
		TH2D tmp("","",n_t,0,n_t,n_t,0,n_t);
		for(int i=0; i<n_t; i++){
			for(int j=0; j<n_t; j++){
				tmp.SetBinContent(i+1,j+1 , U(i,j));
			}
		}

		return tmp;
	}

	TH1D UnfoldAlgoBase::GetErrU(){

		for(int i=0; i<n_t; i++){
			hist_u->SetBinContent(i+1 , std::fabs(std::sqrt(U(i,i))));
		}

		return *hist_u;
	}


	TH1D UnfoldAlgoBase::GetHistT(){

		for(int i=0; i<n_t; i++){
			hist_u->SetBinContent(i+1,t(i));
		}

		return *hist_u;
	}

	TH1D UnfoldAlgoBase::GetHistEff(){

		for(int i=0; i<n_t; i++){
			hist_u->SetBinContent(i+1,ep(i));
		}

		return *hist_u;
	}

	TH2D UnfoldAlgoBase::GetCovEff(){
		TH2D tmp("","",n_t,0,n_t,n_t,0,n_t);
		for(int i=0; i<n_t; i++){
			for(int j=0; j<n_t; j++){
				tmp.SetBinContent(i+1,j+1 , Ep(i,j));
			}
		}

		return tmp;
	}


	TH1D UnfoldAlgoBase::GetHistR(){

		for(int i=0; i<n_r; i++){
			hist_r->SetBinContent(i+1 , r(i));
		}

		return *hist_r;
	}

	TH1D UnfoldAlgoBase::GetHistRefold(){
		TVectorD re = A*u;
		for(int i=0; i<n_r; i++){
			hist_r->SetBinContent(i+1 , re(i));
		}

		return *hist_r;
	}
	TH1D UnfoldAlgoBase::GetHistRandRefold(TRandom3 *rangen){
		//After unfolding, perform a poissonian draw from the unfolded spectra, then refold it.	
		//TRandom3 * rangen = new TRandom3(seed);
		TVectorD pois_u = u;

		for(int a=0; a<n_t; a++){
			pois_u(a) = rangen->Poisson( u(a) ); 
		}


		TVectorD re = A*pois_u;


		for(int i=0; i<n_r; i++){

			hist_r->SetBinContent(i+1 , re(i) );
		}

		return *hist_r;
	}


	TH1D UnfoldAlgoBase::GetHistD(){

		for(int i=0; i<n_r; i++){
			hist_r->SetBinContent(i+1 , d(i));
		}

		return *hist_r;
	}



	/**********************************************************************
	 *
	 *			Other general and Testing functions
	 *
	 * *******************************************************************/


	void UnfoldAlgoBase::GenPoissonNoise(){
		SP_DEBUG() << "start" << std::endl;
		TRandom3 * rangen = new TRandom3(seed);

		SP_DEBUG()<<"Setting data as noisey reco MC with seed "<<rangen->GetSeed()<<std::endl;

		d.ResizeTo(n_r);

		for(int i=0; i<n_r; i++)
			d(i) = rangen->Poisson( r(i) );

		if(this->logger().level() == msg::kDEBUG) {
			SP_DEBUG() << "Gen Poisson d ==> " << std::endl;
			d.Print();
		}

		SP_DEBUG() << "end" << std::endl;
	}



	void UnfoldAlgoBase::TestUnfolding(std::string filename){
		SP_DEBUG() << "start" << std::endl;

		TH1D tmpT = this->GetHistT();	
		TH1D tmpR = this->GetHistR();	

		std::vector<int> tmp_colors ={kOrange-3, kMagenta-3, kGreen-6,kBlue-7,kRed-7};

		TCanvas * c = new TCanvas();
		c->Divide(2,2);
		c->cd(1);	

		tmpR.SetMarkerStyle(21);
		tmpR.SetLineColor(kBlack);
		tmpR.Scale(1,"width");
		tmpR.Draw();

		tmpR.GetXaxis()->SetTitle("E_{QE} [GeV]");
		tmpR.GetYaxis()->SetTitle("Events/MeV");
		tmpR.GetXaxis()->SetTitleSize(0.04);
		tmpR.GetXaxis()->SetLabelSize(0.04);
		tmpR.GetYaxis()->SetLabelSize(0.04);

		c->cd(2);
		tmpT.SetMarkerStyle(21);
		tmpT.SetLineColor(kBlack);
		tmpT.Scale(1,"width");
		tmpT.Draw();

		tmpT.GetXaxis()->SetTitle("E_{#nu} [MeV]");
		tmpT.GetYaxis()->SetTitle("Events/MeV");
		tmpT.GetXaxis()->SetTitleSize(0.04);
		tmpT.GetXaxis()->SetLabelSize(0.04);
		tmpT.GetYaxis()->SetLabelSize(0.04);

		std::vector<TH1D> v_D;
		std::vector<TH1D> v_U;

		int n_noise = 4;

		for(int i=0; i< n_noise; i++){ 
			v_D.push_back(this->GetHistD());
			v_U.push_back(this->GetHistU());
		}

		for(int i=0; i< n_noise; i++){
			SP_DEBUG() << "Generate Poisson noise @ i="<<i<<std::endl;

			this->GenPoissonNoise();
			v_D.at(i) = this->GetHistD();

			v_D.at(i).SetName( ("hope_D_"+std::to_string(i)).c_str());

			SP_DEBUG() << "Unfold @ i="<<i<<std::endl;
			this->Unfold();

			v_U.at(i) = this->GetHistU();
			v_U.at(i).SetName( ("hope_U_"+std::to_string(i)).c_str());

			c->cd(1);
			v_D.at(i).SetLineColor(tmp_colors.at(i));
			v_D.at(i).SetLineWidth(2);
			v_D.at(i).Scale(1,"width");
			v_D.at(i).Draw("same hist");

			c->cd(2);
			v_U.at(i).SetLineWidth(2);
			v_U.at(i).SetLineColor(tmp_colors.at(i));
			v_U.at(i).Scale(1,"width");
			v_U.at(i).Draw("same hist");
		}

		c->cd(3);
		TH2D covU = this->GetCovU();
		double mm=0;
		double mx=0;

		for(int a=1; a<n_t; a++){
			for(int b=1; b<n_t; b++)
			{
				double ank = covU.GetBinContent(a,b) / ( u(a-1) * u(b-1) );
				covU.SetBinContent(a,b ,ank);
				if(ank < mm) mm = ank;
				if(ank > mx) mx = ank;
			}
		}
		SP_DEBUG()<<"Max: "<<mx<<" Min: "<<mm<<std::endl;
		covU.Draw("COLZ3");	

		c->cd(4)->SetLogy();
		TH1D errU = this->GetErrU();
		errU.SetLineColor(kBlack);
		errU.Draw("hist");
		c->SaveAs( (filename+".pdf").c_str(),"pdf");	
		SP_DEBUG() << "end" << std::endl;
	}

	TH1D UnfoldAlgoBase::SampleCovarianceU(TVectorT<double> &result){
		TRandom3 * rangen = new TRandom3(0);
		TDecompChol * chol = new TDecompChol(U,0.0);
		TMatrixT<double> upper_trian(n_t,n_t);
		upper_trian = chol->GetU();

		TVectorT<double> gaus_sample(n_t);
		for(int i=0; i<n_t; i++){
			gaus_sample(i) = rangen->Gaus( u(i), sqrt(U(i,i)));	
		}
		gaus_sample = U*gaus_sample;


		for(int i=0; i<n_t; i++){
			hist_u->SetBinContent(i+1, gaus_sample(i));
			result(i) = gaus_sample(i); 
		}
		return *hist_u;


	}
	void UnfoldAlgoBase::MCbiasCalc(TRandom3 * rangen ){
		//Step 1. Get random (but slightly different "true" spectra)
		TVectorD pois_u = u;
/*		
		b.Zero();
		int num_mc = 1;
		for(int i=0; i< num_mc; i++){

			for(int a=0; a<n_t; a++){
				pois_u(a) = rangen->Poisson( u(a) ); 
			}


			TVectorD re = A*pois_u;
			auto tmp_alg = new UnfoldAlgoBase(*this);
			tmp_alg->Setd(&re);
			tmp_alg->Unfold();

			for(int a=0; a<n_t; a++){
				b(a) += (tmp_alg.u(a)-pois_u(a))/((double)num_mc);
			}			
}
*/
		}



	void UnfoldAlgoBase::TestRegularization(std::string filename, double low_kreg, double high_kreg, int num_kreg){
		TCanvas * c = new TCanvas("","",1200,1200);
		c->Divide(2,3);
		bool logscale = false;

		SP_DEBUG()<<"Starting test from k="<<low_kreg<<" to "<<high_kreg<<std::endl;
		
		if(low_kreg < 0){//must be a log scale
			logscale= true;
		}


		std::vector<double> MSE(num_kreg,0);
		std::vector<double> MSEp(num_kreg,0);
		std::vector<double> MV(num_kreg,0);
		std::vector<double> x(num_kreg,0);    	
		std::vector<double> MB(num_kreg,0);    	
		std::vector<double> MBoE(num_kreg,0);    	
		

		for(double k=0; k<num_kreg; k++){
			double kreg = ceil( (high_kreg-low_kreg)*k/num_kreg+low_kreg );
			if(logscale){
				kreg = pow(10, (high_kreg-low_kreg)*k/num_kreg+low_kreg );
			}

			SP_DEBUG()<<"On k="<<kreg<<std::endl;

			this->SetRegularization(kreg);
			this->Unfold();

			for(int a=0; a<n_t; a++){
				//if( hist_u->GetBinCenter(a) > 800 ) break;
				SP_DEBUG()<<"Bin: "<<a<<" U: "<<U(a,a)<<" b^2 "<<fabs(b(a)*b(a))<<std::endl;
				MSE[k] += 1.0/((double)n_t)*( U(a,a) + b(a)*b(a));
				MSEp[k] += 1.0/((double)n_t)*(  (U(a,a) + b(a)*b(a))/u(a) );
				MV[k] += 1.0/((double)n_t)*U(a,a);
				MB[k] += 1.0/((double)n_t)*fabs(b(a)*b(a));
				MBoE[k] += 1.0/((double)n_t)*fabs(b(a)*b(a))/B(a,a);
					

			}

			SP_DEBUG()<<"RS: "<<" "<<(double)n_t*MV[k]<<" "<<(double)n_t*MB[k]<<std::endl;

			if(logscale){
			x[k]=log10(kreg);
			}else{
			x[k]=kreg;
			}

			this->U.Zero();
			this->b.Zero();


		}


		TGraph *gMV = new TGraph(num_kreg,&x[0],&MV[0]);
		TGraph *gMB = new TGraph(num_kreg,&x[0],&MB[0]);
		TGraph *gMSE = new TGraph(num_kreg,&x[0],&MSE[0]);
		TGraph *gMSEp = new TGraph(num_kreg,&x[0],&MSEp[0]);
		TGraph *gCurve = new TGraph(num_kreg,&MV[0],&MB[0]);
		TGraph *gBiasErr = new TGraph(num_kreg,&x[0],&MBoE[0]);



		c->cd(1)->SetLogy();
		gMV->SetTitle("Minimum Variance");
		gMV->GetXaxis()->SetTitle("Number of Iterations");
		gMV->GetYaxis()->SetTitle("Average Variance");
		
		gMV->SetMarkerStyle(2);
		gMV->SetMarkerSize(3);

		gMV->Draw("ACP");


		c->cd(2)->SetLogy();

		gMB->SetTitle("Minimum Bias");
		gMB->GetXaxis()->SetTitle("Number of Iterations");
		gMB->GetYaxis()->SetTitle("Average Bias");
		gMB->SetMarkerStyle(2);
		gMB->SetMarkerSize(3);
		gMB->Draw("ACP");

		c->cd(3)->SetLogy();
		gMSE->SetTitle("Minimum Square Error MSE");
		gMSE->GetXaxis()->SetTitle("Number of Iterations");
		gMSE->GetYaxis()->SetTitle("Square Error");
		gMSE->SetMarkerStyle(2);
		gMSE->SetMarkerSize(3);

		gMSE->Draw("ACP");

		c->cd(4)->SetLogy();
		gMSEp->SetTitle("Modified Minimum Square Error MSE");
		gMSEp->GetXaxis()->SetTitle("Number of Iterations");
		gMSEp->GetYaxis()->SetTitle("Square Error");
	gMSEp->SetMarkerStyle(2);
		gMSEp->SetMarkerSize(3);

		gMSEp->Draw("ACP");

		TPad *p4 = (TPad*)c->cd(5);
		p4->SetLogy();
		p4->SetLogx();
		gCurve->SetTitle("Variance - Bias curve ");
		gCurve->GetYaxis()->SetTitle("Avg Bias");
		gCurve->GetXaxis()->SetTitle("Avg Variance");
	gCurve->SetMarkerStyle(2);
		gCurve->SetMarkerSize(3);
	
	gCurve->Draw("ACP");

		TPad *p5 = (TPad*)c->cd(6);
		p5->SetLogy();
		gBiasErr->SetTitle("Bias over Bias Error / Ndof ");
		gBiasErr->GetYaxis()->SetTitle("Bias/Err/Ndof");
		gBiasErr->GetXaxis()->SetTitle("Number of Iterations");
		gBiasErr->SetMarkerStyle(2);
		gBiasErr->SetMarkerSize(3);
		gBiasErr->Draw("ACP");




		c->SaveAs( (filename+".pdf").c_str(),"pdf");		
	}
}
#endif
