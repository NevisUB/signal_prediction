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

		std::cout<<"sp:UnfoldAlgoBase::Initialize || Starting on response: "<<response_in->_name<<std::endl;
		n_t = response_in->_true_param->_hist.GetNbinsX()+2;
		n_r = response_in->_reco_param->_hist.GetNbinsX()+2;

		hist_u = (TH1D*)response_in->_true_param->_hist.Clone("unfold");
		hist_r = (TH1D*)response_in->_reco_param->_hist.Clone("reco");

		u.ResizeTo(n_t);
		b.ResizeTo(n_t);
		ep.ResizeTo(n_t);

		d.ResizeTo(n_r);

		D.ResizeTo(n_r,n_r); D.Zero();	
		U.ResizeTo(n_t,n_t); U.Zero();
		UD.ResizeTo(n_t,n_t); UD.Zero();
		UA.ResizeTo(n_t,n_t); UA.Zero();


		std::cout<<"sp::UnfoldAlgoBase::Initialize || Initilizing truth MC."<<std::endl;
		t.ResizeTo(n_t);
		t = hist2TVec(&response_in->_true_param->_hist);
		std::cout<<"sp::UnfoldAlgoBase::Initialize || Initilized truth MC, dimension: "<<t.GetNrows()<<std::endl;

		for(int i =0; i<n_t; i++){

			if(t(i)==0){
				std::cout<<"sp::UnfoldAlgoBase::Initialize || WARNING Truth has a 0 @ position: "<<i<<" of "<<t.GetNrows()<<std::endl;

			}	
		}


		std::cout<<"sp::UnfoldAlgoBase::Initialize || Initilizing reco MC."<<std::endl;
		r.ResizeTo(n_r);
		r = hist2TVec(&response_in->_reco_param->_hist);
		std::cout<<"sp::UnfoldAlgoBase::Initialize || Initilized reco MC, dimension: "<<r.GetNrows()<<std::endl;

		std::cout<<"sp::UnfoldAlgoBase::Initialize || Initilizing Number of Events."<<std::endl;
		N.ResizeTo(n_r,n_t);
		std::cout<<"sp::UnfoldALgoBase::Initialize || input reponse has dimensions: "<<response_in->_response.GetNrows()<<" "<<response_in->_response.GetNcols()<<std::endl;
		N = response_in->_response;
		std::cout<<"sp::UnfoldAlgoBase::Initialize || Initilized Number of events N, dimension: "<<N.GetNrows()<<" "<<N.GetNcols()<<std::endl;

		std::cout<<"sp::UnfoldAlgoBase::Initialize || Initilizing Response."<<std::endl;
		A.ResizeTo(n_r,n_t);

		for(int i=0;i<n_r; i++){
			for(int a=0; a<n_t; a++){
				A(i,a) = N(i,a)/t(a);
				if(A(i,a)>1){
					std::cout<<"sp::UnfoldAlgoBase::Initialize || ERROR! A has entries > 1."<<std::endl;
					std::cout<<i<<" "<<a<<" A:"<<A(i,a)<<" N:"<<N(i,a)<<" t:"<<t(a)<<" Nin:"<<response_in->_response(i,a)<<" tin: "<<response_in->_true_param->_hist.GetBinContent(a)<<" r: "<<r(i)<<" rin: "<<response_in->_reco_param->_hist.GetBinContent(i)<<std::endl;
					exit(EXIT_FAILURE);
				}
			}
		}
		std::cout<<"sp::UnfoldAlgoBase::Initialize || Initilized Response A, dimension: "<<A.GetNrows()<<" "<<A.GetNcols()<<std::endl;


		std::cout<<"sp::UnfoldAlgoBase::Unfold || Calculating a global efficiency vector."<<std::endl;
		ep.Zero();
		for(int b=0; b<n_t; b++){
			for(int l=0; l<n_r; l++){
				ep(b) += A(l,b) ;
			}
			std::cout<<"sp::UnfoldAlgoBase::Unfold || Efficiency is: "<<ep(b)<<" at bin "<<b<<" "<<std::endl;
			if(ep(b)==0){
				std::cout<<"sp::UnfoldAlgoBase::Unfold || ERROR! Efficiency is: "<<ep(b)<<" at bin #"<<b<<" with t: "<<t(b)<<std::endl;
				for(int l=0; l<n_r; l++){
					std::cout<<"ERRORout: bin: "<<l<<" N: "<<N(l,b)<<" r: "<<r(l)<<" A: "<<A(l,b)<<" t: "<<t(b)<<std::endl;
				}
				std::cout<<"sp::UnfoldAlgoBase::Unfold || ERROR! Efficiency cannot be zero! Suggestion to merge last (first) true bins to solve."<<std::endl;
				exit(EXIT_FAILURE);
			}
		}





	};

	/**********************************************************************
	 *
	 *			Unfolding functions, not sure these work
	 *
	 * *******************************************************************/



	void UnfoldAlgoBase::Unfold(const TVectorD * data, const TMatrixD * Din){
		D.ResizeTo(n_r,n_r);
		D=*Din;

		d.ResizeTo(n_r);
		d = *data;

		this->Unfold();
	}

	void UnfoldAlgoBase::Unfold(const TVectorD * data){
		d.Zero();
		d = *data;

		D.ResizeTo(n_r,n_r);
		D.Zero();
		for(int i=0; i<n_r; i++){
			D(i,i)=d(i);
		}

		this->Unfold();

	}


	/**********************************************************************
	 *
	 *			Setting functions
	 *
	 * *******************************************************************/


	void UnfoldAlgoBase::Setd(TVectorD * din){
		if( din->GetNrows() != n_r){
			std::cout<<"sp::UnfoldAlgoBase::Setd || ERROR: Passed a D_in that is size: "<<din->GetNrows()<<" but n_r is "<<n_r<<std::endl;
			exit(EXIT_FAILURE);
		}
		d.ResizeTo(n_r);
		for(int i=0; i<n_r;i++){
			d(i) = (*din)(i);
		}

	}



	void UnfoldAlgoBase::SetD(TMatrixD * din){
		if( din->GetNrows() != n_r){
			std::cout<<"sp::UnfoldAlgoBase::SetD || ERROR: Passed a D_in that is size: "<<din->GetNrows()<<" but n_r is "<<n_r<<std::endl;
			exit(EXIT_FAILURE);
		}
		D.ResizeTo(n_r,n_r);
		for(int i=0; i<n_r;i++){
			for(int j=0; j<n_r;j++){
				D(i,j) = (*din)(i,j);
			}
		}

	}


	void UnfoldAlgoBase::SetRegularization(double reg){

		std::cout<<"sp::UnfoldAlgoBase::SetRegularization() || Setting regularization parameter to "<<reg<<std::endl;
		regularization = reg;


	}


	void UnfoldAlgoBase::SetSeed(int inseed){
		seed=inseed;
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
				tmp.SetBinContent(i ,a, A(i,a));
			}

		}

		return tmp;
	}



	TH1D UnfoldAlgoBase::GetHistU(){

		for(int i=0; i<n_t; i++){
			hist_u->SetBinContent(i , u(i));
		}

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
			hist_u->SetBinContent(i , fabs(sqrt(U(i,i))));
		}

		return *hist_u;
	}


	TH1D UnfoldAlgoBase::GetHistT(){

		for(int i=0; i<n_t; i++){
			hist_u->SetBinContent(i , t(i));
		}

		return *hist_u;
	}

	TH1D UnfoldAlgoBase::GetHistEff(){

		for(int i=0; i<n_t; i++){
			hist_u->SetBinContent(i , ep(i));
		}

		return *hist_u;
	}


	TH1D UnfoldAlgoBase::GetHistR(){

		for(int i=0; i<n_r; i++){
			hist_r->SetBinContent(i , r(i));
		}

		return *hist_r;
	}

	TH1D UnfoldAlgoBase::GetHistRefold(){
		TVectorD re = A*u;
		for(int i=0; i<n_r; i++){
			hist_r->SetBinContent(i , re(i));
		}

		return *hist_r;
	}


	TH1D UnfoldAlgoBase::GetHistD(){

		for(int i=0; i<n_r; i++){
			hist_r->SetBinContent(i , d(i));
		}

		return *hist_r;
	}



	/**********************************************************************
	 *
	 *			Other general and Testing functions
	 *
	 * *******************************************************************/


	void UnfoldAlgoBase::GenPoissonNoise(){

		TRandom3 * rangen = new TRandom3(seed);

		std::cout<<"sp::UnfoldAlgoBase::GenPoissonNoise() || Setting data as noisey reco MC with seed "<<rangen->GetSeed()<<std::endl;

		d.ResizeTo(n_r);
		d.Zero();

		for(int i=0; i<n_r; i++){
			d(i) = rangen->Poisson( r(i) );
		}

	}



	void UnfoldAlgoBase::TestUnfolding(std::string filename){
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
		//tmpR.SetTitle("Reconstructed E_{QE}");
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
		//tmpT.SetTitle("Reconstructed E_{QE}");
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
			std::cout<<"sp::UnfoldAlgoBase::TestUnfolding || on run "<<i<<std::endl;

			this->GenPoissonNoise();
			v_D.at(i) = this->GetHistD();

			v_D.at(i).SetName( ("hope_D_"+std::to_string(i)).c_str());

			this->Unfold();
			v_U.at(i) = this->GetHistU();

			v_U.at(i).SetName( ("hope_U_"+std::to_string(i)).c_str());

			for(int j=0;j<n_r;j++){
				//		std::cout<<tmpD.GetBinContent(j)<<" "<<tmpU.GetBinContent(j)<<std::endl;
			}

		
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
				double ank = covU.GetBinContent(a,b)/(u(a-1)*u(b-1));
				covU.SetBinContent(a,b ,ank);
				//std::cout<<"Frac: "<<a<<" "<<b<<" "<<covU.GetBinContent(a,b)/(u(a-1)*u(b-1))<<std::endl;
				if(ank < mm) mm=ank;
				if(ank>mx) mx=ank;
			}
		}
		std::cout<<"Max: "<<mx<<" Min: "<<mm<<std::endl;
		covU.Draw("COLZ3");	


		c->cd(4)->SetLogy();
		TH1D errU = this->GetErrU();
		errU.SetLineColor(kBlack);
		errU.Draw("hist");



		c->SaveAs( (filename+".pdf").c_str(),"pdf");	


	}




	void UnfoldAlgoBase::TestRegularization(std::string filename, double low_kreg, double high_kreg, int num_kreg){
		TCanvas * c = new TCanvas("","",1200,800);
		c->Divide(2,2);

		std::cout<<"UnfoldAlgoBase::TestRegularization() || Starting test from k="<<low_kreg<<" to "<<high_kreg<<std::endl;


		std::vector<double> MSE(num_kreg,0);
		std::vector<double> MSEp(num_kreg,0);
		std::vector<double> MV(num_kreg,0);
		std::vector<double> x(num_kreg,0);    	
		std::vector<double> MB(num_kreg,0);    	

		for(double k=0; k<num_kreg; k++){
			double kreg = ceil( (high_kreg-low_kreg)*k/num_kreg+low_kreg );

			std::cout<<"UnfoldAlgoBase::TestRegularization() || On k="<<kreg<<std::endl;

			this->SetRegularization(kreg);
			this->Unfold();

			for(int a=0; a<n_t; a++){
				std::cout<<"R: "<<sqrt(U(a,a))<<" "<<fabs(b(a))<<std::endl;
				MSE[k] += 1.0/((double)n_t)*( U(a,a) + b(a)*b(a));
				MSEp[k] += 1.0/((double)n_t)*(  (U(a,a) + b(a)*b(a))/u(a) );
				MV[k] += 1.0/((double)n_t)*U(a,a);
				MB[k] += 1.0/((double)n_t)*fabs(b(a)*b(a));
			}

			std::cout<<"RS: "<<" "<<(double)n_t*MV[k]<<" "<<(double)n_t*MB[k]<<std::endl;

			
			x[k]=kreg;

			this->U.Zero();
			this->b.Zero();


		}


		TGraph *gMV = new TGraph(num_kreg,&x[0],&MV[0]);
		TGraph *gMB = new TGraph(num_kreg,&x[0],&MB[0]);
		TGraph *gMSE = new TGraph(num_kreg,&x[0],&MSE[0]);
		TGraph *gMSEp = new TGraph(num_kreg,&x[0],&MSEp[0]);

		c->cd(1)->SetLogy();
		gMV->SetTitle("Minimum Variance");
		gMV->GetXaxis()->SetTitle("Number of Iterations");
		gMV->GetYaxis()->SetTitle("Average Variance");
		gMV->Draw("ACP");


		c->cd(2)->SetLogy();

		gMB->SetTitle("Minimum Bias");
		gMB->GetXaxis()->SetTitle("Number of Iterations");
		gMB->GetYaxis()->SetTitle("Average Bias");
		gMB->Draw("ACP");

		c->cd(3)->SetLogy();
		gMSE->SetTitle("Minimum Square Error MSE");
		gMSE->GetXaxis()->SetTitle("Number of Iterations");
		gMSE->GetYaxis()->SetTitle("Square Error");
		gMSE->Draw("ACP");

		c->cd(4)->SetLogy();
		gMSEp->SetTitle("Modified Minimum Square Error MSE");
		gMSEp->GetXaxis()->SetTitle("Number of Iterations");
		gMSEp->GetYaxis()->SetTitle("Modified Square Error");
		gMSEp->Draw("ACP");

		c->SaveAs( (filename+".pdf").c_str(),"pdf");		
	}
}
#endif