#include "miniClass.h"


miniClass::miniClass(std::vector<std::string> varin, std::vector<std::string> selin, std::string sigin, std::vector<double> low, std::vector<double> high) : vars(varin), selections(selin), signal_name(sigin), low_bin(low), high_bin(high) {

	int nb =25;
	int nbt = 25;

	for(int i=0; i<selections.size(); i++){
		if(signal_name == selections.at(i)){
			Nsignal = i;
		}
	}

	for(int v =0; v<vars.size(); v++){ 
		std::vector<pass_fail> tmp;
		for(int s =0; s< selections.size(); s++){
			pass_fail tmppf(vars[v], selections[s], TH1D( (vars[v]+"_"+selections[s]+"_pass").c_str() ,"",nb,low.at(v),high.at(v)), TH1D((vars[v]+"_"+selections[s]+"_fail").c_str() ,"",nb,low.at(v), high.at(v)) );
			tmp.push_back(tmppf);
		}
		hists.push_back(tmp);
	}


	pass_colors.resize(selin.size(),kBlack);
	fail_colors.resize(selin.size(),kBlack);

	Ntruth=0;
	truth_var = vars[Ntruth];

	for(int v =0; v<vars.size(); v++){ 
		if(v==Ntruth) continue;
		std::string nam ="Truth_"+vars[Ntruth]+"_Reco_"+vars[v];
		TH2D tp = TH2D( (nam+"_pass").c_str() ,"", nbt, low.at(Ntruth), high.at(Ntruth), nb,low.at(v),high.at(v));
		TH2D tf = TH2D( (nam+"_fail").c_str() ,"", nbt, low.at(Ntruth), high.at(Ntruth), nb,low.at(v),high.at(v));
		TH2D ta = TH2D( (nam+"_all").c_str() ,"", nbt, low.at(Ntruth), high.at(Ntruth), nb,low.at(v),high.at(v));
		TH1D ttruth = TH1D((nam+"_t").c_str(), "", nbt,low.at(Ntruth), high.at(Ntruth));
		TH1D treco = TH1D( (nam+"_r").c_str(), "", nb,low.at(v), high.at(v));
		pass_fail2D tmppf2(vars[v], tp, tf ,ta, ttruth,treco);
		hists2D.push_back(tmppf2);
	}

}

int miniClass::setTruthVar(std::string var){
	for(int i=0; i<vars.size(); i++){
		if(var == vars.at(i)){
			Ntruth = i;
		}
	}

	return Ntruth;

}

int miniClass::setColors(std::vector<int> cols_pass, std::vector<int> cols_fail){
	pass_colors = cols_pass;
	fail_colors = cols_fail;

	for(auto& v1: hists){
		for(int i=0; i< v1.size(); i++){
			v1.at(i).pass.SetLineColor(kBlack);
			v1.at(i).fail.SetLineColor(kBlack);

			v1.at(i).pass.SetFillColor(pass_colors.at(i));
			v1.at(i).fail.SetFillColor(fail_colors.at(i));

		}

	}

	for(int i=0; i< hists.size();i++){
		hists.at(i).at(Nsignal).pass.GetXaxis()->SetTitle("Events");
		hists.at(i).at(Nsignal).fail.GetXaxis()->SetTitle("Events");

		hists.at(i).at(Nsignal).pass.GetXaxis()->SetTitle(vars.at(i).c_str());
		hists.at(i).at(Nsignal).fail.GetXaxis()->SetTitle(vars.at(i).c_str());
	}


}



int miniClass::Fill(std::string which_var, std::string which_sel, bool fPassOsc, double value  , double fWeight){

	for(auto &v1: hists){
		for(auto &v2: v1){
			if(v2.s_name == which_sel && v2.v_name == which_var){

				//std::cout<<"Filling: "<<which_sel<<" "<<which_var<<" pass: "<<fPassOsc<<std::endl;
				if(fPassOsc){
					v2.pass.Fill((double)value,(double)fWeight);

				}else{
					v2.fail.Fill((double)value,(double)fWeight);
				}
			}
		}
	}

	return 0;
}


int miniClass::Fill2D(std::string which_var, std::string which_sel, bool fPassOsc, double valueTruth, double value   , double fWeight){

	if(which_sel == signal_name){
		for(auto &pf2: hists2D){
			if(pf2.v_name == which_var){
				pf2.all.Fill((double)valueTruth,(double)value,(double)fWeight);
				pf2.truth.Fill((double)valueTruth,(double)fWeight);
				if(fPassOsc){
					pf2.pass.Fill((double)valueTruth,(double)value,(double)fWeight);
					pf2.reco.Fill((double)value,(double)fWeight);
				}else{
					pf2.fail.Fill((double)valueTruth,(double)value,(double)fWeight);
				}
			}
		}
	}

	return 0;

}



int miniClass::writeOut(std::string nam){

	TFile *f2 = new TFile(nam.c_str(),"RECREATE" ); 

	for(auto& v1: hists){
		for(auto&v2: v1){
			v2.pass.Write();
			v2.fail.Write();
		}
	}

	for(auto& v2: hists2D){
		v2.pass.Write();
		v2.fail.Write();
	}


	//For each Variable, Make a 2D colz histogram  before and after, aganst defined "Truth" variable
	for(int i=0; i<hists2D.size(); i++){
		TCanvas * temp =  new TCanvas(("c_"+hists2D.at(i).v_name).c_str(),"");


		TH2D h_eff = hists2D.at(i).pass;
		h_eff.SetName( ("mig_"+hists2D.at(i).v_name).c_str());
		TH2D h_eff_log = hists2D.at(i).pass;

		TMatrixT <double> meffInv(h_eff.GetNbinsX()+2 ,h_eff.GetNbinsY()+2);
		TMatrixT <double> meff(h_eff.GetNbinsX()+2 ,h_eff.GetNbinsY()+2);
		TVectorT <double> vt = hists2D.at(i).getTruth();
		TVectorT <double> vr = hists2D.at(i).getReco();
		TVectorT <double> vr2(h_eff.GetNbinsX()+2);


		for(int t=0; t<=h_eff.GetNbinsX()+1; t++){
			for(int r=0; r<=h_eff.GetNbinsY()+1; r++){

				meff(r,t)  = 0;

				//h_eff_log.SetBinContent(t,r, log10(  h_eff.GetBinContent(t,r)/ hists2D.at(i).truth.GetBinContent(t)  ));			
				double d = hists2D.at(i).truth.GetBinContent(t);  
				double n = hists2D.at(i).pass.GetBinContent(t,r);	

				h_eff.SetBinContent(r,t, n/d  );		
				if( d/n != d/n){std::cout<<"ERROR: NAN "<<hists2D.at(i).v_name<<std::endl;}
				meff(r,t)  = n/d;

			}
		}

		std::cout<<"Det: "<<hists2D.at(i).v_name<<" "<<meff.Determinant()<<" X: "<<meff.GetNrows()<<" Y: "<<meff.GetNcols()<<std::endl;	



		if(hists2D.at(i).v_name == "Eqe"){
			hists2D.at(i).truth.Write();
			hists2D.at(i).reco.Write();
			hists2D.at(i).pass.Write();

			meffInv = meff;
			meffInv.Invert();
			std::vector<int> tmp_colors = {kRed,kBlue,kGreen,kOrange};

			meff.Write();
			TRandom3 * rangen= new TRandom3(12128124);
			//Generate 4 new poissonian "reco's"
			std::vector<TH1D> treco;
			std::vector<TH1D> ttruth;
			for(int k=0;k<4;k++){
				TH1D temp = hists2D.at(i).reco;
				for(int p=0;p<=temp.GetNbinsX(); p++){
					double fluc = rangen->Poisson( temp.GetBinContent(p));
					temp.SetBinContent(p,fluc);
				}

				temp.SetLineColor(tmp_colors[k]);
				temp.SetLineWidth(3);

				treco.push_back(temp);
			}

			for(int k=0; k<4; k++){
				int n = treco[k].GetNbinsX()+2;
				TVectorT<double> tmp(n);
				for(int y=0; y< n; y++){
					tmp(y) = treco[k].GetBinContent(y);
				}
				TVectorT<double> temp_truth = meffInv*tmp;

				TH1D temp = hists2D.at(i).reco;
				for(int y=0; y< n; y++){
					temp.SetBinContent(y,temp_truth(y));
				}
				temp.SetLineColor(tmp_colors[k]);
				temp.SetLineWidth(3);

				ttruth.push_back(temp);
			}

			
			TVectorT<double> revR = hists2D.at(i).getReco();
			TVectorT<double> revT = meffInv*revR;
			TH1D h_rev = hists2D.at(i).reco;
			std::cout<<"Start ReEngineer"<<std::endl;
			for(int y=0; y< h_rev.GetNbinsX()+2; y++){
				h_rev.SetBinContent(y, revT(y));
				std::cout<<revT(y)<<" "<<hists2D.at(i).truth.GetBinContent(y)<<std::endl;
			}
	


			TCanvas * test = new TCanvas("Eqe test","");
			test->Divide(2,1);
			test->cd(1);
			hists2D.at(i).truth.SetMaximum(400);
			hists2D.at(i).truth.SetStats(0);
			hists2D.at(i).truth.SetLineWidth(3);
			hists2D.at(i).truth.SetLineColor(kBlack);
			h_rev.SetLineColor(kBlack);
			hists2D.at(i).truth.GetXaxis()->SetTitle("Enu Truth [GeV]");

			hists2D.at(i).truth.SetMarkerStyle(21);
			hists2D.at(i).truth.Draw();
			h_rev.Draw("same hist");
			for(int y=0; y<4; y++){
				ttruth.at(y).Draw("same hist");
			}	


			test->cd(2);
			hists2D.at(i).reco.SetStats(0);
			hists2D.at(i).reco.SetLineWidth(3);
			hists2D.at(i).reco.SetLineColor(kBlack);
			hists2D.at(i).reco.GetXaxis()->SetTitle("Eqe [GeV]");
			hists2D.at(i).reco.SetMarkerStyle(21);
			hists2D.at(i).reco.Draw();
			for(int y=0; y<4; y++){
				treco.at(y).Draw("same hist");
			}	

			test->Write();

			TCanvas * test2 = new TCanvas("Eqe test2","");
			test2->Divide(2,1);
			test2->cd(1);
			hists2D.at(i).truth.Draw();
			h_rev.Draw("same hist");
			test2->cd(2);
			hists2D.at(i).reco.Draw();
			test2->Write();
	

		}

		vr2 = meff*vt;
		std::cout <<"Test of Truth: "<<hists2D.at(i).v_name<<std::endl;
		for(int p=0;p<vt.GetNrows(); p++){
			std::cout<<std::setprecision (10)<< vt(p)<<" ";
		}
		std::cout<<"Test of migration: "<<hists2D.at(i).v_name<<std::endl;
		std::cout<<"Real: ";
		for(int p=0;p<vt.GetNrows(); p++){
			std::cout<<std::setprecision (10)<<  vr(p)<<" ";
		}
		std::cout<<std::endl;
		std::cout<<"Migr: ";
		for(int p=0;p<vt.GetNrows(); p++){
			std::cout<<vr2(p)<<" ";
		}
		std::cout<<std::endl;


		temp->Divide(2,2);
		temp->cd(1);

		hists2D.at(i).pass.GetYaxis()->SetTitle(hists2D.at(i).v_name.c_str());
		hists2D.at(i).pass.GetXaxis()->SetTitle(truth_var.c_str());
		hists2D.at(i).pass.Draw("colz");


		temp->cd(2);
		hists2D.at(i).fail.GetYaxis()->SetTitle(hists2D.at(i).v_name.c_str());
		hists2D.at(i).fail.GetXaxis()->SetTitle(truth_var.c_str());
		hists2D.at(i).fail.Draw("colz");

		temp->cd(3);

		h_eff_log.Sumw2();

		//h_eff_log.Draw("colz");

		h_eff.Write();
		temp->cd(4);
		h_eff.Sumw2();
		h_eff.Draw("colz");

		temp->Write();
	}






	//For each Variable, Make a stacked histogram of all selections, before and after, along with selection "efficiency" and background "efficiency"
	for(int i=0; i<hists.size(); i++){
		std::string vnam = vars.at(i);
		std::cout<<"Making stacked histograms for: "<<vnam<<std::endl;
		THStack * spass = new THStack((vnam+"_pass").c_str(),"");
		THStack * sfail = new THStack((vnam+"_fail").c_str(),"");

		double l =0.6;
		double r = 0.89;

		if(vnam=="CosTheta"){
			l=l-0.49;
			r=r-0.49;
		}



		TLegend * leg_fail = new TLegend(l,0.6,r,0.89);
		TLegend * leg_pass = new TLegend(l,0.6,r,0.89);

		TH1D tempSigPass = hists.at(i).at(Nsignal).pass;
		TH1D tempSigFail = hists.at(i).at(Nsignal).fail;

		int Ngar = Nsignal+1;
		TH1D tempBkgPass = hists.at(i).at(Ngar).pass;
		TH1D tempBkgFail = hists.at(i).at(Ngar).fail;

		for(int j=0; j<hists.at(i).size(); j++){

			std::cout<<"Adding: "<<selections.at(j)<<" to " <<vnam<<" stack."<<std::endl;
			spass->Add(&hists.at(i).at(j).pass);
			sfail->Add(&hists.at(i).at(j).fail);

			leg_fail->AddEntry(&hists.at(i).at(j).fail  ,selections.at(j).c_str(),"f");
			leg_pass->AddEntry(&hists.at(i).at(j).pass  ,selections.at(j).c_str(),"f");

			if(j!= Nsignal && j!= Ngar){
				tempBkgPass.Add(&hists.at(i).at(j).pass);
				tempBkgFail.Add(&hists.at(i).at(j).fail);

			}
		}

		TH1D tempSigAll =  tempSigPass;
		tempSigAll.Add(&tempSigFail);
		TH1D tempBkgAll =  tempBkgPass;
		tempBkgAll.Add(&tempBkgFail);
		tempSigPass.Divide(&tempSigAll);
		tempBkgPass.Divide(&tempBkgAll);


		TCanvas * tc = new TCanvas(  vnam.c_str(), vnam.c_str());
		tc->Divide(2,2);
		tc->cd(1);
		spass->Draw("hist");
		leg_pass->Draw();
		tc->cd(2);
		sfail->Draw("hist");
		leg_fail->Draw();
		tc->cd(3);

		tempSigPass.SetLineColor(kBlack);
		tempSigPass.Sumw2();
		tempSigPass.SetStats(0);      // No statistics on lower plot
		tempSigPass.SetMarkerStyle(21);
		tempSigPass.Draw("ep");       // Draw the ratio plot

		tc->cd(4);
		tempBkgPass.SetLineColor(kBlack);
		tempBkgPass.Sumw2();
		tempBkgPass.SetStats(0);      // No statistics on lower plot
		tempBkgPass.SetMarkerStyle(21);
		tempBkgPass.Draw("ep");       // Draw the ratio plot




		tc->Update();
		tc->Write();
	}




	f2->Close();	
	return 0;
}


/*
   std::vector<TH1D *> miniClass::getVariables(std::string select  ,bool passOrfail){




   }

   std::vector<TH1D *> miniClass::getVariables(std::string select  ,bool passOrfail){




   }
   */
