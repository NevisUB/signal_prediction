#include "miniClass.h"


miniClass::miniClass(std::vector<std::string> varin, std::vector<std::string> selin, std::string sigin, std::vector<double> low, std::vector<double> high) : vars(varin), selections(selin), signal_name(sigin) {

	for(int i=0; i<selections.size(); i++){
		if(signal_name == selections.at(i)){
			Nsignal = i;
		}
	}

	for(int v =0; v<vars.size(); v++){ 
		std::vector<pass_fail> tmp;
		for(int s =0; s< selections.size(); s++){
			pass_fail tmppf(vars[v], selections[s], TH1D( (vars[v]+"_"+selections[s]+"_pass").c_str() ,"",20,low.at(v),high.at(v)), TH1D((vars[v]+"_"+selections[s]+"_fail").c_str() ,"",20,low.at(v), high.at(v)) );
			tmp.push_back(tmppf);
		}
		hists.push_back(tmp);
	}


	pass_colors.resize(selin.size(),kBlack);
	fail_colors.resize(selin.size(),kBlack);


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

}



int miniClass::Fill(std::string which_var, std::string which_sel, bool fPassOsc, double value  , double fWeight){

	for(auto &v1: hists){
		for(auto &v2: v1){
			if(v2.s_name == which_sel && v2.v_name == which_var){

				//std::cout<<"Filling: "<<which_sel<<" "<<which_var<<" pass: "<<fPassOsc<<std::endl;
				if(fPassOsc){
					v2.pass.Fill(value,fWeight);

				}else{
					v2.fail.Fill(value,fWeight);
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
