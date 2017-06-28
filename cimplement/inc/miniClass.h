#ifndef MINICLASS_H_
#define MINICLASS_H_

#include <cmath>
#include <vector>
#include <iostream>
#include <TH1D.h>
#include <TH2D.h>
#include <string>
#include <TF1.h>
#include <THStack.h>
#include <TLegend.h>
#include <TCanvas.h>
//#include <TROOT.h>
#include <TColor.h>
#include <TMatrixT.h>
#include <TVectorT.h>

#include <ctime>
#include <TFile.h>
#include <TRandom3.h>

struct pass_fail{
	std::string name;
	std::string s_name;
	std::string v_name;
	TH1D pass;
	TH1D fail;

	pass_fail(std::string vname, std::string sname, TH1D passin, TH1D failin) : s_name(sname), v_name(vname), pass(passin), fail(failin) {}; 
};

struct pass_fail2D{
	std::string name;
	std::string v_name;
	TH1D truth;
	TH1D reco;

	TH2D pass;
	TH2D fail;
	TH2D all;
	pass_fail2D(std::string vname,  TH2D passin, TH2D failin , TH2D allin, TH1D truthin, TH1D recoin) : v_name(vname), pass(passin), fail(failin), all(allin), truth(truthin),reco(recoin) {}; 

	TVectorT<double> getReco(){
		int n = reco.GetNbinsX()+2;
		TVectorT<double> tmp(n);
		for(int i=0; i< n; i++){
			tmp(i) = reco.GetBinContent(i);
		}
		return tmp;
	};

	TVectorT<double> getTruth(){
		int n = truth.GetNbinsX()+2;
		TVectorT<double> tmp(n);
		for(int i=0; i< n; i++){
			tmp(i) = truth.GetBinContent(i);
		}
		return tmp;
	};
};



class miniClass{
	
	public:
	miniClass(std::vector<std::string>, std::vector<std::string>, std::string,std::vector<double>, std::vector<double>);
	std::vector<std::vector<pass_fail>> hists;
	std::vector<pass_fail2D> hists2D;

	std::vector<double> low_bin;
	std::vector<double> high_bin;


	std::string signal_name;
	std::string truth_var;

	int Nsignal;
	int Ntruth;

	std::vector<int> pass_colors;
	std::vector<int> fail_colors;

	std::vector<std::string> vars;
	std::vector<std::string> selections;

//	std::vector<TH1D *> getVariables(std::string select, bool passOrfail );
//	std::vector<TH1D *> getSelections(std::string var, bool passOrfail );

	int Fill(std::string which_var, std::string which_sel, bool fPassOsc, double value  , double fWeight);
	int Fill2D(std::string which_var, std::string which_sel,  bool fPassOsc, double valueTruth, double value   , double fWeight);
	int writeOut(std::string);
	int setColors(std::vector<int> cols_pass, std::vector<int> cols_fail);

	int setTruthVar(std::string);

};




#endif
