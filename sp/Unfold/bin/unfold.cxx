#include "Unfold/SPIO.h"
#include "Unfold/SPManager.h"
#include "Unfold/UnfoldAlgoDAgnostini.h"
#include "TCanvas.h"

int main(int argc, char** argv) {

  sp::SPIO a;

  a.add_mc_in_file("/rootfiles/filterd_ccqe_nue_nuebar.root");
  a.set_mc_tree_name("MiniBooNE_CCQE");
  a.initialize();
  
  std::vector<std::string> var_v = {"Energy"};
  std::vector<double> bins_lo_v = {200,300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500,1600,1700,1800,1900,2000,2500,3000};
  std::vector<double> bins_lo_v2 = { 200.,300.,375.,475.,550.,675.,800.,950.,1100.,1300.,1500.,3000.};
 // a.add_reco_parameter(var_v,bins_lo_v);

  var_v[0] = "NuMomT";
  a.add_true_parameter(var_v,bins_lo_v2,sp::kOP_GEV2MEV);


  var_v = {"RecoEnuQE"};
  a.add_reco_parameter(var_v,bins_lo_v2,sp::kOP_GEV2MEV);


  std::cout<<"int_response matrix"<<std::endl;
  a.init_response_matrix();


  std::cout<<"fill resp"<<std::endl;
  a.fill_responses();

  std::cout<<"write"<<std::endl;
  a.write_unfold_file();

  std::cout<<"Begininning"<<std::endl;

 
  sp::UnfoldAlgoDAgnostini alg;
  alg.Initialize( &a.Responses().at(0) );
  alg.SetRegularization(5);
//  alg.GenPoissonNoise();
//  alg.Unfold();
/*
  TCanvas * c = new TCanvas();
  c->cd();
  TH1D tmp = alg.GetHistU();
  TH1D tmpT = alg.GetHistT();


  tmp.Scale(1,"width");
  tmp.Draw();

  tmpT.SetLineColor(kRed-6);
  tmpT.Scale(1,"width");
  tmpT.Draw("same");
  
  c->SaveAs("test.pdf","pdf");
*/

  double pot_scale = 6.46/41.10;

  alg.TestUnfolding("hope");

  //Got to tihnk more bout flows
  std::vector<double> miniobs = {0,232,  156,  156,   79,   81,   70,   63,   65,   62,   34,   70, 0};
  std::vector<double> minibkg = {0,180.80171,108.22448,120.03353,63.887782,89.806966,67.249431,69.855878,57.014477,51.846417,38.586738,69.381391,0};
		
  

  TVectorD mini_signal(13);
  for(int i=0; i<13; i++){
		mini_signal(i) = (miniobs.at(i)-minibkg.at(i))+alg.r(i)*6.46/44.0    ;
		//mini_signal(i) = (miniobs.at(i))    ;
  }

  alg.Setd(&mini_signal);

  for(int i=0; i< alg.n_r; i++){
		std::cout<<"d: "<<i<<" "<<alg.d(i)<<std::endl;
  } 
  alg.Unfold();

  TH1D ans =  alg.GetHistU();
  TH1D sig = alg.GetHistD();
  TH1D reco = alg.GetHistR();
  TH1D truth = alg.GetHistT();

  TFile *f = new TFile("ans.root","RECREATE");
  TCanvas *c = new TCanvas();

  
  sig.SetLineColor(kRed); 
  ans.SetLineColor(kBlue); 
  truth.SetLineColor(kBlack);
  reco.SetLineColor(kBlack); 
 
  sig.SetMarkerStyle(21);
  sig.Scale(1,"width");
  ans.Scale(1,"width");
  truth.Scale(6.46/44.0,"width");
  reco.Scale(6.46/44.0,"width");

  sig.Draw();
  reco.Draw("same");

  c->Write();

  TCanvas *c2 =  new TCanvas();
  c2->cd();


  std::vector<int> cols = {kCyan,kGreen+3, kOrange-3, kMagenta-3, kGreen-6,kBlue-7,kRed-7};
  std::vector<TH1D> us(6);



  for(int k=1; k<7; k++){

	  alg.SetRegularization(k);
	  alg.Unfold();

	  us.at(k-1) = alg.GetHistU();
	  us.at(k-1).SetLineColor(cols.at(k));
	  us.at(k-1).Scale(1,"width");
	  if(k==1){ us.at(k-1).Draw();
	  }else{
		  us.at(k-1).Draw(" same");
	  } 


  }
  truth.SetLineColor(kBlack);
  truth.SetMarkerStyle(21);
  truth.Draw("same");

  c2->Write();


  

  TCanvas *cr = new TCanvas();
  TH1D ratio = us.at(0);
  TH1D ratio2 = us.at(5);
  ratio.Divide(&truth);
  ratio2.Divide(&truth);
  ratio.Draw();
  ratio2.Draw("same");
  for(int i=0; i<=ratio.GetNbinsX()+2; i++){
	std::cout<<"Ratio "<<ratio.GetBinContent(i)<<std::endl;
  }
  cr->Write();

  f->Close();
  return 0;
}
