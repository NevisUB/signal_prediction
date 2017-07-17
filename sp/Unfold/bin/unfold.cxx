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
 // a.add_reco_parameter(var_v,bins_lo_v);

  var_v[0] = "NuMomT";
  a.add_true_parameter(var_v,bins_lo_v,sp::kOP_GEV2MEV);


    var_v = {"RecoEnuQE"};
    a.add_reco_parameter(var_v,bins_lo_v,sp::kOP_GEV2MEV);


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

  alg.TestUnfolding("hope");


  return 0;
}
