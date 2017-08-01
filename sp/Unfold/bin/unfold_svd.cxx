#include "Unfold/Core/SPIO.h"
#include "Unfold/Algo/UnfoldAlgoSVD.h"
#include "Unfold/Algo/ModelNueCCQE.h"

int main(int argc, char** argv) {

  sp::SPIO a;
  a.set_verbosity((sp::msg::Level_t)0);
a.add_mc_in_file("/rootfiles/filterd_ccqe_nue_nuebar.root");

//  a.add_mc_in_file("/home/vgenty/signal/simplifyTreeOsc/filtered_ccqe_nue_nuebar/filterd_ccqe_nue_nuebar.root");
  a.set_mc_tree_name("MiniBooNE_CCQE");
  
  sp::ModelNueCCQE model;
  a.set_model(&model);

  a.initialize();

  std::vector<std::string> var_v;
  std::vector<double> bins_lo_v2;

  var_v.resize(1);
  bins_lo_v2 = {200.,300.,375.,475.,550.,675.,800.,950.,1100.,1300.,1500.,3000.};

  var_v[0] = "NuMomT";
  a.add_true_parameter(var_v,bins_lo_v2,sp::kOP_GEV2MEV);

  var_v[0] = "RecoEnuQE";
  a.add_reco_parameter(var_v,bins_lo_v2,sp::kOP_GEV2MEV);

  a.init_response_matrix();
  a.fill_responses();
  a.write_unfold_file();

  sp::UnfoldAlgoSVD alg;
  alg.set_verbosity((sp::msg::Level_t)0);

  std::cout << "Set regularization: 2" << std::endl;
  alg.SetRegularization(2);

  std::cout << "Initialize" << std::endl;
  alg.Initialize(&(a.Responses().front()));

  // std::cout << "Testing unfolding" << std::endl;
  // alg.TestUnfolding("CCQE_data");
  
  //Got to think more about over/under flows
  double eps = 1e-10;
  std::vector<double> miniobs = {eps,232,  156,  156,   79,   81,   70,   63,   65,   62,   34,   70, eps};
  std::vector<double> minibkg = {eps,180.80171,108.22448,120.03353,63.887782,89.806966,67.249431,69.855878,57.014477,51.846417,38.586738,69.381391, eps};
  const size_t n_observed = miniobs.size();
  std::vector<double> sigerr(n_observed,0);
  double Nsignal = 0;

  double pot_scale = 6.46/41.10;
  std::cout << "Setting signal and error" << std::endl;
  TVectorD mini_signal(n_observed);
  for(int i=0; i<n_observed; i++){

    mini_signal(i) = (miniobs[i] - minibkg[i]) + alg.r(i) * pot_scale;
    Nsignal += mini_signal(i);

    sigerr[i] = std::sqrt(miniobs[i] + minibkg[i]);
  }

  std::cout << "Initialize correlation matrix" << std::endl;
  TMatrixD sigcorr(n_observed,n_observed);

  for(int i=0; i<n_observed; i++)
    sigcorr(i,i) = sigerr[i] * sigerr[i];

  std::cout << "Set signal" << std::endl;
  alg.Setd(&mini_signal);

  std::cout << "Set signal correlation" << std::endl;
  alg.SetD(&sigcorr);

  std::cout << "Unfold it" << std::endl;
  alg.Unfold();

  return 0;
}
