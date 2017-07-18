#include "Unfold/Core/SPIO.h"
#include "Unfold/Algo/ModelNueCCQE.h"

int main(int argc, char** argv) {

  sp::SPIO a;

  a.add_mc_in_file("/home/vgenty/signal/simplifyTreeOsc/filtered_ccqe_nue_nuebar/filterd_ccqe_nue_nuebar.root");
  a.set_mc_tree_name("MiniBooNE_CCQE");
  a.initialize();
  
  std::vector<std::string> var_v;
  std::vector<double> bins_lo_v;

  var_v = {"Energy"};
  bins_lo_v = {0,100,200,300,400,500,600,700,800,900,1000};
  a.add_reco_parameter(var_v,bins_lo_v);

  var_v = {"NuMomT"};
  bins_lo_v = {0,100,200,300,400,500,600,700,800,900,1000};
  a.add_true_parameter(var_v,bins_lo_v,sp::kOP_GEV2MEV);

  var_v = {"CosTheta"};
  bins_lo_v = {-1.0,-0.75,-0.5,-0.25,0.0,.25,.5,.75};
  a.add_reco_parameter(var_v,bins_lo_v);
  
  var_v = {"Energy","CosTheta"};
  bins_lo_v = {0,100,200,300,400,500,600,700,800,900,1000};
  a.add_reco_parameter(var_v,bins_lo_v,sp::kOP_EQE);
  
  sp::ModelNueCCQE model;
  a.set_model(&model);

  a.init_response_matrix();
  a.fill_responses();
  a.write_unfold_file();

  return 0;
}
