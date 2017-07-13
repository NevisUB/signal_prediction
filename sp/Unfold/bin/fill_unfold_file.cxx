#include "Unfold/SPIO.h"

int main(int argc, char** argv) {

  sp::SPIO a;

  a.add_mc_in_file("/rootfiles/output_osc_mc_detail_1.root");
  a.set_mc_tree_name("MiniBooNE_CCQE");
  a.initialize();
  
  std::vector<std::string> var_v = {"Energy"};
  std::vector<double> bins_lo_v = {0,100,200,300,400,500,600,700,800,900,1000};
  std::vector<double> bins_hi_v = {0,200,300,400,500,600,700,800,900,1000,1100};
  a.add_reco_parameter(var_v,bins_lo_v,bins_hi_v);

  var_v[0] = "NuMomT";
  bins_lo_v = {0,100,200,300,400,500,600,700,800,900,1000};
  bins_hi_v = {0,200,300,400,500,600,700,800,900,1000,1100};
  a.add_true_parameter(var_v,bins_lo_v,bins_hi_v,sp::kOP_GEV2MEV);

  var_v[0] = "CosTheta";
  bins_lo_v = {-1.0,-0.75,-0.5,-0.25,0.0,.25,.5,.75};
  bins_hi_v = {-0.75,-0.5,-0.25,0.0,.25,.5,.75,1.0};
  a.add_reco_parameter(var_v,bins_lo_v,bins_hi_v);

  if (argc>1) {
    var_v.resize(2);
    var_v = {"Energy","CosTheta"};
    bins_lo_v = {0,100,200,300,400,500,600,700,800,900,1000};
    bins_hi_v = {0,200,300,400,500,600,700,800,900,1000,1100};
    a.add_reco_parameter(var_v,bins_lo_v,bins_hi_v,sp::kOP_EQE);
  }


  a.init_response_matrix();
  a.fill_responses();
  a.write_unfold_file();

  return 0;
}
