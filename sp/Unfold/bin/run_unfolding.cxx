#include "Unfold/SPIO.h"

int main(int argc, char** argv) {

  sp::SPManager() sm;
  
  sp::SPIO& a = sm.get_io();

  a.add_mc_in_file("/home/vgenty/filteredoutput_osc_mc_detail_1.root");
  a.set_mc_tree_name("MiniBooNE_CCQE");
  a.initialize();

  std::vector<std::string> var_v;
  std::vector<double> bins_lo_v;
  std::vector<double> bins_hi_v;

  var_v = {"Energy"};
  bins_lo_v = {0,100,200,300,400,500,600,700,800,900,1000};
  bins_hi_v = {0,200,300,400,500,600,700,800,900,1000,1100};
  a.add_reco_parameter(var_v,bins_lo_v,bins_hi_v);

  var_v = "NuMomT";
  bins_lo_v = {0,100,200,300,400,500,600,700,800,900,1000};
  bins_hi_v = {0,200,300,400,500,600,700,800,900,1000,1100};
  a.add_true_parameter(var_v,bins_lo_v,bins_hi_v,sp::kOP_GEV2MEV);

  var_v = {"Energy","CosTheta"};
  bins_lo_v = {0,100,200,300,400,500,600,700,800,900,1000};
  bins_hi_v = {0,200,300,400,500,600,700,800,900,1000,1100};
  a.add_reco_parameter(var_v,bins_lo_v,bins_hi_v,sp::kOP_EQE);


  a.init_response_matrix();
  a.fill_responses();
  a.write_unfold_file();

  return 0;
}
