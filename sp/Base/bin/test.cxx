#include "Base/SPIO.h"

int main() {

  sp::SPIO a;

  a.add_mc_in_file("/home/vgenty/filteredoutput_osc_mc_detail_1.root");
  a.set_mc_tree_name("MiniBooNE_CCQE");
  a.initialize();

  std::string name = "Energy";
  std::vector<double> bins_lo_v = {0,100,200,300,400,500};
  std::vector<double> bins_hi_v = {0,200,300,400,500,600};
  a.add_reco_parameter(name,bins_lo_v,bins_hi_v);

  name = "NuMomT";
  bins_lo_v = {0,.100,.200,.300,.400,.500};
  bins_hi_v = {0,.200,.300,.400,.500,.600};

  a.add_true_parameter(name,bins_lo_v,bins_hi_v);


  a.init_response_matrix();
  a.fill_responses();
  a.write_unfold_file();


  for(auto& response : a.Responses()) {
    std::cout << std::endl;
    std::cout << "Response @ " << &response << std::endl;
    response.dump();
    std::cout << std::endl;
  }

  return 0;
}
