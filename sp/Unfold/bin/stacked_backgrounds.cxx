#include "Unfold/Core/SPIO.h"

#include "Combined/CombinedTypes.h"
#include "Combined/CombinedUtil.h"

#include "TApplication.h"
#include "TCanvas.h"
#include "THStack.h"
#include "TLegend.h"

#include <iostream>

int main(int argc, char** argv) {

  std::string fname(argv[1]);
  std::string param = "RecoEnuQE";
  std::vector<double> bins_lo_v = {200.,300.,375.,475.,550.,675.,800.,950.,1100.,1300.,1500.,3000.};

  for(auto& v : bins_lo_v) v /= 1000.0;

  sp::SPIO spio;
  spio.set_verbosity((sp::msg::Level_t)0);

  //
  // Get the vector of histograms for each background type
  //
  auto th1d_v = spio.gen_background(fname,param,bins_lo_v);

  TApplication app("app", 0, 0);

  TCanvas c1("c1","test",512,512);
  c1.cd();

  THStack ths("ths","");
  TLegend tl(0.1,0.7,0.48,0.9);

  for(size_t bkgd_id = 0; bkgd_id < (size_t) sp::kBKGD_MAX; ++bkgd_id) {
    if ((sp::StackedBkgdType_t)bkgd_id == sp::kBKGD_INVALID) continue;
    if ((sp::StackedBkgdType_t)bkgd_id == sp::kBKGD_DIRT)    continue;

    auto& th1d = th1d_v[bkgd_id];

    for(size_t bin_id=1; bin_id < bins_lo_v.size(); ++bin_id) {
      auto dx = bins_lo_v.at(bin_id) - bins_lo_v.at(bin_id-1);

      auto modifier = 1.0 / (dx * 1000.0);
      modifier *= 0.157; // POT normalize

      auto bin_content = th1d.GetBinContent(bin_id) * modifier;
      auto bin_error   = th1d.GetBinError(bin_id) * modifier;
      
      th1d.SetBinContent(bin_id, bin_content);
      th1d.SetBinError(bin_id, bin_error);
    }

    th1d.SetFillColor((Color_t)bkgd_id);
    ths.Add(&th1d);
    tl.AddEntry(&th1d,sp::StackedBkgd2String((sp::StackedBkgdType_t)bkgd_id).c_str());
  }
  
  ths.Draw("hist");
  tl.Draw("sames");
  c1.Update();
  c1.Modified();
  app.Run();
  
  return 0;
}
