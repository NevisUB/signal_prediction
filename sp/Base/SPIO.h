#ifndef __SPIO_H__
#define __SPIO_H__

#include <vector>
#include <string>
#include "TChain.h"
#include "TFile.h"


namespace sp {

  class SPIO {

  public:
    SPIO();
    ~SPIO() {}

    void add_mc_in_file(const std::string& filename);
    void set_mc_tree_name(const std::string& treename);


    void set_unfold_in_file(const std::string& filename);
    
    void add_true_var(const std::string& true_var_name,
		      const std::vector<double>& bin_lo_v,
		      const std::vector<double>& bin_hi_v) {}
    
    void add_reco_var(const std::string reco_var_name,
		      const std::vector<double>& bin_lo_v,
		      const std::vector<double>& bin_hi_v) {}
    
    
    bool initialize();

    
    
    
  private:

    std::vector<std::string> _in_mc_file_v;
    std::string _mc_tree_name;

    std::string _unfold_file_name;
    
    TChain *_in_tree;
    TFile  *_unfold_file;
    
    size_t _in_tree_index;
    size_t _in_n_entries;

    std::set<std::string> _in_mc_branch_v;

  public:
    void dump_branches();
    
  };
  
}


#endif
