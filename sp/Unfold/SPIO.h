#ifndef __SPIO_H__
#define __SPIO_H__

#include "Storage.h"

#include <string>
#include "TChain.h"
#include "TFile.h"
#include "UnfoldTypes.h"

namespace sp {
  
  class SPIO {

  public:

    SPIO();
    ~SPIO();

    //
    // Core methods
    //
  public:
    bool initialize();
    bool init_response_matrix();
    bool fill_responses();
    bool write_unfold_file();

    //
    // Public setters 
    //
  public:
    void add_mc_in_file(const std::string& filename);
    void set_mc_tree_name(const std::string& treename);

    void set_unfold_in_file(const std::string& filename);
    
    void add_true_parameter(const std::vector<std::string>& true_var_name,
			    const std::vector<double>& bin_lo_v,
			    const std::vector<double>& bin_hi_v,
			    Operation_t op=kOP_INVALID);
    
    void add_reco_parameter(const std::vector<std::string>& var_v,
			    const std::vector<double>& bin_lo_v,
			    const std::vector<double>& bin_hi_v,
			    Operation_t op=kOP_INVALID);
    
    
    //
    // Private variables
    //

  private:

    std::vector<std::string> _in_mc_file_v;
    std::string _mc_tree_name;

    std::string _unfold_file_name;
    
    TChain *_in_tree;
    TFile  *_unfold_file;
    
    size_t _in_tree_index;
    size_t _in_n_entries;
    
    std::set<std::string> _in_mc_branch_v;

    std::vector<Parameter> _true_parameter_v;
    std::vector<Parameter> _reco_parameter_v;

    std::vector<Response> _response_v;

    std::vector<const Parameter*> _unfold_parameter_ptr_v;
    std::vector<const Response*> _unfold_response_ptr_v;

    //
    // Private functions
    //
  private:
    void prepare_unfold_file_parameters();
    void prepare_unfold_file_responses();
    const Parameter* search_unfold_parameters(const Parameter& in_param);
    const Response*  search_unfold_responses(const Response& in_response);
    
    //
    // Public convenience
    //
  public:
    void dump_branches();
    const std::vector<Response>& Responses() const { return _response_v; }

  };
  
}


#endif
