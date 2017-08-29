#ifndef __SPIO_H__
#define __SPIO_H__

#include "Base/sp_base.h"

#include "Response.h"

#include <string>
#include "TChain.h"
#include "TFile.h"
#include "UnfoldTypes.h"
#include "SPModelBase.h"


namespace sp {
  
  class SPIO : public sp_base{

  public:

    SPIO();
    ~SPIO();

  public:
    //
    // Core methods
    //
    bool initialize();
    bool init_response_matrix();
    bool fill_responses();
    bool write_unfold_file();
    
  public:
    //
    // Public setters
    //

    // set the MC input file
    void add_mc_in_file(const std::string& filename);

    // set the MC tree name (MiniBooNE_CCQE)
    void set_mc_tree_name(const std::string& treename);

    // set the root file which holds parameters and responses
    void set_unfold_in_file(const std::string& filename);
    
    // add a true parameter
    void add_true_parameter(const std::vector<std::string>& true_var_name,
			    const std::vector<double>& bin_lo_v,
			    Operation_t op=kOP_INVALID);

    // add a reco parameter    
    void add_reco_parameter(const std::vector<std::string>& var_v,
			    const std::vector<double>& bin_lo_v,
			    Operation_t op=kOP_INVALID);
    
    void set_model(SPModelBase* model) { _model = model; }

  private:
    //
    // Private variables
    //

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
    std::vector<const Response*>  _unfold_response_ptr_v;

    SPModelBase* _model;

  private:
    //
    // Private functions
    //

    // called @ initialize -- read the unfold file and get pointers to existing parameters
    void prepare_unfold_file_parameters();

    // called @ initialize -- read the unfold file and get pointers to existing responses
    void prepare_unfold_file_responses();

    // given a parameter, search the file to see if it exists
    const Parameter* search_unfold_parameters(const Parameter& in_param);

    // given a response, search the file to check if it exists
    const Response*  search_unfold_responses(const Response& in_response);

  public:
    //
    // Public convenience
    //

    // dump TTree in mc input file
    void dump_branches();

    // return a const vector of response matrices
    const std::vector<Response>& Responses() const { return _response_v; }
    
    // generate vector of backgrounds
    std::vector<TH1D> gen_background(const std::string& filename, 
				     const std::string& filename_dirt,
		    		     const std::string& param,
				     const std::vector<double>& bins_lo_v);
    

  };
  
}


#endif
