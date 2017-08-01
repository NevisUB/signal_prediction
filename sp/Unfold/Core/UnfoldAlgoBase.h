#ifndef __UNFOLDALGOBASE_H__
#define __UNFOLDALGOBASE_H__

#include "Base/sp_base.h"
#include "Response.h"
#include "UnfoldUtil.h"
#include <TRandom3.h>
#include <TCanvas.h>
#include <TColor.h>
#include <string.h>
#include <TGraph.h>
#include "math.h"

namespace sp {

  class UnfoldAlgoBase : public sp_base {
 
  public:

    std::string name;
    int seed;

    TH1D *hist_u;
    TH1D *hist_r;

    double regularization;
    int n_r; // number of reco bins
    int n_t; // number of true bins

    TVectorD r;  // reconstructed 
    TVectorD t;  // true

    TMatrixD R;  // ?
    TMatrixD T;  // ?
 
    TMatrixD N;  // number of events
    TMatrixD A;  // response
    std::vector<std::vector<std::vector<double>>> covA;

    TVectorD d;  // data
    TMatrixD D;  // ?

    TVectorD b;  // bias
    TVectorD ep; // efficiency n_t
    TMatrixD Ep; //Covariance of efficienct

    TVectorD u;  // unfolded spectrum
    TMatrixD U;  // ?
    TMatrixD UA; // ?
    TMatrixD UD; // ?


  public:

    UnfoldAlgoBase(const std::string n = "UnfoldAlgoBase") : 
      sp_base(n), 
      name(n), 
      seed(0) 
    {}

    virtual ~UnfoldAlgoBase() {}
    

    //
    // Initialize from a given Response class
    // 

    void Initialize(const Response *);

    //
    // Pure virtual methods called by Manager
    //

    virtual void Unfold() = 0;

    //
    // Unfold given input spectrum
    //

    void Unfold(const TVectorD* d_in);
    void Unfold(const TVectorD* d_in, const TMatrixD* D_in);

    //
    // 
    //
    void TestRegularization( std::string filename, double low, double high, int num);
    void TestUnfolding(std::string in);

    //
    // Getters
    //

    void GenPoissonNoise();
    double GetRegularization();

    //
    // Setters
    //
    void SetRegularization(double reg);
    void SetSeed(int);
    void Setd(TVectorD* d_in);
    void SetD(TMatrixD * din);

    //
    // Retrieval of algorithm matricies as TH[1,2]D
    //
    TH1D GetHistU();
    TH1D GetHistT();
    TH1D GetHistR();
    TH1D GetHistD();
    TH1D GetHistRefold();

    TH1D GetHistEff();
    TH2D GetCovEff();

    TH2D GetHistA();
    TH2D GetCovU();
    TH1D GetErrU();


  protected:
    virtual void _Initialize_() {};
    
  };


}

#endif
