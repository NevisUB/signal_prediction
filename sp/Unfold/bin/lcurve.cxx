#include "Unfold/SPIO.h"
#include "Unfold/SPManager.h"
#include "Unfold/UnfoldAlgoDAgnostini.h"
#include "Unfold/UnfoldAlgoInverse.h"
#include "TCanvas.h"

int main(int argc, char** argv) {

  sp::SPIO a;

  a.add_mc_in_file("/rootfiles/filterd_ccqe_nue_nuebar.root");
  a.set_mc_tree_name("MiniBooNE_CCQE");
  a.initialize();
  
  std::vector<std::string> var_v = {"Energy"};
  std::vector<double> bins_lo_v = {200,300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500,1600,1700,1800,1900,2000,2500,3000};
  std::vector<double> bins_lo_v2 = { 200.,300.,375.,475.,550.,675.,800.,950.,1100.,1300.,1500.,3000.};
 // a.add_reco_parameter(var_v,bins_lo_v);

  var_v[0] = "NuMomT";
  a.add_true_parameter(var_v,bins_lo_v2,sp::kOP_GEV2MEV);


  var_v = {"RecoEnuQE"};
  a.add_reco_parameter(var_v,bins_lo_v2,sp::kOP_GEV2MEV);


  std::cout<<"int_response matrix"<<std::endl;
  a.init_response_matrix();


  std::cout<<"fill resp"<<std::endl;
  a.fill_responses();

  std::cout<<"write"<<std::endl;
  a.write_unfold_file();

  std::cout<<"Begininning"<<std::endl;

 
  sp::UnfoldAlgoDAgnostini alg; 
//  sp::UnfoldAlgoInverse alg;




  alg.Initialize( &a.Responses().at(0) );
  alg.SetRegularization(5);
//  alg.GenPoissonNoise();
//  alg.Unfold();
/*
  TCanvas * c = new TCanvas();
  c->cd();
  TH1D tmp = alg.GetHistU();
  TH1D tmpT = alg.GetHistT();


  tmp.Scale(1,"width");
  tmp.Draw();

  tmpT.SetLineColor(kRed-6);
  tmpT.Scale(1,"width");
  tmpT.Draw("same");
  
  c->SaveAs("test.pdf","pdf");
*/

	double pot_scale = 6.46/41.10;

  //Got to tihnk more bout flows
  std::vector<double> miniobs = {0,232,  156,  156,   79,   81,   70,   63,   65,   62,   34,   70, 0};
  std::vector<double> minibkg = {0,180.80171,108.22448,120.03353,63.887782,89.806966,67.249431,69.855878,57.014477,51.846417,38.586738,69.381391,0};
		
  	std::vector<double> sigerr(13,0);


	TVectorD mini_signal(13);
	for(int i=0; i<13; i++){
		mini_signal(i) = (miniobs.at(i)-minibkg.at(i))+alg.r(i)*pot_scale    ;
		//mini_signal(i) = alg.r(i)*pot_scale  ;
		//mini_signal(i) = miniobs.at(i)-minibkg.at(i);
		sigerr.at(i) =sqrt(miniobs.at(i)+minibkg.at(i)) ;
		if(sigerr.at(i)!=sigerr.at(i)) {
			std::cout<<"Failure, should fabs"<<std::endl;
			exit(EXIT_FAILURE);

		}
	}

	TMatrixD sigcorr(13,13);
	sigcorr.Zero();
	for(int i=0; i<13;i++){
		sigcorr(i,i)=pow(sigerr.at(i),2.0);
	}
	alg.Setd(&mini_signal);
	alg.SetD(&sigcorr);

  alg.TestRegularization("lcurve", 1,20,19);


  return 0;
}
