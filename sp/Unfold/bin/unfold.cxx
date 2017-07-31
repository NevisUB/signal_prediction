#include "Unfold/Core/SPIO.h"
#include "Unfold/Algo/UnfoldAlgoDAgnostini.h"
#include "Unfold/Algo/UnfoldAlgoInverse.h"
#include "Unfold/Algo/UnfoldAlgoSVD.h"
#include "Unfold/Algo/ModelNueCCQE.h"

#include "TCanvas.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TLine.h"

int main(int argc, char** argv) {

  gStyle->SetOptStat(0);
  sp::SPIO a;
  a.set_verbosity((sp::msg::Level_t)0);

  a.add_mc_in_file("/home/vgenty/signal/simplifyTreeOsc/filtered_ccqe_nue_nuebar/filterd_ccqe_nue_nuebar.root");
  a.set_mc_tree_name("MiniBooNE_CCQE");
  
  sp::ModelNueCCQE model;
  a.set_model(&model);

  a.initialize();

  std::vector<std::string> var_v = {"Energy"};
  std::vector<double> bins_lo_v = {200,300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500,1600,1700,1800,1900,2000,2500,3000};
  std::vector<double> bins_lo_v2 = {200.,300.,375.,475.,550.,675.,800.,950.,1100.,1300.,1500.,3000.};
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
  alg.set_verbosity((sp::msg::Level_t)0);
  //  sp::UnfoldAlgoInverse alg;


  alg.Initialize(&a.Responses().at(0));
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

  alg.TestUnfolding("CCQE_data");

  //Got to tihnk more bout flows
  std::vector<double> miniobs = {0,232,  156,  156,   79,   81,   70,   63,   65,   62,   34,   70, 0};
  std::vector<double> minibkg = {0,180.80171,108.22448,120.03353,63.887782,89.806966,67.249431,69.855878,57.014477,51.846417,38.586738,69.381391,0};
  std::vector<double> sigerr(miniobs.size(),0);
  double Nsignal = 0;


  //changed briefly... 
  TVectorD mini_signal(miniobs.size());
  for(int i=0; i<miniobs.size(); i++){
    mini_signal(i) = (miniobs.at(i)-minibkg.at(i))+alg.r(i)*pot_scale    ;
    Nsignal += mini_signal(i);
    //mini_signal(i) = alg.r(i)*pot_scale  ;
    //mini_signal(i) = miniobs.at(i)-minibkg.at(i);
    sigerr.at(i) =sqrt(miniobs.at(i)+minibkg.at(i)) ;
    if(sigerr.at(i)!=sigerr.at(i)) {
      std::cout<<"Failure, should fabs"<<std::endl;
      exit(EXIT_FAILURE);

    }
  }

  TMatrixD sigcorr(miniobs.size(), miniobs.size());
  sigcorr.Zero();
  for(int i=0; i<miniobs.size();i++){
    sigcorr(i,i)=pow(sigerr.at(i),2.0);
  }
  alg.Setd(&mini_signal);
  alg.SetD(&sigcorr);
  alg.Unfold();
















  TH1D ans =  alg.GetHistU();
  TH1D sig = alg.GetHistD();

  TH1D reco = alg.GetHistR();
  TH1D truth = alg.GetHistT();
  sig.SetLineColor(kRed); 
  ans.SetLineColor(kBlue); 
  truth.SetLineColor(kBlack);
  reco.SetLineColor(kBlack); 

  sig.SetMarkerStyle(21);
  sig.SetError(&sigerr[0]);
  sig.Scale(1,"width");
  ans.Scale(1,"width");
  truth.Scale(pot_scale,"width");
  reco.Scale(pot_scale,"width");

  TFile *f = new TFile("CCQE_ans.root","RECREATE");







  TCanvas *cb = new TCanvas("a","a",1200,1200);
  cb->Divide(2,2);
  cb->cd(1);//->SetLogy();
  TH1D tt = truth;
  TH1D rr = reco;
  tt.SetMarkerStyle(21);
  rr.SetMarkerStyle(21);
  tt.SetMarkerColor(kBlue-6);
  rr.SetMarkerColor(kRed-7);
  tt.SetLineWidth(2);
  rr.SetLineWidth(2);
  tt.SetLineColor(kBlue-6);
  rr.SetLineColor(kRed-7);
  tt.GetXaxis()->SetTitle("Energy (Truth or Reco) [MeV]");
  tt.GetYaxis()->SetTitle("Events/MeV");
  tt.GetXaxis()->SetTitleSize(0.04);
  tt.GetXaxis()->SetLabelSize(0.04);
  tt.GetYaxis()->SetLabelSize(0.04);
  tt.GetYaxis()->SetTitleOffset(1.45);
  tt.SetMinimum(0);
  tt.SetTitle("Intrinsic #nu_{e} CCQE only");
  tt.SetMaximum(5);
  tt.GetXaxis()->SetRange(1,10);

  TLegend * legr = new TLegend(0.58,0.6,0.89,0.89);
  legr->AddEntry(&tt,"True E_{#nu} MC","lep");
  legr->AddEntry(&rr,"Reco E_{QE} MC","lep");
  legr->SetFillStyle(0);
  //	legr->SetBorderSize(0.1);

  tt.Draw();
  rr.Draw("same ap");
  legr->Draw();

  cb->cd(2)->SetRightMargin(0.175);
  TH2D prob_MC = alg.GetHistA();
	
  prob_MC.GetYaxis()->SetTitle("Truth E_{#nu} Bin");
  prob_MC.GetXaxis()->SetTitle("Reco E_{QE} Bin");
  prob_MC.SetTitle("Response Matrix");
  prob_MC.GetXaxis()->SetTitleOffset(1.1);
  prob_MC.GetYaxis()->SetTitleOffset(1.1);
	
  prob_MC.Draw("colz");
  prob_MC.GetXaxis()->SetRange(1,10);
  prob_MC.GetYaxis()->SetRange(1,10);

  cb->cd(3);
  TH1D eff = alg.GetHistEff();
  eff.GetYaxis()->SetTitle("Efficiency");
  eff.GetXaxis()->SetTitle("Truth Bin");
  eff.Draw("hist");
  eff.SetMinimum(0);
  eff.SetMaximum(0.4);


  cb->SaveAs("CCQE_signals_response.pdf","pdf");



  //return 0;












  ///////////////////////////////////////////////////////////////////////

  TCanvas *c = new TCanvas("c","c",1200,800);



  sig.SetTitle("Reconstructed E^{QE}" );
  sig.GetYaxis()->SetTitle("Events/MeV");
  sig.GetXaxis()->SetTitle("Reco E^{QE}_{#nu} [MeV]");
  sig.SetLineColor(kRed-7);
  sig.SetMarkerColor(kRed-7);
  sig.SetMarkerStyle(21);
  sig.SetLineWidth(2);
  sig.Draw("e1");
  sig.GetXaxis()->SetRange(1,10);

  reco.SetFillColor(kGray);
  reco.SetMarkerColor(kBlack);
  reco.SetMarkerStyle(21);
  //	reco.Draw("same e2");
  reco.Draw("same hist");
  sig.Draw("same");

  TLegend * leg1 = new TLegend(0.58,0.6,0.89,0.89);
  leg1->AddEntry(&sig,"6.46e20 POT #nu-mode data","lep");
  //leg1->AddEntry(&reco,"MC intrinsic #nu_{e} CCQE","lef");
  leg1->AddEntry(&reco,"MC intrinsic #nu_{e} CCQE","l");
  leg1->SetFillStyle(0);
  leg1->SetBorderSize(0.0);
  leg1->Draw();







  c->Write();
  c->SaveAs("CCQE_signal.pdf","pdf");


  TCanvas *c2 =  new TCanvas("c1","c1",1200,800);
  TCanvas *cu2 =  new TCanvas("cu2","cu2",1200,800);
  TCanvas *cU =  new TCanvas("c2","c2",3000,2000);
  TCanvas *cR =  new TCanvas("cR","cR",3000,2000);
  c2->Divide(3,2);
  cU->Divide(3,2);
  cR->Divide(3,2);


  std::vector<int> cols = {kBlue-7,kGreen-6,kRed-7, kOrange-3, kMagenta-3, kGreen+3};
  std::vector<TH1D> us(6);
  std::vector<TH2D> US(6);
  std::vector<TH1D> uR(6);
  std::vector<TLegend*> leg(6);
  std::vector<TH1D> us_stat(6);
  std::vector<int> kreg = {1,2,4,6,8,100};
  //std::vector<int> kreg = {1,2,5,8,100};

  for(int k=0; k<6; k++){
    std::string nam = std::to_string(kreg.at(k))+ " Iterations";

    c2->cd(k+1)->SetLeftMargin(0.1);
    alg.SetRegularization(kreg.at(k));
    alg.Unfold();

    //Section to re-fold and calc numbers of events.
    TVectorD refolded = alg.A*alg.u;
    double Nrefold = 0;
    for(int i=0; i<alg.n_r; i++){
      Nrefold+=refolded(i);
    }
    std::cout<<std::setprecision(10)<<"REFOLD || On k="<<kreg.at(k)<<" and refolded-number-of-events is "<<Nrefold<<" . Data is "<<Nsignal<<std::endl;





    std::vector<double> errA(alg.n_t,0.0);
    std::vector<double> errD(alg.n_t,0.0);
    std::vector<double> err(alg.n_t,0.0);
    std::vector<double> errS(alg.n_t,0.0);

    for(int i=0;i<err.size(); i++){
      errD.at(i) = sqrt(fabs(alg.U(i,i)));
      errA.at(i) = sqrt(fabs(alg.UA(i,i)));
      errS.at(i) = sqrt(fabs(alg.u(i)));
      err.at(i) = sqrt(fabs(alg.UA(i,i)+alg.U(i,i)));
      std::cout<<"u= "<<alg.u(i)<<" +/-(stat) "<<errS.at(i)<<" +/-(D) "<<errD.at(i)<<" +/-(A) "<<errA.at(i)<<" +/-(AD) "<<err.at(i)<<" bias "<<alg.b(i)<<std::endl;
    }


    us.at(k) = alg.GetHistU();
    us.at(k).SetTitle(  nam.c_str() );
    us.at(k).GetYaxis()->SetTitle("Events/MeV");
    us.at(k).GetXaxis()->SetTitle("True E_{#nu} [MeV]");
    us.at(k).SetError(&err[0]);
    us.at(k).SetLineColor(kBlack);
    us.at(k).SetLineWidth(2);
    us.at(k).SetFillColor(cols.at(k));
    us.at(k).Scale(1,"width");
    us.at(k).Draw("E2");

    us.at(k).SetMinimum(0);
    us.at(k).GetXaxis()->SetRange(1,10);
    us.at(k).SetMaximum(10);//has to be after!

    us_stat.at(k) = alg.GetHistU();
    us_stat.at(k).SetError(&errS[0]);
    us_stat.at(k).SetLineColor(kBlack);
    us_stat.at(k).SetLineWidth(1);
    us_stat.at(k).Scale(1,"width");
    us_stat.at(k).Draw("same E1");

    truth.SetLineColor(kBlack);
    truth.SetLineWidth(2);
    truth.SetMarkerStyle(21);
    truth.SetMarkerColor(kBlack);
    truth.Draw("same hist");

    leg.at(k) = new TLegend(0.5,0.65,0.89,0.89);
    leg.at(k)->AddEntry(&us.at(k),"Raw excess","lef");
    leg.at(k)->AddEntry(&truth,"MC True #nu_{e} CCQE","l");
    leg.at(k)->SetFillStyle(0);
    leg.at(k)->SetBorderSize(0);
    leg.at(k)->Draw();

    if(k==1){
      cu2->cd();
      us_stat.at(k).SetFillColor(kGreen+2);
      us_stat.at(k).SetLineColor(kGreen+2);
      us_stat.at(k).SetMarkerColor(kBlack);
      us_stat.at(k).SetMarkerStyle(29);
      us.at(k).Draw("E2");
      us_stat.at(k).Draw("same E2");
      us_stat.at(k).Draw("same P");
      truth.Draw("same hist");
      leg.at(k)->Draw();

    }


    cU->cd(k+1)->SetLogz();
    US.at(k) = alg.GetCovU();
    US.at(k).SetTitle(  nam.c_str() );
    US.at(k).GetYaxis()->SetTitle("True E_{#nu} Bin a");
    US.at(k).GetXaxis()->SetTitle("True E_{#nu} Bin b");
    US.at(k).Draw("colz");

    cR->cd(k+1);
    uR.at(k) = alg.GetHistRefold();
    sig.Draw();
    uR.at(k).Scale(1,"width");
    uR.at(k).SetMarkerStyle(29);
    uR.at(k).SetMarkerColor(kBlack);
    uR.at(k).SetLineColor(kBlack);
    uR.at(k).SetMarkerSize(2);
    uR.at(k).Draw("same,hist");

			



  }
  c2->Write();
  cU->Write();
  cu2->SaveAs("CCQE_unfolded.pdf","pdf");
  c2->SaveAs("CCQE_reg_vary.pdf","pdf");
  cU->SaveAs("CCQE_reg_corr.pdf","pdf");
  cR->SaveAs("CCQE_reg_refold.pdf","pdf");



  TCanvas *cr = new TCanvas();
  TH1D ratio = us.at(1);	
  ratio.Divide(&truth);
  ratio.SetFillColor(kBlue-7);
  ratio.GetYaxis()->SetTitle("Ratio to MiniBooNE MC Central Value");
  ratio.SetMarkerStyle(5);
  ratio.SetMarkerSize(1);
  ratio.Draw("e2");
  ratio.SetMaximum(6);


  TLegend * leg2 = new TLegend(0.58,0.6,0.89,0.89);
  leg2->AddEntry(&ratio,"Intrinsic #nu_{e} CCQE Model","lf");
  leg2->SetFillStyle(0);
  leg2->SetBorderSize(0.0);
  leg2->Draw();

  TLine *line = new TLine(200,1,1500,1);
  line->SetLineColor(kBlack);
  line->SetLineStyle(2);
  line->Draw();


  for(int i=0; i<=ratio.GetNbinsX()+2; i++){
    std::cout<<"Ratio "<<ratio.GetBinContent(i)<<std::endl;
  }
  cr->Write();
  cr->SaveAs("CCQE_model_ratio.pdf","pdf");

  f->Close();
  return 0;
}
