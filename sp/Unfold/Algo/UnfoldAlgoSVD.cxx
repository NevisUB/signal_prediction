#ifndef __UNFOLDALGOSVD_CXX__
#define __UNFOLDALGOSVD_CXX__

#include "UnfoldAlgoSVD.h"

namespace sp {

  UnfoldAlgoSVD::UnfoldAlgoSVD() : UnfoldAlgoBase() {
    name = "SVD";

    C.ResizeTo(n_t,n_t);
    invC.ResizeTo(n_t,n_t);

    xi=1e-5;


    for(int i=0; i<n_t;i++){
      for(int j=0; j<n_t;j++){
	if(i==j){
	  if(i==0 || i==n_t-1){
	    C(i,j)=-1.0+xi;
	  }else{
	    C(i,j)=-2.0+xi;
	  }
	}else{
	  C(i,j)=0.0;
	}
      }
    }
    for(int i=0; i<n_t-1;i++){
      C(i,i+1) = 1;
    }
    for(int i=0; i<n_t-1;i++){
      C(i+1,i) = 1;
    }

    TDecompSVD *svd_C = new TDecompSVD(C);

    TVectorT<double> cs = svd_C->GetSig();
    TMatrixT<double> c_inv_s(n_t,n_t);

    for(int i=0; i<n_t; i++){
      for(int j=0; j<n_t; j++){
	if(i==j){
	  c_inv_s(i,j)=1.0/cs(i);
	}else{
	  c_inv_s(i,j)=0.0;
	}
      }
    }

    TMatrixT<double> cU = svd_C->GetU();
    TMatrixT<double> cV = svd_C->GetV();

    invC = cV*c_inv_s*cU.Transpose(cU);




  };



  void UnfoldAlgoSVD::Unfold(){
    std::cout<<"sp::UnfoldAlgoSVD::Unfold || Beginning unfolding with Algorithm: "<<name<<std::endl;


    TMatrixT<double> tilde_A(n_t,n_t);
    TVectorT<double> tilde_d(n_t);

    std::cout<<"sp::UnfoldAlgoSVD::Unfold || Decomposing D "<<std::endl;
    TDecompSVD * svd_D = new TDecompSVD(D);
    TVectorT<double> r = svd_D->GetSig();
    TMatrixT<double> Q = svd_D->GetU();


    std::cout<<"sp::UnfoldAlgoSVD::Unfold || Rotate+Rescale A"<<std::endl;	
    rotateRescale(&tilde_A, &r, &Q, &A);
    std::cout<<"sp::UnfoldAlgoSVD::Unfold || Rotate+Rescale d"<<std::endl;	
    rotateRescale(&tilde_d, &r, &Q, &d);


    std::cout<<"Invering X exactly (no decomp)"<<std::endl;	
	

    //GOTTEN TO HERE
    /*	
	for(int j=0; j< nt; j++){
	for(int k=0; k< nt; k++){
	(*inv_X)(j,k)=0;
	for(int i=0; i<nt; i++){
	double denom = ( (*xini)(j)* (*xini)(k));
	//	std::cout<<j<<" "<<k<<" "<<i<<" "<<denom<<std::endl;
	(*inv_X) (j,k) += 1.0/denom*tilde_A(i,j)*tilde_A(i,k);
	}
	}
	}

	std::cout<<"Tilde_A"<<std::endl;
	tilde_A.Print();
	std::cout<<"Inv C"<<std::endl;
	inv_C->Print();
	//Step 4:
	TMatrixT<double>  * tilde_A_inv_C =  new TMatrixT<double>(nt,nt);
	*tilde_A_inv_C = tilde_A*(*inv_C);


	std::cout<<"Decomp tAC^-1"<<std::endl;	
	tilde_A_inv_C->Print();
	TDecompSVD * svd_taic = new TDecompSVD(*tilde_A_inv_C);
	TVectorT<double> s = svd_taic->GetSig();
	TMatrixT<double> U = svd_taic->GetU();
	TMatrixT<double> V = svd_taic->GetV();


	TVectorT<double> d = (U.Transpose(U))*(tilde_b);
	for(int i=0; i<nt;i++){
	std::cout<<i<<" "<<log10(fabs(d(i)))<<std::endl;
	}	

	Double_t x[nt], y[nt], y2[nt];
	for (int i=0; i<nt; i++) {
	x[i] = i*1.0;
	y[i] = fabs(d(i)) ;
	y2[i] = s(i);
	}

	TGraph *gr1 = new TGraph (nt, x, y);
	gr1->SetTitle("|d_{i}|");
	TGraph *gr2 = new TGraph (nt, x, y2);
	gr2->SetTitle("Single Values");
	TCanvas *c = new TCanvas();
	c->Divide(2,1);
	TPad *p1 = (TPad *)(c->cd(1)); 
	p1->SetLogy();
	gr1->Draw("AB");
	TPad *p2 = (TPad *)(c->cd(2)); 

	p2->SetLogy();
	gr2->SetMinimum(1e-3);
	gr2->Draw("AB");


	c->Write();



	double tau = s(k_reg)*s(k_reg);
	std::cout<<"Tau: "<<tau<<" for a k of "<<k_reg<<std::endl;

	std::cout<<"Calculating z, w Z and W"<<std::endl;
	TVectorT<double> z(nt);
	TVectorT<double> w(nt);
	TMatrixT<double> Z(nt,nt);
	TMatrixT<double> W(nt,nt);

	for(int i=0; i<nt; i++){
	z(i) = d(i)*s(i)/(s(i)*s(i)+tau);	
	Z(i,i) = s(i)*s(i)/pow(s(i)*s(i)+tau,2);
	}

	TMatrixT<double> tmp(nt,nt);
	tmp = (*inv_C)*V;
	w = tmp*z;


	TMatrixT<double> VT(nt,nt);
	TMatrixT<double> CiT(nt,nt);
	CiT = *inv_C;
	VT =V;
	VT.Transpose(VT);
	CiT.Transpose(CiT);
	W = (*inv_C)*V*Z*VT*CiT;

	TVectorT<double> xtau(nt);
	TMatrixT<double> Xtau(nt,nt);


	for(int i=0; i<nt;i++){
	xtau(i) = (*xini)(i)*w(i);
	ans->SetBinContent(i, xtau(i));
	std::cout<<"MC truth "<<(*xini)(i)<<" Ufolded Truth: "<<xtau(i)<<" wi: "<<w(i)<<std::endl;
	for(int j=0; j<nt;j++){
	Xtau(i,j) = (*xini)(i)*W(i,j)*(*xini)(j);

	}
	}

	xtau.Print();

	f->Close();

    */


  }

  int UnfoldAlgoSVD::rotateRescale(TMatrixT<double> * tilde_Ain, TVectorT<double> *rin, TMatrixT<double> * Q, TMatrixT<double> * Ain ){
    //std::cout<<Q->GetNcols()<<" "<<tilde_A->GetNcols()<<" "<<A->GetNcols()<<" "<<r->GetNrows()<<std::endl;

    for(int i=0; i< Ain->GetNrows(); i++){
      for(int j=0; j< Ain->GetNrows(); j++){
	double denom = (*rin)(i);
	double num = ((*Q)*(*Ain))(i,j);
	//std::cout<<"denom "<<i<<" "<<j<<" "<<denom<<std::endl;
	//std::cout<<"nom "<<i<<" "<<j<<" "<<num<<" "<<std::endl;
	(*tilde_Ain)(i,j) = num/denom;

	if(num/denom != num/denom){exit(EXIT_FAILURE);}

      }
    }

    return 0;

  }

  int UnfoldAlgoSVD::rotateRescale(TVectorT<double> * tilde_b, TVectorT<double> *rin, TMatrixT<double> * Q, TVectorT<double> * bin ){
	
    *tilde_b = *(Q)*(*bin);

    for(int i=0; i< rin->GetNrows(); i++){
      (*tilde_b)(i) = 1.0/ (*rin)(i) * (*tilde_b)(i);

    }

    return 0;

  }

}

#endif
