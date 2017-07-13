#ifndef MINISVD_H_
#define MINISVD_H_

#include <cmath>
#include <vector>
#include <iostream>
#include <TH1D.h>
#include <TH2D.h>
#include <string>
#include <TF1.h>
#include <THStack.h>
#include <TLegend.h>
#include <TCanvas.h>
//#include <TROOT.h>
#include <TColor.h>
#include <TMatrixT.h>
#include <TVectorT.h>
#include <TGraph.h>


#include <ctime>
#include <TFile.h>
#include <TRandom3.h>
#include <TDecompSVD.h>

class miniSVD{
	
	public:
	int nr;
	int nt;

	TMatrixT<double> * A;
	
	TMatrixT<double> * B;
	TVectorT<double> * b;
	TVectorT<double> * xini;

	miniSVD(TMatrixT<double>*, TVectorT<double>*, TVectorT<double>*, int);
	miniSVD(TH2D *, TH1D *, TH1D *,int);
	int k_reg;
	TH1D *ans;

	TMatrixT<double> * C;
	TMatrixT<double> * inv_C;
	
	TMatrixT<double> * X;
	TMatrixT<double> * inv_X;


	int init(double, int); //initilises C and inv_C and B ;
	TH1D * unfold();

	int rotateRescale(TMatrixT<double> * tilde_A, TVectorT<double> *r, TMatrixT<double> * Q, TMatrixT<double> * A );
	int rotateRescale(TVectorT<double> * tilde_b, TVectorT<double> *r, TMatrixT<double> * Q, TVectorT<double> * b );

};

#endif
