#ifndef MINIBBB_H_
#define MINIBBB_H_

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

class miniBBB{
	
	public:
	int nr;
	int nt;

	TMatrixT<double> * A;//This A is number of events
	TMatrixT<double> * K;// This is proper repsonce
	
	TMatrixT<double> * B;
	TVectorT<double> * b;
	TVectorT<double> * xini;
	TVectorT<double> * truthpass;
	TVectorT<double> * recomc;

	miniBBB(TMatrixT<double>*, TVectorT<double>*, TVectorT<double>*);
	miniBBB(TH2D *, TH1D *, TH1D *, TH1D *, TH1D *);

	TH1D *ans;

	TMatrixT<double> * C;
	TMatrixT<double> * inv_C;
	
	TMatrixT<double> * X;
	TMatrixT<double> * inv_X;


	int init(double); //initilises C and inv_C and B ;
	TH1D * unfold();
	TH1D * unfold(int);

	int rotateRescale(TMatrixT<double> * tilde_A, TVectorT<double> *r, TMatrixT<double> * Q, TMatrixT<double> * A );
	int rotateRescale(TVectorT<double> * tilde_b, TVectorT<double> *r, TMatrixT<double> * Q, TVectorT<double> * b );

};

#endif
