#ifndef __UNFOLDALGOBASE_H__
#define __UNFOLDALGOBASE_H__

#include "Response.h"
#include "UnfoldUtil.h"
#include <TRandom3.h>
#include <TCanvas.h>
#include <TColor.h>
#include <string.h>
#include <TGraph.h>
#include "math.h"

namespace sp {

	class UnfoldAlgoBase {
		public:
			std::string name;
			int seed;

			TH1D *hist_u;
			TH1D *hist_r;

			double regularization;
			int n_r;
			int n_t;

			TVectorD r;
			TVectorD t;

			TMatrixD N; //number of events
			TMatrixD A;// response

			TVectorD d;
			TMatrixD D;

			TVectorD b;//bias
			TVectorD ep;//efficiency n_t

			TVectorD u;
			TMatrixD U;
			TMatrixD UA;
			TMatrixD UD;

		public:
			UnfoldAlgoBase() { 
				name="BASE";
				seed = 0;
			}
			virtual ~UnfoldAlgoBase() {}

			void Initialize(const Response *);

			//
			// Pure virtual moth-ods called by Manager
			//

			virtual void Unfold() = 0;



			void Unfold(const TVectorD* d_in);
			void Unfold(const TVectorD* d_in, const TMatrixD* D_in);


			void TestRegularization( std::string filename, double low, double high, int num);
			void TestUnfolding(std::string in);
//
			// Some generic useful functions on an alorithm
			//

			void GenPoissonNoise();
			void SetRegularization(double reg);
			double GetRegularization();
			void SetSeed(int);

			void Setd(TVectorD* d_in);
			void SetD(TMatrixD * din);


			TH1D GetHistU();
			TH1D GetHistT();
			TH1D GetHistR();
			TH1D GetHistD();

			TH1D GetHistEff();

			TH2D GetHistA();
			TH2D GetCovU();
			TH1D GetErrU();


	};


}

#endif
