#ifndef __UNFOLDALGOBASE_H__
#define __UNFOLDALGOBASE_H__

#include "Response.h"
#include "UnfoldUtil.h"
#include <TRandom3.h>

namespace sp {

	class UnfoldAlgoBase {
		public:
			std::string name;
			int seed;

			TH1D *hist_u;

			double regularization;
			int n_r;
			int n_t;

			TVectorD r;
			TVectorD t;

			TMatrixD N;
			TMatrixD A;

			TVectorD d;
			TMatrixD D;

			TVectorD u;
			TMatrixD U;
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
			virtual void Unfold(const TVectorD* d_in) = 0;
			virtual void Unfold(const TVectorD* d_in, const TMatrixD* D_in) = 0;


//
			// Some generic useful functions on an alorithm
			//

			void GenPoissonNoise();
			void SetRegularization(double reg);
			double GetRegularization();
			void SetSeed(int);

			TH1D GetHistU();

	};


}

#endif
