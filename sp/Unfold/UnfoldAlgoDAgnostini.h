#ifndef __UNFOLDALGODAGNOSTINI_H__
#define __UNFOLDALGODAGNOSTINI_H__

#include "UnfoldAlgoBase.h"

namespace sp {

	class UnfoldAlgoDAgnostini : public UnfoldAlgoBase {

		private:
		public:
			UnfoldAlgoDAgnostini() : UnfoldAlgoBase() {
				name = "DAgnostini";
			}

			void Unfold();
			void Unfold(const TVectorD* d) ;
			void Unfold(const TVectorD* d, const TMatrixD * D);

	};


}

#endif
