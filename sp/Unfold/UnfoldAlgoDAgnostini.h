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

	};


}

#endif
