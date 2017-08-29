#ifndef __SPMODELBASE_H__
#define __SPMODELBASE_H__
#include <vector>
#include <string>
#include <iostream>

namespace sp {
	class SPModelBase {
		public:
			SPModelBase() {}
			virtual ~SPModelBase() {}
			virtual bool Valid() = 0;
			virtual float Operate(std::string) = 0;

			int NuType;
			int NUANCEChan;
			int NFSP;

			float RecoEnuQE;
			float NuMomT;

			std::vector<int> *FSPType ;
			std::vector<float>* Time ;
			std::vector<float>* MomX ;
			std::vector<float>* MomY ;
			std::vector<float>* MomZ ;
			std::vector<float>* MomT ;
			std::vector<float>* TimePMT ;



		protected:

	};
}

#endif
