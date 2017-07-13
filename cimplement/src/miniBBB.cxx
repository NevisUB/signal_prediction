#include "miniBBB.h"

int miniBBB::rotateRescale(TMatrixT<double> * tilde_A, TVectorT<double> *r, TMatrixT<double> * Q, TMatrixT<double> * A ){
	std::cout<<Q->GetNcols()<<" "<<tilde_A->GetNcols()<<" "<<A->GetNcols()<<" "<<r->GetNrows()<<std::endl;

	for(int i=0; i< A->GetNrows(); i++){
		for(int j=0; j< A->GetNrows(); j++){
			double denom = (*r)(i);
			double num = ((*Q)*(*A))(i,j);
			std::cout<<"denom "<<i<<" "<<j<<" "<<denom<<std::endl;
			std::cout<<"nom "<<i<<" "<<j<<" "<<num<<" "<<std::endl;
			(*tilde_A)(i,j) = num/denom;
			if(num/denom != num/denom){exit(EXIT_FAILURE);}

		}
	}

	return 0;

}

int miniBBB::rotateRescale(TVectorT<double> * tilde_b, TVectorT<double> *r, TMatrixT<double> * Q, TVectorT<double> * b ){
	std::cout<<Q->GetNrows()<<" "<<tilde_b->GetNrows()<<" "<<b->GetNrows()<<" "<<r->GetNrows()<<std::endl;
	*tilde_b = *(Q)*(*b);

	for(int i=0; i< A->GetNrows(); i++){
		(*tilde_b)(i) = 1.0/ (*r)(i) * (*tilde_b)(i);

	}

	return 0;

}

miniBBB::miniBBB(TMatrixT<double>* Ain, TVectorT<double>* bin, TVectorT<double> *xiniin){

	A = Ain;
	b=bin;
	xini=xiniin;

	miniBBB::init(0.001);
}


miniBBB::miniBBB(TH2D *Ain, TH1D * bin, TH1D * xiniin, TH1D * truthpassin, TH1D * recomcin){
	int nn = Ain->GetNbinsX()+2;
	A= new TMatrixT<double>(nn,nn);
	K = new TMatrixT<double>(nn,nn);
	b = new TVectorT<double>(nn);
	xini = new TVectorT<double>(nn);
	truthpass =  new TVectorT<double>(nn);
	recomc =  new TVectorT<double>(nn);

	ans = (TH1D*)xiniin->Clone("ans");

	for(int i=0; i<nn; i++){
		(*b)(i) = bin->GetBinContent(i);
		(*xini)(i) = xiniin->GetBinContent(i);
		(*truthpass)(i) = truthpassin->GetBinContent(i);
		(*recomc)(i) = recomcin->GetBinContent(i);
		for(int j=0; j<nn; j++){
			(*A)(i,j) = Ain->GetBinContent(j,i);
			(*K)(i,j) = Ain->GetBinContent(j,i)/xiniin->GetBinContent(j);
		}
	}

	nt=nn;


	std::cout<<"GENERAL tests"<<std::endl;
	for(int i=0; i<nn;i++){
		std::cout<<i<<" "<<(*recomc)(i)<<" "<<((*K)*(*xini))(i)<<std::endl;

	}

}





int miniBBB::init(double xi){
	return 0;
}

TH1D* miniBBB::unfold(){
	for(int j=0; j< nt; j++){
		double tml =0;
		for(int i=0; i< nt; i++){
			tml += (*A)(i,j)*     (*b)(i)/( (*recomc)(i) )  *(*xini)(j)/ ((*truthpass)(j)) ;
			std::cout<<j<<" "<<i<<" "<<tml<<" truthpas "<<(*truthpass)(j)<<" A "<<(*A)(i,j)<<" recomc "<<(*recomc)(i)<<" xini "<<(*xini)(j)<<std::endl;
		}
		ans->SetBinContent(j,tml);
	}

	return ans;
}

TH1D* miniBBB::unfold(int k){

	//this->unfold();
	for(int j=0; j<nt; j++){
		ans->SetBinContent(j, (*xini)(j));
	}

	for(int kk=0; kk<k; kk++){
		TH1D * temp_answer =(TH1D*)ans->Clone( ("temp"+std::to_string(kk)).c_str() );
		std::cout<<"ITERATION: "<<kk<<std::endl;
		
		for(int j=0; j< nt; j++){
			double tml =0;
			double ep = 0;

			for(int i=0; i<nt; i++){
				ep += (*K)(i,j) ;
			}
		//	std::cout<<"EP: "<<1./ep<<" "<<(*xini)(j)/((*truthpass)(j))<<std::endl;

			for(int i=0; i< nt; i++){
				double de =0;
				for(int l=0;l<nt; l++){
					de += (*K)(i,l)*temp_answer->GetBinContent(l);
			//		std::cout<<l<<" "<<temp_answer->GetBinContent(l)<<" k "<<(*K)(i,l)<<" a "<<(*A)(i,l)<<std::endl;
				}
			//	std::cout<<"CC: i"<<i<<" "<<(*recomc)(i)<<" de "<<de<<std::endl;
				tml += (*K)(i,j)*(*b)(i)/de ;
			//	std::cout<<" K "<<(*K)(i,j)<<" b "<<(*b)(i)<<" de "<<de<<std::endl;
				if((*K)(i,j)!=(*K)(i,j) || (*b)(i)!= (*b)(i) || de == 0) exit(EXIT_FAILURE);
			}
			ans->SetBinContent(j,    tml*temp_answer->GetBinContent(j)/ep);
		}
	}

	return ans;

}
