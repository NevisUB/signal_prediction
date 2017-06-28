#include "miniSVD.h"

int miniSVD::rotateRescale(TMatrixT<double> * tilde_A, TVectorT<double> *r, TMatrixT<double> * Q, TMatrixT<double> * A ){
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

int miniSVD::rotateRescale(TVectorT<double> * tilde_b, TVectorT<double> *r, TMatrixT<double> * Q, TVectorT<double> * b ){
	std::cout<<Q->GetNrows()<<" "<<tilde_b->GetNrows()<<" "<<b->GetNrows()<<" "<<r->GetNrows()<<std::endl;
	*tilde_b = *(Q)*(*b);

	for(int i=0; i< A->GetNrows(); i++){
		(*tilde_b)(i) = 1.0/ (*r)(i) * (*tilde_b)(i);

	}

	return 0;

}

miniSVD::miniSVD(TMatrixT<double>* Ain, TVectorT<double>* bin, TVectorT<double> *xiniin){

	A = Ain;
	b=bin;
	xini=xiniin;

	miniSVD::init(0.001);
}


miniSVD::miniSVD(TH2D *Ain, TH1D * bin, TH1D * xiniin){
	int nn = Ain->GetNbinsX()+2;
	A= new TMatrixT<double>(nn,nn);
	b = new TVectorT<double>(nn);
	xini = new TVectorT<double>(nn);

	ans = (TH1D*)xiniin->Clone("ans");

	for(int i=0; i<nn; i++){
		(*b)(i) = bin->GetBinContent(i);
		(*xini)(i) = xiniin->GetBinContent(i);
	for(int j=0; j<nn; j++){
		(*A)(i,j) = Ain->GetBinContent(i,j);
	}
	}


	miniSVD::init(0.001);
}





int miniSVD::init(double xi){

	//A b and xini will be loaded above here
	b->Print();
	xini->Print();


	nt = A->GetNrows();
	nr = A->GetNrows();

	C = new TMatrixT<double>(nt,nt);
	B = new TMatrixT<double>(nt,nt);
	X = new TMatrixT<double>(nt,nt);
	inv_C = new TMatrixT<double>(nt,nt);
	inv_X = new TMatrixT<double>(nt,nt);

	std::cout<<"Filling C"<<std::endl;	
	for(int i=0; i<nt;i++){
		for(int j=0; j<nt;j++){
			if(i==j){
				if(i==0 || i==nt-1){
					(*C)(i,j)=-1.0+xi;
				}else{
					(*C)(i,j)=-2.0+xi;
				}
			}else{
				(*C)(i,j)=0.0;
			}
		}
	}
	for(int i=0; i<nt-1;i++){
		(*C)(i,i+1) = 1;
	}
	for(int i=0; i<nt-1;i++){
		(*C)(i+1,i) = 1;
	}

	std::cout<<"SVD Decompose C"<<std::endl;	
	TDecompSVD *svd_inv_C = new TDecompSVD(*C);

	TVectorT<double> cs = svd_inv_C->GetSig();
	TMatrixT<double> c_inv_s(nt,nt);

	std::cout<<"Inverting SVD C"<<std::endl;	
	for(int i=0; i<nt; i++){
		for(int j=0; j<nt; j++){
			if(i==j){
				c_inv_s(i,j)=1.0/cs(i);
			}else{
				c_inv_s(i,j)=0.0;
			}
		}
	}

	TMatrixT<double> cU = svd_inv_C->GetU();
	TMatrixT<double> cV = svd_inv_C->GetV();


	*inv_C = cV*c_inv_s*cU.Transpose(cU);



	std::cout<<"Filling X and B, purely stat"<<std::endl;	
	//for(int i=0; i<nt; i++){
	for(int j=0; j<nt; j++){
		//std::cout<<sqrt((*b)(j))<<" "<<sqrt((*xini)(j))<<std::endl;
		(*B)(j,j)=std::max(1.0, sqrt((*b)(j)));
		(*X)(j,j)=std::max(1.0, sqrt((*xini)(j)));
	}


	//}
	return 0;
}


TH1D* miniSVD::unfold(){

	std::cout<<"Unfolding"<<std::endl;	
	TMatrixT<double> tilde_A(nt,nt);
	TVectorT<double> tilde_b(nt);

	std::cout<<"Decomposing B"<<std::endl;	
	TDecompSVD * svd_B = new TDecompSVD(*B);
	TVectorT<double> r = svd_B->GetSig();
	TMatrixT<double> Q = svd_B->GetU();


	std::cout<<"Rotate+Rescale A"<<std::endl;	
	rotateRescale(&tilde_A, &r, &Q, A);
	std::cout<<"Rotate+Rescale b"<<std::endl;	
	rotateRescale(&tilde_b, &r, &Q, b);


	std::cout<<"Invering X exactly (no decomp)"<<std::endl;	
	//Calculate the invers of the covariance matrix of the unfolded truth vector (x)
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

	int k=2;

	double tau = s(k)*s(k);
	std::cout<<"Tau: "<<tau<<" for a k of "<<k<<std::endl;

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



	return ans;

}
