#ifndef __UNFOLDALGODAGNOSTINI_CXX__
#define __UNFOLDALGODAGNOSTINI_CXX__

#include "UnfoldAlgoDAgnostini.h"
#include <sstream>

namespace sp {


  void UnfoldAlgoDAgnostini::Unfold(){
    SP_DEBUG()<<"Beginning unfolding with Algorithm: "<<name<<std::endl;

    TMatrixT<double> dudd(n_t, n_r);
    TVectorT<double> u_last(n_t);

    std::vector<std::vector<std::vector<double>>> dudA_last(n_t, std::vector<std::vector<double>>(n_r, std::vector<double>(n_t,0)));
  


    SP_DEBUG()<<"Setting initial guess."<<std::endl;
    u.ResizeTo(n_t);
    //what is the initial guess for bayes theorem;
    u=t;

    SP_DEBUG()<<"Truth is : ";
    double ntruth =0;
    std::stringstream ss;
    for(int a=0; a<n_t; a++){
      ntruth+=t(a);
      ss<<t(a)<<" ";
    }
    SP_DEBUG()<<ss.str()<<std::endl;

    SP_DEBUG()<<"Reco is : ";
    std::stringstream sr;
    for(int i=0; i<n_r; i++){
      sr<<r(i)<<" ";
    }
    SP_DEBUG()<<sr.str()<<std::endl;



    u_last = u;
    for(int a=0; a<n_t; a++){
      // UNCOMMENT  this for flat prior. Flat prior seems very poor.  
      //u(a) = ntruth/((double)n_t); 
    }

	

    for(int k=0; k<regularization; k++){
      SP_DEBUG()<<"On Iteration "<<k<<" of "<<regularization<<std::endl;


      TMatrixT<double> Munfold(n_t, n_r);	

      for(int a=0; a< n_t; a++){ //loop over the u elemements

	//Calculate "unfolding matrix"
	for(int i=0; i< n_r; i++){
	  double de =0;

	  for(int b=0;b<n_t; b++){
	    de += A(i,b)*u(b);
	  }
	  
	  Munfold(a,i) += A(i,a)*u(a)/(ep(a)*de) ;

				
	  if(de==0){
	    SP_ERROR()<<"ERROR denominator is: "<<de <<" leading to NAN's"<<std::endl;
	    SP_ERROR()<<"ERROR This means the current reconstructed variable in this iteration is 0 for bin "<<i<<std::endl;
	    SP_ERROR()<<"ERROR suggestion is to merge the last or first two reco bins, and repeat"<<std::endl;
	    exit(EXIT_FAILURE);
	  }
	}

      }//end of u loop

      //Unfold!
      u = Munfold*d;

      /***************************************
       *	For each iteration, need to keep track of covariance and error matrix better
       * ************************************/
      TMatrixT<double> dudd_cur(n_t, n_r);
      dudd_cur = Munfold;

      SP_DEBUG()<<"Begininning calculation of covariance matrix due to measurement uncertainty D."<<std::endl;
      if(k !=0){
	for(int a=0; a<n_t; a++){
	  for(int i =0; i<n_r; i++){
	    //already have munfold in

	    if(u_last(a)==0){
	      SP_DEBUG()<<"WARNING: u_last(a) is zero, leading to infinity in dudd_cur: "<<dudd_cur(a,i)<<" @ a = "<<a<<". Previous dudd(a,i)= "<<dudd(a,i)<<" and u(i): "<<u(a)<<std::endl;
	    }	


	    dudd_cur(a,i) += u(a)/u_last(a)*dudd(a,i);


	    for(int j=0; j< n_r; j++){
	      for(int b=0; b< n_t; b++){


		dudd_cur(a,i)-= d(j)*ep(b)/u_last(b)*Munfold(a,j)*Munfold(b,j)*dudd(b,i);


	      }
	    }


	  }
	}


      }

      /***********************************
       *	Then calculate U due to D propagation
       * ********************************/

      TMatrixT<double> tt(n_t,n_r);
      TMatrixT<double> tr(n_r,n_t);
      tr.Transpose(dudd_cur);

      for(int a =0; a<n_t; a++){
	for(int i=0; i<n_r; i++){
	  //		SP_DEBUG()<<"M "<<a<<" "<<i<<" "<<dudd_cur(a,i)<<" "<<D(a,i)<<" "<<tt(a,i)<<std::endl;

	}
      }


      U.Zero();
      U = dudd_cur*D*tr;




      /***********************************
       *	Then calculate U due to Underlying A
       * ********************************/
      std::vector<std::vector<std::vector<double>>> dudA(n_t, std::vector<std::vector<double>>(n_r, std::vector<double>(n_t,0)));

      TVectorD f(n_r);
      f = A*u_last;

		


      SP_DEBUG()<<"Begininning calculation of covariance matrix due to unertainty in response matrix A."<<std::endl;

      for(int a =0; a<n_t; a++){
	for(int j =0; j<n_r; j++){
	  for(int b=0; b<n_t; b++){
	    double dab=0;
	    if(a==b) dab = 1.0;

	    double t1 =1.0/ep(a)*(u_last(a)*d(j)/f(j)-u(a))*dab;			
	    double t2 = -u_last(b)*d(j)/f(j)*Munfold(a,j);
	    double t3 = u(a)/u_last(a)*dudA_last.at(a).at(j).at(b);

	    double t4 =0;
	    for(int l=0; l<n_r; l++){
	      for(int g=0;g<n_t; g++){
		t4+=d(l)*Munfold(a,l)*Munfold(g,l)*dudA_last.at(g).at(j).at(b);

	      }
	    }
	    t4 = -t4*ep(a)/u_last(a);



	    dudA.at(a).at(j).at(b) = t1+t2+t3+t4;

	  }
	}
      }


      //these two are outer indicied
      for(int a=0; a<n_t; a++){			
	for(int b=0; b<n_t; b++){			
	  UA(a,b) = 0;


	  for(int j=0; j<n_r; j++){
	    for(int i=0; i<n_r; i++){
	      for(int g=0; g<n_t; g++){
		UA(a,b)+= dudA.at(a).at(j).at(g)*covA.at(j).at(g).at(i)*dudA.at(b).at(i).at(g);
	      }
	    }


	  }//end of interior loop
	}
      }



      SP_DEBUG()<<"Begininning calculation of bias, using derivative + taylor approx."<<std::endl;
      // And bias! This may not be taking it into account properly as of yet. Using bias approximation from cowans book. should be sufficient. 
      for(int a =0; a<n_t;a++){

	b(a) = 0;
	for(int j=0; j<n_r; j++){
	  double vj = 0;
	  for(int b=0; b<n_t; b++){
	    vj+= A(j,b)*u(b);
	  }
	  b(a) += dudd_cur(a,j)*(vj-d(j));
	  //b(a) += Munfold(a,j)*(vj-d(j));
	}
      }




      dudA_last = dudA;
      dudd = dudd_cur;
      u_last=u;



    }//end of iterative loop


  }

}

#endif
