#include <iostream>
#include <Eigen/Dense>

using Eigen::MatrixXf;
using Eigen::Matrix2f;
using Eigen::VectorXf;

int main()
{
  MatrixXf m(2,2);
  m(0,0) = 3;
  m(1,0) = 2.5;
  m(0,1) = -1;
  m(1,1) = m(1,0) + m(0,1);
  std::cout << m << std::endl;


  VectorXf a(10);
  VectorXf b(9);
  for(size_t i=0; i<a.size(); ++i) a[i] = 1.0;
  for(size_t i=0; i<b.size(); ++i) b[i] = 1.0;
  
  Eigen::Matrix<float,10,10> mat1 = a.asDiagonal();
  mat1.diagonal<1>() = b;
  std::cout << mat1 << std::endl;
}


