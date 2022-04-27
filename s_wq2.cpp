#include <iostream>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/LU>

using namespace Eigen;

void solve_linear_system()
{
  std::cout << "[[[[[solve_linear_system]]]]]" << std::endl;
  MatrixXd A = MatrixXd::Random(100,100);
  FullPivLU< MatrixXd > lu(A);
 
  VectorXd b = VectorXd::Random(100);
  VectorXd x;
 
  std::cout << "Matrix A" << std::endl << A << std::endl;
  std::cout << "Vector b" << std::endl << b << std::endl << std::endl;
   
  x=lu.solve(b);
  std::cout << "Solve Ax=b:" << std::endl
        << x << std::endl << std::endl;
  std::cout << "Check solution by Ax-b: " << std::endl
        << A*x-b << std::endl << std::endl;
}
 
int main()
{
  solve_linear_system();
}