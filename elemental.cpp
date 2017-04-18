#include <El.hpp>

int main (int argc, char *argv[]) {
  El::Environment env( argc, argv );
  El::DistMatrix<double> A(1,1);
  El::DistMatrix<double> B(1,1);
  El::DistMatrix<double> S(1,1);
  El::DistMatrix<double> X(1,1);
  const std::string filenameA = "1ststageM.dmp";
  const std::string filenameB = "1ststageRHS.dmp";
  const std::string filenameS = "1ststageSol.dmp";
  El::Read( A, filenameA, El::BINARY, false);
  El::Read( B, filenameB, El::BINARY, false);
  El::Read( S, filenameS, El::BINARY, false);
  El::Output("Size: ",A.Height(),"x",A.Width());
  El::DistPermutation P;
  El::Timer timer;
 // El::Output("Condition number: ", El::TwoCondition(A));
  timer.Start();
  El::LU( A, P );
  El::Output("Factorized in ", timer.Stop()," [sec]");
  timer.Start();
  El::lu::SolveAfter(El::OrientationNS::NORMAL, A, P, B);
  El::Output("Solved in ", timer.Stop()," [sec]");
  El::DistMatrix<double> S_T(1,1);
  El::Transpose(S,S_T);
  B-=S_T;
  El::Output("Difference to stored solution: ", El::FrobeniusNorm(B));
}
