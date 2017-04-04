#include <El.hpp>


template<typename Field>
void FormAndUpdateLU( El::Int n, El::Int rank, bool print )
{
    El::Output("Testing LU updates with ",El::TypeName<Field>());
    typedef El::Base<Field> Real;
    El::Timer timer;

    // Generate a random n x n matrix A.
    El::Matrix<Field> A;
    El::Uniform( A, n, n );
    if( print )
        El::Print( A, "A" );

    // Overwrite A with its LU factorization (with partial pivoting).
    timer.Start();
    El::Permutation P;
    El::LU( A, P );
    El::Output("Initial LU: ",timer.Stop()," [sec]");
    if( print )
        El::Print( A, "In-place LU of A" );

    // Update the LU factors of A to those of A + 2 X Y^H.
    El::Matrix<Field> X, Y;
    El::Uniform( X, n, rank );
    El::Uniform( Y, n, rank );
    if( print )
    {
        El::Print( X, "X" );
        El::Print( Y, "Y" );
    }
    timer.Start();
    const bool conjugate = true;
    El::LUMod( A, P, X, Y, conjugate );
    El::Output("LU update of rank ",rank,": ",timer.Stop()," [sec]");
    if( print )
        El::Print( A, "In-place LU of A + 2 X Y^H" );
    El::Output("");
}

int main (int argc, char *argv[]) {
  El::Environment env( argc, argv );
  El::DistMatrix<double> A(1,1);
  const std::string filename = "1ststageM.dmp";
  El::Read( A, filename);
  //FormAndUpdateLU<double>(n, rank, print );
    
}
