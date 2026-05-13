//------------------------------------------
//
//          Simple matrix class
//
//-------------------------------------------

//--- required libraries ---
#include "input.hpp"

//------- structure --------
class Matrix
{
public:

  //--- fields ---
  int row, col;
  double * coord;

  //--- constructors, destructors ---
  Matrix();
  ~Matrix();
  Matrix(int,int);
  double & operator() (int,int);
  void set_size(int,int);
  void zeros();
  
  // overloaded operations
  Matrix & operator= (const Matrix &);
  Matrix operator+ (const Matrix &);
  Matrix operator- (const Matrix &);
  Matrix operator* (const Matrix &);
  Matrix operator* (const double); 
  
};


