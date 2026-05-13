//------------------------------------------
//
//          Simple vector class
//
//-------------------------------------------

//--- required libraries ---
#include "input.hpp"

//--- definition of the structure ---
class Vector
{
public:

  //--- fields ---
  int longueur;
  double * coord;

  //--- constructors, destructors, operations ---
  Vector();
  ~Vector();
  Vector(int);
  double & operator() (int);
  void set_size(int);
  void zeros();
  Vector & operator= (const Vector &);
 
};


