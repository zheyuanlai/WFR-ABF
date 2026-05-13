//------------------------------------------
//
//          Simple vector class
//
//-------------------------------------------

#include "vector.hpp"

//------------ creator ------------
Vector::Vector()
{
  longueur = 0;
  coord = NULL;
}

//------------ destructor ------------
Vector::~Vector()
{
  if (coord!=NULL) 
    delete[] coord;
  coord = NULL;
  longueur = 0;
}

//------------ initialization ------------
void Vector::set_size(int Longueur)
{
  longueur = Longueur;
  if (coord != NULL) 
    delete[] coord;
  coord = new double[longueur];
}

Vector::Vector(int Longueur)
{
  longueur = Longueur;
  coord = new double[longueur];
}

//-------- creation of the same vector ------------
Vector & Vector::operator= (const Vector & a)
{
  if (this!=&a) 
    {
      longueur = a.longueur;
      if (coord != NULL) 
	delete[] coord;
      if (a.coord != NULL) 
	{
	  coord = new double[longueur];
	  for (int i=0;i<longueur;i++) 
	    coord[i] = a.coord[i];
	} else {
	  coord = NULL;
	}
    }
  return *this;
}

//------------ set to zero ------------
void Vector::zeros()
{
  if (coord!=NULL) 
    {
      for (int i=0; i<longueur; i++) 
	coord[i] = 0.0;
    }
}

//------------ access an element ------------
double & Vector::operator() (int a)
{
  if (coord != NULL) 
    {
      return coord[a];
    } else {
      cerr << " ### ERROR ### requiring the component of an empty vector! " << endl;
      exit(1);
    }
}
