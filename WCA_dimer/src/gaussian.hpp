//------------------------------------------
//
//  Functions for sampling random numbers
//
//-------------------------------------------

//---- required libraries ------
#include <iostream>
#include <math.h>
#include <stdlib.h>

//----- definition of the structure
class gaussian
{
  int random_flag;

public :
  //-- sampling a gaussian random number --
  double rand_number(double,double);
  //-- sampling uniform law --
  double random_number();

};




