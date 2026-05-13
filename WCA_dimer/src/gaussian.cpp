//------------------------------------------
//
//  Functions for sampling random numbers
//  (See 'main.cpp' for the initial seed)
//
//-------------------------------------------

#include "gaussian.hpp"

//--- global fields ---
double r1 = 0;
double r2 = 0;
 
//-- uniform random number generator --
double gaussian::random_number()
{
  return ((double)rand()/RAND_MAX);
}
 
//-- gaussian random number, mean m, variance \sigma^2 --
// (method of change of variables to polar coordinates)
double gaussian::rand_number(double m, double sigma)
{
  double temp1, temp2;
  if (random_flag)
    {
      temp1 = random_number();
      temp2 = random_number();
      if (temp1 == 0) 
	temp1 = 1e-9;
      r1 = m + sigma*(sqrt(-2.*log(temp1))*cos(2.*M_PI*temp2));
      r2 = m + sigma*(sqrt(-2.*log(temp1))*sin(2.*M_PI*temp2));
      random_flag = 0;
      return r1;
    }
  else
    {
      random_flag = 1;
      return r2;
    }
}
 


