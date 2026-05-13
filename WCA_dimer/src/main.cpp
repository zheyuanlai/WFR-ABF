#include"algorithm.hpp"

int main(void)
{

  //---------- initialization of random number generator ----
  srand(time(NULL));
  
  //---------- reading data from "input_file" --------------
  input I;
  I.load();
  
  //----------- creating structures ------------------------
  algorithm* A;
  A = new algorithm(I);
  A->H = new hamiltonian(I);
  
  //----------- effective computation ------------------------
  A->load();
  
  return EXIT_SUCCESS;
}
