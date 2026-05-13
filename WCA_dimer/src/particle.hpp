#include "input.hpp"
#include "vector.hpp"
#include "matrix.hpp"

class particle
{
public:
  
  //---------------------------------
  //  All the masses are set to 1 !
  //---------------------------------
  
  //-------------------------------------------------
  //  General structure of the fields :  A(k,i)
  //  with k = index of the replica (0 -> nb_simu-1)
  //       i = particle number for a given replica
  //       (0,1 = dimer, 2 -> N-1 = solvent)
  //  For vectors : V(k) -> k = replica index
  //-------------------------------------------------

  Matrix qx, qy;     // positions 
  Matrix px, py;     // momenta
  Matrix RFx, RFy;   // noise term
  
  Vector energy;           // total energy  
  Vector potential_energy; // potential energy 

  Vector Birth, Death;       // birth and death times for selection process
  double nb_birth, nb_death; // total number of deaths/births

  Vector Work, Work2;    // for nonequilibrium switching 

  Matrix old_qx, old_qy, old_px, old_py; // previous configurations when Metropolizing
  Vector old_energy, old_potential_energy;

  particle(){};
  particle(const input& I);

};


