//-----------------------------------------------------
//
//     Reading parameters from "input_file_..."
//
//-----------------------------------------------------


#ifndef INPUT_HPP
#define INPUT_HPP

#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <string>
#include <cstring>
#include <cstdlib>

using namespace std;

//--------------- definition of the structure --------
class input
{
public:
 
  //------------------------ FIELDS --------------------------
  
  //--- Simulation parameters ---
  int    method;         // chosen method 
  double t_step;         // time step for the time integrator algorithm
  int    nb_iter;        // number of iterations
  int    nb_iter_thm;    // preliminary thermalization (when needed)
  int    freq_xmakemol;  // XMakeMol output frequency
  int    freq;           // outputs written every 'freq' steps (energy, reaction coordinate, ...)
  double xi;             // friction coefficient for Langevin dynamics
  int    nb_replicas;    // number of replicas
  double beta;        // inverse temperature

  //--- Thermodynamic conditions ---
  int    Ndim;        // Ndim^2 particles -> 2 for dimer, Ndim^2 - 2 for solvent  
  int    Nsolvent;    // number of solvent particles = Ndim^2 - 2
  double a;           // domain [0,Ndim*a]^2, volume = (Ndim*a)^2

  //--- Potential parameters ---
  double h, w, sig, eps;

  //--- M-BAR method ---
  double K;               // restraining potential
  double zmin, zmax;      // upper and lower bounds of the restraining pot. centers
  int    Nz;              // number of restraining centers
  double xi_min, xi_max;  // range of values of reaction coordinate
  int    Nxi;             // number of points where PMF computed  
  double tol;             // tolerance for self-consistent convergence

  //--- ABF method ---
  int    freq_histo;      // frequency for bias / mean force / values of RC
  int    selection;       // selection to enhance diffusion
  double c_selection;     // selection intensity

  //---- Jarzynski -----
  double Time;            // switching time

  //---------------------- READING FUNCTIONS --------------------------
  
  template <class T>
  //--- elementary reading function -------
  void read_item(istream& ,char *, T *); 
  //--- copying datas from input_file to the parameters class ---
  void load(void);
  //--- reading the input file for "sampling" ------
  void read_sampling(int&,double&,int&,int&,int&,double&,double&,int&,double&,double&,double&,double&,double&);
   //--- reading the input file for "mbar" ------
  void read_mbar(double&,int&,int&,int&,double&,double&,int&,double&,double&,double&,double&,double&,double&,double&,double&,int&,double&,double&,int&,double&);
  //--- reading the input file for "ABF" ------
  void read_ABF(int&,int&,double&,int&,int&,int&,double&,double&,int&,double&,double&,double&,double&,double&,double&,double&,int&,int&,double&);
  //--- reading the input file for "TI" ------
  void read_TI(int&,double&,int&,int&,int&,double&,double&,int&,double&,double&,double&,double&,double&,double&,double&,int&);
  //--- reading the input file for "JARZ OVD" ------
  void read_jarz(int&,double&,int&,int&,int&,double&,double&,int&,double&,double&,double&,double&,double&,double&,double&,double&,int&);
  
};

#endif

