#include"particle.hpp"

class hamiltonian
{
public:
  
  //---------------------------------------------------------------
  //                          Fields 
  //---------------------------------------------------------------

  //--- constants for force computations ---
  int    np;         // total number of particle = I.Ndim^2
  int    nb_sys;     // number of replicas simulated in parallel
  double dcut;       // cut-off radius for WCA potential
  double sig, eps;   // parameters for WCA potential
  double h, w;       // parameters for dimer potential
  double boxlg;      // periodicity length
  double K;          // restraining force (M-BAR)

  //--- Force matrices ---
  Matrix fx, fy;
  double act, dx, dy, dist;  // useful variables for force computations

  //--- Creation: initialize some fields ---
  hamiltonian(const input I)
  { 
    // constants
    np     = (int)pow(I.Ndim,2.);         
    nb_sys = I.nb_replicas;        
    dcut   = pow(2.,1./6)*I.sig;
    eps    = I.eps;
    sig    = I.sig;
    w      = I.w;
    h      = I.h;
    boxlg  = I.Ndim*I.a;
    K      = I.K;
    // matrices
    fx.set_size(nb_sys,np);
    fx.zeros();
    fy.set_size(nb_sys,np);
    fy.zeros();
  }
  
  //---------------------------------------------------------------
  //                        Functions 
  //---------------------------------------------------------------

  //--- Kinetic energy ---
  void add_k_energy(particle&);

  //--- Elementary forces ---
  void   lengthBC (double,double,double,double);
  double WCA      (double);
  double fWCA     (double);
  double DW       (double);
  double fDW      (double);
  
  //--- global force computations ---
  void   compute_force(particle&);
  double RC           (double);

  //--- restraint force ---
  double restraining_potential (double,double);
  double restraining_force     (double,double);
  void   add_restraint_force   (particle&,double);

};




