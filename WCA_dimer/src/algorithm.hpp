#ifndef ALGORITH_HPP
#define ALGORITH_HPP

#include"gaussian.hpp"
#include"hamiltonian.hpp"

class algorithm
{
public:

  //---------------------------------------------------------------
  //                          Fields 
  //---------------------------------------------------------------
  
  //--- structures ---
  gaussian g;
  particle X;
  hamiltonian* H;
  
  //--- useful constants ---
  double dt, beta, a, boxlg, xi, sig, w, dcut;
  int nb_replicas, np, Ndim;
  int n_iter, n_iter_thm, freq, freq_xmakemol;
  int dynamics;
  double dx, dy, dist;
  
  //--- for reaction coordinate ---
  double xi_min, xi_max, current_xi;
  int    Nxi;

  //--- specific M-BAR ---
  double K, zmin, zmax, tol;
  int    Nz;
  
  //--- specific ABF ---
  Matrix AverageForce;  
  Vector Marginal, MarginalUpdate;
  Vector Bias, CurrentMeanForce;
  int    freq_histo;
  double DeltaXi;
  int    selection;
  double c_selection;

  //--- specific TI ---
  Vector MeanForce, e_x, e_y;
  Matrix old_qx, old_qy;
  Matrix LagrangeMeanForce;
  double SumMeanForces, SumLagrange, SumLagrangeReduced, SumLagrangeTR;
  double SumLagrangeAveraged, SumLagrangePosition, SumLagrangeMomenta;
  double rejection;

  //--- specific Jarzynski -----
  double current_vxi;
  double Time;
  Vector Xi, Lagrange;

  //--- Creation: initialize some fields ---
  algorithm(const input I)
  { 
    nb_replicas   = I.nb_replicas;
    np            = I.Ndim*I.Ndim;
    beta          = I.beta;
    a             = I.a;
    Ndim          = I.Ndim;
    dt            = I.t_step;
    n_iter        = I.nb_iter;
    n_iter_thm    = I.nb_iter_thm;
    boxlg         = I.Ndim*I.a;
    xi            = I.xi;
    freq          = I.freq;
    freq_xmakemol = I.freq_xmakemol;
    sig           = I.sig;
    w             = I.w;
    dynamics      = I.method;
    dcut          = pow(2,1./6)*I.sig;
    Nz            = I.Nz;
    Nxi           = I.Nxi;
    K             = I.K;
    zmin          = I.zmin;
    zmax          = I.zmax;
    xi_min        = I.xi_min;
    xi_max        = I.xi_max;
    tol           = I.tol;
    freq_histo    = I.freq_histo;
    selection     = I.selection;
    c_selection   = I.c_selection;
    Time          = I.Time;
  }

  //---------------------------------------------------------------
  //                        Functions 
  //---------------------------------------------------------------
  
  //---- Initialization ---
  void Initialize(); 
  void Periodic  ();
  
  //---- Integrators (sampling) ---
  void Overdamped();
  void Langevin  ();

  //--- Functions for M-BAR ---
  void Langevin_restraint (double);
  
  //--- Functions for ABF ---
  int  index_RC       (double);
  void mean_force     ();  
  void ApplyBias      ();
  void Langevin_ABF   ();
  void Overdamped_ABF ();
  void Select         ();

  //--- Functions for TI ---
  void Langevin_TI         ();
  void Overdamped_TI       ();
  void mean_force_Langevin ();  

  //------- Functions for Jarzynski -------
  double schedule               (double);
  void   Overdamped_constrained ();
  void   Langevin_constrained   ();

  //--- Main function ---
  void load ();

  //--- Auxiliary functions ---
  void   length   (double,double,double,double);
  double RC       (int);
   
};
#endif
