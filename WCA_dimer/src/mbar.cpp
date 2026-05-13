#include "algorithm.hpp"

//-----------------------------
//      Initialization
//----------------------------

void algorithm::Initialize()
{
  //--- Initialize the fields ---
  X.qx.set_size(nb_replicas,np);
  X.qx.zeros();
  X.qy.set_size(nb_replicas,np);
  X.qy.zeros();
  X.px.set_size(nb_replicas,np);
  X.px.zeros();
  X.py.set_size(nb_replicas,np);
  X.py.zeros();
  X.RFx.set_size(nb_replicas,np);
  X.RFx.zeros();
  X.RFy.set_size(nb_replicas,np);
  X.RFy.zeros();
  X.energy.set_size(nb_replicas);
  X.energy.zeros();
  X.potential_energy.set_size(nb_replicas);
  X.potential_energy.zeros();

  //--------- Loop over the replicas ---------
  for (int k = 0; k < nb_replicas; k++)
    {
      //----  Initial positions: on a lattice ---
      for (int i = 0; i < Ndim; ++i)  
	{
	  for (int j = 0; j < Ndim; ++j)  
	    {
	      X.qx(k,i*Ndim+j) = (0.5 + i)*a;
	      X.qy(k,i*Ndim+j) = (0.5 + j)*a;
	    }
	}
      //--- Set dimer bond length to equilibrium ---
      X.qy(k,1) = X.qy(k,0) + H->dcut;
      //--- Momenta and random forces ---
      for (int i = 0; i < np; ++i)
	{
      	  X.px(k,i)  = sqrt(1./beta)*g.rand_number(0,1);
      	  X.py(k,i)  = sqrt(1./beta)*g.rand_number(0,1);
      	  X.RFx(k,i) = g.rand_number(0,1);
      	  X.RFy(k,i) = g.rand_number(0,1);
      	}
    }	 
  
  //------------------- Screen output --------------------------
  cout << endl;
  cout << "------------------------------------------" << endl;
  cout << "       DIMER IN WCA SOLVENT : M-BAR       " << endl;
  cout << "------------------------------------------" << endl;
  cout << endl;
  cout << "        Initialization performed          " << endl;
  cout << "                                          " << endl;
  cout << "-----  Parameters of the computation -----" << endl;
  cout << endl;
  cout << " Dynamics                       : LANGEVIN "      << endl;
  cout << " Time step                 (dt) : " << dt         << endl;
  cout << " Number of configs / restraint  : " << n_iter     << endl;
  cout << " Thermalization steps           : " << n_iter_thm << endl;
  cout << " Subsampling rate               : " << freq       << endl; 
  cout << " Friction coefficient      (xi) : " << xi         << endl; 
  cout << " Inverse temperature     (beta) : " << beta       << endl; 
  cout << " Total number of particles (np) : " << np         << endl; 
  cout << " Elementary cell size       (a) : " << a          << endl;
  cout << " Particle density               : " << 1./(a*a)   << endl;
  cout << " WCA equilibrium distance (sig) : " << sig        << endl; 
  cout << " WCA energy epsilon       (eps) : " << H->eps     << endl; 
  cout << " Dimer barrier height       (h) : " << H->h       << endl; 
  cout << " Dimer bond length          (w) : " << w          << endl; 
  cout << " Restraining potential      (K) : " << K          << endl; 
  cout << " Lower restraint center  (zmin) : " << zmin       << endl; 
  cout << " Upper restraint center  (zmax) : " << zmax       << endl; 
  cout << " Number of restraint cent. (Nz) : " << Nz         << endl; 
  cout << " Lower value of RC     (xi_min) : " << xi_min     << endl; 
  cout << " Upper value of RC     (xi_max) : " << xi_max     << endl;
  cout << " Number of PMF RC values  (Nxi) : " << Nxi        << endl; 
  cout << " Tolerance for SC convergence   : " << tol        << endl; 
  cout << endl;
  cout << "----------------------------------------- " << endl;

  cout << endl;
      
}

//---------------------------------------------------------
//   Integration of the dynamics and auxiliary functions
//---------------------------------------------------------

//--- value of the reaction coordinate, given the bond length ---
double algorithm::RC(int dumb)   // argument not used, only one replica
{
  double value;
  value = (dist-dcut)/(2*w);  
  return value;
}

//-------- Return length of the vector, up to periodic BC ----------
void algorithm::length(double x1, double y1, double x2, double y2)
{
  dx = x1 - x2;
  dx -= boxlg*rint(dx*1./boxlg);
  dy = y1 - y2;
  dy -= boxlg*rint(dy*1./boxlg);
  dist = dx*dx + dy*dy; 
  dist = sqrt(dist);
}

//-------- Apply periodic boundary conditions ---------
void algorithm::Periodic()
{
  Matrix temp;

  //-- x direction ---
  temp = X.qx;
  for (int k = 0; k < nb_replicas; k++)
    for (int i = 0; i < np; i++)
      temp(k,i) -= boxlg*floor( temp(k,i)*1./boxlg );
  X.qx = temp;

  //--- y direction
  temp = X.qy;
  for (int k = 0; k < nb_replicas; k++)
    for (int i = 0; i < np; i++) 
      temp(k,i) = temp(k,i) - boxlg*floor( temp(k,i)*1./boxlg );
  X.qy = temp;

}

//------------ Langevin dynamics -----------------------
void algorithm::Langevin_restraint(double center)
{
  double sigma = sqrt(2*xi/beta);
  
  //--- BBK like algorithm ---
  for (int k = 0; k < nb_replicas; k++)
    {
      //--- Random terms
      for (int i = 0; i < np; ++i)
	{
	  X.RFx(k,i) = g.rand_number(0,1);
	  X.RFy(k,i) = g.rand_number(0,1);
	} 
      //--- random part of the force ---
      for (int i = 0; i < np; ++i)
	{
	  H->fx(k,i) += - xi*X.px(k,i) + sqrt(1/dt)*sigma*X.RFx(k,i);
	  H->fy(k,i) += - xi*X.py(k,i) + sqrt(1/dt)*sigma*X.RFy(k,i);
	} 
      for (int i = 0; i < np; ++i)
	{
	  //--- momenta update ---
	  X.px(k,i) += 0.5*dt*H->fx(k,i);
	  X.py(k,i) += 0.5*dt*H->fy(k,i);
	  //--- position update
	  X.qx(k,i) += dt*X.px(k,i);
	  X.qy(k,i) += dt*X.py(k,i);
	} 
    }
  //--- Apply periodic BC on positions --- 
  Periodic();
  //--- Update the forces ---
  H->compute_force(X);
  H->add_restraint_force(X,center);
  //--- Final update of the momenta ---
  for (int k = 0; k < nb_replicas; k++)
    {
      for(int i = 0 ; i < np; i++)
	{
	  X.px(k,i) = 1/( 1 + xi*dt*0.5 )
	    *( X.px(k,i) + H->fx(k,i)*dt*0.5 + 0.5*sqrt(dt)*sigma*X.RFx(k,i) );
	  X.py(k,i) = 1/( 1 + xi*dt*0.5 )
	    *( X.py(k,i) + H->fy(k,i)*dt*0.5 + 0.5*sqrt(dt)*sigma*X.RFy(k,i) );
	} 
    }
  //--- Add kinetic energy to the total energy ---
  H->add_k_energy(X);
  
  //--- End of Langevin ---
}

//-----------------------------------------------------------------
//
//   Compute the trajectory + M-BAR SCF loop + PMF approximation
//
//-----------------------------------------------------------------

void algorithm::load()
{
  
  //--------------- Useful fields ---------------------
  ofstream ENERGY ("../output/energy");   // energy as a function of time
  ofstream COORD  ("../output/coord");    // reaction coordinate vs. time
  ofstream RATIO  ("../output/ratio");    // ratios of partition functions in M-BAR loop
  ofstream PMF    ("../output/PMF");      // PMF approximation at chosen points
  //--- Fields for sampling procedure ---
  double dz;
  int Nsample = Nz*n_iter;
  Vector Coord, Energy;
  Coord.set_size(Nsample);
  Coord.zeros();
  Energy.set_size(Nsample);
  Energy.zeros();
  //--- Fields for SCF loop in M-BAR ---
  double Dist = 1;
  Vector denom; 
  denom.set_size(Nsample);
  denom.zeros();
  Vector z, z_old;
  z.set_size(Nz);
  z_old.set_size(Nz);
  Matrix Q; 
  Q.set_size(Nz,Nsample);
  Q.zeros();
  //---- Fields for PMF computation ---
  double DeltaXi = (xi_max-xi_min)/Nxi;
  double Xi, value;
  double top, bottom;

  //------------- Initialization ----------------------
  Initialize();
  H->compute_force(X);
  cout << "------- Creation of statistical sample -------" << endl;
  cout << endl;
  
  //---- center positions ----
  Vector centers;
  centers.set_size(Nz);
  dz = (zmax - zmin) / (Nz-1);
  for (int i = 0; i < Nz; i++)
    centers(i) = zmin + i*dz;
  
  //-----------------------------------------------------
  //            CREATION OF STATISTICAL SAMPLE
  //----------- Loop over the restraining centers -------
  for (int i = 0; i < Nz; i++)
    {
      cout << "  Restraint potential centered on " << centers(i) << endl;
      
      //--- thermalization ---
      H->compute_force(X);
      H->add_restraint_force(X,centers(i));
      for(int k = 0; k < n_iter_thm; k++)
	Langevin_restraint(centers(i));
      
      //--- Loop over time ---
      for(int Niter = 0; Niter < n_iter; Niter++)
	{
	  //---- Subsamling at rate 'fres' ----
	  for(int k = 0; k < freq; k++)
	    Langevin_restraint(centers(i));
	  //---- keep the corresponding configuration ----
	  length(X.qx(0,0),X.qy(0,0),X.qx(0,1),X.qy(0,1));
	  Energy(i*n_iter+Niter) = X.energy(0);
	  Coord (i*n_iter+Niter) = RC(0);
	  ENERGY << Energy(i*n_iter+Niter) << " ";
	  COORD  << Coord (i*n_iter+Niter) << " ";
	}
      
      //--- end of outputs ---
      ENERGY << endl;
      COORD << endl;
    }
  
  //-----------------------------------------------------
  //         ESTIMATION OF FREE ENERGY DIFFERENCES
  //----- Self-consistent procedure based on data -------
  cout << endl;
  cout << "----- Computation of free energy differences ------" << endl;
  cout << endl;
  
  //--- Gather the weights f_k(x^i) in a matrix ---
  for (int i = 0; i < Nsample; i++)
    for(int k = 0; k < Nz; k++)
      Q(k,i) = exp( -beta*( Energy(i) + H->restraining_potential(Coord(i),centers(k)) ) );
  //--- Initialization of the ratios of partition functions ---
  // could start from more educated initial guess, see e.g. [Shirts/Chodera,2008], Eq. (C4) 
  for (int k = 0; k < Nz; k++)
    z(k) = 1;   
  //--- Self consistent loop ---
  while (Dist > tol)
    {
      //--- Precomputations ---
      for (int i = 0; i < Nsample; i++)
	{
	  denom(i) = 0;
	  for (int k = 0; k < Nz; k++)
	    denom(i) += Q(k,i) / z(k);  
	  // CAUTION -- This last update should be modified accordingly 
	  // if the samples for a given restraining center have different sizes
	}
      //--- Obtain new estimate for each normalization factor ---
      for (int k = 0; k < Nz; k++)
	{
	  z_old(k) = z(k);
	  z(k) = 0;
	  for (int i = 0; i < Nsample; i++)
	    z(k) += Q(k,i)/denom(i);
	}
      //--- Renormalization procedure ---
      double zref = z(Nz/2);
      for (int k = 0; k < Nz; k++)
	z(k) = z(k)/zref;
      //--- Update the convergence criterion ---
      Dist = 0;
      for (int k = 0; k < Nz; k++)
	{
	  if (fabs( (z(k)-z_old(k))/z_old(k) ) > Dist) 
	    Dist = fabs( (z(k)-z_old(k))/z_old(k) );
	}
      //--- Monitor the convergence of the ratios --- 
      for (int k = 0; k < Nz; k++)
	RATIO << z(k) << "  ";
      RATIO << endl;
      //--- end of SCF loop ---
    }
  //--- Final output of ratios of partition functions ----
  for (int k = 0; k < Nz; k++)
    RATIO << z(k) << "  ";
  RATIO << endl;


  //-----------------------------------------------------
  //                 COMPUTATION OF PMF
  //----- Using the computed free energy differences ----
  cout << "------- Computation of PMF -------" << endl;
  cout << endl;
  //---  Some precomputations ---
  for (int i = 0; i < Nsample; i++)
    {
      denom(i) = 0;
      for (int k = 0; k < Nz; k++)
	denom(i) += Q(k,i) / z(k);  
    }
  //--- Loop over the required values of PMF ---
  for (int s = 0; s < Nxi+1; s++)
    {
      Xi = xi_min + s*DeltaXi;
      top = 0;
      bottom = 0;
      //--- Compute an approximation of the average of indicator function ---
      for (int i = 0; i < Nsample; i++)
	{
	  bottom += exp(-beta*Energy(i))/denom(i);
	  if ( fabs(Coord(i)-Xi) < 0.5*DeltaXi )
	    top += exp(-beta*Energy(i))/denom(i);
	}
      value = top/bottom;
      value = -log(value/DeltaXi);
      PMF << Xi << " " << value/beta << endl;
    }

  //------------- End of computation ---------------
  cout << endl;
  cout << "------- End of computation ------" << endl;
  cout << endl;
  cout << endl;
}









