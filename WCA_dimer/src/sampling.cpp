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
  cout << "     DIMER IN WCA SOLVENT : SAMPLING      " << endl;
  cout << "------------------------------------------" << endl;
  cout << endl;
  cout << "        Initialization performed          " << endl;
  cout << "                                          " << endl;
  cout << "-----  Parameters of the computation -----" << endl;
  cout << endl;
  cout << " Dynamics                       : ";
  switch (dynamics) {
  case 0:
    cout << "OVERDAMPED " << endl;
    break;
  case 1:
    cout << "LANGEVIN " << endl;
    break;
  }
  cout << " Time step                 (dt) : " << dt            << endl;
  cout << " Number of iterations           : " << n_iter        << endl;
  cout << " Xmakemol output every...steps  : " << freq_xmakemol << endl; 
  cout << " Other outputs   every...steps  : " << freq          << endl; 
  cout << " Friction coefficient      (xi) : " << xi            << endl; 
  cout << " Inverse temperature     (beta) : " << beta          << endl; 
  cout << " Total number of particles (np) : " << np            << endl; 
  cout << " Elementary cell size       (a) : " << a             << endl;
  cout << " Particle density               : " << 1./(a*a)      << endl;
  cout << " WCA equilibrium distance (sig) : " << sig           << endl; 
  cout << " WCA energy epsilon       (eps) : " << H->eps        << endl; 
  cout << " Dimer barrier height       (h) : " << H->h          << endl; 
  cout << " Dimer bond length          (w) : " << w             << endl; 
  cout << endl;
  cout << "----------------------------------------- " << endl;

  cout << endl;
      
}

//---------------------------------------------------------
//   Integration of the dynamics and auxiliary functions
//---------------------------------------------------------

//--- value of the reaction coordinate, given the bond length ---
double algorithm::RC(int dumb)  // argument not used, only one replica
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
void algorithm::Langevin()
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

//---------------- Overdamped Langevin dynamics --------------
void algorithm::Overdamped()
{
  double sigma = sqrt(2*dt/beta);
  
  //--- Force computation ---
  H->compute_force(X);
  //--- Loop over the replicas ---
  for (int k = 0; k < nb_replicas; k++)
    {
      for (int i = 0; i < np; i++)
      	{
	  X.RFx(k,i) = g.rand_number(0,1);
	  X.RFy(k,i) = g.rand_number(0,1);
	  X.qx(k,i) += dt*H->fx(k,i) + sigma*X.RFx(k,i);
	  X.qy(k,i) += dt*H->fy(k,i) + sigma*X.RFy(k,i);
	} 
    }
  //--- Apply periodic BC ----
  Periodic();
  
  //--- End of overdamped Langevin ---
  
}


//----------------------------------------------------------
//
//     Sampling procedure = compute the trajectory
//
//----------------------------------------------------------

void algorithm::load()
{
  
  //--------------- Useful fields ---------------------
  ofstream ENERGY   ("../output/energy");        // energy as a function of time
  ofstream XMAKEMOL ("../output/xmakemol.xyz");  // XMakeMol output 
  ofstream COORD    ("../output/coord");         // reaction coordinate vs. time
  double time;
  
  //------------- Initialization ----------------------
  Initialize();
  H->compute_force(X);
  int acquisition = -1;
  int acquisition_xmakemol = -1;
  cout << "------- Computation started -------" << endl;
  cout << endl;
  
  //-------------- Loop over time ---------------------
  for(int Niter = 0; Niter < n_iter+1; Niter++)
    {
      //---- Chosen dynamics ----
      switch (dynamics) {
      case 0: 
	Overdamped();
      case 1:
	Langevin();
	break;
      }
      
      //--- Energy, reaction coordinate ---
      acquisition += 1;
      if (acquisition == freq)
	{
	  acquisition = 0;
	  time = Niter*dt;
	  length(X.qx(0,0),X.qy(0,0),X.qx(0,1),X.qy(0,1));
	  ENERGY << time << " " << X.energy(0) << " " << X.potential_energy(0) << " " << endl;
	  COORD  << time << " " << RC(0) << endl;
	}
      
      //---- Positions (read with XMakeMol software) ------
      acquisition_xmakemol += 1;
      if (acquisition_xmakemol == freq_xmakemol)
	{
	  cout <<  "  Output at step " << Niter << endl;
	  acquisition_xmakemol = 0;
	  XMAKEMOL << np << endl;
	  XMAKEMOL << "Time = " << time << endl;
	  //----- Solvent molecules, domain [-Na/2,Na/2]^2 ----
	  for (int i = 2; i < np; ++i)
	    {
	      if ( X.qx(0,i) > boxlg*0.5) 
		X.qx(0,i) += -boxlg;
	      if ( X.qy(0,i) > boxlg*0.5) 
		X.qy(0,i) += -boxlg;
	      XMAKEMOL << "H " << X.qx(0,i) << "  " << X.qy(0,i) << "  " << 0 << endl;
	    }
	  //---- Dimer: color depends on the bond length ----
	  if (RC(0) < 0.5)
	    {
	      for (int i = 0; i < 2; ++i)
		{
		  if ( X.qx(0,i) > boxlg*0.5) 
		    X.qx(0,i) += -boxlg;
		  if ( X.qy(0,i) > boxlg*0.5) 
		    X.qy(0,i) += -boxlg;
		  XMAKEMOL << "O " << X.qx(0,i) << "  " << X.qy(0,i) << "  " << 0 << endl;
		}
	    } else {
	      for (int i = 0; i < 2; ++i)
		{
		  if ( X.qx(0,i) > boxlg*0.5) 
		    X.qx(0,i) += -boxlg;
		  if ( X.qy(0,i) > boxlg*0.5) 
		    X.qy(0,i) += -boxlg;
		  XMAKEMOL << "B " << X.qx(0,i) << "  " << X.qy(0,i) << "  " << 0 << endl;
		}
	    }
	}     
    }
  
  //------------- End of computation ---------------
  cout << endl;
  cout << "------- End of computation ------" << endl;
  cout << endl;
  cout << endl;
}









