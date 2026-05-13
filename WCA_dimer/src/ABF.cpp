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

  //------ Bias ------
  DeltaXi = (xi_max-xi_min)/Nxi;
  AverageForce.set_size(Nxi,2);
  AverageForce.zeros();
  Bias.set_size(Nxi);
  Bias.zeros();
  Marginal.set_size(Nxi);
  Marginal.zeros();
  MarginalUpdate.set_size(Nxi);
  MarginalUpdate.zeros();
  CurrentMeanForce.set_size(nb_replicas);
  CurrentMeanForce.zeros();

  //------ Selection -------
  X.Birth.set_size(nb_replicas);
  X.Birth.zeros();
  X.Death.set_size(nb_replicas);
  X.Death.zeros();
  X.nb_death = 0.;
  X.nb_birth = 0.;

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
      //--- birth and death times ---
      if (selection == 1)
	{
	  X.Birth(k) = -log( g.random_number() );
	  X.Death(k) = -log( g.random_number() );
	}
    }	 
  
  //------------------- Screen output --------------------------
  cout << endl;
  cout << "------------------------------------------" << endl;
  cout << "       DIMER IN WCA SOLVENT : ABF         " << endl;
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
  cout << " Number of replicas             : " << nb_replicas   << endl;
  cout << " Number of iterations           : " << n_iter        << endl;
  cout << " Histrograms   every...steps    : " << freq_histo    << endl; 
  cout << " Other outputs every...steps    : " << freq          << endl; 
  cout << " Friction coefficient      (xi) : " << xi            << endl; 
  cout << " Inverse temperature     (beta) : " << beta          << endl; 
  cout << " Total number of particles (np) : " << np            << endl; 
  cout << " Elementary cell size       (a) : " << a             << endl;
  cout << " Particle density               : " << 1./(a*a)      << endl;
  cout << " WCA equilibrium distance (sig) : " << sig           << endl; 
  cout << " WCA energy epsilon       (eps) : " << H->eps        << endl; 
  cout << " Dimer barrier height       (h) : " << H->h          << endl; 
  cout << " Dimer bond length          (w) : " << w             << endl; 
  cout << " Lower value of RC     (xi_min) : " << xi_min     << endl; 
  cout << " Upper value of RC     (xi_max) : " << xi_max     << endl;
  cout << " Number of PMF RC values  (Nxi) : " << Nxi        << endl; 
  if (selection == 1)
    {
      cout << " SELECTION ON   " << endl;
      cout << " Selection intensity            : " << c_selection << endl; 
    }
  cout << endl;
  cout << "----------------------------------------- " << endl;

  cout << endl;
      
}

//---------------------------------------------------------
//   Integration of the dynamics and auxiliary functions
//---------------------------------------------------------

//--- value of the reaction coordinate, given the bond length ---
double algorithm::RC(int k)
{
  length(X.qx(k,0),X.qy(k,0),X.qx(k,1),X.qy(k,1));
  double value = (dist-dcut)/(2*w);  
  return value;
}

//------- Index of the RC in the discretization ----
int algorithm::index_RC(double RCvalue)
{
  int index;
  index = (int)floor( (RCvalue-xi_min)/DeltaXi );
  if ( (index < 0) || (index > Nxi-1) )
    index = -1;
  return index;
}

//--------- Mean force computation -------
void algorithm::mean_force()
{
  for (int k = 0; k < nb_replicas; k++)
    {
      //--- recompute the length and unit vectors ---
      length(X.qx(k,0),X.qy(k,0),X.qx(k,1),X.qy(k,1));
      //--- formula for the mean force (force already up-to-date) ---
      CurrentMeanForce(k) = (H->fx(k,1)-H->fx(k,0))*dx 
	+ (H->fy(k,1)-H->fy(k,0))*dy - 2./beta;
      CurrentMeanForce(k) *= w/dist;
    }
}

//------ Apply the biasing force -----
void algorithm::ApplyBias()
{
  int index; 
  double BiasingForce, BiasingForce_x, BiasingForce_y;
  for (int k = 0; k < nb_replicas; k++)
    {
      index = index_RC( RC(k) );
      if (index > -1)
	{
	  if (AverageForce(index,0) > 0)
	    BiasingForce = AverageForce(index,1)/AverageForce(index,0)/(2*w);
	  else 
	    BiasingForce = 0;
	  BiasingForce_x = BiasingForce * dx/dist;
	  BiasingForce_y = BiasingForce * dy/dist;
	  H->fx(k,0) += BiasingForce_x;
	  H->fx(k,1) -= BiasingForce_x;
	  H->fy(k,0) += BiasingForce_y;
	  H->fy(k,1) -= BiasingForce_y;
	}
    }
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
void algorithm::Langevin_ABF()
{
  double sigma = sqrt(2*xi/beta);
  int index;
  //--- BBK like algorithm ---
  for (int k = 0; k < nb_replicas; k++)
    {
      //--- Random terms
      for (int i = 0; i < np; ++i)
	{
	  X.RFx(k,i) = g.rand_number(0,1);
	  X.RFy(k,i) = g.rand_number(0,1);
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
  //------ Update the mean force (Force already up-to-date -------
  mean_force();
  MarginalUpdate.zeros();
  for (int k = 0; k < nb_replicas; k++)
    {
      index = index_RC( RC(k) );
      AverageForce(index,0) += 1;
      AverageForce(index,1) += CurrentMeanForce(k);
      MarginalUpdate(index) += 1;
    }
  //--- Apply now the bias (Beware: it modifies the force!) ---
  ApplyBias(); 
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
  //--- End of ABF Langevin ---
}

//---------------- Overdamped Langevin dynamics --------------
void algorithm::Overdamped_ABF()
{
  double sigma = sqrt(2*dt/beta);
  int index;

  //--- Force computation ---
  H->compute_force(X);
  //------ Update the mean force (Force already up-to-date -------
  mean_force();
  MarginalUpdate.zeros();
  for (int k = 0; k < nb_replicas; k++)
    {
      index = index_RC( RC(k) );
      AverageForce(index,0) += 1;
      AverageForce(index,1) += CurrentMeanForce(k);
      MarginalUpdate(index) += 1;
    }
  //--- Apply now the bias (Beware: it modifies the force!) ---
  ApplyBias();
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

//---------------- Selection procedure ------------
void algorithm::Select()
{
  double diff = 0;
  int    ind;
  //------ Update the marginal ---
  for (int i = 0; i < Nxi; i++)
    Marginal(i) = MarginalUpdate(i);
  //----- Update the birth/death times -----
  for (int k = 0; k < nb_replicas; k++)
    {
      ind = index_RC( RC(k) );
      if ( (ind > 0) && (ind < Nxi-2) )
	//diff = -(AverageForce(ind+1,0)-2*AverageForce(ind,0)+AverageForce(ind-1,0))
	//    / AverageForce(ind,0) / pow(DeltaXi,2);
	diff = -(Marginal(ind+1)-2*Marginal(ind)+Marginal(ind-1))
	  / Marginal(ind) / pow(DeltaXi,2);
      else 
	diff = 0;
      //-- when locally too much exploration
      //   decrease the death time, otherwise, decrease 
      //   the birth time (favor this configuration)
      if ( diff > 0)
	X.Death(k) -= c_selection*diff*dt; 
      else
	X.Birth(k) += c_selection*diff*dt; 
    }
  //----- Resampling step when some time < 0 ----
  for(int k = 0 ; k < nb_replicas; k++)
    {
      //--- Check death ---
      if (X.Death(k) < 0)  
	{
	  //--- Update total number of deaths ---
	  X.nb_death += 1;
	  //--- New exponential clock ---
	  X.Death(k) = -log( g.random_number() );
	  //--- Choose the new position, and replace ---
	  ind = (int) floor(g.random_number()*nb_replicas);
	  for (int i = 0; i < np; i++)
	    {
	      X.qx(k,i) = X.qx(ind,i);
	      X.qy(k,i) = X.qy(ind,i);
	      // p exchange useless if overdamped dynamics...
	      if (dynamics == 1)
		{
		  X.px(k,i) = X.px(ind,i);
		  X.py(k,i) = X.py(ind,i);
		}
	    }
	  X.energy(k)           = X.energy(ind);
	  X.potential_energy(k) = X.potential_energy(ind);
	}
      //--- Check birth ---
      if (X.Birth(k) < 0)  
	{
	  //--- Update total number of deaths ---
	  X.nb_birth += 1;
	  //--- New exponential clock ---
	  X.Birth(k) = -log( g.random_number() );
	  //--- Choose the new position, and replace ---
	  ind = (int) floor(g.random_number()*nb_replicas);
	  for (int i = 0; i < np; i++)
	    {
	      X.qx(ind,i) = X.qx(k,i);
	      X.qy(ind,i) = X.qy(k,i);
	      // p exchange useless if overdamped dynamics...
	      if (dynamics == 1)
		{
		  X.px(ind,i) = X.px(k,i);
		  X.py(ind,i) = X.py(k,i);
		}
	    }
	  X.energy(ind)           = X.energy(k);
	  X.potential_energy(ind) = X.potential_energy(k);
	}
    }
  //------------ End of selection procedure -----------
}

//----------------------------------------------------------
//
//     Sampling procedure = compute the trajectory
//
//----------------------------------------------------------

void algorithm::load()
{
  
  //--------------- Useful fields ---------------------
  ofstream ENERGY     ("../output/energy");        // energy as a function of time
  ofstream COORD      ("../output/coord");         // reaction coordinate at given time
  ofstream COORD1     ("../output/coord_first");   // reaction coordinate of first system
  ofstream MEAN_FORCE ("../output/mean_force");    // mean force as a function of time
  ofstream BIAS       ("../output/PMF");           // associated PMF as a function of time
  ofstream HISTO      ("../output/histo");         // cumulated histogram of visits
  double time, CurrentForce, normalization;
  int    Middle  = (int)rint(Nxi/2);

  //------------- Initialization ----------------------
  Initialize();
  H->compute_force(X);
  int acquisition = -1;
  int acquisition_histo = -1;
  for (int i = 0; i < Nxi; i++)
    {
      MEAN_FORCE << xi_min + i*DeltaXi << "  ";
      BIAS       << xi_min + i*DeltaXi << "  ";
    }
  MEAN_FORCE << endl;
  BIAS       << endl;
  cout << "------- Computation started -------" << endl;
  cout << endl;
  
  //-------------- Loop over time ---------------------
  for(int Niter = 0; Niter < n_iter+1; Niter++)
    {
      //---- Chosen dynamics ----
      switch (dynamics) {
      case 0: 
	Overdamped_ABF();
      case 1:
	Langevin_ABF();
	break;
      }

      //--- Selection procedure if requested ---
      if (selection == 1)
	Select();
      
      //--- Energy and RC of the first system only ---
      acquisition += 1;
      if (acquisition == freq)
	{
	  acquisition = 0;
	  time = Niter*dt;
	  ENERGY << time << " " << X.energy(0) << " " << X.potential_energy(0) << " " << endl;
	  COORD1 << time << " " << RC(0) << endl;
	}
      
      //--- Histograms: mean force, bias, reaction coordinates ---
      acquisition_histo += 1;
      if (acquisition_histo == freq_histo)
	{
	  cout <<  "  Output at step " << Niter << endl;
	  //---- Rate of selection -------
	  if (selection == 1)
	    cout << "         Birth rate : " << X.nb_birth*1./(Niter*dt)/nb_replicas 
		 << ", death rate : " << X.nb_death*1./(Niter*dt)/nb_replicas << endl;
	  acquisition_histo = 0;
	  //--- Values of the reaction coordinate ---
	  for (int k = 0; k < nb_replicas; k++)
	    COORD  << RC(k) << " ";
	  COORD << endl;
	  //--- histogram of the marginal ---
	  normalization = 0;
	  for (int k = 0; k < Nxi; k++)
	    normalization += AverageForce(k,0);
	  for (int k = 0; k < Nxi; k++)
	    HISTO << AverageForce(k,0)/normalization << " ";
	  HISTO << endl;
	  //--- Mean force and bias: loop over the bins ---
	  for (int k = 0; k < Nxi; k++)
	    {
	      if (AverageForce(k,0) > 0)
		CurrentForce = AverageForce(k,1)/AverageForce(k,0);
	      else 
		CurrentForce = 0;
	      MEAN_FORCE << CurrentForce << " ";
	      if (k == 0)
		Bias(k) = 0;
	      else 
		Bias(k) = Bias(k-1) + CurrentForce * DeltaXi; 
	    }
	  for (int k = 0; k < Nxi; k++)
	    BIAS << Bias(k) - Bias(Middle) << " ";
	  MEAN_FORCE << endl;
	  BIAS       << endl;
	}
    }
  
  //------------- End of computation ---------------
  cout << endl;
  cout << "------- End of computation ------" << endl;
  cout << endl;
  cout << endl;
}









