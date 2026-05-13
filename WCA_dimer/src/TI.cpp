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

  //------ TI ------
  DeltaXi = (xi_max-xi_min)/Nxi;
  current_xi = xi_min-DeltaXi; // since value is increased by dXi at first iteration
  SumMeanForces = 0.;
  MeanForce.set_size(Nxi+1);
  MeanForce.zeros();
  Bias.set_size(Nxi);
  Bias.zeros();
  CurrentMeanForce.set_size(nb_replicas);
  CurrentMeanForce.zeros();
  LagrangeMeanForce.set_size(nb_replicas,4);
  LagrangeMeanForce.zeros();
  old_qx.set_size(nb_replicas,np);
  old_qy.set_size(nb_replicas,np);
  e_x.set_size(nb_replicas); 
  e_y.set_size(nb_replicas);
  //-- for Metropolization --
  X.old_qx.set_size(nb_replicas,np);
  X.old_qx.zeros();
  X.old_qy.set_size(nb_replicas,np);
  X.old_qy.zeros();
  X.old_px.set_size(nb_replicas,np);
  X.old_px.zeros();
  X.old_py.set_size(nb_replicas,np);
  X.old_py.zeros();
  X.old_energy.set_size(nb_replicas);
  X.old_energy.zeros();
  X.old_potential_energy.set_size(nb_replicas);
  X.old_potential_energy.zeros();
  rejection = 0;
  
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
  cout << "       DIMER IN WCA SOLVENT : TI          " << endl;
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
  cout << " Outputs every...steps          : " << freq          << endl; 
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

//--------- Mean force computations -------
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

void algorithm::mean_force_Langevin()
{
  double Dx, Dy;
  for (int k = 0; k < nb_replicas; k++)
    {
      //--- recompute the length and unit vectors ---
      length(X.qx(k,0),X.qy(k,0),X.qx(k,1),X.qy(k,1));
      //--- formula for the mean force (force already up-to-date) ---
      CurrentMeanForce(k) = (H->fx(k,1)-H->fx(k,0))*dx 
	+ (H->fy(k,1)-H->fy(k,0))*dy;
      Dx = dx/dist;
      Dy = dy/dist;
      CurrentMeanForce(k) -= pow(X.px(k,0)-X.px(k,1),2) + pow(X.py(k,0)-X.py(k,1),2)
	- pow( Dx*(X.px(k,0)-X.px(k,1)) + Dy*(X.py(k,0)-X.py(k,1)) ,2);
      CurrentMeanForce(k) *= w/dist;
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

//------------------------------------
//         Langevin dynamics 
//------------------------------------

void algorithm::Langevin_TI()
{

  double b, c, Delta, lagrange, new_energy, ratio;
  int    nbx, nby;
  double coeff = exp(-xi*dt/2);
  double sigma = sqrt( (1-pow(coeff,2))/beta );

  //--- First Ornstein-Uhlenbeck part (dt/2) ---
  for (int k = 0; k < nb_replicas; k++)
    {
      //--- Random terms ---
      for (int i = 0; i < np; ++i)
	{
	  X.RFx(k,i) = g.rand_number(0,1);
	  X.RFy(k,i) = g.rand_number(0,1);
	  X.px(k,i) = coeff*X.px(k,i) + sigma*X.RFx(k,i);
	  X.py(k,i) = coeff*X.py(k,i) + sigma*X.RFy(k,i);
	} 
      //--- satisfaction of the constraint ---
      length(X.qx(k,0),X.qy(k,0),X.qx(k,1),X.qy(k,1));
      e_x(k) = dx/dist;
      e_y(k) = dy/dist;
      lagrange = -e_x(k)*(X.px(k,0)-X.px(k,1)) - e_y(k)*(X.py(k,0)-X.py(k,1));
      lagrange *= w;
      X.px(k,0) += lagrange/2/w*e_x(k);
      X.px(k,1) -= lagrange/2/w*e_x(k);
      X.py(k,0) += lagrange/2/w*e_y(k);
      X.py(k,1) -= lagrange/2/w*e_y(k);
    }
  H->add_k_energy(X);
  
  //-------------- RATTLE part ---------------
  for (int k = 0; k < nb_replicas; k++)
    {
      // keep old configuration
      X.old_energy(k) = X.energy(k);
      X.old_potential_energy(k) = X.potential_energy(k);
      for (int i = 0; i < np; ++i)
	{
	  X.old_qx(k,i) = X.qx(k,i);
	  X.old_qy(k,i) = X.qy(k,i);
	  X.old_px(k,i) = X.px(k,i);
	  X.old_py(k,i) = X.py(k,i);
	}
      // update without constraints
      length(X.qx(k,0),X.qy(k,0),X.qx(k,1),X.qy(k,1));
      e_x(k) = dx/dist;
      e_y(k) = dy/dist;
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
  //---- computation of Lagrange multipliers ----
  for (int k = 0; k < nb_replicas; k++)
    {
      // first, translate q_1 so that the vector q_1-q_2 
      // is not defined up to a translation by the periodic box 
      nbx = rint((X.qx(k,0)-X.qx(k,1))/boxlg);
      nby = rint((X.qy(k,0)-X.qy(k,1))/boxlg);
      X.qx(k,0) -= nbx*boxlg;
      X.qy(k,0) -= nby*boxlg;
      // then, use second order equation
      //cout << sqrt( pow(X.qx(k,0)-X.qx(k,1),2) + pow(X.qy(k,0)-X.qy(k,1),2) ) 
      //   <<  " " << 2*w*current_xi+dcut << endl;
      b = w*( e_x(k)*(X.qx(k,0)-X.qx(k,1)) + e_y(k)*(X.qy(k,0)-X.qy(k,1)) );
      c = pow(X.qx(k,0)-X.qx(k,1),2) + pow(X.qy(k,0)-X.qy(k,1),2) - pow(2*w*current_xi+dcut,2);
      c *= pow(w,2);
      Delta = pow(b,2)-c;
      if ( Delta > 0)
	{
	  if (b > 0)
	    LagrangeMeanForce(k,0) = -b + sqrt(Delta);
	  else 
	    LagrangeMeanForce(k,0) = -b - sqrt(Delta);
	  LagrangeMeanForce(k,0) /= dt;
	  //--- update of positions and momenta ---
	  X.px(k,0) += LagrangeMeanForce(k,0)/2/w*e_x(k);
	  X.px(k,1) -= LagrangeMeanForce(k,0)/2/w*e_x(k);
	  X.py(k,0) += LagrangeMeanForce(k,0)/2/w*e_y(k);
	  X.py(k,1) -= LagrangeMeanForce(k,0)/2/w*e_y(k);
	  X.qx(k,0) += LagrangeMeanForce(k,0)/2/w*e_x(k)*dt;
	  X.qx(k,1) -= LagrangeMeanForce(k,0)/2/w*e_x(k)*dt;
	  X.qy(k,0) += LagrangeMeanForce(k,0)/2/w*e_y(k)*dt;
	  X.qy(k,1) -= LagrangeMeanForce(k,0)/2/w*e_y(k)*dt;
	  LagrangeMeanForce(k,0) /= dt*0.5;
	  X.qx(k,0) += nbx*boxlg;
	  X.qy(k,0) += nby*boxlg;
	}
      else
	cout << " Error in the projection step, for z = " << current_xi 
	     << ", cf. discriminant in the 2nd order equation  = " << Delta << endl;
    }
  //--- Apply periodic BC on positions --- 
  Periodic();
  //--- Update the forces ---
  H->compute_force(X);
  //--- Final update of the momenta ---
  for (int k = 0; k < nb_replicas; k++)
    {
      length(X.qx(k,0),X.qy(k,0),X.qx(k,1),X.qy(k,1));
      e_x(k) = dx/dist;
      e_y(k) = dy/dist;
      for(int i = 0 ; i < np; i++)
	{
	  X.px(k,i) += 0.5*dt*H->fx(k,i);
	  X.py(k,i) += 0.5*dt*H->fy(k,i);
	} 
      LagrangeMeanForce(k,1) = -e_x(k)*(X.px(k,0)-X.px(k,1)) - e_y(k)*(X.py(k,0)-X.py(k,1));
      LagrangeMeanForce(k,1) *= w;
      X.px(k,0) += LagrangeMeanForce(k,1)/2/w*e_x(k);
      X.px(k,1) -= LagrangeMeanForce(k,1)/2/w*e_x(k);
      X.py(k,0) += LagrangeMeanForce(k,1)/2/w*e_y(k);
      X.py(k,1) -= LagrangeMeanForce(k,1)/2/w*e_y(k);
      LagrangeMeanForce(k,1) /= 0.5*dt;
      // check satisfaction of constraints (optional): positions then momenta
      //cout << (dist-dcut)/2/w << "  " 
      //<< e_x(k)*(X.px(k,0)-X.px(k,1)) + e_y(k)*(X.py(k,0)-X.py(k,1)) << endl;
    }
  //--- Add kinetic energy to the potential energy ---
  H->add_k_energy(X);
  //--- acceptance or rejection step : generalized Metropolis Hastings ---
  for (int k = 0; k < nb_replicas; k++)
    {
      new_energy = X.energy(k);
      ratio = exp( -beta*(new_energy-X.old_energy(k)) );
      if (ratio < g.random_number() ) 
	{
	  rejection += 1;
	  // rejection with momentum reversal
	  X.energy(k) = X.old_energy(k);
	  X.potential_energy(k) = X.old_potential_energy(k);
	  for (int i = 0; i < np; ++i)
	    {
	      X.qx(k,i) = X.qx(k,i);
	      X.qy(k,i) = X.qy(k,i);
	      X.px(k,i) = -X.old_px(k,i);
	      X.py(k,i) = -X.old_py(k,i);
	    }
	}
    }
  //--- Second Ornstein-Uhlenbeck part (dt/2) ---
  for (int k = 0; k < nb_replicas; k++)
    {
      //--- Random terms ---
      for (int i = 0; i < np; ++i)
	{
	  X.RFx(k,i) = g.rand_number(0,1);
	  X.RFy(k,i) = g.rand_number(0,1);
	  X.px(k,i) = coeff*X.px(k,i) + sigma*X.RFx(k,i);
	  X.py(k,i) = coeff*X.py(k,i) + sigma*X.RFy(k,i);
	} 
      //--- satisfaction of the constraint ---
      length(X.qx(k,0),X.qy(k,0),X.qx(k,1),X.qy(k,1));
      e_x(k) = dx/dist;
      e_y(k) = dy/dist;
      lagrange = -e_x(k)*(X.px(k,0)-X.px(k,1)) - e_y(k)*(X.py(k,0)-X.py(k,1));
      lagrange *= w;
      X.px(k,0) += lagrange/2/w*e_x(k);
      X.px(k,1) -= lagrange/2/w*e_x(k);
      X.py(k,0) += lagrange/2/w*e_y(k);
      X.py(k,1) -= lagrange/2/w*e_y(k);
    }
  //------- end of Langevin TI -------
}

//---------------------------------------- 
//          Overdamped Langevin 
//----------------------------------------

void algorithm::Overdamped_TI()
{

  double sigma = sqrt(2*dt/beta);
  double lagrange, lagrange_TR, Dist, factor, q1x, q1y, q2x, q2y;
  int nbx, nby;
  double old_dx, old_dy, old_dist;
  
  //--- Force computation ---
  H->compute_force(X);
  //--- Loop over the replicas ---
  for (int k = 0; k < nb_replicas; k++)
    {
      for (int i = 0; i < np; i++)
      	{
	  X.RFx(k,i) = g.rand_number(0,1);
	  X.RFy(k,i) = g.rand_number(0,1);
	}
      // computation of the variance reduction term in the lagrange mult.
      // important to do it before, for the correction to be non anticipating
      length(X.qx(k,0),X.qy(k,0),X.qx(k,1),X.qy(k,1));
      LagrangeMeanForce(k,1) = sigma*w*( dx*(X.RFx(k,0)-X.RFx(k,1))
					 +dy*(X.RFy(k,0)-X.RFy(k,1)) )/dist/dt;
      for (int i = 0; i < np; i++)
      	{
	  old_qx(k,i) = X.qx(k,i);
	  old_qy(k,i) = X.qy(k,i);
	}
      old_dx = dx;
      old_dy = dy;
      old_dist = dist;
      // time reversed update
      for (int i = 0; i < np; i++)
      	{
	  X.qx(k,i) = old_qx(k,i) + dt*H->fx(k,i) - sigma*X.RFx(k,i);
	  X.qy(k,i) = old_qy(k,i) + dt*H->fy(k,i) - sigma*X.RFy(k,i);
	} 
      // computation of the lagrange multiplier
      Dist = RC(k);
      lagrange_TR = (current_xi-Dist)*2*pow(w,2);  // length = dcut+2*w*RC
      // real update the positions
      for (int i = 0; i < np; i++)
      	{
	  X.qx(k,i) = old_qx(k,i) + dt*H->fx(k,i) + sigma*X.RFx(k,i);
	  X.qy(k,i) = old_qy(k,i) + dt*H->fy(k,i) + sigma*X.RFy(k,i);
	} 
      // computation of the lagrange multiplier
      Dist = RC(k);
      lagrange = (current_xi-Dist)*2*pow(w,2);  // length = dcut+2*w*RC
      lagrange_TR += lagrange;
      // projection step
      factor = lagrange/w/(dcut+2*w*current_xi);
      if ( fabs(factor-1) > tol ) // avoid singularities
	factor = 1/(1-factor);
      else 
	factor = 0;
      // first, translate q_1 so that the vector q_1-q_2 
      // is not defined up to a translation by the periodic box 
      nbx = rint((X.qx(k,0)-X.qx(k,1))/boxlg);
      nby = rint((X.qy(k,0)-X.qy(k,1))/boxlg);
      X.qx(k,0) -= nbx*boxlg;
      X.qy(k,0) -= nby*boxlg;
      // then compute the new positions
      q1x = X.qx(k,0);
      q1y = X.qy(k,0);
      q2x = X.qx(k,1);
      q2y = X.qy(k,1);
      X.qx(k,0) = 0.5*(1+factor)*q1x + 0.5*(1-factor)*q2x;
      X.qy(k,0) = 0.5*(1+factor)*q1y + 0.5*(1-factor)*q2y;
      X.qx(k,1) = 0.5*(1-factor)*q1x + 0.5*(1+factor)*q2x;
      X.qy(k,1) = 0.5*(1-factor)*q1y + 0.5*(1+factor)*q2y;
      // come back to the original representation
      X.qx(k,0) += nbx*boxlg;
      X.qy(k,0) += nby*boxlg;
      // computation of mean force from lagrange multipliers
      LagrangeMeanForce(k,0)  = lagrange/dt;  
      LagrangeMeanForce(k,1) += lagrange/dt;  
      LagrangeMeanForce(k,2)  = 0.5*lagrange_TR/dt;  
    }
  //--- for the sake of security : reperiodize ---
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
  ofstream ENERGY     ("../output/energy");             // energy as a function of time
  ofstream COORD      ("../output/coord");              // reaction coordinate at given time
  ofstream CURRENT_MF ("../output/current_mean_force"); // current estimate of mean force
  ofstream MEAN_FORCE ("../output/mean_force");         // mean force as a function of time
  ofstream BIAS       ("../output/PMF");                // associated PMF as a function of time
  ofstream COMPARE    ("../output/Lagrange");           // compare mean force and Lagrange mult.
  double time;
  int    Middle  = rint(Nxi/2);

  //------------- Initialization ----------------------
  Initialize();
  //--- Set dimer bond length to first value ---
  for (int k = 0; k < nb_replicas; k++)
    X.qy(k,1) = X.qy(k,0) + dcut + xi_min*2*w;
  H->compute_force(X);
  int acquisition = -1;
  for (int i = 0; i < Nxi+1; i++)
    {
      MEAN_FORCE << xi_min + i*DeltaXi << "  ";
      BIAS       << xi_min + i*DeltaXi << "  ";
    }
  MEAN_FORCE << endl;
  BIAS       << endl;
  cout << "------- Computation started -------" << endl;
  cout << endl;
  
  //-------- Loop over the values of the reaction coordinate ---------
  for(int i = 0; i < Nxi+1; i++)
    {
      //--- update value of the RC ---
      current_xi += DeltaXi;
      SumMeanForces       = 0.;
      SumLagrange         = 0.;
      SumLagrangeAveraged = 0.;
      SumLagrangePosition = 0.;
      SumLagrangeMomenta  = 0.;
      SumLagrangeReduced  = 0.;
      SumLagrangeTR       = 0.;
      cout << "   Sampling at z = " << current_xi << endl;
      
      //--------------- Preliminary thermalization --------
      cout << "           --- thermalization " << endl;
      for(int Niter = 0; Niter < n_iter_thm; Niter++)
	{
	  //---- Projected dynamics ----
	  switch (dynamics) {
	  case 0: 
	    Overdamped_TI();
	    break;
	  case 1:
	    Langevin_TI();
	    break;
	  }
	}

      //-------------- Loop over time ---------------------
      cout << "           --- production " << endl;
      rejection = 0;
      acquisition = -1;
      for(int Niter = 0; Niter < n_iter+1; Niter++)
	{
	  //---- Projected dynamics ----
	  switch (dynamics) {
	  case 0:  
	    //--------- Overdamped case --------
	    Overdamped_TI();
	    //--- computation of local mean force ----
	    mean_force();
	    //--- computation of mean force ----
	    for (int k = 0; k < nb_replicas; k++)
	      {
		SumMeanForces      += CurrentMeanForce(k);
		SumLagrange        += LagrangeMeanForce(k,0);
		SumLagrangeReduced += LagrangeMeanForce(k,1);
		SumLagrangeTR      += LagrangeMeanForce(k,2);
	      }
	    break;
	  case 1:
	     //--------- Langevin ---------
	    Langevin_TI();
	    //--- computation of averaged local mean force ----
	    mean_force();
	    for (int k = 0; k < nb_replicas; k++)
	      {
		LagrangeMeanForce(k,3) = CurrentMeanForce(k);
		SumLagrangeAveraged += CurrentMeanForce(k);
	      }
	    //--- computation of local mean force ----
	    mean_force_Langevin();
	    //--- computation of mean force ----
	    for (int k = 0; k < nb_replicas; k++)
	      {
		LagrangeMeanForce(k,2) = 0.5*(LagrangeMeanForce(k,0)+LagrangeMeanForce(k,1));
		SumMeanForces       += CurrentMeanForce(k);
		SumLagrange         += LagrangeMeanForce(k,2);
		SumLagrangePosition += LagrangeMeanForce(k,0);
		SumLagrangeMomenta  += LagrangeMeanForce(k,1);
	      }
	    break;
	  }

	  //--- Energy and RC ---
	  acquisition += 1;
	  if (acquisition == freq)
	    {
	      cout << "       Output for xi = " << current_xi << ", at step " << Niter << endl; 
	      acquisition = 0;
	      time = Niter*dt;
	      ENERGY     << X.energy(0) << " " 
			 << X.potential_energy(0) << " " << endl;
	      COORD      << RC(0) << endl;
	      switch (dynamics) {
	      case 0: 
		COMPARE << CurrentMeanForce(0) << "  " << LagrangeMeanForce(0,2)
			<< "  " << LagrangeMeanForce(0,0) 
			<< "  " << LagrangeMeanForce(0,1) << endl;
		CURRENT_MF << SumMeanForces/(nb_replicas*(Niter+1)) << " " 
			   << SumLagrange/(nb_replicas*(Niter+1)) << "  "
			   << SumLagrangeReduced/(nb_replicas*(Niter+1)) << " " 
			   << SumLagrangeTR/(nb_replicas*(Niter+1)) << endl;
		break;
	      case 1:
		cout << "            -- rejection rate : " << rejection*1./(Niter+1) << endl;
		COMPARE << CurrentMeanForce(0) << "  " << LagrangeMeanForce(0,2)
			<< "  " << LagrangeMeanForce(0,3) 
			<< "  " << LagrangeMeanForce(0,0) 
			<< "  " << LagrangeMeanForce(0,1) << endl;
		CURRENT_MF << SumMeanForces/(nb_replicas*(Niter+1)) << " " 
			   << SumLagrange/(nb_replicas*(Niter+1)) << "  "
			   << SumLagrangeAveraged/(nb_replicas*(Niter+1)) << "  "
			   << SumLagrangePosition/(nb_replicas*(Niter+1)) << "  "
			   << SumLagrangeMomenta/(nb_replicas*(Niter+1)) << endl;
		break;
	      }
	    }
	  
	}
      
      //--- keep track of the estimated mean force ---
      MeanForce(i) = SumMeanForces/(nb_replicas*(n_iter+1));
      MEAN_FORCE << MeanForce(i) << " ";
    }

  //--- Compute finally approximation of PMF (simple quadrature) ---
  MEAN_FORCE << endl;
  for (int k = 1; k < Nxi+1; k++)
    Bias(k) = Bias(k-1) + MeanForce(k) * DeltaXi;  
  for (int k = 0; k < Nxi+1; k++)
    BIAS << Bias(k) - Bias(Middle) << " ";
  BIAS << endl;
  
  
  //-------------- End of computation ---------------
  cout << endl;
  cout << "------- End of computation ------" << endl;
  cout << endl;
  cout << endl;
}









