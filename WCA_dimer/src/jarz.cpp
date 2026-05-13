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

  //------ jarzynski ------
  Nxi = (int)floor(Time/dt);
  Xi.set_size(Nxi+1);
  Xi.zeros();  
  X.Work.set_size(nb_replicas);
  X.Work.zeros();
  X.Work2.set_size(nb_replicas);
  X.Work2.zeros();
  Lagrange.set_size(nb_replicas);
  Lagrange.zeros();
  CurrentMeanForce.set_size(nb_replicas);
  CurrentMeanForce.zeros();
  LagrangeMeanForce.set_size(nb_replicas,2);
  LagrangeMeanForce.zeros();
  e_x.set_size(nb_replicas); 
  e_y.set_size(nb_replicas);
  
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
  cout << "---------------------------------------------------" << endl;
  cout << "   DIMER IN WCA SOLVENT : JARZYNSKI / OVERDAMPED   " << endl;
  cout << "---------------------------------------------------" << endl;
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
  cout << " Switching time                 : " << Time          << endl;
  cout << " Number of replicas             : " << nb_replicas   << endl;
  cout << " Lower value of RC     (xi_min) : " << xi_min     << endl; 
  cout << " Upper value of RC     (xi_max) : " << xi_max     << endl;
  cout << " Time step                 (dt) : " << dt            << endl;
  cout << " Number of iterations           : " << Nxi           << endl;
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
  cout << endl;
  cout << "----------------------------------------- " << endl;
  cout << endl;
      
}

//---------------------------------------------------------
//   Integration of the dynamics and auxiliary functions
//---------------------------------------------------------

//------- schedule ----------
double algorithm::schedule(double x)
{
  double val = x;
  return val;
}

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
      //       (remember that force = -gradient of V...)
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

//--------- Overdamped Langevin ---------
void algorithm::Overdamped_constrained()
{

  double sigma = sqrt(2*dt/beta);
  double lagrange, Dist, factor, q1x, q1y, q2x, q2y;
  int nbx, nby;
  Lagrange.zeros();

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
      // important to do it before for the correction to be non anticipating
      length(X.qx(k,0),X.qy(k,0),X.qx(k,1),X.qy(k,1));
      Lagrange(k) = sigma*w*( dx*(X.RFx(k,0)-X.RFx(k,1))
			      +dy*(X.RFy(k,0)-X.RFy(k,1)) )/dist;
      // update the positions
      for (int i = 0; i < np; i++)
      	{
	  X.qx(k,i) += dt*H->fx(k,i) + sigma*X.RFx(k,i);
	  X.qy(k,i) += dt*H->fy(k,i) + sigma*X.RFy(k,i);
	} 
      // computation of the lagrange multiplier
      Dist = RC(k);                             // distance = 2*w*Dist + dcut
      lagrange = (current_xi-Dist)*2*pow(w,2);  // lagr. mult. = w*(target lenght-current length)
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
      Lagrange(k) += lagrange; 
    }
  //--- for the sake of security : reperiodize ---
  Periodic();
  //--- End of overdamped Langevin ---
}

//------------------------------------
//         Langevin dynamics 
//------------------------------------

void algorithm::Langevin_constrained()
{

  double b, c, Delta, lagrange;
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
      lagrange += current_vxi;
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
      LagrangeMeanForce(k,1) += current_vxi;
      LagrangeMeanForce(k,1) *= w;
      X.px(k,0) += LagrangeMeanForce(k,1)/2/w*e_x(k);
      X.px(k,1) -= LagrangeMeanForce(k,1)/2/w*e_x(k);
      X.py(k,0) += LagrangeMeanForce(k,1)/2/w*e_y(k);
      X.py(k,1) -= LagrangeMeanForce(k,1)/2/w*e_y(k);
      LagrangeMeanForce(k,1) /= 0.5*dt;
//       cout << CurrentMeanForce(k) << " : " 
// 	   << 0.5*(LagrangeMeanForce(k,0)+LagrangeMeanForce(k,1)) << "  |  "  
// 	   << LagrangeMeanForce(k,0) << "  " 
// 	   << LagrangeMeanForce(k,1) << endl;
      //       cout << (dist-dcut)/2/w << " |  " 
// 	   << current_vxi << "  " 
// 	   << e_x(k)*(X.px(k,0)-X.px(k,1)) + e_y(k)*(X.py(k,0)-X.py(k,1)) << endl;      
    }
  //--- Add kinetic energy to the potential energy ---
  H->add_k_energy(X);
  
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
      lagrange += current_vxi;
      lagrange *= w;
      X.px(k,0) += lagrange/2/w*e_x(k);
      X.px(k,1) -= lagrange/2/w*e_x(k);
      X.py(k,0) += lagrange/2/w*e_y(k);
      X.py(k,1) -= lagrange/2/w*e_y(k);
    }
  //------- end of Langevin -------
}


//----------------------------------------------------------
//
//     Sampling procedure = compute the trajectory
//
//----------------------------------------------------------

void algorithm::load()
{
  
  //--------------- Useful fields ---------------------
  ofstream ENERGY        ("../output/energy");        // energy as a function of time
  ofstream COORD         ("../output/coord");         // reaction coordinate at given time
  ofstream PMF           ("../output/PMF");           // associated PMF as a function of time
  ofstream WORK          ("../output/Work");          // work distributions (local mean force)      
  ofstream WORK_LAGRANGE ("../output/Work_lagrange"); // work distributions (lagrange mult.)     
  ofstream CURRENT_WORK  ("../output/current_work");  // current works for the first system
  double   time;
  double   SumWork, SumWorkLagrange;

  //------------- Initialization ----------------------
  Initialize();
  H->compute_force(X);
  int acquisition       = 0;
  int acquisition_histo = 0;
  
  //---- computation of the switching schedule ----
  for (int i = 0; i < Nxi+1; i++)
    Xi(i) = xi_min + (xi_max-xi_min)*schedule(i*1./(Nxi));
  PMF << Xi(0) << "  " << (Xi(1)-Xi(0))/dt << "  " << 0 << "  " << 0 << endl;
  cout << "------- Computation started -------" << endl;
  cout << endl;
  
  //----- creation of initial conditions ------
  current_xi  = Xi(0);
  current_vxi = 0;
  cout << "   --- thermalization " << endl;
  for(int i = 0; i < n_iter_thm; i++)
    {
      switch (dynamics) {
      case 0: 
	Overdamped_constrained();
	break;
      case 1:
	Langevin_constrained();
	break;
      }
    }
  // for Langevin dynamics: add the effective velocity
  current_vxi = (Xi(1)-Xi(0))/dt;
  if (dynamics == 1)
    {
      for (int k = 0; k < nb_replicas; k++)
	{
	  length(X.qx(k,0),X.qy(k,0),X.qx(k,1),X.qy(k,1));
	  e_x(k) = dx/dist;
	  e_y(k) = dy/dist;
	  X.px(k,0) += w*current_vxi*e_x(k);
	  X.px(k,1) -= w*current_vxi*e_x(k);
	  X.py(k,0) += w*current_vxi*e_y(k);
	  X.py(k,1) -= w*current_vxi*e_y(k);
	} 
    }
  
  //-------- Switching procedure ---------
  cout << "   --- production " << endl;
  for(int Niter = 1; Niter < Nxi+1; Niter++)
    {
      //--- update value of the RC ---
      current_xi = Xi(Niter);
      current_vxi = (Xi(Niter)-Xi(Niter-1))/dt;

      //---- Projected dynamics ----
      switch (dynamics) {
      case 0: 
	Overdamped_constrained();
	break;
      case 1:
	Langevin_constrained();
	break;
      }
      
      //--- computation of mean force and work updates ------
      for (int k = 0; k < nb_replicas; k++)
	{
	  switch (dynamics) {
	  case 0: 
	    mean_force();
	    Lagrange(k) += -(Xi(Niter)-Xi(Niter-1))*2*pow(w,2);
	    Lagrange(k) /= dt;
	    break;
	  case 1:
	    mean_force_Langevin();
	    Lagrange(k) = (LagrangeMeanForce(k,0)+LagrangeMeanForce(k,1))*0.5;
	    break;
	  }
	  X.Work(k) += CurrentMeanForce(k)*(Xi(Niter)-Xi(Niter-1));
	  X.Work2(k)+= Lagrange(k)*(Xi(Niter)-Xi(Niter-1));
	}
      CURRENT_WORK << CurrentMeanForce(0) << "  " << Lagrange(0) << endl; 

      //---- Computation of free energy difference -----
      SumWork         = 0;
      SumWorkLagrange = 0;
      for (int k = 0; k < nb_replicas; k++)
	{
	  SumWork         += exp(-beta*X.Work (k));
	  SumWorkLagrange += exp(-beta*X.Work2(k));
	}
      PMF  << Xi(Niter) << "  " << (Xi(Niter)-Xi(Niter-1))/dt << "  " 
	   << -log(SumWork /nb_replicas)/beta << "  " 
	   << -log(SumWorkLagrange/nb_replicas)/beta << endl;
      

      //--- Energy and RC ---
      acquisition += 1;
      if (acquisition == freq)
	{
	  acquisition = 0;
	  time = Niter*dt;
	  ENERGY     << time << " " << X.energy(0) << " " 
		     << X.potential_energy(0) << " " << endl;
	  COORD      << time << " " << RC(0) << endl;
	}
      
      //--- Work distributions ----
      acquisition_histo += 1;
      if (acquisition_histo == freq_histo)
	{
	  acquisition_histo = 0;
	  cout << "       output (Work) at time " << time << endl;
	  for (int k = 0; k < nb_replicas; k++)
	    WORK << X.Work(k) << " ";
	  WORK << endl;
	  for (int k = 0; k < nb_replicas; k++)
	    WORK_LAGRANGE << X.Work2(k) << " ";
	  WORK_LAGRANGE << endl;
	}
      
    }
  PMF  << endl;
  
  //-------------- End of computation ---------------
  cout << endl;
  cout << "------- End of computation ------" << endl;
  cout << endl;
  cout << endl;
}









