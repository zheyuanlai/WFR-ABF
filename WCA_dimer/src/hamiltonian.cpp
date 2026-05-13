#include"hamiltonian.hpp"

//--- Compute kinetic energy ---
void hamiltonian::add_k_energy(particle & X)
{
  for (int k = 0; k < nb_sys; k++)
    {
      X.energy(k) = X.potential_energy(k);
      for (int i = 0; i < np; i++)
	X.energy(k) += 0.5*X.px(k,i)*X.px(k,i) + 0.5*X.py(k,i)*X.py(k,i);
    }
}


///-----------------------------------------------
//          Elementary interactions
//----------------------------------------------

//--- value of the reaction coordinate, given the bond length ---
double hamiltonian::RC(double d)
{
  double value;
  value = (d-dcut)/(2*w);  
  return value;
}

//------ Compute length of the vector, up to periodic BC --------
void hamiltonian::lengthBC(double x1, double y1, double x2, double y2)
{
  dx = x1 - x2;
  dx -= boxlg*rint(dx*1./boxlg);
  dy = y1 - y2;
  dy -= boxlg*rint(dy*1./boxlg);
  dist = sqrt(dx*dx + dy*dy);    // distance up to BC 
  dx = dx/dist;                  // unit vector (dx,dy)
  dy = dy/dist;
}

//------------ WCA potential energy ----------
double hamiltonian::WCA(double d)  
{
  double pot;
  if (d > dcut )
    pot = 0;
  else pot = eps + 4*eps*( pow(sig/d,12.) - pow(sig/d,6.) );
  return pot; 
}

//------------ WCA force -------------
double hamiltonian::fWCA(double d)   
{
  double ff;
  if (d > dcut)
    ff = 0;
  else ff = -24.*eps*( 2*pow(sig,12.)*pow(1./d,13.) - pow(sig,6)*pow(1./d,7.) );
  return ff;
}

//----------- Dimer potential -------------
double hamiltonian::DW(double d)  
{
  double pot;
  double dd = (d-dcut-w*sig)/(w*sig);
  pot = h*pow( 1-pow(dd,2.), 2.);
  return pot; 
}

//----------- Dimer force -----------------
double hamiltonian::fDW(double d)  
{
  double ff;
  double dd = (d-dcut-w*sig)/(w*sig);
  ff = -4*h*( 1-pow(dd,2.) )*dd/(w*sig);
  return ff;
}

//-------- Restraining potential for M-BAR ---------
double hamiltonian::restraining_potential(double d, double d0)
{
  double pot;
  pot = 0.5*K*pow(d-d0,2.);
  return pot;
}

//-------- Restraining force for M-BAR ---------
double hamiltonian::restraining_force(double d, double d0)
{
  double ff;
  ff = K*(d-d0)/(2*w);  // cf. expression of RC...
  return ff;
}

//----------------------------------------------------
//    Global force and potential energy computation
//----------------------------------------------------

void hamiltonian::compute_force(particle& X)
{
  
  //--- reset forces to 0 ---
  fx.zeros();
  fy.zeros();

  //---------- Loop over the replicas --------------
  for (int k = 0; k < nb_sys; k++)
    {
      //--- double loop over the particles, depending on their types ---
      double enpot = 0.;
      
      // 1 - solvent/solvent interactions
      for (int i = 2; i < np; ++i)
	{
	  for (int j = i+1; j < np; ++j)
	    { 
	      lengthBC(X.qx(k,j),X.qy(k,j),X.qx(k,i),X.qy(k,i));
	      if (dist < dcut)
		{
		  enpot = enpot + WCA(dist);
		  act = fWCA(dist);
		  fx(k,i) += act*dx;  
		  fy(k,i) += act*dy;
		  fx(k,j) -= act*dx;
		  fy(k,j) -= act*dy;
		}
	    }
	}
      
      // 2 - solvent/dimer interactions
      for (int i=0; i<2; ++i)
	{
	  for (int j=2; j<np; ++j)
	    { 
	      lengthBC(X.qx(k,j),X.qy(k,j),X.qx(k,i),X.qy(k,i));
	      if (dist < dcut)
		{
		  enpot += WCA(dist);
		  act = fWCA(dist);
		  fx(k,i) += act*dx;  
		  fy(k,i) += act*dy;
		  fx(k,j) -= act*dx;
		  fy(k,j) -= act*dy;
		}
	    }
	}

      // 3 - dimer/dimer interactions
      lengthBC(X.qx(k,1),X.qy(k,1),X.qx(k,0),X.qy(k,0));
      enpot += DW(dist);
      act = fDW(dist);
      fx(k,0) += act*dx;  
      fy(k,0) += act*dy;
      fx(k,1) -= act*dx;
      fy(k,1) -= act*dy;
      
      // potential energy of replica 'k'
      X.energy(k) = enpot;
      X.potential_energy(k) = enpot;
      
    }
  
  //--------- end of force computations ----
  
}

//----------------------------------------------------
//           Restraint force
//----------------------------------------------------

void hamiltonian::add_restraint_force(particle& X, double center)
{
  double enpot_restraint, RClength;

  //------ Loop over the replicas, even if only 1 in M-BAR case... ------
  for (int k = 0; k < nb_sys; k++)
    {
      // restraint acts only on the dimer
      lengthBC(X.qx(k,1),X.qy(k,1),X.qx(k,0),X.qy(k,0));
      RClength = RC(dist);
      enpot_restraint = restraining_potential(RClength,center);
      X.energy(k) += enpot_restraint;
      X.potential_energy(k) += enpot_restraint;
      act = restraining_force(RClength,center);
      fx(k,0) += act*dx;  
      fy(k,0) += act*dy;
      fx(k,1) -= act*dx;
      fy(k,1) -= act*dy;
    }
}


