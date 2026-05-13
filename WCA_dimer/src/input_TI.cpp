//--------------------------------------
//
//   Functions required in 'input.hpp'
//
//--------------------------------------

#include"input.hpp" 

//--------- start reading of the input file -------------
void input::load(void)
{
  read_TI(method,t_step,nb_iter,nb_iter_thm,freq,xi,beta,Ndim,a,sig,eps,h,w,xi_min,xi_max,Nxi); 
  nb_replicas = 1;
  Nsolvent = Ndim*Ndim - 2;
}

//------------- Reading the file line by line --------------------
void input::read_TI(int& method_, double& t_step_, int& nb_iter_, int& nb_iter_thm_, int& freq_,double& xi_, double& beta_, int& Ndim_, double& a_, double& sig_,double& eps_, double& h_, double& w_, double& xi_min_, double& xi_max_, int& Nxi_)
{
  ifstream from("../input/input_file_TI");
  read_item<int>   (from,"Dynamics (0=Ovrdmpd, 1=Lngvn)  ",&method_);  
  read_item<double>(from,"Time step                 (dt) ",&t_step_);
  read_item<int>   (from,"Number of iterations           ",&nb_iter_);
  read_item<int>   (from,"Number of thermalization steps ",&nb_iter_thm_);
  read_item<int>   (from,"Output frequency               ",&freq_);
  read_item<double>(from,"Friction coefficient      (xi) ",&xi_);
  read_item<double>(from,"Inverse temperature     (beta) ",&beta_);
  read_item<int>   (from,"Total number of particles (^2) ",&Ndim_);
  read_item<double>(from,"Elementary cell size       (a) ",&a_);
  read_item<double>(from,"WCA equilibrium distance (sig) ",&sig_);
  read_item<double>(from,"WCA energy epsilon       (eps) ",&eps_);
  read_item<double>(from,"Dimer barrier height       (h) ",&h_);
  read_item<double>(from,"Dimer bond length          (w) ",&w_);
  read_item<double>(from,"Lower value of RC     (xi_min) ",&xi_min_);
  read_item<double>(from,"Upper value of RC     (xi_max) ",&xi_max_);
  read_item<int>   (from,"Number of PMF RC values  (Nxi) ",&Nxi_);
  
}

template <class T>

//---------- Reading a line in the input file -----------------
void input::read_item(istream& from,char *_title, T * val)
{
  char title[80];
  char c;

  from.get(title,80,':');
  from.get(c);
  if(strcmp(title,_title))
    {
      cerr << "The item reads '" << title << "' different from expected '"
	   << _title << "'." << endl;
      exit(1);
    }
  from >> *val;  
  from.get(title,80,'\n');
  from.get(c);
}




