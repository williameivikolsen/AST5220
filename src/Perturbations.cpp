#include"Perturbations.h"

//====================================================
// Constructors
//====================================================

Perturbations::Perturbations(
    BackgroundCosmology *cosmo, 
    RecombinationHistory *rec) : 
  cosmo(cosmo), 
  rec(rec)
{}

//====================================================
// Do all the solving
//====================================================

void Perturbations::solve(){

  // Integrate all the perturbation equation and spline the result
  integrate_perturbations();

  // Compute source functions and spline the result
  compute_source_functions();
}

//====================================================
// The main work: integrate all the perturbations
// and spline the results
//====================================================

void Perturbations::integrate_perturbations(){
  Utils::StartTiming("integrateperturbation");

  //===================================================================
  // TODO: Set up the k-array for the k's we are going to integrate over
  // Start at k_min end at k_max with n_k points with either a
  // quadratic or a logarithmic spacing
  //===================================================================
  int n_ell      = Constants.n_ell_theta;
  double H0      = Constants.H0_over_h*cosmo->get_h();
  double Omega_R = cosmo->get_OmegaR();
  double c       = Constants.c;

  Vector k_array = exp(Utils::linspace(log(k_min), log(k_max), n_k));
  Vector x_full  = Utils::linspace(x_start, x_end, n_x);
  Vector delta_cdm(n_x * n_k);
  Vector delta_b(n_x * n_k);
  Vector v_cdm(n_x * n_k);
  Vector v_b(n_x * n_k);
  Vector Phi(n_x * n_k);
  Vector Psi(n_x * n_k);
  Vector Theta(n_ell * n_x * n_k);

  // Loop over all wavenumbers
  for(int ik = 0; ik < n_k; ik++){

    // Progress bar...
    if( (10*ik) / n_k != (10*ik+10) / n_k ) {
      std::cout << (100*ik+100)/n_k << "% " << std::flush;
      if(ik == n_k-1) std::cout << std::endl;
    }

    // Current value of k
    double k = k_array[ik];

    // Find value to integrate to
    double idx_end_tight = get_tight_coupling_time(k);

    //===================================================================
    // TODO: Tight coupling integration
    // Remember to implement the routines:
    // set_ic : The IC at the start
    // rhs_tight_coupling_ode : The dydx for our coupled ODE system
    //===================================================================

    // Set up initial conditions in the tight coupling regime
    auto y_tight_coupling_ini = set_ic(x_start, k);

    // The tight coupling ODE system
    ODEFunction dydx_tight_coupling = [&](double x, const double *y, double *dydx){
      return rhs_tight_coupling_ode(x, k, y, dydx);
    };

    // Integrate from x_start -> x_end_tight
    // Her har vi endret get_tight_coupling time til å gi ut indeksen istedenfor x
    double x_end_tight = x_full[idx_end_tight];
    Vector x_tc        = Utils::linspace(x_start, x_end_tight, idx_end_tight + 1);
    ODESolver tc_ode;
    tc_ode.solve(dydx_tight_coupling, x_tc, y_tight_coupling_ini);
    auto y_tc          = tc_ode.get_final_data();

    //===================================================================
    // TODO: Full equation integration
    // Remember to implement the routines:
    // set_ic_after_tight_coupling : The IC after tight coupling ends
    // rhs_full_ode : The dydx for our coupled ODE system
    //===================================================================

    // Set up initial conditions (y_tight_coupling is the solution at the end of tight coupling)
    auto y_full_ini = set_ic_after_tight_coupling(y_tc, x_end_tight, k);

    // The full ODE system
    ODEFunction dydx_full = [&](double x, const double *y, double *dydx){
      return rhs_full_ode(x, k, y, dydx);
    };

    // Integrate from x_end_tight -> x_end
    Vector x_after_tc = Utils::linspace(x_end_tight, x_end, n_x - idx_end_tight);
    ODESolver after_tc_ode;
    after_tc_ode.solve(dydx_full, x_after_tc, y_full_ini);

    //===================================================================
    // TODO: remember to store the data found from integrating so we can
    // spline it below
    //
    // To compute a 2D spline of a function f(x,k) the data must be given 
    // to the spline routine as a 1D array f_array with the points f(ix, ik) 
    // stored as f_array[ix + n_x * ik]
    // Example:
    // Vector x_array(n_x);
    // Vector k_array(n_k);
    // Vector f(n_x * n_k);
    // Spline2D y_spline;
    // f_spline.create(x_array, k_array, f_array);
    // We can now use the spline as f_spline(x, k)
    //
    //===================================================================

    for(int ix = 0; ix < idx_end_tight; ix++){
      delta_cdm[ix + n_x*ik] = tc_ode.get_data_by_component(Constants.ind_deltacdm_tc)[ix];
      delta_b[ix + n_x*ik]   = tc_ode.get_data_by_component(Constants.ind_deltab_tc)[ix];
      v_cdm[ix + n_x*ik]     = tc_ode.get_data_by_component(Constants.ind_vcdm_tc)[ix];
      v_b[ix + n_x*ik]       = tc_ode.get_data_by_component(Constants.ind_vb_tc)[ix];
      Phi[ix + n_x*ik]       = tc_ode.get_data_by_component(Constants.ind_Phi_tc)[ix];
    }
    for(int ix = idx_end_tight; ix < n_x; ix++){
      delta_cdm[ix + n_x*ik] = after_tc_ode.get_data_by_component(Constants.ind_deltacdm)[ix - idx_end_tight];
      delta_b[ix + n_x*ik]   = after_tc_ode.get_data_by_component(Constants.ind_deltab)[ix - idx_end_tight];
      v_cdm[ix + n_x*ik]     = after_tc_ode.get_data_by_component(Constants.ind_vcdm)[ix - idx_end_tight];
      v_b[ix + n_x*ik]       = after_tc_ode.get_data_by_component(Constants.ind_vb)[ix - idx_end_tight];
      Phi[ix + n_x*ik]       = after_tc_ode.get_data_by_component(Constants.ind_Phi)[ix - idx_end_tight];
    }

    // Multipole quantities
    for(int ell = 0; ell < n_ell; ell++){
      // Tight coupling
      for(int ix = 0; ix < idx_end_tight; ix++){
        if(ell < 2){
          Theta[ix + n_x*(ell + n_ell*ik)] = tc_ode.get_data_by_component(Constants.ind_start_theta_tc + ell)[ix];
        }
        else{
          Theta[ix + n_x*(ell + n_ell*ik)] = 0;
        }
      }
      // After tight coupling
      for(int ix = idx_end_tight; ix < n_x; ix++){
        Theta[ix + n_x*(ell + n_ell*ik)] = after_tc_ode.get_data_by_component(Constants.ind_start_theta + ell)[ix - idx_end_tight];
      }
    }
  }
  Utils::EndTiming("integrateperturbation");

  //=============================================================================
  // TODO: Make all splines needed: Theta0,Theta1,Theta2,Phi,Psi,...
  // =============================================================================
  delta_cdm_spline.create(x_full, k_array, delta_cdm);
  delta_b_spline.create(x_full, k_array, delta_b);
  v_cdm_spline.create(x_full, k_array, v_cdm);
  v_b_spline.create(x_full, k_array, v_b);
  Phi_spline.create(x_full, k_array, Phi);

  Theta_spline = std::vector<Spline2D>(n_ell);
  Vector Theta_xk(n_x * n_k);
  for(int ell = 0; ell < n_ell; ell++){
    for(int ik = 0; ik < n_k; ik++){
      for(int ix = 0; ix < n_x; ix++){
        Theta_xk[ix + n_x*ik] = Theta[ix + n_x*(ell + n_ell*ik)];
      }
    }
    std::string name = "theta" + std::to_string(ell) + "_x_k";
    Theta_spline[ell].create(x_full, k_array, Theta_xk, name);
  }
  for(int ik = 0; ik < n_k; ik++){
      for(int ix = 0; ix < n_x; ix++){
        Psi[ix + n_x*ik] = - Phi[ix + n_x*ik] - 12*pow(H0/(c*k_array[ik]*exp(x_full[ix])),2)*Omega_R*Theta[ix + n_x*(2 + n_ell*ik)];
      }
  }
  Psi_spline.create(x_full, k_array, Psi);
}

//====================================================
// Set IC at the start of the run (this is in the
// tight coupling regime)
//====================================================
Vector Perturbations::set_ic(const double x, const double k) const{
  // Konstanter
  double c           = Constants.c;
  double Hp          = cosmo->Hp_of_x(x);
  double dtaudx      = rec->dtaudx_of_x(x);

  // The vector we are going to fill
  Vector y_tc(Constants.n_ell_tot_tc);

  //=============================================================================
  // Compute where in the y_tc array each component belongs
  // This is just an example of how to do it to make it easier
  // Feel free to organize the component any way you like
  //=============================================================================
  
  // For integration of perturbations in tight coupling regime (Only 2 photon multipoles + neutrinos needed)
  const int n_ell_theta_tc      = Constants.n_ell_theta_tc;
  const int n_ell_neutrinos_tc  = Constants.n_ell_neutrinos_tc;
  const int n_ell_tot_tc        = Constants.n_ell_tot_tc;
  const bool polarization       = Constants.polarization;
  const bool neutrinos          = Constants.neutrinos;

  // References to the tight coupling quantities
  double &delta_cdm    =  y_tc[Constants.ind_deltacdm_tc];
  double &delta_b      =  y_tc[Constants.ind_deltab_tc];
  double &v_cdm        =  y_tc[Constants.ind_vcdm_tc];
  double &v_b          =  y_tc[Constants.ind_vb_tc];
  double &Phi          =  y_tc[Constants.ind_Phi_tc];
  double *Theta        = &y_tc[Constants.ind_start_theta_tc];
  double *Nu           = &y_tc[Constants.ind_start_nu_tc];

  //=============================================================================
  // TODO: Set the initial conditions in the tight coupling regime
  //=============================================================================

  // SET: Scalar quantities (Gravitational potential, baryons and CDM)
  delta_cdm       = 1.0;
  delta_b         = delta_cdm;
  v_cdm           = c*k/(3.0*Hp);
  v_b             = v_cdm;
  Phi             = 2.0/3.0;

  // SET: Photon temperature perturbations (Theta_ell)
  Theta[0]        = 1.0/3.0;
  Theta[1]        = -c*k/(9.0*Hp);

  // SET: Neutrino perturbations (N_ell)
  if(neutrinos){
    // ...
    // ...
  }

  return y_tc;
}

//====================================================
// Set IC for the full ODE system after tight coupling 
// regime ends
//====================================================

Vector Perturbations::set_ic_after_tight_coupling(
    const Vector &y_tc, 
    const double x, 
    const double k) const{

  // Make the vector we are going to fill
  Vector y(Constants.n_ell_tot_full);
  
  //=============================================================================
  // Compute where in the y array each component belongs and where corresponding
  // components are located in the y_tc array
  // This is just an example of how to do it to make it easier
  // Feel free to organize the component any way you like
  //=============================================================================

  // Number of multipoles we have in the full regime
  const int n_ell_theta         = Constants.n_ell_theta;
  const int n_ell_thetap        = Constants.n_ell_thetap;
  const int n_ell_neutrinos     = Constants.n_ell_neutrinos;
  const bool polarization       = Constants.polarization;
  const bool neutrinos          = Constants.neutrinos;

  // Number of multipoles we have in the tight coupling regime
  const int n_ell_theta_tc      = Constants.n_ell_theta_tc;
  const int n_ell_neutrinos_tc  = Constants.n_ell_neutrinos_tc;

  // References to the tight coupling quantities
  const double &delta_cdm_tc    =  y_tc[Constants.ind_deltacdm_tc];
  const double &delta_b_tc      =  y_tc[Constants.ind_deltab_tc];
  const double &v_cdm_tc        =  y_tc[Constants.ind_vcdm_tc];
  const double &v_b_tc          =  y_tc[Constants.ind_vb_tc];
  const double &Phi_tc          =  y_tc[Constants.ind_Phi_tc];
  const double *Theta_tc        = &y_tc[Constants.ind_start_theta_tc];
  const double *Nu_tc           = &y_tc[Constants.ind_start_nu_tc];

  // References to the quantities we are going to set
  double &delta_cdm       =  y[Constants.ind_deltacdm_tc];
  double &delta_b         =  y[Constants.ind_deltab_tc];
  double &v_cdm           =  y[Constants.ind_vcdm_tc];
  double &v_b             =  y[Constants.ind_vb_tc];
  double &Phi             =  y[Constants.ind_Phi_tc];
  double *Theta           = &y[Constants.ind_start_theta_tc];
  double *Theta_p         = &y[Constants.ind_start_thetap_tc];
  double *Nu              = &y[Constants.ind_start_nu_tc];

  //=============================================================================
  // TODO: fill in the initial conditions for the full equation system below
  // NB: remember that we have different number of multipoles in the two
  // regimes so be careful when assigning from the tc array
  //=============================================================================

  double Hp     = cosmo->Hp_of_x(x);
  double dtaudx = rec->dtaudx_of_x(x);
  double c      = Constants.c;

  // SET: Scalar quantities (Gravitational potental, baryons and CDM)
  delta_cdm       = delta_cdm_tc;
  delta_b         = delta_b_tc;
  v_cdm           = v_cdm_tc;
  v_b             = v_b_tc;
  Phi             = Phi_tc;

  // SET: Photon temperature perturbations (Theta_ell)
  Theta[0]        = Theta_tc[0];
  Theta[1]        = Theta_tc[1];
  Theta[2]        = -4.0*c*k/(9.0*Hp*dtaudx)*Theta[1];
  for(int ell = 2; ell < n_ell_theta; ell++){
    Theta[ell] = -ell/(2*ell+1)*c*k/(Hp*dtaudx)*Theta[ell-1];
  }

  // SET: Photon polarization perturbations (Theta_p_ell)
  if(polarization){
    // ...
    // ...
  }

  // SET: Neutrino perturbations (N_ell)
  if(neutrinos){
    // ...
    // ...
  }

  return y;
}

//====================================================
// The time when tight coupling end
//====================================================

int Perturbations::get_tight_coupling_time(const double k) const{
  // Merk: Jeg har endret funksjonen til å returnere indeks istedenfor x-verdi
  int idx_tight_coupling_end;
  double x_rec = -7.6;        // Time of recombination
  //=============================================================================
  // TODO: compute and return x for when tight coupling ends
  // Remember all the three conditions in Callin
  //=============================================================================
  Vector x_array = Utils::linspace(x_start, x_end, n_x);
  double c       = Constants.c;
  for(int i = 0; i < n_x; i++){
    double x          = x_array[i];
    double Hp         = cosmo->Hp_of_x(x);
    double dtaudx     = rec->dtaudx_of_x(x);
    double max = 1.0;
    if(c*k/Hp > max){
      max = c*k/Hp;
    }
    if(abs(dtaudx) < 10.0*max){
      idx_tight_coupling_end = i;
      break;
    }
    if(x > x_rec){
      idx_tight_coupling_end = i;
      break;
    }
  }

  return idx_tight_coupling_end;
}

//====================================================
// After integrsating the perturbation compute the
// source function(s)
//====================================================
void Perturbations::compute_source_functions(){
  Utils::StartTiming("source");

  //=============================================================================
  // TODO: Make the x and k arrays to evaluate over and use to make the splines
  //=============================================================================
  Vector k_array = exp(Utils::linspace(log(k_min), log(k_max), n_k));
  Vector x_array = Utils::linspace(x_start, x_end, n_x);

  // Make storage for the source functions (in 1D array to be able to pass it to the spline)
  Vector ST_array(k_array.size() * x_array.size());
  Vector SE_array(k_array.size() * x_array.size());

  // Compute source functions
  for(auto ix = 0; ix < x_array.size(); ix++){
    const double x = x_array[ix];
    for(auto ik = 0; ik < k_array.size(); ik++){
      const double k = k_array[ik];

      // NB: This is the format the data needs to be stored 
      // in a 1D array for the 2D spline routine source(ix,ik) -> S_array[ix + nx * ik]
      const int index = ix + n_x * ik;

      //=============================================================================
      // TODO: Compute the source functions
      //=============================================================================
      // Fetch all the things we need...
      const double Hp       = cosmo->Hp_of_x(x);
      const double tau      = rec->tau_of_x(x);
      const double g        = rec->g_tilde_of_x(x);
      const double c        = Constants.c;
      double Theta0         = Theta_spline[0](x,k);
      double Theta2         = Theta_spline[2](x,k);
      double dPsidx         = Psi_spline.deriv_x(x,k);
      double dPhidx         = Phi_spline.deriv_x(x,k);
      double dHpdx          = cosmo->dHpdx_of_x(x);
      double ddHpddx        = cosmo->ddHpddx_of_x(x);
      double dgdx           = rec->dgdx_tilde_of_x(x);
      double ddgddx         = rec->ddgddx_tilde_of_x(x);
      double v_b            = v_b_spline(x,k);
      double dv_bdx         = v_b_spline.deriv_x(x,k);
      double ddv_bddx       = v_b_spline.deriv_xx(x,k);

      // Temperature source
      // ST_array[index] = g*(Theta0 + Psi_spline(x,k) + 0.25*Theta2) + exp(-tau)*(dPsidx - dPhidx)
      //                 - 1.0/(c*k)*(dHpdx*g*v_b + Hp*dgdx*v_b + Hp*g*dv_bdx) + 3.0/(4.0*pow(c*k,2))*(
      //                   pow(dHpdx,2)*g*v_b + Hp*ddHpddx*g*v_b + Hp*dHpdx*dgdx*v_b + Hp*dHpdx*g*dv_bdx
      //                   + 2*dHpdx*Hp*dgdx*v_b + pow(Hp,2)*ddgddx*v_b + pow(Hp,2)*dgdx*dv_bdx
      //                   + 2*dHpdx*Hp*g*dv_bdx + pow(Hp,2)*dgdx*dv_bdx + pow(Hp,2)*g*ddv_bddx);
      ST_array[index] = 1.0;

      // Polarization source
      if(Constants.polarization){
        SE_array[index] = 0.0;
      }
    }
  }

  // Spline the source functions
  ST_spline.create(x_array, k_array, ST_array, "Source_Temp_x_k");
  if(Constants.polarization){
    SE_spline.create(x_array, k_array, SE_array, "Source_Pol_x_k");
  }

  Utils::EndTiming("source");
}

//====================================================
// The right hand side of the perturbations ODE
// in the tight coupling regime
//====================================================

// Derivatives in the tight coupling regime
int Perturbations::rhs_tight_coupling_ode(double x, double k, const double *y, double *dydx){
  double c         = Constants.c;
  double H0        = Constants.H0_over_h*cosmo->get_h();
  double Hp        = cosmo->Hp_of_x(x);
  double dHpdx     = cosmo->dHpdx_of_x(x);
  double Omega_R   = cosmo->get_OmegaR();
  double Omega_B   = cosmo->get_OmegaB();
  double Omega_CDM = cosmo->get_OmegaCDM();
  double dtaudx    = rec->dtaudx_of_x(x);
  double ddtauddx  = rec->ddtauddx_of_x(x);

  //=============================================================================
  // Compute where in the y / dydx array each component belongs
  // This is just an example of how to do it to make it easier
  // Feel free to organize the component any way you like
  //=============================================================================
  
  // For integration of perturbations in tight coupling regime (Only 2 photon multipoles + neutrinos needed)
  const int n_ell_theta_tc      = Constants.n_ell_theta_tc;
  const int n_ell_neutrinos_tc  = Constants.n_ell_neutrinos_tc;
  const bool neutrinos          = Constants.neutrinos;

  // The different quantities in the y array
  const double &delta_cdm       =  y[Constants.ind_deltacdm_tc];
  const double &delta_b         =  y[Constants.ind_deltab_tc];
  const double &v_cdm           =  y[Constants.ind_vcdm_tc];
  const double &v_b             =  y[Constants.ind_vb_tc];
  const double &Phi             =  y[Constants.ind_Phi_tc];
  const double *Theta           = &y[Constants.ind_start_theta_tc];
  const double *Nu              = &y[Constants.ind_start_nu_tc];

  // References to the quantities we are going to set in the dydx array
  double &ddelta_cdmdx    =  dydx[Constants.ind_deltacdm_tc];
  double &ddelta_bdx      =  dydx[Constants.ind_deltab_tc];
  double &dv_cdmdx        =  dydx[Constants.ind_vcdm_tc];
  double &dv_bdx          =  dydx[Constants.ind_vb_tc];
  double &dPhidx          =  dydx[Constants.ind_Phi_tc];
  double *dThetadx        = &dydx[Constants.ind_start_theta_tc];
  double *dNudx           = &dydx[Constants.ind_start_nu_tc];

  //=============================================================================
  // TODO: fill in the expressions for all the derivatives
  //=============================================================================
  double Theta2 = -4.0/9.0*c*k/(Hp*dtaudx)*Theta[1];
  double Psi    = -Phi - 12.0*pow(H0/(c*k*exp(x)),2)*Omega_R*Theta2;
  double R      = 4.0*Omega_R/(3.0*Omega_B*exp(x));

  // SET: Scalar quantities (Phi, delta, v, ...)
  dPhidx        = Psi - pow(c*k/Hp,2)/3.0*Phi + 0.5*pow(H0/Hp,2)*(Omega_CDM*exp(-x)*delta_cdm + Omega_B*exp(-x)*delta_b + 4*Omega_R*exp(-2.0*x)*Theta[0]);
  ddelta_cdmdx  = c*k/Hp*v_cdm - 3.0*dPhidx;
  dv_cdmdx      = -v_cdm - c*k/Hp*Psi;
  ddelta_bdx    = c*k/Hp*v_b - 3.0*dPhidx;
  dThetadx[0]   = -c*k/Hp*Theta[1] - dPhidx;
  double q      = (-((1.0-R)*dtaudx + (1.0+R)*ddtauddx)*(3.0*Theta[1] + v_b) - c*k/Hp*Psi + (1.0 - dHpdx/Hp)*c*k/Hp*(-Theta[0] + 2.0*Theta2) - c*k/Hp*dThetadx[0])/((1.0+R)*dtaudx + dHpdx/Hp - 1.0);
  dv_bdx        = 1.0/(1.0+R)*(-v_b - c*k/Hp*Psi + R*(q + c*k/Hp*(-Theta[0] + 2.0*Theta2) - c*k/Hp*Psi));
  dThetadx[1]   = 1.0/3.0*(q - dv_bdx);

  // SET: Neutrino mutlipoles (Nu_ell)
  if(neutrinos){
    // ...
    // ...
    // ...
  }

  return GSL_SUCCESS;
}

//====================================================
// The right hand side of the full ODE
//====================================================

int Perturbations::rhs_full_ode(double x, double k, const double *y, double *dydx){
  
  //=============================================================================
  // Compute where in the y / dydx array each component belongs
  // This is just an example of how to do it to make it easier
  // Feel free to organize the component any way you like
  //=============================================================================

  // Index and number of the different quantities
  const int n_ell_theta         = Constants.n_ell_theta;
  const int n_ell_thetap        = Constants.n_ell_thetap;
  const int n_ell_neutrinos     = Constants.n_ell_neutrinos;
  const bool polarization       = Constants.polarization;
  const bool neutrinos          = Constants.neutrinos;

  // The different quantities in the y array
  const double &delta_cdm       =  y[Constants.ind_deltacdm];
  const double &delta_b         =  y[Constants.ind_deltab];
  const double &v_cdm           =  y[Constants.ind_vcdm];
  const double &v_b             =  y[Constants.ind_vb];
  const double &Phi             =  y[Constants.ind_Phi];
  const double *Theta           = &y[Constants.ind_start_theta];
  const double *Theta_p         = &y[Constants.ind_start_thetap];
  const double *Nu              = &y[Constants.ind_start_nu];

  // References to the quantities we are going to set in the dydx array
  double &ddelta_cdmdx    =  dydx[Constants.ind_deltacdm];
  double &ddelta_bdx      =  dydx[Constants.ind_deltab];
  double &dv_cdmdx        =  dydx[Constants.ind_vcdm];
  double &dv_bdx          =  dydx[Constants.ind_vb];
  double &dPhidx          =  dydx[Constants.ind_Phi];
  double *dThetadx        = &dydx[Constants.ind_start_theta];
  double *dTheta_pdx      = &dydx[Constants.ind_start_thetap];
  double *dNudx           = &dydx[Constants.ind_start_nu];

  // Cosmological parameters and variables
  double c         = Constants.c;
  double Hp        = cosmo->Hp_of_x(x);
  double H0        = Constants.H0_over_h*cosmo->get_h();
  double Omega_R   = cosmo->get_OmegaR();
  double Omega_B   = cosmo->get_OmegaB();
  double Omega_CDM = cosmo->get_OmegaCDM();
  double eta       = cosmo->eta_of_x(x);

  // Recombination variables
  double dtaudx    = rec->dtaudx_of_x(x);

  //=============================================================================
  // TODO: fill in the expressions for all the derivatives
  //=============================================================================

  // SET: Scalar quantities (Phi, delta, v, ...)
  double Psi    = -Phi - 12.0*pow(H0/(c*k*exp(x)),2)*Omega_R*Theta[2];
  double R      = 4.0*Omega_R/(3.0*Omega_B*exp(x));
  dPhidx        = Psi - pow(c*k/Hp,2)/3.0*Phi + 0.5*pow(H0/Hp,2)*(Omega_CDM*exp(-x)*delta_cdm + Omega_B*exp(-x)*delta_b + 4.0*Omega_R*exp(-2.0*x)*Theta[0]);
  ddelta_cdmdx  = c*k/Hp*v_cdm - 3.0*dPhidx;
  dv_cdmdx      = -v_cdm - c*k/Hp*Psi;
  ddelta_bdx    = c*k/Hp*v_b - 3.0*dPhidx;
  dv_bdx        = -v_b -c*k/Hp*Psi + dtaudx*R*(3.0*Theta[1] + v_b);

  // SET: Photon multipoles (Theta_ell)
  dThetadx[0] = -c*k/Hp*Theta[1] - dPhidx;
  dThetadx[1] = c*k/(3.0*Hp)*Theta[0] - 2.0*c*k/(3.0*Hp)*Theta[2] + c*k/(3.0*Hp)*Psi + dtaudx*(Theta[1] + v_b/3.0);
  for(int ell = 2; ell < n_ell_theta - 1; ell++){
    dThetadx[ell] = ell*c*k/((2.0*ell+1.0)*Hp)*Theta[ell-1] - (ell+1.0)*c*k/((2.0*ell+1.0)*Hp)*Theta[ell+1] + dtaudx*Theta[ell];
    if(ell == 2){
      dThetadx[ell] = dThetadx[ell] - dtaudx*0.1*Theta[2];
    }
  }
  int ellmax       = n_ell_theta - 1;
  dThetadx[ellmax] = c*k/Hp*Theta[ellmax-1] - c*(ellmax+1.0)/(Hp*eta)*Theta[ellmax] + dtaudx*Theta[ellmax];

  // SET: Photon polarization multipoles (Theta_p_ell)
  if(polarization){
    // ...
    // ...
    // ...
  }

  // SET: Neutrino mutlipoles (Nu_ell)
  if(neutrinos){
    // ...
    // ...
    // ...
  }

  return GSL_SUCCESS;
}

//====================================================
// Get methods
//====================================================

double Perturbations::get_delta_cdm(const double x, const double k) const{
  return delta_cdm_spline(x,k);
}
double Perturbations::get_delta_b(const double x, const double k) const{
  return delta_b_spline(x,k);
}
double Perturbations::get_v_cdm(const double x, const double k) const{
  return v_cdm_spline(x,k);
}
double Perturbations::get_v_b(const double x, const double k) const{
  return v_b_spline(x,k);
}
double Perturbations::get_Phi(const double x, const double k) const{
  return Phi_spline(x,k);
}
double Perturbations::get_Psi(const double x, const double k) const{
  return Psi_spline(x,k);
}
double Perturbations::get_Pi(const double x, const double k) const{
  return Pi_spline(x,k);
}
double Perturbations::get_Source_T(const double x, const double k) const{
  return ST_spline(x,k);
}
double Perturbations::get_Source_E(const double x, const double k) const{
  return SE_spline(x,k);
}
double Perturbations::get_Theta(const double x, const double k, const int ell) const{
  return Theta_spline[ell](x,k);
}
double Perturbations::get_Theta_p(const double x, const double k, const int ell) const{
  return Theta_p_spline[ell](x,k);
}
double Perturbations::get_Nu(const double x, const double k, const int ell) const{
  return Nu_spline[ell](x,k);
}

//====================================================
// Print some useful info about the class
//====================================================

void Perturbations::info() const{
  std::cout << "\n";
  std::cout << "Info about perturbations class:\n";
  std::cout << "x_start:       " << x_start                << "\n";
  std::cout << "x_end:         " << x_end                  << "\n";
  std::cout << "n_x:     " << n_x              << "\n";
  std::cout << "k_min (1/Mpc): " << k_min * Constants.Mpc  << "\n";
  std::cout << "k_max (1/Mpc): " << k_max * Constants.Mpc  << "\n";
  std::cout << "n_k:     " << n_k              << "\n";
  if(Constants.polarization)
    std::cout << "We include polarization\n";
  else
    std::cout << "We do not include polarization\n";
  if(Constants.neutrinos)
    std::cout << "We include neutrinos\n";
  else
    std::cout << "We do not include neutrinos\n";

  std::cout << "Information about the perturbation system:\n";
  std::cout << "ind_deltacdm:       " << Constants.ind_deltacdm         << "\n";
  std::cout << "ind_deltab:         " << Constants.ind_deltab           << "\n";
  std::cout << "ind_v_cdm:          " << Constants.ind_vcdm             << "\n";
  std::cout << "ind_v_b:            " << Constants.ind_vb               << "\n";
  std::cout << "ind_Phi:            " << Constants.ind_Phi              << "\n";
  std::cout << "ind_start_theta:    " << Constants.ind_start_theta      << "\n";
  std::cout << "n_ell_theta:        " << Constants.n_ell_theta          << "\n";
  if(Constants.polarization){
    std::cout << "ind_start_thetap:   " << Constants.ind_start_thetap   << "\n";
    std::cout << "n_ell_thetap:       " << Constants.n_ell_thetap       << "\n";
  }
  if(Constants.neutrinos){
    std::cout << "ind_start_nu:       " << Constants.ind_start_nu       << "\n";
    std::cout << "n_ell_neutrinos     " << Constants.n_ell_neutrinos    << "\n";
  }
  std::cout << "n_ell_tot_full:     " << Constants.n_ell_tot_full       << "\n";

  std::cout << "Information about the perturbation system in tight coupling:\n";
  std::cout << "ind_deltacdm:       " << Constants.ind_deltacdm_tc      << "\n";
  std::cout << "ind_deltab:         " << Constants.ind_deltab_tc        << "\n";
  std::cout << "ind_v_cdm:          " << Constants.ind_vcdm_tc          << "\n";
  std::cout << "ind_v_b:            " << Constants.ind_vb_tc            << "\n";
  std::cout << "ind_Phi:            " << Constants.ind_Phi_tc           << "\n";
  std::cout << "ind_start_theta:    " << Constants.ind_start_theta_tc   << "\n";
  std::cout << "n_ell_theta:        " << Constants.n_ell_theta_tc       << "\n";
  if(Constants.neutrinos){
    std::cout << "ind_start_nu:       " << Constants.ind_start_nu_tc    << "\n";
    std::cout << "n_ell_neutrinos     " << Constants.n_ell_neutrinos_tc << "\n";
  }
  std::cout << "n_ell_tot_tc:       " << Constants.n_ell_tot_tc         << "\n";
  std::cout << std::endl;
}

//====================================================
// Output some results to file for a given value of k
//====================================================

void Perturbations::output(const double k, const std::string filename) const{
  std::ofstream fp(filename.c_str());
  const int npts = 5000;
  auto x_array = Utils::linspace(x_start, x_end, npts);
  auto print_data = [&] (const double x) {
    double arg = k * (cosmo->eta_of_x(0.0) - cosmo->eta_of_x(x));
    fp << x                  << " ";
    fp << get_Theta(x,k,0)   << " ";
    fp << get_Theta(x,k,1)   << " ";
    fp << get_Theta(x,k,2)   << " ";
    fp << get_Phi(x,k)       << " ";
    fp << get_Psi(x,k)       << " ";
    fp << get_delta_b(x,k)   << " ";
    fp << get_delta_cdm(x,k) << " ";
    fp << get_v_b(x,k)       << " ";
    fp << get_v_cdm(x,k)     << " ";
    fp << get_Source_T(x,k)  << " ";
    fp << get_Source_T(x,k) * Utils::j_ell(5,   arg)           << " ";
    fp << get_Source_T(x,k) * Utils::j_ell(50,  arg)           << " ";
    fp << get_Source_T(x,k) * Utils::j_ell(500, arg)           << " ";
    fp << "\n";
  };
  std::for_each(x_array.begin(), x_array.end(), print_data);
}

