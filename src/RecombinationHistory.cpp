#include"RecombinationHistory.h"

//====================================================
// Constructors
//====================================================
   
RecombinationHistory::RecombinationHistory(
    BackgroundCosmology *cosmo, 
    double Yp) :
  cosmo(cosmo),
  Yp(Yp)
{}

//====================================================
// Do all the solving we need to do
//====================================================

void RecombinationHistory::solve(){
    
  // Compute and spline Xe, ne
  solve_number_density_electrons();
   
  // Compute and spline tau, dtaudx, ddtauddx, g, dgdx, ddgddx, ...
  solve_for_optical_depth_tau();
}

//====================================================
// Solve for X_e and n_e using Saha and Peebles and spline the result
//====================================================

void RecombinationHistory::solve_number_density_electrons(){
  Utils::StartTiming("Xe");
  
  //=============================================================================
  // TODO: Set up x-array and make arrays to store X_e(x) and n_e(x) on
  //=============================================================================
  Vector x_array = Utils::linspace(x_start, x_end, npts_rec_arrays);
  Vector Xe_arr(npts_rec_arrays);
  Vector ne_arr(npts_rec_arrays);

  // Physical constants and cosmological parameters
  const double G           = Constants.G;
  const double OmegaB      = cosmo->get_OmegaB();
  const double H0          = cosmo->get_H0();
  const double m_H         = Constants.m_H;
  double rho_c0            = 3*pow(H0,2)/(8*M_PI*G);

  // Calculate recombination history
  bool saha_regime = true;
  for(int i = 0; i < npts_rec_arrays; i++){

    //==============================================================
    // TODO: Get X_e from solving the Saha equation so
    // implement the function electron_fraction_from_saha_equation
    //==============================================================
    auto Xe_ne_data = electron_fraction_from_saha_equation(x_array[i]);

    // Electron fraction and number density
    const double Xe_current = Xe_ne_data.first;
    const double ne_current = Xe_ne_data.second;

    // Are we still in the Saha regime?
    if(Xe_current < Xe_saha_limit)
      saha_regime = false;

    if(saha_regime){
      
      //=============================================================================
      // TODO: Store the result we got from the Saha equation
      //=============================================================================
      //...
      //...
      Xe_arr[i] = Xe_current;
      ne_arr[i] = ne_current;

    } else {

      //==============================================================
      // TODO: Compute X_e from current time til today by solving 
      // the Peebles equation (NB: if you solve all in one go remember to
      // exit the for-loop!)
      // Implement rhs_peebles_ode
      //==============================================================
      //...
      //...

      // The Peebles ODE equation
      ODESolver peebles_Xe_ode;
      ODEFunction dXedx = [&](double x, const double *Xe, double *dXedx){
        return rhs_peebles_ode(x, Xe, dXedx);
      };
      
      //=============================================================================
      // TODO: Set up IC, solve the ODE and fetch the result 
      //=============================================================================
      //...
      //...

      Vector x_peebles{x_array[i-1], x_array[i]};
      Vector Xe_ic{Xe_arr[i-1]};
      peebles_Xe_ode.solve(dXedx, x_peebles, Xe_ic);
      Xe_arr[i] = peebles_Xe_ode.get_data_by_component(0)[1];
      double a = exp(x_array[i]);
      double nb = OmegaB*rho_c0/(m_H*pow(a,3));
      ne_arr[i] = nb*Xe_arr[i];
    }
  }

  //=============================================================================
  // TODO: Spline the result. Implement and make sure the Xe_of_x, ne_of_x 
  // functions are working
  //=============================================================================
  //...
  //...
  auto log_Xe_arr = log(Xe_arr);
  log_Xe_of_x_spline.create(x_array, log_Xe_arr, "Xe");
  auto log_ne_arr = log(ne_arr);
  log_ne_of_x_spline.create(x_array, log_ne_arr, "ne");

  Utils::EndTiming("Xe");
}

//====================================================
// Solve the Saha equation to get ne and Xe
//====================================================
std::pair<double,double> RecombinationHistory::electron_fraction_from_saha_equation(double x) const{
  const double a           = exp(x);
 
  // Physical constants
  const double k_b         = Constants.k_b;
  const double G           = Constants.G;
  const double m_e         = Constants.m_e;
  const double hbar        = Constants.hbar;
  const double m_H         = Constants.m_H;
  const double epsilon_0   = Constants.epsilon_0;
  const double H0_over_h   = Constants.H0_over_h;

  // Fetch cosmological parameters
  //const double OmegaB      = cosmo->get_OmegaB();
  //...
  //...
  const double OmegaB      = cosmo->get_OmegaB();
  const double OmegaR      = cosmo->get_OmegaR();
  const double OmegaNu     = cosmo->get_OmegaNu();
  const double OmegaCDM    = cosmo->get_OmegaCDM();
  const double OmegaLambda = cosmo->get_OmegaLambda();
  const double OmegaK      = cosmo->get_OmegaK();
  const double H0          = cosmo->get_H0();
  const double h           = cosmo->get_h();
  const double Neff        = cosmo->get_Neff();
  const double TCMB        = cosmo->get_TCMB();

  // Electron fraction and number density
  double Xe = 0.0;
  double ne = 0.0;
  
  //=============================================================================
  // TODO: Compute Xe and ne from the Saha equation
  //=============================================================================
  //...
  //...
  double rho_c0 = 3*pow(H0,2)/(8*M_PI*G);
  double nb = OmegaB*rho_c0/(m_H*pow(a,3));
  double Tb = TCMB/a;
  double C = 1/nb*pow(m_e*k_b*Tb/(2*M_PI*pow(hbar,2)),1.5)*exp(-epsilon_0/(k_b*Tb));
  if(C > 1e10){
    Xe = 1;
  }
  else{
    Xe = -C/2 + sqrt(pow(C, 2) + 4*C)/2;    // Egentlig +- ???
  }
  ne = Xe*nb;

  return std::pair<double,double>(Xe, ne);
}

//====================================================
// The right hand side of the dXedx Peebles ODE
//====================================================
int RecombinationHistory::rhs_peebles_ode(double x, const double *Xe, double *dXedx){

  // Current value of a and X_e
  const double X_e         = Xe[0];
  const double a           = exp(x);

  // Physical constants in SI units
  const double k_b         = Constants.k_b;
  const double G           = Constants.G;
  const double c           = Constants.c;
  const double m_e         = Constants.m_e;
  const double hbar        = Constants.hbar;
  const double m_H         = Constants.m_H;
  const double sigma_T     = Constants.sigma_T;
  const double lambda_2s1s = Constants.lambda_2s1s;
  const double epsilon_0   = Constants.epsilon_0;
  const double alpha       = Constants.alpha;

  // Cosmological parameters
  // const double OmegaB      = cosmo->get_OmegaB();
  // ...
  // ...
  const double OmegaB      = cosmo->get_OmegaB();
  const double OmegaR      = cosmo->get_OmegaR();
  const double OmegaNu     = cosmo->get_OmegaNu();
  const double OmegaCDM    = cosmo->get_OmegaCDM();
  const double OmegaLambda = cosmo->get_OmegaLambda();
  const double OmegaK      = cosmo->get_OmegaK();
  const double H0          = cosmo->get_H0();
  const double h           = cosmo->get_h();
  const double Neff        = cosmo->get_Neff();
  const double TCMB        = cosmo->get_TCMB();
  const double H           = cosmo->H_of_x(x);

  //=============================================================================
  // TODO: Write the expression for dXedx
  //=============================================================================
  //...
  //...
  double rho_c0          = 3*pow(H0,2)/(8*M_PI*G);
  double nb              = OmegaB*rho_c0/(m_H*pow(a,3));
  double Tb              = TCMB/a;
  double n1s             = (1 - X_e)*nb;
  double lambda_alpha    = H*pow(3*epsilon_0/(hbar*c),3)/(pow(8*M_PI,2)*n1s);
  double phi2            = 0.448*log(epsilon_0/(k_b*Tb));
  double alpha2          = 64*M_PI*pow(alpha/m_e,2)*sqrt(epsilon_0/(27*M_PI*k_b*Tb))*phi2*pow(hbar,2)/c;
  double beta            = alpha2*pow(m_e*k_b*Tb/(2*M_PI*pow(hbar,2)),1.5)*exp(-epsilon_0/(k_b*Tb));
  double beta2           = alpha2*pow(m_e*k_b*Tb/(2*M_PI*pow(hbar,2)),1.5)*exp(-epsilon_0/(4*k_b*Tb));
  double Cr              = (lambda_2s1s + lambda_alpha)/(lambda_2s1s + lambda_alpha + beta2);
  
  dXedx[0] = Cr/H*(beta*(1-X_e) - nb*alpha2*pow(X_e,2));

  return GSL_SUCCESS;
}

//====================================================
// Solve for the optical depth tau, compute the 
// visibility function and spline the result
//====================================================

void RecombinationHistory::solve_for_optical_depth_tau(){
  Utils::StartTiming("opticaldepth");

  // Set up x-arrays to integrate over. We split into three regions as we need extra points in reionisation
  const int npts = npts_rec_arrays;
  Vector x_array = Utils::linspace(x_start, x_end, npts);

  // Find index of x value closest to zero:
  int i_zero = 0;
  for(int i = 1; i < npts; i++){
    if(abs(x_array[i]) < abs(x_array[i_zero])){
      i_zero = i;
    }
  }

  // Physical constants
  const double c       = Constants.c;
  const double sigma_T = Constants.sigma_T;

  // The ODE system dtau/dx, dtau_noreion/dx and dtau_baryon/dx
  ODEFunction dtaudx = [&](double x, const double *tau, double *dtaudx){

    //=============================================================================
    // TODO: Write the expression for dtaudx
    //=============================================================================
    //...
    //...
    double H       = cosmo->H_of_x(x);

    // Set the derivative for photon optical depth
    dtaudx[0] = -c*exp(log_ne_of_x(x))*sigma_T/H;

    return GSL_SUCCESS;
  };

  //=============================================================================
  // TODO: Set up and solve the ODE and make tau splines
  //=============================================================================
  //...
  //...
  Vector tau_ic{0};                                                                   // Sjekk at initialbetingelsene her er riktig!!!
  ODESolver tau_ode;
  tau_ode.solve(dtaudx, x_array, tau_ic);
  auto tau_array = tau_ode.get_data_by_component(0);
  Vector dtaudx_array(npts);

  // Subtract tau(x = 0) to make sure tau(x = 0) = 0 and store the value of the derivative
  double tau0 = tau_array[i_zero];
  for(int i = 0; i < npts; i++){
    tau_array[i] = tau_array[i] - tau0;
    dtaudx_array[i] = -c*exp(log_ne_of_x(x_array[i]))*sigma_T/cosmo->H_of_x(x_array[i]);
  }
  tau_of_x_spline.create(x_array, tau_array, "tau");
  dtaudx_spline.create(x_array, dtaudx_array, "dtaudx");

  //=============================================================================
  // TODO: Compute visibility functions and spline everything
  //=============================================================================
  //...
  //...
  Vector g_tilde_arr(npts);
  double x;
  double H;
  for(int i = 0; i < npts; i++){
    x = x_array[i];
    H = cosmo->H_of_x(x);
    g_tilde_arr[i] = c*exp(log_ne_of_x(x))*sigma_T/H*exp(-tau_of_x_spline(x));
  }
  g_tilde_of_x_spline.create(x_array, g_tilde_arr, "g_tilde_of_x");

  Utils::EndTiming("opticaldepth");
}

//====================================================
// Get methods
//====================================================

double RecombinationHistory::tau_of_x(double x) const{
  return tau_of_x_spline(x);
}

double RecombinationHistory::dtaudx_of_x(double x) const{

  //=============================================================================
  // TODO: Implement. Either from the tau-spline tau_of_x_spline.deriv_(x) or 
  // from a separate spline if you choose to do this
  //=============================================================================
  //...
  //...

  return dtaudx_spline(x);
}

double RecombinationHistory::ddtauddx_of_x(double x) const{

  //=============================================================================
  // TODO: Implement
  //=============================================================================
  //...
  //...

  return dtaudx_spline.deriv_x(x);
}

double RecombinationHistory::g_tilde_of_x(double x) const{
  return g_tilde_of_x_spline(x);
}

double RecombinationHistory::dgdx_tilde_of_x(double x) const{

  //=============================================================================
  // TODO: Implement
  //=============================================================================
  //...
  //...

  return g_tilde_of_x_spline.deriv_x(x);
}

double RecombinationHistory::ddgddx_tilde_of_x(double x) const{

  //=============================================================================
  // TODO: Implement
  //=============================================================================
  //...
  //...

  return g_tilde_of_x_spline.deriv_xx(x);
}

double RecombinationHistory::Xe_of_x(double x) const{

  //=============================================================================
  // TODO: Implement
  //=============================================================================
  //...
  //...

  return exp(log_Xe_of_x_spline(x));
}

double RecombinationHistory::log_ne_of_x(double x) const{

  //=============================================================================
  // TODO: Implement
  //=============================================================================
  //...
  //...

  return log_ne_of_x_spline(x);
}

double RecombinationHistory::get_Yp() const{
  return Yp;
}

double RecombinationHistory::get_x_dec() const{
  // Compute time of decoupling
  return Utils::binary_search_for_value(tau_of_x_spline, 1);
}

//====================================================
// Print some useful info about the class
//====================================================
void RecombinationHistory::info() const{
  std::cout << "\n";
  std::cout << "Info about recombination/reionization history class:\n";
  std::cout << "Yp:          " << Yp          << "\n";
  std::cout << std::endl;
} 

//====================================================
// Output the data computed to file
//====================================================
void RecombinationHistory::output(const std::string filename) const{
  std::ofstream fp(filename.c_str());
  const int npts       = 2000;
  const double x_min   = x_start;
  const double x_max   = x_end;

  Vector x_array = Utils::linspace(x_min, x_max, npts);
  auto print_data = [&] (const double x) {
    fp << x                    << " ";
    fp << Xe_of_x(x)           << " ";
    fp << log_ne_of_x(x)       << " ";
    fp << tau_of_x(x)          << " ";
    fp << dtaudx_of_x(x)       << " ";
    fp << ddtauddx_of_x(x)     << " ";
    fp << g_tilde_of_x(x)      << " ";
    fp << dgdx_tilde_of_x(x)   << " ";
    fp << ddgddx_tilde_of_x(x) << " ";
    fp << "\n";
  };
  std::for_each(x_array.begin(), x_array.end(), print_data);
  double x_decoupling = get_x_dec();
  double a_decoupling = exp(x_decoupling);
  double z_decoupling = 1/a_decoupling - 1;
  std::cout << "Time of decoupling: x = " << x_decoupling << std::endl;
  std::cout << "a = " << a_decoupling << std::endl;
  std::cout << "z = " << z_decoupling << std::endl;
  double x_rec = Utils::binary_search_for_value(log_Xe_of_x_spline, log(0.5));
  double a_rec = exp(x_rec);
  double z_rec = 1/a_rec - 1;
  std::cout << "Halfway in recombination: x = " << x_rec << std::endl;
  std::cout << "a = " << a_rec << std::endl;
  std::cout << "z = " << z_rec << std::endl;
  Vector Xe_saha(npts);
  for(int i = 0; i < npts; i++){
    Xe_saha[i] = electron_fraction_from_saha_equation(x_array[i]).first;
  }
  Spline Xe_saha_spline;
  Xe_saha_spline.create(x_array, Xe_saha, "Xe_saha");
  double x_saha = Utils::binary_search_for_value(Xe_saha_spline, 0.5);
  double a_saha = exp(x_saha);
  double z_saha = 1/a_saha - 1;
  std::cout << "Expected halfway by Saha: x = " << x_saha << std::endl;
  std::cout << "a = " << a_saha << std::endl;
  std::cout << "z = " << z_saha << std::endl;
  std::cout << "Freeze-out abundance of electrons: Xe = " << Xe_of_x(0) << std::endl;
}