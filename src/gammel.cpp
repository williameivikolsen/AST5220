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
  int n = 1000;
  Vector x_array = Utils::linspace(x_start, x_end, n);
  Vector Xe_arr(npts_rec_arrays);     // ER DETTE RIKTIG???
  Vector ne_arr(npts_rec_arrays);     // ------------------
  int i_peebles = -1;                      // Index i when exiting Saha regime

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
      Xe_arr[i] = Xe_current;
      ne_arr[i] = ne_current;

    } else {
      i_peebles = i;

      //==============================================================
      // TODO: Compute X_e from current time til today by solving 
      // the Peebles equation (NB: if you solve all in one go remember to
      // exit the for-loop!)
      // Implement rhs_peebles_ode
      //==============================================================
      //...
      //...
      Vector x_peebles = Utils::linspace(x_array[i], x_end, n-i)

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
      Vector Xe_ic{Xe_saha_limit};

      // Solve the ODE
      ODESolver ode;
      ode.solve(dXedx, x_peebles, Xe_ic);

      // Get the solution
      auto Xe_peebles = ode.get_data_by_component(0);
      break;
    }
  }

  if(i_peebles != -1){
    for(int i = i_peebles; i < npts_rec_arrays; i++){
      Xe_arr[i] = Xe_peebles[i];
      ne_arr[i] = Xe_peebles[i];        // MÅ ENDRES!!!
    }
  }
  //=============================================================================
  // TODO: Spline the result. Implement and make sure the Xe_of_x, ne_of_x 
  // functions are working
  //=============================================================================
  //...
  //...

  Xe_of_x_spline.create(x_array, Xe_arr, "Xe_of_x");
  ne_of_x_spline.create(x_array, ne_arr, "ne_of_x");

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
  const double OmegaM      = cosmo->get_OmegaM();
  const double OmegaR      = cosmo->get_OmegaR();
  const double OmegaRtot   = cosmo->get_OmegaRtot();
  const double OmegaNu     = cosmo->get_OmegaNu();
  const double OmegaCDM    = cosmo->get_OmegaCDM();
  const double OmegaLambda = cosmo->get_OmegaLambda();
  const double OmegaK      = cosmo->get_OmegaK();
  const double OmegaMnu    = cosmo->get_OmegaMnu();
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
  double nb = 3*pow(H0,2)*OmegaB(0)/(8*m_PI*G*m_H*pow(a,3));
  double Tb = TCMB/a;
  double C = 1/nb*pow(m_e*k_b*Tb/(2*m_PI*pow(hbar,2)),1.5)*exp(-epsilon_0/(k_b*Tb));
  Xe = -C/2 + sqrt(pow(C, 2) - 4*C);    // Egentlig +- ???
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
  const double OmegaM      = cosmo->get_OmegaM();
  const double OmegaR      = cosmo->get_OmegaR();
  const double OmegaRtot   = cosmo->get_OmegaRtot();
  const double OmegaNu     = cosmo->get_OmegaNu();
  const double OmegaCDM    = cosmo->get_OmegaCDM();
  const double OmegaLambda = cosmo->get_OmegaLambda();
  const double OmegaK      = cosmo->get_OmegaK();
  const double OmegaMnu    = cosmo->get_OmegaMnu();
  const double H0          = cosmo->get_H0();
  const double h           = cosmo->get_h();
  const double Neff        = cosmo->get_Neff();
  const double TCMB        = cosmo->get_TCMB();
  const double H           = cosmo->H_of_x();

  //=============================================================================
  // TODO: Write the expression for dXedx
  //=============================================================================
  //...
  //...
  double Tb = TCMB/a;
  double nb = 3*pow(H0,2)*OmegaB(0)/(8*m_PI*G*m_H*pow(a,3));
  double n1s = (1 - Xe)*n_H;
  double lambda_alpha = H*pow(3*epsilon_0,3)/(pow(8*m_PI,2)*n1s);
  double phi2 = 0.448*log(epsilon_0/Tb);
  double alpha2 = 64*m_PI*pow(alpha,2)/pow(m_e,2)*sqrt(epsilon_0/(27*m_PI*Tb))*phi2;
  double beta = alpha2*pow(m_e*Tb/(2*m_PI),1.5)*exp(-epsilon_0/Tb);
  double beta2 = beta*exp(3*epsilon_0/(4*Tb));
  double Cr = (lambda_2s1s + lambda_alpha)/(lambda_2s1s + lambda_alpha + beta2);
  
  // BØR H VÆRE EN VARIABEL???

  dXedx[0] = Cr/H(x)*(beta*(1-Xe) - nb*alpha2*pow(Xe,2));

  return GSL_SUCCESS;
}

//====================================================
// Solve for the optical depth tau, compute the 
// visibility function and spline the result
//====================================================

void RecombinationHistory::solve_for_optical_depth_tau(){
  Utils::StartTiming("opticaldepth");

  // Set up x-arrays to integrate over. We split into three regions as we need extra points in reionisation
  const int npts = 1000;
  Vector x_array = Utils::linspace(x_start, x_end, npts);

  const double H = cosmo->H_of_x();

  // The ODE system dtau/dx, dtau_noreion/dx and dtau_baryon/dx
  ODEFunction dtaudx = [&](double x, const double *tau, double *dtaudx){

    //=============================================================================
    // TODO: Write the expression for dtaudx
    //=============================================================================
    //...
    //...

    // Set the derivative for photon optical depth
    dtaudx[0] = -Constants.c*ne_of_x(x)*Constants.sigma_T/H(x);

    return GSL_SUCCESS;
  };

  //=============================================================================
  // TODO: Set up and solve the ODE and make tau splines
  //=============================================================================
  //...
  //...

  // Initial condition
  Vector tau_ic{0};

  // Solve the ODE
  ODESolver ode;
  ode.solve(dtaudx, x_array, tau_ic);
  auto tau_array = ode.get_data_by_component(0);

  // Make spline
  tau_of_x_spline.create(x_array, tau_array, "tau_of_x");

  //=============================================================================
  // TODO: Compute visibility functions and spline everything
  //=============================================================================
  //...
  //...

  Vector g_tilde_arr(npts_rec_arrays); 
  for(int i = 0; i < npts_rec_arrays; i++){
    g_tilde_arr[i] = Constants.c*ne_of_x(x_array[i])*Constants.sigma_T/cosmo->H_of_x(x_array[i])*exp(-tau_of_x_spline(x_array[i]));
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

  return tau_of_x_spline.deriv_x;
}

double RecombinationHistory::ddtauddx_of_x(double x) const{

  //=============================================================================
  // TODO: Implement
  //=============================================================================
  //...
  //...

  return tau_of_x_spline.deriv_xx;
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

  return g_tilde_of_x_spline.deriv_x;
}

double RecombinationHistory::ddgddx_tilde_of_x(double x) const{

  //=============================================================================
  // TODO: Implement
  //=============================================================================
  //...
  //...

  return g_tilde_of_x_spline.deriv_xx;
}

double RecombinationHistory::Xe_of_x(double x) const{

  //=============================================================================
  // TODO: Implement
  //=============================================================================
  //...
  //...

  return Xe_of_x_spline(x);
}

double RecombinationHistory::ne_of_x(double x) const{

  //=============================================================================
  // TODO: Implement
  //=============================================================================
  //...
  //...

  return ne_of_x_spline(x);
}

double RecombinationHistory::get_Yp() const{
  return Yp;
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
  const int npts       = 5000;
  const double x_min   = x_start;
  const double x_max   = x_end;

  Vector x_array = Utils::linspace(x_min, x_max, npts);
  auto print_data = [&] (const double x) {
    fp << x                    << " ";
    fp << Xe_of_x(x)           << " ";
    fp << ne_of_x(x)           << " ";
    fp << tau_of_x(x)          << " ";
    fp << dtaudx_of_x(x)       << " ";
    fp << ddtauddx_of_x(x)     << " ";
    fp << g_tilde_of_x(x)      << " ";
    fp << dgdx_tilde_of_x(x)   << " ";
    fp << ddgddx_tilde_of_x(x) << " ";
    fp << "\n";
  };
  std::for_each(x_array.begin(), x_array.end(), print_data);
}

