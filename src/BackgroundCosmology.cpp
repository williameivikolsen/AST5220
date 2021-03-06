#include "BackgroundCosmology.h"

//====================================================
// Constructors
//====================================================
    
BackgroundCosmology::BackgroundCosmology(
    double h, 
    double OmegaB, 
    double OmegaCDM, 
    double OmegaK,
    double Neff, 
    double TCMB) :
  h(h),
  OmegaB(OmegaB),
  OmegaCDM(OmegaCDM),
  OmegaK(OmegaK),
  Neff(Neff), 
  TCMB(TCMB)
{

  //=============================================================================
  // TODO: Compute OmegaR, OmegaNu, OmegaK, H0, ...
  //=============================================================================
  H0 = Constants.H0_over_h*h;
  OmegaNu = Neff*7/8*pow(4.0/11.0, 4.0/3.0);
  OmegaR = 2*pow(M_PI,2)/30*pow(Constants.k_b*TCMB, 4)/(pow(Constants.hbar,3)*pow(Constants.c, 5))*8*M_PI*Constants.G/(3*pow(H0, 2));
  OmegaLambda = 1 - OmegaB - OmegaCDM - OmegaR - OmegaNu - OmegaK;
}

//====================================================
// Do all the solving. Compute eta(x)
//====================================================

// Solve the background
void BackgroundCosmology::solve(){
  Utils::StartTiming("Eta");
    
  //=============================================================================
  // TODO: Set the range of x and the number of points for the splines
  // For this Utils::linspace(x_start, x_end, npts) is useful
  //=============================================================================
  Vector x_array = Utils::linspace(x_start, x_end, 1000);

  // The ODE for deta/dx
  ODEFunction detadx = [&](double x, const double *eta, double *detadx){

    //=============================================================================
    // TODO: Set the rhs of the detadx ODE
    //=============================================================================

    detadx[0] = Constants.c/Hp_of_x(x);

    return GSL_SUCCESS;
  };

  //=============================================================================
  // TODO: Set the initial condition, set up the ODE system, solve and make
  // the spline eta_of_x_spline 
  //=============================================================================

  double etaini = Constants.c/Hp_of_x(x_start);
  Vector eta_ic{etaini};

  // Solve the ODE
  ODESolver ode;
  ode.solve(detadx, x_array, eta_ic);

  // Get the solution (we only have one component so index 0 holds the solution)
  auto eta_array = ode.get_data_by_component(0);
  
  // Make spline
  eta_of_x_spline.create(x_array, eta_array, "eta_of_x");

  Utils::EndTiming("Eta");

  // The ODE for dt/dx
  ODEFunction dtdx = [&](double x, const double *t, double *dtdx){

    dtdx[0] = 1/H_of_x(x);

    return GSL_SUCCESS;
  };

  double t_ini = 1/(2*H_of_x(x_start));
  Vector t_ic{t_ini};

  // Solve the ODE
  ode.solve(dtdx, x_array, t_ic);

  // Get the solution (we only have one component so index 0 holds the solution)
  auto t_array = ode.get_data_by_component(0);
  
  // Make spline
  t_of_x_spline.create(x_array, t_array, "t_of_x");
}

//====================================================
// Get methods
//====================================================

double BackgroundCosmology::H_of_x(double x) const{

  //=============================================================================
  // TODO: Implement...
  //=============================================================================

  return H0*sqrt((OmegaB + OmegaCDM)*exp(-3*x) + (OmegaR + OmegaNu)*exp(-4*x) + OmegaK*exp(-2*x) + OmegaLambda);
}

double BackgroundCosmology::Hp_of_x(double x) const{

  //=============================================================================
  // TODO: Implement...
  //=============================================================================

  return exp(x)*H_of_x(x);
}

double BackgroundCosmology::dHpdx_of_x(double x) const{

  //=============================================================================
  // TODO: Implement...
  //=============================================================================

  return pow(H0, 2)/(2*Hp_of_x(x))*(-(OmegaB + OmegaCDM)*exp(-x) -2*(OmegaR + OmegaNu)*exp(-2*x) + 2*OmegaLambda*exp(2*x));
}

double BackgroundCosmology::ddHpddx_of_x(double x) const{

  //=============================================================================
  // TODO: Implement...
  //=============================================================================
  //...
  //...

  return pow(H0,2)/(2*Hp_of_x(x))*((OmegaB + OmegaCDM)*exp(-x) + 4*(OmegaR + OmegaNu)*exp(-2*x) + 4*OmegaLambda*exp(2*x))
  - pow(H0,4)/(4*pow(Hp_of_x(x),3))*pow(-(OmegaB + OmegaCDM)*exp(-x) - 2*(OmegaR + OmegaNu)*exp(-2*x) + 2*OmegaLambda*exp(2*x), 2);
}

double BackgroundCosmology::get_OmegaB(double x) const{ 
  if(x == 0.0) return OmegaB;

  //=============================================================================
  // TODO: Implement...
  //=============================================================================

  return OmegaB/(exp(3*x)*pow((H_of_x(x)/H0),2));
}

double BackgroundCosmology::get_OmegaR(double x) const{ 
  if(x == 0.0) return OmegaR;

  //=============================================================================
  // TODO: Implement...
  //=============================================================================

  return OmegaR/(exp(4*x)*pow((H_of_x(x)/H0),2));
}

double BackgroundCosmology::get_OmegaNu(double x) const{ 
  if(x == 0.0) return OmegaNu;

  //=============================================================================
  // TODO: Implement...
  //=============================================================================

  return OmegaNu/(exp(4*x)*pow((H_of_x(x)/H0),2));
}

double BackgroundCosmology::get_OmegaCDM(double x) const{ 
  if(x == 0.0) return OmegaCDM;

  //=============================================================================
  // TODO: Implement...
  //=============================================================================

  return OmegaCDM/(exp(3*x)*pow((H_of_x(x)/H0),2));
}

double BackgroundCosmology::get_OmegaLambda(double x) const{ 
  if(x == 0.0) return OmegaLambda;

  //=============================================================================
  // TODO: Implement...
  //=============================================================================

  return OmegaLambda/pow((H_of_x(x)/H0),2);
}

double BackgroundCosmology::get_OmegaK(double x) const{ 
  if(x == 0.0) return OmegaK;

  //=============================================================================
  // TODO: Implement...
  //=============================================================================

  return OmegaK/(exp(2*x)*pow((H_of_x(x)/H0),2));
}

double BackgroundCosmology::eta_of_x(double x) const{
  return eta_of_x_spline(x);
}

double BackgroundCosmology::t_of_x(double x) const{
  return t_of_x_spline(x);
}

double BackgroundCosmology::get_H0() const{ 
  return H0; 
}

double BackgroundCosmology::get_h() const{ 
  return h; 
}

double BackgroundCosmology::get_Neff() const{ 
  return Neff; 
}

double BackgroundCosmology::get_TCMB() const{ 
  return TCMB; 
}

//====================================================
// Print out info about the class
//====================================================
void BackgroundCosmology::info() const{ 
  std::cout << "\n";
  std::cout << "Info about cosmology class:\n";
  std::cout << "OmegaB:      " << OmegaB      << "\n";
  std::cout << "OmegaCDM:    " << OmegaCDM    << "\n";
  std::cout << "OmegaLambda: " << OmegaLambda << "\n";
  std::cout << "OmegaK:      " << OmegaK      << "\n";
  std::cout << "OmegaNu:     " << OmegaNu     << "\n";
  std::cout << "OmegaR:      " << OmegaR      << "\n";
  std::cout << "Neff:        " << Neff        << "\n";
  std::cout << "h:           " << h           << "\n";
  std::cout << "TCMB:        " << TCMB        << "\n";
  std::cout << std::endl;
} 

//====================================================
// Output some data to file
//====================================================
void BackgroundCosmology::output(const std::string filename) const{
  const double x_min = x_start;
  const double x_max =  x_end;
  const int    n_pts =  1000;
  
  Vector x_array = Utils::linspace(x_min, x_max, n_pts);

  std::ofstream fp(filename.c_str());
  auto print_data = [&] (const double x) {
    fp << x                  << " ";
    fp << eta_of_x(x)        << " ";
    fp << Hp_of_x(x)         << " ";
    fp << dHpdx_of_x(x)      << " ";
    fp << ddHpddx_of_x(x)    << " ";
    fp << get_OmegaB(x)      << " ";
    fp << get_OmegaCDM(x)    << " ";
    fp << get_OmegaLambda(x) << " ";
    fp << get_OmegaR(x)      << " ";
    fp << get_OmegaNu(x)     << " ";
    fp << get_OmegaK(x)      << " ";
    fp <<"\n";
  };
  std::for_each(x_array.begin(), x_array.end(), print_data);
  std::cout << "Age of Universe:" << t_of_x(0.0)/(365.2425*24*60*60) << std::endl;
  std::cout << "eta_0:" << eta_of_x(0.0) << std::endl;
  double a_r  = OmegaR/(OmegaB + OmegaCDM);
  double k_eq = Hp_of_x(log(a_r))*Constants.Mpc/Constants.c;
  std::cout << "k_eq: " << k_eq/h << " h/Mpc" << std::endl;
}
