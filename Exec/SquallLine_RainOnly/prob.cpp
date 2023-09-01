#include "prob.H"
#include "prob_common.H"

#include "EOS.H"
#include "AMReX_ParmParse.H"
#include "AMReX_MultiFab.H"
#include "IndexDefines.H"
#include "TileNoZ.H"

using namespace amrex;

ProbParm parms;

double p0 = 1e5;
double T0 = 300.0;
double Cp = 1004.0;
double Rd = 287.05;
double Rv = 461.50;

double zmin = 0.0;
double zmax = 12800.0;
int nz = 64;

double z_tr = 12000.0;
double g = 9.81;

double height = 1200.0;

double delz = (zmax-zmin)/nz;
double Rv_by_Rd = Rv/Rd;

double compute_theta(double z)
{
    double theta_0=300.0, theta_tr=343.0, T_tr=213.0;
    if(z <= z_tr){
        return theta_0 + (theta_tr - theta_0)*std::pow(z/z_tr,1.25);
    }
    else{
        return theta_tr*exp(g/(Cp*T_tr)*(z - z_tr));
    }

}

double compute_saturation_pressure(const double T_b)
{

    double p_s = exp(34.494 - 4924.99/(T_b - 273.15 + 237.1))/std::pow(T_b - 273.15 + 105.0,1.57);

    double T = T_b - 273.15;

    //double p_s = 0.61121e3*exp((18.678 - T/234.5)*(T/(257.14 + T)));

    return p_s;
}


double compute_relative_humidity(double z, double p_b, double T_b)
{
    double p_s = compute_saturation_pressure(T_b);

    double q_s = 0.622*p_s/(p_b - p_s);

    if(z <= height){
        return 0.014/q_s;
    }
    else if(z <= z_tr){
        return 1.0 - 0.75*std::pow(z/z_tr,1.25);
    }
    else{
        return 0.25;
    }
}

double compute_vapor_pressure(double p_s, double RH)
{

    return p_s*RH;
}


double vapor_mixing_ratio(const double z, const double p_b, const double T_b, const double RH){

    double p_s = compute_saturation_pressure(T_b);
    double p_v = compute_vapor_pressure(p_s, RH);

    //std::cout << "Vapor pressure is " << p_s << " " << p_v << "\n";

    double q_v = 0.622*p_v/(p_b - p_v);

        if(z < height){
            return 0.014;
        }
        else{
            return q_v;
        }
}

double compute_temperature(const double p_b, const double theta_b)
{
    return theta_b*std::pow(p_b/p0,Rd/Cp);
}

double compute_dewpoint_temperature(const double z, const double p_b, const double T_b, const double RH)
{

    double T_dp, gamma, T;
    T = T_b - 273.15;

    double a = 6.11e-3*p0, b = 18.678, c = 257.14, d = 234.5;
	std::cout << "Value is " << RH << " " << (b - T/d)*T/(c + T) << "\n";
    gamma = log(RH*exp((b - T/d)*T/(c + T)));
	
    T_dp = c*gamma/(b - gamma);

	std::cout << "T_dp value is " << T_dp << "\n";

    return T_dp;

}

void
init_isentropic_hse_no_terrain(Real *theta, Real* r, Real* p, 
                               const Real& dz, const Real&  prob_lo_z,
                               const int& khi)
{
    double z, p_b, theta_b, T_b, rho_b, RH, T_dp, rho0, theta0, q_v=0.0;
    p_b = p0;
    T_b = T0;

    FILE *file_IC, *file_parcel;
    file_IC = fopen("input_sounding.txt","w");

	// Do the first integration from z = 0 to z = 0.5*dz
	z = 0;
	theta0 = compute_theta(z);
	T_b    = compute_temperature(p0, theta0);
	RH     = compute_relative_humidity(z, p0, T0);
    q_v    = 0.0;//vapor_mixing_ratio(z, p[k-1], T_b, RH);
	rho0   = p0/(Rd*T_b*(1.0 + Rv_by_Rd*q_v));

	p[0]   = p0 - rho0*g*(1.0 + q_v)*delz/2.0;

	// Compute the quantities at z = 0.5*dz

	z = prob_lo_z + 0.5*dz;
	theta[0] = compute_theta(z);
	T_b      = compute_temperature(p[0], theta[0]);
	RH       = compute_relative_humidity(z, p[0], T_b);
	q_v      = 0.0;//vapor_mixing_ratio(z, p[0], T_b, RH);
    T_dp     = compute_dewpoint_temperature(z, p[0], T_b, RH);
	r[0]     = p[0]/(Rd*T_b*(1.0 + Rv_by_Rd*q_v));
    fprintf(file_IC, "%0.15g %0.15g %0.15g %0.15g %0.15g\n " , z, T_b-273.15, T_dp, p[0], r[0]);
	

    for (int k=1;k<=khi;k++){
		std::cout << "k val is " << k << "\n";
        z = prob_lo_z + (k+0.5)*dz;

        p[k] = p[k-1] - r[k-1]*g*(1.0 + q_v)*delz;

        theta[k] = compute_theta(z);
		T_b      = compute_temperature(p[k], theta[k]);
		RH       = compute_relative_humidity(z, p[k], T_b);
		q_v      = 0.0;//vapor_mixing_ratio(z, p[k-1], T_b, RH);
        r[k]     = p[k]/(Rd*T_b*(1.0 + Rv_by_Rd*q_v));


    	fprintf(file_IC, "%0.15g %0.15g %0.15g %0.15g %0.15g\n " , z, T_b-273.15, T_dp, p[k], r[k]);
        std::cout << "Temperature is " << z << " " << T_b  - 273.15 << " " << T_dp << " " << rho_b  << " " << p[k] << "\n";

    }
    fclose(file_IC);

  r[khi+1] = r[khi];
}

void
erf_init_dens_hse(MultiFab& rho_hse,
                  std::unique_ptr<MultiFab>& z_phys_nd,
                  std::unique_ptr<MultiFab>& z_phys_cc,
                  Geometry const& geom)
{
    const Real prob_lo_z = geom.ProbLo()[2];
    const Real dz        = geom.CellSize()[2];
    const int khi        = geom.Domain().bigEnd()[2];

	std::cout << "Value of khi is " << khi  << "\n";
	exit(0);

    const Real T_sfc    = 300.;
    const Real rho_sfc  = p_0 / (R_d*T_sfc);
    const Real Thetabar = T_sfc;

    // use_terrain = 1
    if (z_phys_nd) {

        if (khi > 255) amrex::Abort("1D Arrays are hard-wired to only 256 high");

        for ( MFIter mfi(rho_hse, TileNoZ()); mfi.isValid(); ++mfi )
        {
            Array4<Real      > rho_arr  = rho_hse.array(mfi);
            Array4<Real const> z_cc_arr = z_phys_cc->const_array(mfi);

            // Create a flat box with same horizontal extent but only one cell in vertical
            const Box& tbz = mfi.nodaltilebox(2);
            Box b2d = tbz; // Copy constructor
            b2d.grow(0,1); b2d.grow(1,1); // Grow by one in the lateral directions
            b2d.setRange(2,0);

            ParallelFor(b2d, [=] AMREX_GPU_DEVICE (int i, int j, int) {
              Array1D<Real,0,255> r;;
              Array1D<Real,0,255> p;;

              //init_isentropic_hse_terrain(i,j,rho_sfc,Thetabar,&(r(0)),&(p(0)),z_cc_arr,khi);

              for (int k = 0; k <= khi; k++) {
                 rho_arr(i,j,k) = r(k);
              }
              rho_arr(i,j,   -1) = rho_arr(i,j,0);
              rho_arr(i,j,khi+1) = rho_arr(i,j,khi);
            });
        } // mfi
    } else { // use_terrain = 0

        // These are at cell centers (unstaggered)
        Vector<Real> h_r(khi+2);
        Vector<Real> h_p(khi+2);
        Vector<Real> h_t(khi+2);

        amrex::Gpu::DeviceVector<Real> d_r(khi+2);
        amrex::Gpu::DeviceVector<Real> d_p(khi+2);
        amrex::Gpu::DeviceVector<Real> d_t(khi+2);

        init_isentropic_hse_no_terrain(h_t.data(), h_r.data(),h_p.data(),dz,prob_lo_z,khi);

        amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, h_r.begin(), h_r.end(), d_r.begin());
        amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, h_p.begin(), h_p.end(), d_p.begin());
        amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, h_t.begin(), h_t.end(), d_t.begin());

        Real* r = d_r.data();

#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
          for ( MFIter mfi(rho_hse,TilingIfNotGPU()); mfi.isValid(); ++mfi)
          {
              const Box& bx = mfi.growntilebox(1);
              const Array4<Real> rho_hse_arr = rho_hse[mfi].array();
              ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
              {
                  int kk = std::max(k,0);
                  rho_hse_arr(i,j,k) = r[kk];
              });
          } // mfi
    } // no terrain
}



void
init_custom_prob(
    const Box& bx,
    const Box& xbx,
    const Box& ybx,
    const Box& zbx,
    Array4<Real      > const& state,
    Array4<Real      > const& x_vel,
    Array4<Real      > const& y_vel,
    Array4<Real      > const& z_vel,
    Array4<Real      > const& /*r_hse*/,
    Array4<Real      > const& /*p_hse*/,
    Array4<Real const> const& /*z_nd*/,
    Array4<Real const> const& /*z_cc*/,
#if defined(ERF_USE_MOISTURE)
    Array4<Real      > const& qv,
    Array4<Real      > const& qc,
    Array4<Real      > const& qi,
#elif defined(ERF_USE_WARM_NO_PRECIP)
    Array4<Real      > const&   ,
    Array4<Real      > const&   ,
#endif
    GeometryData const& geomdata,
    Array4<Real const> const& /*mf_m*/,
    Array4<Real const> const& /*mf_u*/,
    Array4<Real const> const& /*mf_v*/,
    const SolverChoice& sc)
{
  const int khi = geomdata.Domain().bigEnd()[2];

  AMREX_ALWAYS_ASSERT(bx.length()[2] == khi+1);

  // This is what we do at k = 0 -- note we assume p = p_0 and T = T_0 at z=0
  const amrex::Real& dz        = geomdata.CellSize()[2];
  const amrex::Real& prob_lo_z = geomdata.ProbLo()[2];
  const amrex::Real& prob_hi_z = geomdata.ProbHi()[2];

  const amrex::Real rdOcp   = sc.rdOcp;

  // const amrex::Real thetatr = 343.0;
  // const amrex::Real theta0  = 300.0;
  const amrex::Real ztr     = 12000.0;
  const amrex::Real Ttr     = 213.0;
  const amrex::Real Ttop    = 213.0;
  const amrex::Real deltaz  = 1000.*0.0;
  const amrex::Real zs      = 5000.;
  const amrex::Real us      = 30.;
  const amrex::Real uc      = 15.;

  // Call the routine to calculate the 1d background condition

   Vector<Real> h_r(khi+2);
   Vector<Real> h_p(khi+2);
   Vector<Real> h_t(khi+2);

   amrex::Gpu::DeviceVector<Real> d_r(khi+2);
   amrex::Gpu::DeviceVector<Real> d_p(khi+2);
   amrex::Gpu::DeviceVector<Real> d_t(khi+2);

	std::cout << "Reached here 0" << "\n";
   init_isentropic_hse_no_terrain(h_t.data(), h_r.data(),h_p.data(),dz,prob_lo_z,khi);

	std::cout << "Reached here 1" << "\n";
   amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, h_r.begin(), h_r.end(), d_r.begin());
   amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, h_p.begin(), h_p.end(), d_p.begin());
   amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, h_t.begin(), h_t.end(), d_t.begin());

	std::cout << "Reached here 2" << "\n";

    Real* t = d_t.data();
    Real* r = d_r.data();

  amrex::ParallelForRNG(bx, [=, parms=parms] AMREX_GPU_DEVICE(int i, int j, int k, const amrex::RandomEngine& engine) noexcept
  {
    // Geometry (note we must include these here to get the data on device)
    const auto prob_lo         = geomdata.ProbLo();
    const auto dx              = geomdata.CellSize();
    const amrex::Real z        = prob_lo[2] + (k + 0.5) * dx[2];


    // This version perturbs rho but not p
    state(i, j, k, RhoTheta_comp) = d_r[k]*d_t[k];
    state(i, j, k, Rho_comp)      = d_r[k];

    // Set scalar = 0 everywhere
    state(i, j, k, RhoScalar_comp) = 0.0;

    // mean states
#if defined(ERF_USE_MOISTURE)
    state(i, j, k, RhoQt_comp) = 0.0;//rho*qvapor;
    state(i, j, k, RhoQp_comp) = 0.0;
    qv(i, j, k) = 0.0;//qvapor;
    qc(i, j, k) = 0.0;
    qi(i, j, k) = 0.0;
#elif defined(ERF_USE_WARM_NO_PRECIP)
    state(i, j, k, RhoQv_comp) = 0.0;//rho*qvapor;
    state(i, j, k, RhoQc_comp) = 0.0;
#endif
  });

  // Set the x-velocity
  amrex::ParallelFor(xbx, [=, parms=parms] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    const amrex::Real z = prob_lo_z + (k+0.5) * dz;
    if (z < zs-deltaz) {
      x_vel(i, j, k) = us*(z/zs) - uc;
    } else if (std::abs(z-zs) < deltaz) {
      x_vel(i, j, k) = (-0.8+3.*(z/zs)-1.25*(z/zs)*(z/zs))*us-uc;
    } else {
      x_vel(i, j, k) = us-uc;
    }
	x_vel(i,j,k) = 0.0;
  });

  // Set the y-velocity
  amrex::ParallelFor(ybx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    y_vel(i, j, k) = 0.0;
  });

  // Set the z-velocity
  amrex::ParallelFor(zbx, [=, parms=parms] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    z_vel(i, j, k) = 0.0;
  });

  amrex::Gpu::streamSynchronize();
}

void
init_custom_terrain (const Geometry& /*geom*/,
                           MultiFab& z_phys_nd,
                     const Real& /*time*/)
{
    // Number of ghost cells
    int ngrow = z_phys_nd.nGrow();

    // Bottom of domain
    int k0 = 0;

    for ( MFIter mfi(z_phys_nd, TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
        // Grown box with no z range
        amrex::Box xybx = mfi.growntilebox(ngrow);
        xybx.setRange(2,0);

        Array4<Real> const& z_arr = z_phys_nd.array(mfi);

        ParallelFor(xybx, [=] AMREX_GPU_DEVICE (int i, int j, int) {

            // Flat terrain with z = 0 at k = 0
            z_arr(i,j,k0) = 0.0;
        });
    }
}

void
erf_init_rayleigh(amrex::Vector<amrex::Real>& tau,
                  amrex::Vector<amrex::Real>& ubar,
                  amrex::Vector<amrex::Real>& vbar,
                  amrex::Vector<amrex::Real>& wbar,
                  amrex::Vector<amrex::Real>& thetabar,
                  amrex::Geometry      const& geom)
{
  const int khi = geom.Domain().bigEnd()[2];

  // We just use these values to test the Rayleigh damping
  for (int k = 0; k <= khi; k++)
  {
      tau[k]  = 1.0;
      ubar[k] = 2.0;
      vbar[k] = 1.0;
      wbar[k] = 0.0;
      thetabar[k] = parms.Theta_0;
  }
}

void
amrex_probinit(
  const amrex_real* /*problo*/,
  const amrex_real* /*probhi*/)
{
  // Parse params
  amrex::ParmParse pp("prob");
  pp.query("T_0", parms.T_0);
  pp.query("U_0", parms.U_0);
  pp.query("x_c", parms.x_c);
  pp.query("z_c", parms.z_c);
  pp.query("x_r", parms.x_r);
  pp.query("z_r", parms.z_r);
  pp.query("T_pert", parms.T_pert);
}
