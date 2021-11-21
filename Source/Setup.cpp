#include <cstdio>

#include <AMReX_LevelBld.H>
#include <AMReX_ParmParse.H>
#include <AMReX_buildInfo.H>

#include "ERF.H"
#include "Derive.H"
#include "IndexDefines.H"
#include "prob.H"

static amrex::Box
the_same_box(const amrex::Box& b)
{
  return b;
}
static amrex::Box
grow_box_by_one(const amrex::Box& b)
{
  return amrex::grow(b, 1);
}

void
ERF::variableSetUp()
{
  // ERF::variableSetUp is called in the constructor of Amr.cpp, so
  // it should get called every time we start or restart a job

  // initialize the start time for our CPU-time tracker
  startCPUTime = amrex::ParallelDescriptor::second();

  // Output the git commit hashes used to build the executable.

  if (amrex::ParallelDescriptor::IOProcessor()) {
    const char* erf_hash = amrex::buildInfoGetGitHash(1);
    const char* amrex_hash = amrex::buildInfoGetGitHash(2);
    const char* buildgithash = amrex::buildInfoGetBuildGitHash();
    const char* buildgitname = amrex::buildInfoGetBuildGitName();

    if (strlen(erf_hash) > 0) {
      amrex::Print() << "\n"
                     << "ERF git hash: " << erf_hash << "\n";
    }
    if (strlen(amrex_hash) > 0) {
      amrex::Print() << "AMReX git hash: " << amrex_hash << "\n";
    }
    if (strlen(buildgithash) > 0) {
      amrex::Print() << buildgitname << " git hash: " << buildgithash << "\n";
    }

    amrex::Print() << "\n";
  }

  AMREX_ASSERT(desc_lst.size() == 0);

  // Get options, set phys_bc
  read_params();

  //
  // Set number of state variables and pointers to components
  //

#ifdef NUM_ADV
  NumAdv = NUM_ADV;
#else
  NumAdv = 0;
#endif

  if (NumAdv > 0) {
    FirstAdv = RhoTheta_comp+1;
  }

  amrex::Interpolater* interp;

  // At this cell_cons_interp is the default, but we might at some point use pc_interp for debugging
  int state_interp_order = 1;
  if (state_interp_order == 0) {
    interp = &amrex::pc_interp;
  } else {
    interp = &amrex::cell_cons_interp;
  }

  // Note that the default is state_data_extrap = false,
  // store_in_checkpoint = true.  We only need to put these in
  // explicitly if we want to do something different,
  // like not store the state data in a checkpoint directory
  bool state_data_extrap = false;
  bool store_in_checkpoint;

  int ngrow_state = 1;

  // NVAR is currently set to 2 in IndexDefines.H
  store_in_checkpoint = true;
  desc_lst.addDescriptor(
    State_Type, amrex::IndexType::TheCellType(), amrex::StateDescriptor::Point,
    ngrow_state, NVAR, interp, state_data_extrap, store_in_checkpoint);

  amrex::Vector<amrex::BCRec> bcs(NVAR);
  amrex::Vector<std::string> name(NVAR);

  amrex::BCRec bc;
  int cnt = 0;
  bcs[cnt] = bc;
  name[cnt] = "density";
  cnt++;
  bcs[cnt] = bc;
  name[cnt] = "rhotheta";

  for (int i = 0; i < NumAdv; ++i) {
    char buf[64];
    sprintf(buf, "adv_%d", i);
    cnt++;
    bcs[cnt] = bc;
    name[cnt] = std::string(buf);
  }

  amrex::StateDescriptor::BndryFunc bndryfunc1(erf_bcfill_hyp);
  bndryfunc1.setRunOnGPU(true);

  desc_lst.setComponent(State_Type, Rho_comp, name, bcs, bndryfunc1);

  //
  // Create face-based StateData for velocities on each face,
  //     and for all conserved quantities as well
  //     (unnecessary right now but convenient)
  //
  amrex::StateDescriptor::BndryFunc bndryfunc_null(erf_nullfill);
  bndryfunc_null.setRunOnGPU(true);
  amrex::Vector<amrex::BCRec> vel_bcs(1);
  amrex::Vector<std::string> vel_name(1);
  vel_bcs[0] = bc;

  store_in_checkpoint = true;
  amrex::IndexType xface(amrex::IntVect(1,0,0));
  desc_lst.addDescriptor(X_Vel_Type, xface,
                         amrex::StateDescriptor::Point, 1, 1,
                         interp, state_data_extrap,
                         store_in_checkpoint);

  vel_name[0] = "x_velocity";
  desc_lst.setComponent(X_Vel_Type, 0, vel_name, vel_bcs, bndryfunc_null);

  amrex::IndexType yface(amrex::IntVect(0,1,0));
  desc_lst.addDescriptor(Y_Vel_Type, yface,
                         amrex::StateDescriptor::Point, 1, 1,
                         interp, state_data_extrap,
                         store_in_checkpoint);

  vel_name[0] = "y_velocity";
  desc_lst.setComponent(Y_Vel_Type, 0, vel_name, vel_bcs, bndryfunc_null);

  amrex::IndexType zface(amrex::IntVect(0,0,1));
  desc_lst.addDescriptor(Z_Vel_Type, zface,
                         amrex::StateDescriptor::Point, 1, 1,
                         interp, state_data_extrap,
                         store_in_checkpoint);

  vel_name[0] = "z_velocity";
  desc_lst.setComponent(Z_Vel_Type, 0, vel_name, vel_bcs, bndryfunc_null);

  num_state_type = desc_lst.size();

  //
  // DEFINE DERIVED QUANTITIES
  //
  // Pressure
  //
  derive_lst.add(
    "pressure", amrex::IndexType::TheCellType(), 1, erf_derpres, the_same_box);
  derive_lst.addComponent("pressure", desc_lst, State_Type, Rho_comp, NVAR);

  //
  // Temperature
  //
  derive_lst.add(
    "temp", amrex::IndexType::TheCellType(), 1, erf_dertemp, the_same_box);
  derive_lst.addComponent("temp", desc_lst, State_Type, Rho_comp, NVAR);

  //
  // Potential Temperature
  //
  derive_lst.add(
    "theta", amrex::IndexType::TheCellType(), 1, erf_dertheta, the_same_box);
  derive_lst.addComponent("theta", desc_lst, State_Type, Rho_comp, NVAR);

  //
  // Velocities
  //

  // This calculcates cell-centered x-velocity from x-face-centered values
  // Note that we don't need to use "addComponent" because we are are computing
  //      these in ERF itself, not calling the AMReX derive routines
  derive_lst.add(
    "x_velocity", amrex::IndexType::TheCellType(), 1, erf_dernull, the_same_box);

  // This calculcates cell-centered y-velocity from y-face-centered values
  derive_lst.add(
    "y_velocity", amrex::IndexType::TheCellType(), 1, erf_dernull, the_same_box);

  // This calculcates cell-centered z-velocity from z-face-centered values
  derive_lst.add(
    "z_velocity", amrex::IndexType::TheCellType(), 1, erf_dernull, the_same_box);

  // Problem-specific derives
  add_problem_derives<ProblemDerives>(derive_lst, desc_lst);
}

void
ERF::initDataProb (amrex::MultiFab& S_new,
                   amrex::MultiFab& U_new,
                   amrex::MultiFab& V_new,
                   amrex::MultiFab& W_new)
{

  // Initialize to zero (though we sholdn't actually need to do this)
  S_new.setVal(0.0);
  U_new.setVal(0.0);
  V_new.setVal(0.0);
  W_new.setVal(0.0);

  if (verbose) {
    amrex::Print() << "Initializing the data at level " << level << std::endl;
  }

#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif

  const auto geomdata = geom.data();

  for (amrex::MFIter mfi(S_new, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
  {
    const amrex::Box& bx = mfi.tilebox();
    auto sfab  = S_new.array(mfi);
    auto ufab  = U_new.array(mfi);
    auto vfab  = V_new.array(mfi);
    auto wfab  = W_new.array(mfi);
    erf_init_prob(bx, sfab, ufab, vfab, wfab, geomdata);
  }

  // Ensure that the face-based data are the same on both sides of a periodic domain.
  // The data associated with the lower grid ID is considered the correct value.
  U_new.OverrideSync(geom.periodicity());
  V_new.OverrideSync(geom.periodicity());
  W_new.OverrideSync(geom.periodicity());

  if (verbose) {
    amrex::Print() << "Done initializing level " << level << " data "
                   << std::endl;
  }

}

void
ERF::initHSE()
{
    erf_init_dens_hse(dens_hse[level]);
    erf_enforce_hse(dens_hse[level],pres_hse[level]);
}

void
ERF::erf_enforce_hse(amrex::Vector<amrex::Real>& dens,
                     amrex::Vector<amrex::Real>& pres)
{
    AMREX_ALWAYS_ASSERT(dens.size() == pres.size());

    const auto geomdata = geom.data();
    const Real dz = geomdata.CellSize(2);
    int nz = dens.size();
    
    // We start by assuming pressure on the ground is p_0 (in ERF_Constants.H)
    // Note that gravity is positive

    int l_spatial_order = solverChoice.spatial_order;
    int l_gravity       = solverChoice.gravity;

    Real dens_interp;

    // We integrate to the first cell by using rho in this cell
    pres[0] = p_0 - (0.5*dz) * dens[0] * l_gravity;

    for (int k = 1; k < nz; k++)
    {
        switch (l_spatial_order) {
            case 2:
               dens_interp = 0.5*(dens[k] + dens[k-1]);
               break;
            case 4:
               if (k == 1)
                   dens_interp =  0.5*(dens[k] + dens[k-1]);
               else
                   dens_interp = (7./12.)*(dens[k] + dens[k-1]) - 1./12.*(dens[k+1]+dens[k-2]);
               break;
            case 6:
               if (k == 1)
                   dens_interp =  0.5*(dens[k] + dens[k-1]);
               else if (k == 2)
                   dens_interp = (7./12.)*(dens[k] + dens[k-1]) - 1./12.*(dens[k+1]+dens[k-2]);
               else 
                   dens_interp = (37./60.)*(dens[k  ]+dens[k-1]) 
                               - ( 8./60.)*(dens[k+1]+dens[k-2])
                                +( 1./60.)*(dens[k+2]+dens[k-3]) ;
               break;
        }
        pres[k] = pres[k-1] - (0.5*dz) * dens_interp * l_gravity;
    }
}
