#include <EOS.H>
#include <ERF.H>

// get plotfile name
std::string
ERF::PlotFileName (int lev) const
{
    return amrex::Concatenate(plot_file, lev, 5);
}

void
ERF::setPlotVariables ()
{
    ParmParse pp("amr");

    if (pp.contains("plot_vars"))
    {
        std::string nm;

        int nPltVars = pp.countval("plot_vars");

        for (int i = 0; i < nPltVars; i++)
        {
            pp.get("plot_vars", nm, i);

            if (nm == "ALL") {
                // put all conserved state variables in the plot variables list
                plot_state_names = cons_names;
                // put all velocity components in the plot variables list
                plot_state_names.insert(plot_state_names.end(), velocity_names.begin(), velocity_names.end());
            } else if (nm == "NONE" || nm == "None") {
                // put no state variables in the plot variable list
                plot_state_names.clear();
            } else {
                // add the named variable to our list of plot variables
                // if it is not already in the list
                if (!containerHasElement(plot_state_names, nm)) {
                    plot_state_names.push_back(nm);
                }
            }
        }
    }
    else
    {
        //
        // The default is to add all state conserved and velocity variables to the plot variables list
        //
        plot_state_names = cons_names;
        plot_state_names.insert(plot_state_names.end(), velocity_names.begin(), velocity_names.end());
    }

    // Get state variables in the same order as we define them,
    // since they may be in any order in the input list
    Vector<std::string> tmp_plot_names;
    for (int i = 0; i < Cons::NumVars; ++i) {
        if (containerHasElement(plot_state_names, cons_names[i])) {
            tmp_plot_names.push_back(cons_names[i]);
        }
    }
    // check for velocity since it's not in cons_names
    // if we are asked for any velocity component, we will need them all
    if (containerHasElement(plot_state_names, "x_velocity") ||
        containerHasElement(plot_state_names, "y_velocity") ||
        containerHasElement(plot_state_names, "z_velocity")) {
        tmp_plot_names.push_back("x_velocity");
        tmp_plot_names.push_back("y_velocity");
        tmp_plot_names.push_back("z_velocity");
    }
    plot_state_names = tmp_plot_names;

    if (pp.contains("derive_plot_vars"))
    {
        std::string nm;

        int nDrvPltVars = pp.countval("derive_plot_vars");

        for (int i = 0; i < nDrvPltVars; i++)
        {
            pp.get("derive_plot_vars", nm, i);

            if (nm == "ALL") {
                // put all diagnostic variables in the plot variables list
                plot_deriv_names = derived_names;
            } else if (nm == "NONE" || nm == "None") {
                // put no diagnostic variables in the plot variable list
                plot_deriv_names.clear();
            } else {
                // add the named variable to our list of plot variables
                // if it is not already in the list
                if (!containerHasElement(plot_deriv_names, nm)) {
                    plot_deriv_names.push_back(nm);
                }
            }
        }
    }
    else
    {
        //
        // The default is to add none of the diagnostic variables to the plot variables list
        //
        plot_deriv_names.clear();
    }

    // Get derived variables in the same order as we define them,
    // since they may be in any order in the input list
    Vector<std::string> tmp_deriv_names;
    for (int i = 0; i < derived_names.size(); ++i) {
        if (containerHasElement(plot_deriv_names, derived_names[i])) {
            tmp_deriv_names.push_back(derived_names[i]);
        }
    }
    plot_deriv_names = tmp_deriv_names;
}

// set plotfile variable names
Vector<std::string>
ERF::PlotFileVarNames () const
{
    Vector<std::string> names;

    names.insert(names.end(), plot_state_names.begin(), plot_state_names.end());
    names.insert(names.end(), plot_deriv_names.begin(), plot_deriv_names.end());

    return names;

}

// write plotfile to disk
void
ERF::WritePlotFile () const
{
    const Vector<std::string> varnames = PlotFileVarNames();
    const int ncomp_mf = varnames.size();

    if (ncomp_mf == 0)
        return;

    Vector<MultiFab> mf(finest_level+1);
    for (int lev = 0; lev <= finest_level; ++lev) {
        mf[lev].define(grids[lev], dmap[lev], ncomp_mf, 0);
    }

    for (int lev = 0; lev <= finest_level; ++lev) {
        int mf_comp = 0;

        // First, copy any of the conserved state variables into the output plotfile
        for (int i = 0; i < Cons::NumVars; ++i) {
            if (containerHasElement(plot_state_names, cons_names[i])) {
                MultiFab::Copy(mf[lev],vars_new[lev][Vars::cons],i,mf_comp,1,0);
                mf_comp++;
            }
        }

        // Next, check for velocities and if desired, output them
        if (containerHasElement(plot_state_names, "x_velocity") ||
            containerHasElement(plot_state_names, "y_velocity") ||
            containerHasElement(plot_state_names, "z_velocity")) {

            average_face_to_cellcenter(mf[lev],mf_comp,
                Array<const MultiFab*,3>{&vars_new[lev][Vars::xvel],&vars_new[lev][Vars::yvel],&vars_new[lev][Vars::zvel]});
            mf_comp += AMREX_SPACEDIM;
        }

        // Finally, check for any derived quantities and compute them, inserting
        // them into our output multifab
        auto calculate_derived = [&](const std::string der_name,
                                     const decltype(derived::erf_dernull)& der_function)
        {
            if (containerHasElement(plot_deriv_names, der_name)) {
                MultiFab dmf(mf[lev], amrex::make_alias, mf_comp, 1);

                for (MFIter mfi(dmf, TilingIfNotGPU()); mfi.isValid(); ++mfi)
                {
                    const Box& bx = mfi.tilebox();
                    auto& dfab = dmf[mfi];
                    auto& sfab = vars_new[lev][Vars::cons][mfi];
                    der_function(bx, dfab, 0, 1, sfab, Geom(lev), t_new[0], nullptr, lev);
                }

                mf_comp++;
            }
        };

        calculate_derived("pressure",    derived::erf_derpres);
        calculate_derived("soundspeed",  derived::erf_dersoundspeed);
        calculate_derived("temp",        derived::erf_dertemp);
        calculate_derived("theta",       derived::erf_dertheta);
        calculate_derived("scalar",      derived::erf_derscalar);

        if (containerHasElement(plot_deriv_names, "pres_hse"))
        {
            auto d_pres_hse_lev = d_pres_hse[lev].dataPtr();
            for ( amrex::MFIter mfi(mf[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                const Box& bx = mfi.tilebox();
                const Array4<Real>& derdat = mf[lev].array(mfi);
                amrex::ParallelFor(bx, [=, ng_pres_hse=ng_pres_hse] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    derdat(i, j, k, mf_comp) = d_pres_hse_lev[k+ng_pres_hse];
                });
            }
            mf_comp += 1;
        }

        if (containerHasElement(plot_deriv_names, "dens_hse"))
        {
            auto d_dens_hse_lev = d_dens_hse[lev].dataPtr();
            for ( amrex::MFIter mfi(mf[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                const Box& bx = mfi.tilebox();
                const Array4<Real>& derdat = mf[lev].array(mfi);
                amrex::ParallelFor(bx, [=, ng_dens_hse=ng_dens_hse] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    derdat(i, j, k, mf_comp) = d_dens_hse_lev[k+ng_dens_hse];
                });
            }
            mf_comp ++;
        }

        if (containerHasElement(plot_deriv_names, "pert_pres"))
        {
            auto d_pres_hse_lev = d_pres_hse[lev].dataPtr();
            for ( amrex::MFIter mfi(mf[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                const Box& bx = mfi.tilebox();
                const Array4<Real>& derdat = mf[lev].array(mfi);
                const Array4<Real const>& S_arr = vars_new[lev][Vars::cons].const_array(mfi);
                amrex::ParallelFor(bx, [=, ng_pres_hse=ng_pres_hse] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    const Real rhotheta = S_arr(i,j,k,RhoTheta_comp);
                    derdat(i, j, k, mf_comp) = getPgivenRTh(rhotheta) - d_pres_hse_lev[k+ng_pres_hse];
                });
            }
            mf_comp ++;
        }

        if (containerHasElement(plot_deriv_names, "pert_dens"))
        {
            auto d_dens_hse_lev = d_dens_hse[lev].dataPtr();
            for ( amrex::MFIter mfi(mf[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                const Box& bx = mfi.tilebox();
                const Array4<Real>& derdat  = mf[lev].array(mfi);
                const Array4<Real const>& S_arr = vars_new[lev][Vars::cons].const_array(mfi);
                amrex::ParallelFor(bx, [=, ng_dens_hse=ng_dens_hse] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    derdat(i, j, k, mf_comp) = S_arr(i,j,k,Rho_comp) - d_dens_hse_lev[k+ng_dens_hse];
                });
            }
            mf_comp ++;
        }
    }

    const std::string& plotfilename = PlotFileName(istep[0]);
    amrex::Print() << "Writing plotfile " << plotfilename << "\n";

    if (plotfile_type == "amrex") {
        amrex::WriteMultiLevelPlotfile(plotfilename, finest_level+1, GetVecOfConstPtrs(mf), varnames,
                                       Geom(), t_new[0], istep, refRatio());
        writeJobInfo(plotfilename);
#ifdef ERF_USE_NETCDF
    } else {
         writeNCPlotFile(plotfilename, GetVecOfConstPtrs(mf), varnames, istep, t_new[0]);
#endif
    }
}