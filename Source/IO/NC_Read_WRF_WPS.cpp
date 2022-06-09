#include "ERF.H"
#include "AMReX_FArrayBox.H"

using namespace amrex;

#ifdef ERF_USE_NETCDF
void
ERF::read_from_wrfinput(int lev)
{
    Box input_box;

    NC_xvel_fab[lev].resize(num_boxes_at_level[lev]);
    NC_yvel_fab[lev].resize(num_boxes_at_level[lev]);
    NC_zvel_fab[lev].resize(num_boxes_at_level[lev]);
    NC_rho_fab[lev].resize(num_boxes_at_level[lev]);
    NC_rhotheta_fab[lev].resize(num_boxes_at_level[lev]);

#ifdef ERF_USE_TERRAIN
    NC_PH_fab[lev].resize(num_boxes_at_level[lev]);
    NC_PHB_fab[lev].resize(num_boxes_at_level[lev]);
#endif

    for (int idx = 0; idx < num_boxes_at_level[lev]; idx++)
    {
        if (lev == 0) {
            input_box = geom[0].Domain();
        } else {
            input_box = boxes_at_level[lev][idx];
        }

        // We allocate these here so they exist on all ranks
        Box ubx(input_box); ubx.surroundingNodes(0);
        Box vbx(input_box); vbx.surroundingNodes(1);
        Box wbx(input_box); wbx.surroundingNodes(2);

        NC_xvel_fab[lev][idx].resize(ubx,1);
        NC_yvel_fab[lev][idx].resize(vbx,1);
        NC_zvel_fab[lev][idx].resize(wbx,1);
        NC_rho_fab[lev][idx].resize(input_box,1);
        NC_rhotheta_fab[lev][idx].resize(input_box,1);

#ifdef ERF_USE_TERRAIN
        NC_PH_fab[lev][idx].resize(wbx,1);
        NC_PHB_fab[lev][idx].resize(wbx,1);
#endif

#ifdef AMREX_USE_GPU
        FArrayBox host_NC_xvel_fab   (NC_xvel_fab[lev][idx].box(),     NC_xvel_fab[lev][idx].nComp(),amrex::The_Pinned_Arena());
        FArrayBox host_NC_yvel_fab   (NC_yvel_fab[lev][idx].box(),     NC_yvel_fab[lev][idx].nComp(),amrex::The_Pinned_Arena());
        FArrayBox host_NC_zvel_fab   (NC_zvel_fab[lev][idx].box(),     NC_zvel_fab[lev][idx].nComp(),amrex::The_Pinned_Arena());
        FArrayBox host_NC_rho_fab    (NC_rho_fab[lev][idx].box(),      NC_rho_fab[lev][idx].nComp(),amrex::The_Pinned_Arena());
       FArrayBox host_NC_rhotheta_fab(NC_rhotheta_fab[lev][idx].box(), NC_rhotheta_fab[lev][idx].nComp(),amrex::The_Pinned_Arena());
#ifdef ERF_USE_TERRAIN
        FArrayBox host_NC_PH_fab (NC_PH_fab[lev][idx].box(),  NC_PH_fab[lev][idx].nComp(),amrex::The_Pinned_Arena());
        FArrayBox host_NC_PHB_fab(NC_PHB_fab[lev][idx].box(), NC_PHB_fab[lev][idx].nComp(),amrex::The_Pinned_Arena());
#endif
#else
        FArrayBox host_NC_xvel_fab    (NC_xvel_fab[lev][idx]    , amrex::make_alias, 0, NC_xvel_fab[lev][idx].nComp());
        FArrayBox host_NC_yvel_fab    (NC_yvel_fab[lev][idx]    , amrex::make_alias, 0, NC_yvel_fab[lev][idx].nComp());
        FArrayBox host_NC_zvel_fab    (NC_zvel_fab[lev][idx]    , amrex::make_alias, 0, NC_zvel_fab[lev][idx].nComp());
        FArrayBox host_NC_rho_fab     (NC_rho_fab[lev][idx]     , amrex::make_alias, 0, NC_rho_fab[lev][idx].nComp());
        FArrayBox host_NC_rhotheta_fab(NC_rhotheta_fab[lev][idx], amrex::make_alias, 0, NC_rhotheta_fab[lev][idx].nComp());
#ifdef ERF_USE_TERRAIN
        FArrayBox host_NC_PH_fab      (NC_PH_fab[lev][idx]      , amrex::make_alias, 0, NC_PH_fab[lev][idx].nComp());
        FArrayBox host_NC_PHB_fab     (NC_PHB_fab[lev][idx]     , amrex::make_alias, 0, NC_PHB_fab[lev][idx].nComp());
#endif
#endif

        if (ParallelDescriptor::IOProcessor())
        {
            Vector<FArrayBox*> NC_fabs;
            Vector<std::string> NC_names;
            Vector<enum NC_Data_Dims_Type> NC_dim_types;

            NC_fabs.push_back(&host_NC_xvel_fab);     NC_names.push_back("U");
            NC_fabs.push_back(&host_NC_yvel_fab);     NC_names.push_back("V");
            NC_fabs.push_back(&host_NC_zvel_fab);     NC_names.push_back("W");
            NC_fabs.push_back(&host_NC_rho_fab);      NC_names.push_back("ALB");
            NC_fabs.push_back(&host_NC_rhotheta_fab); NC_names.push_back("T_INIT");

            NC_dim_types.push_back(NC_Data_Dims_Type::Time_BT_SN_WE);
            NC_dim_types.push_back(NC_Data_Dims_Type::Time_BT_SN_WE);
            NC_dim_types.push_back(NC_Data_Dims_Type::Time_BT_SN_WE);
            NC_dim_types.push_back(NC_Data_Dims_Type::Time_BT_SN_WE);
            NC_dim_types.push_back(NC_Data_Dims_Type::Time_BT_SN_WE);

#ifdef ERF_USE_TERRAIN
            NC_fabs.push_back(&host_NC_PH_fab);  NC_names.push_back("PH");
            NC_fabs.push_back(&host_NC_PHB_fab); NC_names.push_back("PHB");

            NC_dim_types.push_back(NC_Data_Dims_Type::Time_BT_SN_WE);
            NC_dim_types.push_back(NC_Data_Dims_Type::Time_BT_SN_WE);
#endif

            // Read the netcdf file and fill these FABs
            // NOTE: right now we are hard-wired to one "domain" per level -- but that can be generalized
            //       once we know how to determine the level for each input file
            BuildFABsFromWRFInputFile(nc_init_file[lev][idx], NC_names, NC_fabs, NC_dim_types);

        } // if ParalleDescriptor::IOProcessor()

        // We put a barrier here so the rest of the processors wait to do anything until they have the data
        amrex::ParallelDescriptor::Barrier();

        // When an FArrayBox is built, space is allocated on every rank.  However, we only
        //    filled the data in these FABs on the IOProcessor.  So here we broadcast
        //    the data to every rank.

        int ioproc = ParallelDescriptor::IOProcessorNumber();  // I/O rank
        ParallelDescriptor::Bcast(host_NC_xvel_fab.dataPtr(),NC_xvel_fab[lev][idx].box().numPts(),ioproc);
        ParallelDescriptor::Bcast(host_NC_yvel_fab.dataPtr(),NC_yvel_fab[lev][idx].box().numPts(),ioproc);
        ParallelDescriptor::Bcast(host_NC_zvel_fab.dataPtr(),NC_zvel_fab[lev][idx].box().numPts(),ioproc);
        ParallelDescriptor::Bcast(host_NC_rho_fab.dataPtr(),NC_rho_fab[lev][idx].box().numPts(),ioproc);
        ParallelDescriptor::Bcast(host_NC_rhotheta_fab.dataPtr(),NC_rhotheta_fab[lev][idx].box().numPts(),ioproc);
#ifdef ERF_USE_TERRAIN
        ParallelDescriptor::Bcast(host_NC_PHB_fab.dataPtr(),NC_PHB_fab[lev][idx].box().numPts(),ioproc);
        ParallelDescriptor::Bcast(host_NC_PH_fab.dataPtr() ,NC_PH_fab[lev][idx].box().numPts() ,ioproc);
#endif

#ifdef AMREX_USE_GPU
         Gpu::copy(Gpu::hostToDevice, host_NC_xvel_fab.dataPtr(), host_NC_xvel_fab.dataPtr()+host_NC_xvel_fab.size(),
                                           NC_xvel_fab[lev][idx].dataPtr());
         Gpu::copy(Gpu::hostToDevice, host_NC_yvel_fab.dataPtr(), host_NC_yvel_fab.dataPtr()+host_NC_yvel_fab.size(),
                                           NC_yvel_fab[lev][idx].dataPtr());
         Gpu::copy(Gpu::hostToDevice, host_NC_zvel_fab.dataPtr(), host_NC_zvel_fab.dataPtr()+host_NC_zvel_fab.size(),
                                           NC_zvel_fab[lev][idx].dataPtr());
         Gpu::copy(Gpu::hostToDevice, host_NC_rho_fab.dataPtr(), host_NC_rho_fab.dataPtr()+host_NC_rho_fab.size(),
                                           NC_rho_fab[lev][idx].dataPtr());
         Gpu::copy(Gpu::hostToDevice, host_NC_rhotheta_fab.dataPtr(), host_NC_rhotheta_fab.dataPtr()+host_NC_rhotheta_fab.size(),
                                           NC_rhotheta_fab[lev][idx].dataPtr());
#ifdef ERF_USE_TERRAIN
         Gpu::copy(Gpu::hostToDevice, host_NC_PH_fab.dataPtr(), host_NC_PH_fab.dataPtr()+host_NC_PH_fab.size(),
                                           NC_PH_fab[lev][idx].dataPtr());
         Gpu::copy(Gpu::hostToDevice, host_NC_PHB_fab.dataPtr(), host_NC_PHB_fab.dataPtr()+host_NC_PHB_fab.size(),
                                           NC_PHB_fab[lev][idx].dataPtr());
#endif
#endif

        // Convert to rho by inverting
        NC_rho_fab[lev][idx].template invert<RunOn::Device>(1.0);

        // The ideal.exe NetCDF file has this ref value subtracted from theta or T_INIT. Need to add in ERF.
        const Real theta_ref = 300.0;
        NC_rhotheta_fab[lev][idx].template plus<RunOn::Device>(theta_ref);

        // Now multiply by rho to get (rho theta) instead of theta
        NC_rhotheta_fab[lev][idx].template mult<RunOn::Device>(NC_rho_fab[lev][idx],0,0,1);

        amrex::Print() <<
          "Successfully loaded data from the wrfinput (output of 'ideal.exe' / 'real.exe') NetCDF file at level " << lev << std::endl;
    } // idx
}
#endif

#ifdef ERF_USE_NETCDF
void
ERF::read_from_wrfbdy()
{
    // *********************************************************
    // Allocate space for all of the boundary planes we may need
    // Here we make only one enough space for one time -- we will
    //    add space for the later time slices later
    // *********************************************************

    int ncomp = 4; // for right now just the velocities + temperature
    const amrex::Box& domain = geom[0].Domain();
    for (amrex::OrientationIter oit; oit != nullptr; ++oit) {
        auto ori = oit();
        if (ori.coordDir() < 2) {
            const auto& lo = domain.loVect();
            const auto& hi = domain.hiVect();

            amrex::IntVect plo(lo);
            amrex::IntVect phi(hi);
            const int normal = ori.coordDir();
            plo[normal] = ori.isHigh() ? hi[normal] + 1 : -1;
            phi[normal] = ori.isHigh() ? hi[normal] + 1 : -1;
            const amrex::Box pbx(plo, phi);
            /*
             TODO: Discuss
             bdy_data_xlo contains 4 different variables U_BXS, V_BXS, W_BXS, T_BXS.
             The dimensions of all these variables are not same as "pbx" or 1  X  NY    X  NZ
             U_BXS -> dim -> 1  X  NY    X  NZ
             V_BXS -> dim -> 1  X (NY+1) X  NZ
             W_BXS -> dim -> 1  X  NY    X (NZ+1)
             T_BXS -> dim -> 1  X  NY    X  NZ
             */
            /*
             Instead of bdy_data_xlo with 4 components, should we use U_xlo, V_xlo, W_xlo, T_xLo separately?
             */
            if (ori.coordDir() == 0 && !ori.isHigh())
                bdy_data_xlo.push_back(amrex::FArrayBox(pbx, ncomp));
            if (ori.coordDir() == 0 && ori.isHigh())
                bdy_data_xhi.push_back(amrex::FArrayBox(pbx, ncomp));
            if (ori.coordDir() == 1 && !ori.isHigh())
                bdy_data_ylo.push_back(amrex::FArrayBox(pbx, ncomp));
            if (ori.coordDir() == 1 && ori.isHigh())
                bdy_data_yhi.push_back(amrex::FArrayBox(pbx, ncomp));

            /*
             * Alternate data structure considering different dimensions of U, V, W, T on the planes
            */
            if (ori.coordDir() == 0) {
                Box x_plane_no_stag(pbx);
                Box x_plane_y_stag = convert(pbx, {0, 1, 0});
                Box x_plane_z_stag = convert(pbx, {0, 0, 1});

                Vector<FArrayBox> x_plane_data;
                x_plane_data.push_back(FArrayBox(x_plane_no_stag, 1)); //U
                x_plane_data.push_back(FArrayBox(x_plane_y_stag, 1)); //V
                x_plane_data.push_back(FArrayBox(x_plane_z_stag, 1)); //W
                x_plane_data.push_back(FArrayBox(x_plane_no_stag, 1)); // T

//                if(!ori.isHigh())
//                    bdy_x_lo.push_back(x_plane_data);
//                if(ori.isHigh())
//                    bdy_x_hi.push_back(x_plane_data);
            }
            if (ori.coordDir() == 1) {
                Box y_plane_no_stag(pbx);
                Box y_plane_x_stag = convert(pbx, {1, 0, 0});
                Box y_plane_z_stag = convert(pbx, {0, 0, 1});

                Vector<FArrayBox> y_plane_data;
                y_plane_data.push_back(FArrayBox(y_plane_x_stag, 1)); //U
                y_plane_data.push_back(FArrayBox(y_plane_no_stag, 1)); //V
                y_plane_data.push_back(FArrayBox(y_plane_z_stag, 1)); //W
                y_plane_data.push_back(FArrayBox(y_plane_no_stag, 1)); // T

//                if(!ori.isHigh())
//                    bdy_y_lo.push_back(y_plane_data);
//                if(ori.isHigh())
//                    bdy_y_hi.push_back(y_plane_data);
            }
        }
    }

    int ntimes;
    if (ParallelDescriptor::IOProcessor())
    {
        // Read the netcdf file and fill these FABs
        ntimes = BuildFABsFromWRFBdyFile(nc_bdy_file, bdy_data_xlo, bdy_data_xhi, bdy_data_ylo, bdy_data_yhi);

    } // if ParalleDescriptor::IOProcessor()

    // We put a barrier here so the rest of the processors wait to do anything until they have the data
    amrex::ParallelDescriptor::Barrier();

    int ioproc = ParallelDescriptor::IOProcessorNumber();  // I/O rank

    // Make sure all processors know how many times are stored
    ParallelDescriptor::Bcast(&ntimes,1,ioproc);

    // When an FArrayBox is built, space is allocated on every rank.  However, we only
    //    filled the data in these FABs on the IOProcessor.  So here we broadcast
    //    the data to every rank.
    for (int nt = 0; nt < ntimes; nt++)
    {
        ParallelDescriptor::Bcast(bdy_data_xlo[nt].dataPtr(),bdy_data_xlo[nt].box().numPts(),ioproc);
        ParallelDescriptor::Bcast(bdy_data_xhi[nt].dataPtr(),bdy_data_xhi[nt].box().numPts(),ioproc);
        ParallelDescriptor::Bcast(bdy_data_ylo[nt].dataPtr(),bdy_data_ylo[nt].box().numPts(),ioproc);
        ParallelDescriptor::Bcast(bdy_data_yhi[nt].dataPtr(),bdy_data_yhi[nt].box().numPts(),ioproc);
    }

    amrex::Print() << "Successfully loaded data from the wrfbdy (output of 'real.exe') NetCDF file" << std::endl << std::endl;
}
#endif // ERF_USE_NETCDF
