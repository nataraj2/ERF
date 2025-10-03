#ifndef ERF_SURFACEDATAINTERPOLATION_H_
#define ERF_SURFACEDATAINTERPOLATION_H_
/**
 * Trilinear interpolation of weather surface data onto the simulation mesh
 * The coarse weather surface data is interpolated in time first to get the surface
 * at the current time, and then spatially interpolated onto the simulation mesh
 */

#include <filesystem>
#include <stdexcept>
#include "ERF.H"
#include "ERF_ReadCustomBinaryIC.H"
#include "ERF_Interpolation_Bilinear.H"

using namespace amrex;
namespace fs = std::filesystem;

void
ERF::InterpSurfaceDataOntoMesh(const int lev,
                               const std::string& filename,
                               Vector<Vector<MultiFab>>& surface_state)
{

    // Open the binary file in input mode
    std::ifstream infile(filename, std::ios::binary);
    if (!infile) {
        std::cerr << "Error: Could not open file " << filename << std::endl;
    }
    Vector<Real> xvec_h, yvec_h, zvec_h;
    Vector<Real> sst_h, q_star_h, t_star_h, u_star_h, ls_mask_h;

    int nx, ny, nz, ndata;
    float value;

    // Read the four integers
    infile.read(reinterpret_cast<char*>(&nx), sizeof(int));
    infile.read(reinterpret_cast<char*>(&ny), sizeof(int));
    infile.read(reinterpret_cast<char*>(&nz), sizeof(int));
    infile.read(reinterpret_cast<char*>(&ndata), sizeof(int));

    amrex::Gpu::DeviceVector<Real> xvec_d(nx*ny*nz), yvec_d(nx*ny*nz), zvec_d(nx*ny*nz);
    for(int i=0; i<nx; i++) {
        infile.read(reinterpret_cast<char*>(&value), sizeof(float));
        xvec_h.emplace_back(value);
    }
    amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, xvec_h.begin(), xvec_h.end(), xvec_d.begin());

    for(int j=0; j<ny; j++) {
        infile.read(reinterpret_cast<char*>(&value), sizeof(float));
        yvec_h.emplace_back(value);
    }
    amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, yvec_h.begin(), yvec_h.end(), yvec_d.begin());

    for(int k=0; k<nz; k++) {
        infile.read(reinterpret_cast<char*>(&value), sizeof(float));
        zvec_h.emplace_back(value);
    }
    amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, zvec_h.begin(), zvec_h.end(), zvec_d.begin());

         // Vector to store the data

    Vector<Real>* data_h = nullptr; // Declare pointer outside the loop

    Real* xvec_d_ptr = xvec_d.data();
    Real* yvec_d_ptr = yvec_d.data();
    Real* zvec_d_ptr = zvec_d.data();

    Real dxvec = (xvec_h[nx-1]-xvec_h[0])/(nx-1);
    Real dyvec = (yvec_h[ny-1]-yvec_h[0])/(ny-1);

    // Read the file
    for(int idx=0; idx<ndata; idx++){
        if(idx == 0){
            data_h = &sst_h;
        } else if (idx==1) {
            data_h = &q_star_h;
        } else if (idx==2) {
            data_h = &t_star_h;
        } else if (idx==3) {
            data_h = &u_star_h;
        } else if(idx==4) {
            data_h = &ls_mask_h;
        }
        for(int k=0; k<nz; k++) {
            for(int j=0; j<ny; j++) {
                for(int i=0; i<nx; i++) {
                    infile.read(reinterpret_cast<char*>(&value), sizeof(float));
                    //if(idx == 3) {
                        //printf("theta is %0.15g, %0.15g, %0.15g %0.15g\n", xvec_h[i], yvec_h[j], zvec_h[k], value);
                    //}
                    data_h->emplace_back(value);
                }
            }
        }
    }

    infile.close();

    amrex::Gpu::DeviceVector<Real> sst_d(nx*ny*nz), q_star_d(nx*ny*nz), t_star_d(nx*ny*nz), u_star_d(nx*ny*nz), ls_mask_d(nx*ny*nz);

    amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, sst_h.begin(), sst_h.end(), sst_d.begin());
    amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, q_star_h.begin(), q_star_h.end(), q_star_d.begin());
    amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, t_star_h.begin(), t_star_h.end(), t_star_d.begin());
    amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, u_star_h.begin(), u_star_h.end(), u_star_d.begin());
    amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, ls_mask_h.begin(), ls_mask_h.end(), ls_mask_d.begin());

    Real* sst_d_ptr   = sst_d.data();
    Real* q_star_d_ptr  = q_star_d.data();
    Real* t_star_d_ptr  = t_star_d.data();
    Real* u_star_d_ptr  = u_star_d.data();
    Real* ls_mask_d_ptr = ls_mask_d.data();

    const auto prob_lo  = geom[lev].ProbLo();
    const auto dx       = geom[lev].CellSize();

        for (amrex::MFIter mfi(surface_state[lev][0]); mfi.isValid(); ++mfi) {
            const Box gtbx = mfi.growntilebox();

            const Array4<Real>& surf_arr = surface_state[lev][0].array(mfi);

            ParallelFor(gtbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {

                const Real x        = prob_lo[0] + (i + 0.5) * dx[0];
                const Real y        = prob_lo[1] + (j + 0.5) * dx[1];

                if(j<0) {
                    //printf("j is negative. value of y is %d %0.15g\n", j, y);
                }


                // First interpolate where the weather data is available from
                Real tmp_u_star, tmp_q_star, tmp_t_star;
                bilinear_interpolation_2d(xvec_d_ptr, yvec_d_ptr,
                                       dxvec, dyvec,
                                       nx, ny,
                                          x, y,
                                          u_star_d_ptr, tmp_u_star);

                bilinear_interpolation_2d(xvec_d_ptr, yvec_d_ptr,
                                       dxvec, dyvec,
                                       nx, ny,
                                          x, y,
                                          q_star_d_ptr, tmp_q_star);

                bilinear_interpolation_2d(xvec_d_ptr, yvec_d_ptr,
                                       dxvec, dyvec,
                                       nx, ny,
                                          x, y,
                                          t_star_d_ptr, tmp_t_star);



                surf_arr(i, j, 0) = tmp_u_star;
                surf_arr(i, j, 1) = -tmp_t_star/(1.226*1004.5*3600.0);
                surf_arr(i, j, 2) = -tmp_q_star/(1.226*2260000.0*3600.0);
                //printf("ustar is (%0.15g,%0.15g) is %0.15g %0.15g %0.15g\n", x, y, tmp_u_star, t_star_arr(i,j,0), q_star_arr(i,j,0));
            });
        }
}

void
ERF::CreateSurfaceStateMultiFabs(Vector<Vector<MultiFab>>& surface_state)
{
    surface_state.resize(max_level+1);
    for (int lev = 0; lev < max_level+1; ++lev) {
        surface_state[lev].resize(1);
        const MultiFab& src = vars_new[lev][0];

        amrex::BoxArray ba = src.boxArray();
        amrex::BoxList bl2d = ba.boxList();
        for (auto& b : bl2d) {
            b.setRange(2, 0);
        }
        amrex::BoxArray ba2d(std::move(bl2d));
        const amrex::DistributionMapping& dm = src.DistributionMap();

        surface_state[lev][0].define(ba2d, dm, 3, src.nGrow());
     }
}

void
ERF::SurfaceDataInterpolation(const Real time)
{
    static Real next_read_surface_time = -1.0;

    if (next_read_surface_time < 0.0) {
        int next_multiple = static_cast<int>(time / 10800);
        next_read_surface_time = next_multiple * 10800;
    }
    if (time >= next_read_surface_time) {

        std::string folder = solverChoice.hindcast_surface_data_dir;

        // Check if folder exists and is a directory
        if (!fs::exists(folder) || !fs::is_directory(folder)) {
            throw std::runtime_error("Error: Folder '" + folder + "' does not exist or is not a directory.");
        }

        std::vector<std::string> bin_files;

        for (const auto& entry : fs::directory_iterator(folder)) {
            if (!entry.is_regular_file()) continue;

            std::string fname = entry.path().filename().string();
            if (fname.size() >= 4 && fname.substr(fname.size() - 4) == ".bin") {
                bin_files.push_back(entry.path().string());
            }
        }
        std::sort(bin_files.begin(), bin_files.end());

    // Check if no .bin files were found
        if (bin_files.empty()) {
            throw std::runtime_error("Error: No .bin files found in folder '" + folder + "'.");
        }

        std::string filename1, filename2;

        int idx1 = static_cast<int>(time / 10800);
        int idx2 = static_cast<int>(time / 10800)+1;
        std::cout << "Reading surface data " << time << " " << idx1 << " " << idx2 <<" " << bin_files.size() << std::endl;

        if (idx2 >= static_cast<int>(bin_files.size())) {
            throw std::runtime_error("Error: Not enough .bin files to cover time " + std::to_string(time));
        }

        filename1 = bin_files[idx1];
        filename2 = bin_files[idx2];

        CreateSurfaceStateMultiFabs(surface_state_1);
        for(int lev = 0; lev< max_level+1; ++lev) {
            InterpSurfaceDataOntoMesh(lev, filename1, surface_state_1);
        }

        CreateSurfaceStateMultiFabs(surface_state_2);
        for(int lev = 0; lev< max_level+1; ++lev) {
            InterpSurfaceDataOntoMesh(lev, filename2, surface_state_2);
        }

        CreateSurfaceStateMultiFabs(surface_state_interp);

        next_read_surface_time += 10800;
    }
    Real alpha1 = 1.0 - (time - next_read_surface_time)/10800;
    Real alpha2 = 1.0 - alpha1;

    MultiFab& erf_mf_cons   = surface_state_interp[0][0];

    MultiFab::LinComb(surface_state_interp[0][0],
                      alpha1, surface_state_1[0][0], 0,
                      alpha2, surface_state_2[0][0], 0,
                      0, erf_mf_cons.nComp(), surface_state_interp[0][0].nGrow());

    /*Vector<std::string> varnames_plot_mf = {
    "rho", "rhotheta", "rhoqv", "rhoqc", "rhoqr", "xvel", "yvel", "zvel", "latitude", "longitude"
    }; // Customize variable names

    std::string pltname = "plt_interp";

    MultiFab plot_mf(erf_mf_cons.boxArray(), erf_mf_cons.DistributionMap(),
                     10, 0);

    plot_mf.setVal(0.0);

    for (MFIter mfi(plot_mf); mfi.isValid(); ++mfi) {
        const Array4<Real> &plot_mf_arr = plot_mf.array(mfi);
        const Array4<Real> &erf_mf_cons_arr = erf_mf_cons.array(mfi);
        const Array4<Real> &erf_mf_xvel_arr = erf_mf_xvel.array(mfi);
        const Array4<Real> &erf_mf_yvel_arr = erf_mf_yvel.array(mfi);
        const Array4<Real> &erf_mf_zvel_arr = erf_mf_zvel.array(mfi);
        const Array4<Real> &erf_mf_latlon_arr = erf_mf_latlon.array(mfi);

        const Box& bx = mfi.validbox();

        ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            plot_mf_arr(i,j,k,0) = erf_mf_cons_arr(i,j,k,Rho_comp);
            plot_mf_arr(i,j,k,1) = erf_mf_cons_arr(i,j,k,RhoTheta_comp);
            plot_mf_arr(i,j,k,2) = erf_mf_cons_arr(i,j,k,RhoQ1_comp);
            plot_mf_arr(i,j,k,3) = erf_mf_cons_arr(i,j,k,RhoQ2_comp);
            plot_mf_arr(i,j,k,4) = erf_mf_cons_arr(i,j,k,RhoQ3_comp);

            plot_mf_arr(i,j,k,5) = (erf_mf_xvel_arr(i,j,k,0) + erf_mf_xvel_arr(i+1,j,k,0))/2.0;
            plot_mf_arr(i,j,k,6) = (erf_mf_yvel_arr(i,j,k,0) + erf_mf_yvel_arr(i,j+1,k,0))/2.0;
            plot_mf_arr(i,j,k,7) = (erf_mf_zvel_arr(i,j,k,0) + erf_mf_zvel_arr(i,j,k+1,0))/2.0;

            plot_mf_arr(i,j,k,8) = erf_mf_latlon_arr(i,j,k,0);
            plot_mf_arr(i,j,k,9) = erf_mf_latlon_arr(i,j,k,1);
        });
    }


    WriteSingleLevelPlotfile(
            pltname,
            plot_mf,
            varnames_plot_mf,
            geom[0],
            time,
            0 // level
        );*/



}

#endif
