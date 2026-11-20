#include <ERF.H>
#include <ERF_DataStruct.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_Utility.H>
#include "ERF_DA_EnKFSRF.H"

#include <vector>
#include <stdexcept>
#include <cassert>
#include <filesystem>
#include <string>

using namespace amrex;
namespace fs = std::filesystem;

// Simple matrix multiplication
Matrix
matrix_multiply(const Matrix& A, const Matrix& B)
{
    int n = A.size();
    Matrix C(n);

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {

            double sum = 0.0;
            for (int k = 0; k < n; ++k) {
                sum += A(i,k) * B(k,j);
            }

            C(i,j) = sum;
        }
    }

    return C;
}

// Print matrix
void matrix_print(const Matrix& A)
{
    int n = A.size();

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            std::cout << A(i,j) << " ";
        }
        std::cout << "\n";
    }
}

void
lapack_testing()
{
 Matrix S(3);
    S(0,0) = 4.0;  S(0,1) = 1.0;  S(0,2) = 1.0;
    S(1,0) = 1.0;  S(1,1) = 3.0;  S(1,2) = 0.5;
    S(2,0) = 1.0;  S(2,1) = 0.5;  S(2,2) = 2.5;

    Matrix Sinv = S.inverse();
    // fill Sinv ...

    Matrix T = Sinv.cholesky_lower();
    Matrix T_trans = T.transpose();

    Matrix T_T_trans = matrix_multiply(T, T_trans);
    if (ParallelDescriptor::IOProcessor()) {
        matrix_print(T_T_trans);
        std::cout << "T*T_trans done" << std::endl;
    }

    if (ParallelDescriptor::IOProcessor()) {
        std::cout  << "Sinv is " << std::endl;
        matrix_print(Sinv);
    }

    Matrix I_mat = matrix_multiply(S, Sinv);
    if (ParallelDescriptor::IOProcessor()) {
        std::cout << "Checking S*Sinv " << std::endl;
        matrix_print(I_mat);
    }
}

std::vector<std::string>
get_plotfile_list(const int da_iter)
{
    const std::string da_dir = "DA_Cycle_" + std::to_string(da_iter) + "/";

    std::vector<std::string> pltfiles;
    const std::string member_prefix = da_dir + "member_";

    std::string pf_dir = member_prefix + "00/plotfiles";

    if (!fs::exists(pf_dir)) {
        Abort("Plotfile directory not found: " + pf_dir);
    }

    for (const auto& entry : fs::directory_iterator(pf_dir)) {
        if (entry.is_directory()) {
            std::string name = entry.path().filename().string();
            if (name.rfind("plt", 0) == 0) {  // starts with "plt"
                pltfiles.push_back(name);
            }
        }
    }

    std::sort(pltfiles.begin(), pltfiles.end());
    return pltfiles;
}

// Reads the plotfile data into cell cenetred multifab
// Does not fill ghost cells
void
read_plot_file(PlotFileData& pf,
               const std::vector<std::string> varnames,
               MultiFab& mf)
{
    // ------------------------------------------------------------
    // Open plotfile
    // ------------------------------------------------------------
    const std::vector<std::string>& var_names_pf = pf.varNames();

    // ------------------------------------------------------------
    // Validate requested variables
    // ------------------------------------------------------------
    for (auto const& v : varnames) {
        bool found = false;
        for (auto const& vpf : var_names_pf) {
            if (v == vpf) {
                found = true;
                break;
            }
        }
        if (!found) {
            Abort("read_plot_file: invalid variable name: " + v);
        }
    }

    // ------------------------------------------------------------
    // Define destination MultiFab (single level only)
    // ------------------------------------------------------------
    const int level = 0;

    BoxArray ba = pf.boxArray(level);
    DistributionMapping dm(ba);

    int ncomp = varnames.size();

    mf.define(ba, dm, ncomp, 0);

    // ------------------------------------------------------------
    // Copy plotfile data → mf
    // ------------------------------------------------------------
    for (int comp = 0; comp < ncomp; ++comp)
    {
        const MultiFab& src = pf.get(level, varnames[comp]);
        MultiFab::Copy(mf, src, 0, comp, 1, 0);
    }
}

MultiFab
read_member_multifab(const int da_iter,
                     const int n,
                     const std::string& pf_name,
                     const Vector<std::string>& varnames)
{

    const std::string da_dir = "DA_Cycle_" + std::to_string(da_iter) + "/";

    const std::string member_prefix = da_dir + "member_";
    std::string member_dir = member_prefix + Concatenate("", n, 2);
    std::string pf_path = member_dir + "/plotfiles/" + pf_name;

    PlotFileData pf(pf_path);

    const BoxArray& ba = pf.boxArray(0);
    const DistributionMapping& dm = pf.DistributionMap(0);
    int ncomp = varnames.size();

    MultiFab mf(ba, dm, ncomp, 0);

    read_plot_file(pf, varnames, mf);

    return mf;
}

MultiFab
compute_ensemble_mean(const int da_iter,
                      const int Nens,
                      const std::string& pf_name,
                      const Vector<std::string>& varnames)
{
    MultiFab mf_mean;
    bool initialized = false;

    for (int n = 0; n < Nens; ++n)
    {
        std::cout << "Reading ensembles for compute_ensemble_mean = " << n << std::endl;
        MultiFab mf_tmp = read_member_multifab(da_iter, n, pf_name, varnames);

        if (!initialized) {
            mf_mean.define(mf_tmp.boxArray(),
                           mf_tmp.DistributionMap(),
                           mf_tmp.nComp(),
                           mf_tmp.nGrow());
            mf_mean.setVal(0.0);
            initialized = true;
        }

        MultiFab::Add(mf_mean, mf_tmp, 0, 0, mf_tmp.nComp(), mf_tmp.nGrow());
    }

    mf_mean.mult(1.0 / Real(Nens));
    return mf_mean;
}

void
compute_d_vec(const MultiFab& mf1,
              const MultiFab& mf2,
              MultiFab& d_vec)
{
    AMREX_ALWAYS_ASSERT(mf1.boxArray() == mf2.boxArray());
    AMREX_ALWAYS_ASSERT(mf1.DistributionMap() == mf2.DistributionMap());
    AMREX_ALWAYS_ASSERT(mf1.nComp() == mf2.nComp());
    AMREX_ALWAYS_ASSERT(mf1.nGrowVect() == mf2.nGrowVect());

    d_vec.define(mf1.boxArray(),
                 mf1.DistributionMap(),
                 mf1.nComp(),
                 mf1.nGrowVect());

    MultiFab::Copy(d_vec,
                   mf1,
                   0, 0,
                   mf1.nComp(),
                   mf1.nGrowVect());

    MultiFab::Subtract(d_vec,
                       mf2,
                       0, 0,
                       mf1.nComp(),
                       mf1.nGrowVect());
}

IntVect
get_obs_thinning ()
{
    int sx = 8, sy = 1, sz = 8;
    ParmParse pp("da");
    pp.queryAdd("obs_thin_x", sx);
    pp.queryAdd("obs_thin_y", sy);
    pp.queryAdd("obs_thin_z", sz);
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(sx >= 1 && sy >= 1 && sz >= 1,
                                     "da.obs_thin_x/y/z must be >= 1");
    return IntVect(sx, sy, sz);
}

void
compute_d_prime_vec(MultiFab& d_prime_vec,
                    const MultiFab& d_vec,
                    const Vector<Real>& R_diag)
{
    // Observations at cells not selected by the thinning get zero weight (d' = 0)
    const IntVect thin = get_obs_thinning();
    const int sx = thin[0], sy = thin[1], sz = thin[2];

    AMREX_ALWAYS_ASSERT(d_vec.nComp() == R_diag.size());

    const int ncomp = d_vec.nComp();

    d_prime_vec.define(d_vec.boxArray(),
                       d_vec.DistributionMap(),
                       ncomp,
                       d_vec.nGrowVect());

    Gpu::DeviceVector<Real> R_diag_d(R_diag.size());
    Gpu::copy(Gpu::hostToDevice, R_diag.begin(), R_diag.end(), R_diag_d.begin());
    const Real* R_diag_d_ptr = R_diag_d.data();

    for (MFIter mfi(d_vec, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        const Array4<const Real> d_arr = d_vec.const_array(mfi);
        const Array4<Real> dp_arr = d_prime_vec.array(mfi);

        ParallelFor(bx, ncomp,
                            [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
        {
            const bool selected = (i % sx == 0) && (j % sy == 0) && (k % sz == 0);
            dp_arr(i,j,k,n) = selected ? d_arr(i,j,k,n) / R_diag_d_ptr[n] : 0.0;
        });
    }
}

void
compute_mean_H_xf(const int da_iter,
                  MultiFab& mean_H_xf,
                  const int Nens,
                  const std::string& last_pf_name,
                  const Vector<std::string>& varnames)
{

    for (int n = 0; n < Nens; ++n)
    {
        MultiFab xf_i = read_member_multifab(da_iter, n, last_pf_name, varnames);
        if (n==0) {
            mean_H_xf.define(xf_i.boxArray(),
                             xf_i.DistributionMap(),
                             8,
                             xf_i.nGrow());
            mean_H_xf.setVal(0.0);
        }

        MultiFab H_xf_i;

        Apply_H(xf_i, H_xf_i);

        if(n==0){
            MultiFab sum_H_xf_i(H_xf_i.boxArray(), H_xf_i.DistributionMap(), H_xf_i.nComp(), H_xf_i.nGrow());
        }

        MultiFab::Add(mean_H_xf, H_xf_i, 0, 0, H_xf_i.nComp(), H_xf_i.nGrow());
    }

    mean_H_xf.mult(1.0 / Nens);
}

Real
Compute_yf_prime_i_T_Rinv_yf_prime_j(const MultiFab& Yi,
                                   const MultiFab& Yj,
                                   const Vector<Real>& R_diag)
{
    const int ncomp = Yi.nComp();

    // Only the cells selected by the observation thinning contribute
    const IntVect thin = get_obs_thinning();
    const int sx = thin[0], sy = thin[1], sz = thin[2];

    ReduceOps<ReduceOpSum> reduce_op;
    ReduceData<Real> reduce_data(reduce_op);

    using ReduceTuple = typename decltype(reduce_data)::Type;

    Gpu::DeviceVector<Real> R_diag_d(R_diag.size());
    Gpu::copy(Gpu::hostToDevice, R_diag.begin(), R_diag.end(), R_diag_d.begin());
    const Real* R_diag_d_ptr = R_diag_d.data();

    for (MFIter mfi(Yi, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();

        auto const& yi = Yi.const_array(mfi);
        auto const& yj = Yj.const_array(mfi);

        reduce_op.eval(bx, reduce_data,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
        {
            Real local = 0.0;

            if ((i % sx == 0) && (j % sy == 0) && (k % sz == 0)) {
                for (int n = 0; n < ncomp; ++n) {
                    local += yi(i,j,k,n) * yj(i,j,k,n)/R_diag_d_ptr[n];
                }
            }

            return {local};
        });
    }

    ReduceTuple hv = reduce_data.value();
    Real result = get<0>(hv);

    ParallelDescriptor::ReduceRealSum(result);

    return result;
}

void
compute_yf_prime (const int da_iter,
                  const int i,
                  const std::string& last_pf_name,
                  const Vector<std::string>& varnames,
                  const MultiFab& mean_H_xf,
                  MultiFab& yf_prime_i)
{
    // Read ensemble member
    MultiFab xf_i = read_member_multifab(da_iter, i, last_pf_name, varnames);

    // Apply observation operator
    MultiFab H_xf_i;
    Apply_H(xf_i, H_xf_i);

    // Allocate output
    yf_prime_i.define(H_xf_i.boxArray(),
                      H_xf_i.DistributionMap(),
                      H_xf_i.nComp(),
                      H_xf_i.nGrow());

    // y'_i = H(xf_i) - mean(H(xf))
    MultiFab::Copy(yf_prime_i,
                   H_xf_i,
                   0, 0,
                   H_xf_i.nComp(),
                   H_xf_i.nGrow());

    MultiFab::Subtract(yf_prime_i,
                       mean_H_xf,
                       0, 0,
                       H_xf_i.nComp(),
                       H_xf_i.nGrow());
}

void
compute_S_matrix(const int da_iter,
                 Matrix& S,
                 const int& Nens,
                 const MultiFab& mean_H_xf,
                 const Vector<Real>& R_diag,
                 const std::string& last_pf_name,
                 const Vector<std::string>& varnames)
{
    for (int i = 0; i < Nens; ++i) {
        MultiFab yf_prime_i;
        compute_yf_prime(da_iter, i, last_pf_name, varnames, mean_H_xf, yf_prime_i);

        for (int j = 0; j < Nens; ++j) {
            MultiFab yf_prime_j;
            compute_yf_prime(da_iter, j, last_pf_name, varnames, mean_H_xf, yf_prime_j);
            Real val = Compute_yf_prime_i_T_Rinv_yf_prime_j(yf_prime_i, yf_prime_j, R_diag);
            S(i,j) = val/(Nens-1);
            if(i==j) {
                S(i,j) = 1.0 + S(i,j);
            }
        }
    }
}

Real
compute_yf_prime_T_d_prime_vec (const MultiFab& yf_prime,
                                const MultiFab& d_prime_vec)
{
    AMREX_ALWAYS_ASSERT(yf_prime.nComp() ==
                        d_prime_vec.nComp());

    const int ncomp = yf_prime.nComp();

    ReduceOps<ReduceOpSum> reduce_op;
    ReduceData<Real> reduce_data(reduce_op);

    using ReduceTuple =
        typename decltype(reduce_data)::Type;

    for (MFIter mfi(yf_prime, TilingIfNotGPU());
         mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();

        auto const& yf_arr =
            yf_prime.const_array(mfi);

        auto const& d_arr =
            d_prime_vec.const_array(mfi);

        reduce_op.eval(
        bx,
        reduce_data,
        [=] AMREX_GPU_DEVICE (int i, int j, int k)
            -> ReduceTuple
        {
            Real local = 0.0;

            for (int n = 0; n < ncomp; ++n)
            {
                local +=
                    yf_arr(i,j,k,n) *
                    d_arr(i,j,k,n);
            }

            return {local};
        });
    }

    ReduceTuple hv = reduce_data.value();

    Real r = get<0>(hv);

    // Global MPI reduction
    ParallelDescriptor::ReduceRealSum(r);

    return r;
}

void
compute_r_vec (const int da_iter,
               int Nens,
               const std::string& last_pf_name,
               const Vector<std::string>& varnames,
               const MultiFab& mean_H_xf,
               const MultiFab& d_prime_vec,
               Vector<Real>& r_vec)
{
    r_vec.resize(Nens);

    for (int i = 0; i < Nens; ++i)
    {
        MultiFab yf_prime_i;

        compute_yf_prime(da_iter, i, last_pf_name, varnames, mean_H_xf, yf_prime_i);

        r_vec[i] = compute_yf_prime_T_d_prime_vec(yf_prime_i, d_prime_vec);
    }
}


// Eigenvalues of S = I + Y'^T R^-1 Y'/(N-1) (symmetric positive definite, S >= I).
// T = S^{-1/2} scales the analysis perturbations along eigenvector k by lambda_k^{-1/2},
// so large eigenvalues mean strong spread collapse. One eigenvalue should be exactly 1
// (the eigenvector 1/sqrt(N), since S*1 = 1); the other N-1 are >= 1.
void
print_S_eigenvalues (const Matrix& S_mat)
{
    int n = S_mat.size();   // non-const: LAPACK takes int*

    // Symmetrize (dsyev reads the lower triangle) into a work copy
    Matrix A(n);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            A(i,j) = 0.5 * (S_mat(i,j) + S_mat(j,i));
        }
    }

    std::vector<double> w(n);
    int lda = n;
    int info;
    char jobz = 'N';   // eigenvalues only
    char uplo = 'L';

    // workspace query
    int lwork = -1;
    double work_query;
    dsyev_(&jobz, &uplo, &n, A.data(), &lda, w.data(), &work_query, &lwork, &info);

    lwork = static_cast<int>(work_query);
    std::vector<double> work(lwork);

    dsyev_(&jobz, &uplo, &n, A.data(), &lda, w.data(), work.data(), &lwork, &info);

    if (info != 0) {
        Print() << "print_S_eigenvalues: dsyev failed, info = " << info << std::endl;
        return;
    }

    // Eigenvalues are returned in ascending order
    double trace_minus_I = 0.0;
    for (int k = 0; k < n; ++k) { trace_minus_I += (w[k] - 1.0); }

    Print() << "---- Eigenvalues of S (N = " << n << ") ----" << std::endl;
    Print() << std::setprecision(8) << std::scientific;
    for (int k = 0; k < n; ++k) {
        Print() << "  lambda[" << k << "] = " << w[k]
                << "   spread factor 1/sqrt(lambda) = " << 1.0 / std::sqrt(w[k])
                << std::endl;
    }
    Print() << "  min = " << w[0] << ", max = " << w[n-1]
            << ", trace(S-I) = " << trace_minus_I << std::endl;
    Print() << "  (expect one eigenvalue = 1, others >= 1; "
            << "analysis spread along eigenvector k is scaled by 1/sqrt(lambda_k))"
            << std::endl;
    Print() << std::defaultfloat;
}

void
compute_alpha_vec (const int& Nens,
                   const Matrix& S,
                   const Vector<Real>& r_vec,
                   Vector<Real>& alpha_vec)
{
    AMREX_ALWAYS_ASSERT(r_vec.size() == Nens);

    Matrix Sinv = S.inverse();

    alpha_vec.resize(Nens, 0.0);

    const Real fac = 1.0_rt / static_cast<Real>(Nens - 1);

    for (int i = 0; i < Nens; ++i)
    {
        Real sum = 0.0_rt;

        for (int j = 0; j < Nens; ++j)
        {
            sum += Sinv(i,j) * r_vec[j];
        }

        alpha_vec[i] = fac * sum;
    }
}

// Perform a matrix time a small vector multiply as a summation of
// column-vector element multiply

void
compute_Xf_prime_times_vector (const int da_iter,
                               const int Nens,
                               const std::string& last_pf_name,
                               const Vector<std::string>& varnames,
                               const MultiFab& xf_bar,
                               const Vector<Real>& vec_in,
                               MultiFab& result)
{
    AMREX_ALWAYS_ASSERT(vec_in.size() == Nens);

    // Read first member to define result
    MultiFab xf_0 =
        read_member_multifab(da_iter,
                             0,
                             last_pf_name,
                             varnames);

    result.define(xf_0.boxArray(),
                  xf_0.DistributionMap(),
                  xf_0.nComp(),
                  xf_0.nGrow());

    result.setVal(0.0);

    for (int n = 0; n < Nens; ++n)
    {
        MultiFab xf_n = read_member_multifab(da_iter, n, last_pf_name, varnames);

        AMREX_ALWAYS_ASSERT(xf_n.boxArray() == result.boxArray());
        AMREX_ALWAYS_ASSERT(xf_n.DistributionMap() == result.DistributionMap());
        AMREX_ALWAYS_ASSERT(xf_n.ixType() == result.ixType());

        AMREX_ALWAYS_ASSERT(xf_bar.boxArray() == result.boxArray());
        AMREX_ALWAYS_ASSERT(xf_bar.DistributionMap() == result.DistributionMap());

        // xf_n_prime = xf_n - xf_bar
        MultiFab xf_n_prime(xf_n.boxArray(),
                            xf_n.DistributionMap(),
                            xf_n.nComp(),
                            xf_n.nGrow());

        MultiFab::Copy(xf_n_prime,
                       xf_n,
                       0, 0,
                       xf_n.nComp(),
                       xf_n.nGrow());

        MultiFab::Subtract(xf_n_prime,
                           xf_bar,
                           0, 0,
                           xf_n.nComp(),
                           xf_n.nGrow());

        // result += vec_in[n] * xf_n_prime
        MultiFab::Saxpy(result,
                        vec_in[n],
                        xf_n_prime,
                        0, 0,
                        xf_n.nComp(),
                        xf_n.nGrow());
    }
}

void
add_multifabs (const MultiFab& xf_bar,
              const MultiFab& Xf_prime_alpha,
              MultiFab& result)
{
    result.define(
        xf_bar.boxArray(),
        xf_bar.DistributionMap(),
        xf_bar.nComp(),
        xf_bar.nGrowVect()
    );
    result.setVal(0.0);
    MultiFab::Copy(result, xf_bar, 0, 0, xf_bar.nComp(), 0);

    MultiFab::Add(result, Xf_prime_alpha,
                  0, 0, Xf_prime_alpha.nComp(), 0);
}

void
compute_T_matrix (const Matrix& S_mat,
                  Matrix& T_mat)
{
    // T = S^{-1/2} (symmetric square root), so that T*T^T = S^{-1} and T*1 = 1.
    // A Cholesky factor of S^{-1} also gives T*T^T = S^{-1} but does not
    // preserve the zero-mean property of the analysis perturbations.
    T_mat = S_mat.inverse_sqrt_symmetric();

    // Sanity check: T*1 should equal 1 (S*1 = 1 because Y'*1 = 0)
    const int n = T_mat.size();
    Real max_dev = 0.0;
    for (int i = 0; i < n; ++i) {
        Real row_sum = 0.0;
        for (int j = 0; j < n; ++j) {
            row_sum += T_mat(i,j);
        }
        max_dev = std::max(max_dev, std::abs(row_sum - 1.0));
    }
    Print() << "compute_T_matrix: max |T*1 - 1| = " << max_dev << std::endl;
}


void
update_ensemble (const int da_iter,
                 const int Nens,
                 const std::string& last_pf_name,
                 const Vector<std::string>& varnames,
                 const MultiFab& xf_bar,
                 const Matrix& T,
                 const int& n,
                 MultiFab& result)
{
    Vector<Real> T_colvec(Nens);
    for (int i = 0; i < Nens; ++i) {
        T_colvec[i] = T(i,n);
    }
    compute_Xf_prime_times_vector(da_iter, Nens, last_pf_name, varnames, xf_bar, T_colvec, result);
}

// Averages DA_Cycle_<da_iter>/member_XX/plt_ens_update over all members and writes
// DA_Cycle_<da_iter>/plt_ens_update_mean. Each plt_ens_update holds X'T[:,n], so the
// mean over n equals X'(T*1)/N = X'*1/N = 0 (to round-off) for the symmetric square root.
void
ERF::ComputeEnsembleUpdateMean (const int da_iter)
{
    const int Nens = solverChoice.n_ensemble;
    const Vector<std::string> varnames = {"density","theta", "x_velocity","y_velocity",
                                          "z_velocity", "qv", "qc", "qrain"};
    const int ncomp = varnames.size();

    MultiFab mf_sum;

    for (int n = 0; n < Nens; ++n)
    {
        const std::string member_dir = "member_" + amrex::Concatenate("", n, 2);
        const std::string pf_path    = MakeFullPath(member_dir + "/plt_ens_update", da_iter);

        PlotFileData pf(pf_path);
        MultiFab mf_n;
        read_plot_file(pf, varnames, mf_n);

        if (n == 0) {
            mf_sum.define(mf_n.boxArray(), mf_n.DistributionMap(), ncomp, 0);
            mf_sum.setVal(0.0);
        }
        MultiFab::Add(mf_sum, mf_n, 0, 0, ncomp, 0);
    }

    mf_sum.mult(1.0 / Real(Nens));

    // Report max |mean update| per variable (should be ~0, round-off level)
    for (int c = 0; c < ncomp; ++c) {
        Print() << "ensemble update mean, max |" << varnames[c] << "| = "
                << mf_sum.norm0(c) << std::endl;
    }

    const std::string pltname = MakeFullPath("plt_ens_update_mean", da_iter);
    WriteSingleLevelPlotfile(pltname, mf_sum, varnames, geom[0], 0.0, 0);
}

void
ERF::SetDirsForPlotfilesAndCheckpointsForDA (const int da_iter,
                                             const int ens_no)
{
    // --------------------------------------------------------
    // Set up cycle- and member-specific output directories
    // --------------------------------------------------------
    std::stringstream cycle_ss;
    cycle_ss << "DA_Cycle_" << da_iter;
    std::string cycle_dir = cycle_ss.str();

    std::stringstream member_ss;
    member_ss << "member_" << std::setw(2) << std::setfill('0') << ens_no;
    std::string member_dir = cycle_dir + "/" + member_ss.str();

    if (ParallelDescriptor::IOProcessor())
    {
        fs::create_directories(member_dir + "/plotfiles");
        fs::create_directories(member_dir + "/chkfiles");
        fs::create_directories(member_dir + "/pertfiles");
    }

    ParallelDescriptor::Barrier();

    // Set file prefixes under the specific cycle and member path
    check_file   = member_dir + "/chkfiles/chk";
    plot3d_file_1 = member_dir + "/plotfiles/plt";
}
// da_iter : DA iteration number
// ens_no  : Ensemble number

std::string
MakeEnsembleCheckpointName (int da_iter, int ens_no)
{
    // Create CheckpointAfterDA if needed
    const std::string da_dir = "DA_Cycle_" + std::to_string(da_iter) + "/";
    fs::create_directories(da_dir);

    const std::string chk_after_da_dir = da_dir + "CheckpointAfterDA";

    // Create DA_Cycle_<da_iter>/CheckpointAfterDA
    fs::create_directories(chk_after_da_dir);

    std::ostringstream chk_name;
    chk_name << chk_after_da_dir
             << "/chk_"
             << std::setw(2) << std::setfill('0') << ens_no
             << "_";

    return chk_name.str();
}

void
ERF::GetEnsembleCheckpointName (int da_iter, int ens_no)
{
    // Create CheckpointAfterDA if needed
    const std::string da_dir = "DA_Cycle_" + std::to_string(da_iter) + "/";
    fs::create_directories(da_dir);

    // Create CheckpointAfterDA/DA_Cycle_<da_iter>
    const std::string chk_after_da_dir = da_dir + "CheckpointAfterDA";
    fs::create_directories(chk_after_da_dir);

    std::ostringstream chk_name;
    chk_name << chk_after_da_dir
             << "/chk_"
             << std::setw(2) << std::setfill('0') << ens_no
             << "_" <<  "00000";

    restart_chkfile = chk_name.str();
}

std::string
MakeFullPath (const std::string& filename, const int da_iter)
{
    std::stringstream ss;
    ss << "DA_Cycle_" << da_iter << "/" << filename;
    return ss.str();
}

void WriteInnovationRMS(const MultiFab& dvec,
                        const std::string& filename,
                        int da_iter)
{
    Real sum_sq = 0.0;
    long count = 0;

    for (MFIter mfi(dvec); mfi.isValid(); ++mfi)
    {
        const auto& arr = dvec.const_array(mfi);
        const auto& box = mfi.validbox();

        for (int n = 0; n < dvec.nComp(); ++n)
        {
            ReduceOps<ReduceOpSum> reduce_op;
            ReduceData<Real> reduce_data(reduce_op);
            using ReduceTuple = typename decltype(reduce_data)::Type;

            reduce_op.eval(box, reduce_data,
                [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                -> ReduceTuple
                {
                    const Real val = arr(i, j, k, n);
                    return {val * val};
                });

            sum_sq += get<0>(reduce_data.value());
        }

        count += static_cast<long>(box.numPts()) * dvec.nComp();
    }

    // Global reductions across MPI ranks
    ParallelDescriptor::ReduceRealSum(sum_sq);
    ParallelDescriptor::ReduceLongSum(count);

    if (ParallelDescriptor::IOProcessor())
    {
        const Real rms =
            (count > 0) ? std::sqrt(sum_sq / count) : 0.0;

        std::ofstream out(filename, std::ios::app);
        out << da_iter << " "
            << std::setprecision(16) << rms << '\n';
    }
}

