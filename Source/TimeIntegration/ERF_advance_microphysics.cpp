#include <ERF.H>

using namespace amrex;

void ERF::advance_microphysics (int lev,
                                MultiFab& cons,
                                const Real& dt_advance,
								const Real& time)
{
    if (solverChoice.moisture_type != MoistureType::None) {
        micro.Update_Micro_Vars_Lev(lev, cons);
        micro.Advance(lev, dt_advance, time);
        micro.Update_State_Vars_Lev(lev, cons);
    }
}
