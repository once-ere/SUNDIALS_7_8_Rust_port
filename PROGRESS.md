# PROGRESS — per-file port checklist

Status legend: todo | ported | building | committed
(impl headers and public include/ headers port together with their module and share its line)

## Phase 1 — sundials_core

- [x] src/sundials/sundials_math.c — committed
- [x] src/sundials/sundials_errors.c — committed
- [x] src/sundials/sundials_context.c — committed
- [x] src/sundials/sundials_nvector.c — committed
- [x] src/sundials/sundials_matrix.c — committed
- [x] src/sundials/sundials_direct.c — committed
- [x] src/sundials/sundials_band.c — committed
- [x] src/sundials/sundials_dense.c — committed
- [x] src/sundials/sundials_iterative.c — committed
- [x] src/sundials/sundials_linearsolver.c — committed
- [x] src/sundials/sundials_nonlinearsolver.c — committed
- [x] src/sundials/sundials_nvector_senswrapper.c — committed
- [x] src/sundials/sundials_memory.c — committed
- [x] src/sundials/sundials_logger.c — committed
- [x] src/sundials/sundials_profiler.c — committed
- [x] src/sundials/sundials_futils.c — committed
- [x] src/sundials/sundials_hashmap.c — committed
- [x] src/sundials/sundials_version.c — committed
- [x] src/sundials/sundials_cli.c — committed
- [x] src/sundials/sundials_adaptcontroller.c — committed
- [x] src/sundials/sundials_stepper.c — committed
- [x] src/sundials/sundials_adjointstepper.c — committed
- [x] src/sundials/sundials_adjointcheckpointscheme.c — committed
- [x] src/sundials/sundials_datanode.c — committed
- [x] src/sundials/sundials_domeigestimator.c — committed
- [x] src/sundials/sundatanode/sundatanode_inmem.c — committed
- [x] src/sundials/stl/sunstl_vector.h — committed (pulled forward: hashmap needs it)
- [x] src/nvector/serial/nvector_serial.c — committed
- [x] src/sunmatrix/band/sunmatrix_band.c — committed
- [x] src/sunmatrix/dense/sunmatrix_dense.c — committed
- [x] src/sunmatrix/sparse/sunmatrix_sparse.c — committed
- [x] src/sunlinsol/band/sunlinsol_band.c — committed
- [x] src/sunlinsol/dense/sunlinsol_dense.c — committed
- [x] src/sunlinsol/pcg/sunlinsol_pcg.c — committed
- [x] src/sunlinsol/spbcgs/sunlinsol_spbcgs.c — committed
- [x] src/sunlinsol/spfgmr/sunlinsol_spfgmr.c — committed
- [x] src/sunlinsol/spgmr/sunlinsol_spgmr.c — committed
- [x] src/sunlinsol/sptfqmr/sunlinsol_sptfqmr.c — committed
- [x] src/sunnonlinsol/newton/sunnonlinsol_newton.c — committed
- [x] src/sunnonlinsol/fixedpoint/sunnonlinsol_fixedpoint.c — committed
- [x] src/sunnonlinsol/auto/sunnonlinsol_auto.c — committed
- [x] src/sunadaptcontroller/soderlind/sunadaptcontroller_soderlind.c — committed
- [x] src/sunadaptcontroller/imexgus/sunadaptcontroller_imexgus.c — committed
- [x] src/sunadaptcontroller/mrihtol/sunadaptcontroller_mrihtol.c — committed
- [x] src/sundomeigest/power/sundomeigest_power.c — committed
- [x] src/sundomeigest/arnoldi/sundomeigest_arnoldi.c — committed
- [x] src/sunadjointcheckpointscheme/fixed/sunadjointcheckpointscheme_fixed.c — committed
- [x] src/sunmemory/system/sundials_system_memory.c — committed

## Phase 2 — cvode

- [x] src/cvode/cvode.c — committed
- [x] src/cvode/cvode_bandpre.c — committed
- [x] src/cvode/cvode_bbdpre.c — committed
- [x] src/cvode/cvode_cli.c — committed
- [x] src/cvode/cvode_diag.c — committed
- [x] src/cvode/cvode_fused_stubs.c — committed
- [x] src/cvode/cvode_io.c — committed
- [x] src/cvode/cvode_ls.c — committed
- [x] src/cvode/cvode_nls.c — committed
- [x] src/cvode/cvode_proj.c — committed
- [x] src/cvode/cvode_resize.c — committed

## Phase 3 — cvodes

- [x] src/cvodes/cvodea.c — building
- [x] src/cvodes/cvodea_io.c — building
- [x] src/cvodes/cvodes.c — building
- [x] src/cvodes/cvodes_bandpre.c — building
- [x] src/cvodes/cvodes_bbdpre.c — building
- [x] src/cvodes/cvodes_cli.c — building
- [x] src/cvodes/cvodes_diag.c — building
- [x] src/cvodes/cvodes_io.c — building
- [x] src/cvodes/cvodes_ls.c — building
- [x] src/cvodes/cvodes_nls.c — building
- [x] src/cvodes/cvodes_nls_sim.c — building
- [x] src/cvodes/cvodes_nls_stg.c — building
- [x] src/cvodes/cvodes_nls_stg1.c — building
- [x] src/cvodes/cvodes_proj.c — building
- [x] src/cvodes/cvodes_resize.c — building

## Phase 4 — kinsol

- [x] src/kinsol/kinsol.c — building
- [x] src/kinsol/kinsol_aa.c — building
- [x] src/kinsol/kinsol_bbdpre.c — building
- [x] src/kinsol/kinsol_cli.c — building
- [x] src/kinsol/kinsol_io.c — building
- [x] src/kinsol/kinsol_ls.c — building
- [x] src/kinsol/kinsol_orth.c — building

## Phase 5 — ida

- [x] src/ida/ida.c — building
- [x] src/ida/ida_bbdpre.c — building
- [x] src/ida/ida_cli.c — building
- [x] src/ida/ida_ic.c — building
- [x] src/ida/ida_io.c — building
- [x] src/ida/ida_ls.c — building
- [x] src/ida/ida_nls.c — building

## Phase 6 — idas

- [ ] src/idas/idaa.c — todo
- [ ] src/idas/idaa_io.c — todo
- [ ] src/idas/idas.c — todo
- [ ] src/idas/idas_bbdpre.c — todo
- [ ] src/idas/idas_cli.c — todo
- [ ] src/idas/idas_ic.c — todo
- [ ] src/idas/idas_io.c — todo
- [ ] src/idas/idas_ls.c — todo
- [ ] src/idas/idas_nls.c — todo
- [ ] src/idas/idas_nls_sim.c — todo
- [ ] src/idas/idas_nls_stg.c — todo

## Phase 7 — arkode

- [ ] src/arkode/arkode.c — todo
- [ ] src/arkode/arkode_adapt.c — todo
- [ ] src/arkode/arkode_arkstep.c — todo
- [ ] src/arkode/arkode_arkstep_io.c — todo
- [ ] src/arkode/arkode_arkstep_nls.c — todo
- [ ] src/arkode/arkode_bandpre.c — todo
- [ ] src/arkode/arkode_bbdpre.c — todo
- [ ] src/arkode/arkode_butcher.c — todo
- [ ] src/arkode/arkode_butcher_dirk.c — todo
- [ ] src/arkode/arkode_butcher_erk.c — todo
- [ ] src/arkode/arkode_cli.c — todo
- [ ] src/arkode/arkode_erkstep.c — todo
- [ ] src/arkode/arkode_erkstep_io.c — todo
- [ ] src/arkode/arkode_forcingstep.c — todo
- [ ] src/arkode/arkode_interp.c — todo
- [ ] src/arkode/arkode_io.c — todo
- [ ] src/arkode/arkode_ls.c — todo
- [ ] src/arkode/arkode_lsrkstep.c — todo
- [ ] src/arkode/arkode_lsrkstep_io.c — todo
- [ ] src/arkode/arkode_mri_tables.c — todo
- [ ] src/arkode/arkode_mristep.c — todo
- [ ] src/arkode/arkode_mristep_controller.c — todo
- [ ] src/arkode/arkode_mristep_io.c — todo
- [ ] src/arkode/arkode_mristep_nls.c — todo
- [ ] src/arkode/arkode_relaxation.c — todo
- [ ] src/arkode/arkode_root.c — todo
- [ ] src/arkode/arkode_splittingstep.c — todo
- [ ] src/arkode/arkode_splittingstep_coefficients.c — todo
- [ ] src/arkode/arkode_sprk.c — todo
- [ ] src/arkode/arkode_sprkstep.c — todo
- [ ] src/arkode/arkode_sprkstep_io.c — todo
- [ ] src/arkode/arkode_sunstepper.c — todo
- [ ] src/arkode/arkode_user_controller.c — todo
- [ ] src/arkode/arkode_butcher_dirk.def — todo
- [ ] src/arkode/arkode_butcher_erk.def — todo
- [ ] src/arkode/arkode_mri_tables.def — todo
- [ ] src/arkode/arkode_splittingstep_coefficients.def — todo

## Example programs (one line per ported program; variants tracked in VERIFICATION.md)

- [ ] arkode_rs example ark_KrylovDemo_prec — todo
- [ ] arkode_rs example ark_advection_diffusion_reaction_splitting — todo
- [ ] arkode_rs example ark_analytic — todo
- [ ] arkode_rs example ark_analytic_lsrk — todo
- [ ] arkode_rs example ark_analytic_lsrk_domeigest — todo
- [ ] arkode_rs example ark_analytic_lsrk_varjac — todo
- [ ] arkode_rs example ark_analytic_mels — todo
- [ ] arkode_rs example ark_analytic_nonlin — todo
- [ ] arkode_rs example ark_analytic_partitioned — todo
- [ ] arkode_rs example ark_analytic_ssprk — todo
- [ ] arkode_rs example ark_brusselator — todo
- [ ] arkode_rs example ark_brusselator1D — todo
- [ ] arkode_rs example ark_brusselator1D_imexmri — todo
- [ ] arkode_rs example ark_brusselator_1D_mri — todo
- [ ] arkode_rs example ark_brusselator_fp — todo
- [ ] arkode_rs example ark_brusselator_lsrk_domeigest — todo
- [ ] arkode_rs example ark_brusselator_lsrk_externaldomeigest — todo
- [ ] arkode_rs example ark_brusselator_mri — todo
- [ ] arkode_rs example ark_conserved_exp_entropy_ark — todo
- [ ] arkode_rs example ark_conserved_exp_entropy_erk — todo
- [ ] arkode_rs example ark_damped_harmonic_symplectic — todo
- [ ] arkode_rs example ark_dissipated_exp_entropy — todo
- [ ] arkode_rs example ark_harmonic_symplectic — todo
- [ ] arkode_rs example ark_heat1D — todo
- [ ] arkode_rs example ark_heat1D_adapt — todo
- [ ] arkode_rs example ark_kepler — todo
- [ ] arkode_rs example ark_kpr_mri — todo
- [ ] arkode_rs example ark_lotka_volterra_ASA — todo
- [ ] arkode_rs example ark_onewaycouple_mri — todo
- [ ] arkode_rs example ark_reaction_diffusion_mri — todo
- [ ] arkode_rs example ark_robertson — todo
- [ ] arkode_rs example ark_robertson_constraints — todo
- [ ] arkode_rs example ark_robertson_root — todo
- [ ] arkode_rs example ark_twowaycouple_mri — todo
- [ ] cvode_rs example cvAdvDiff_bnd — todo
- [ ] cvode_rs example cvAdvDiff_bndL — todo
- [ ] cvode_rs example cvAnalytic_mels — todo
- [ ] cvode_rs example cvDirectDemo_ls — todo
- [ ] cvode_rs example cvDisc_dns — todo
- [ ] cvode_rs example cvDiurnal_kry — todo
- [ ] cvode_rs example cvDiurnal_kry_bp — todo
- [ ] cvode_rs example cvKrylovDemo_ls — todo
- [ ] cvode_rs example cvKrylovDemo_prec — todo
- [ ] cvode_rs example cvParticle_dns — todo
- [ ] cvode_rs example cvPendulum_dns — todo
- [ ] cvode_rs example cvRoberts_dns — todo
- [ ] cvode_rs example cvRoberts_dnsL — todo
- [ ] cvode_rs example cvRoberts_dns_constraints — todo
- [ ] cvode_rs example cvRoberts_dns_negsol — todo
- [ ] cvode_rs example cvRoberts_dns_uw — todo
- [ ] cvode_rs example cvRocket_dns — todo
- [ ] cvode_rs example cvVdp_auto_nls — todo
- [ ] cvodes_rs example cvsAdvDiff_ASAi_bnd — todo
- [ ] cvodes_rs example cvsAdvDiff_FSA_non — todo
- [ ] cvodes_rs example cvsAdvDiff_bnd — todo
- [ ] cvodes_rs example cvsAdvDiff_bndL — todo
- [ ] cvodes_rs example cvsAnalytic_mels — todo
- [ ] cvodes_rs example cvsDirectDemo_ls — todo
- [ ] cvodes_rs example cvsDiurnal_FSA_kry — todo
- [ ] cvodes_rs example cvsDiurnal_kry — todo
- [ ] cvodes_rs example cvsDiurnal_kry_bp — todo
- [ ] cvodes_rs example cvsFoodWeb_ASAi_kry — todo
- [ ] cvodes_rs example cvsFoodWeb_ASAp_kry — todo
- [ ] cvodes_rs example cvsHessian_ASA_FSA — todo
- [ ] cvodes_rs example cvsKrylovDemo_ls — todo
- [ ] cvodes_rs example cvsKrylovDemo_prec — todo
- [ ] cvodes_rs example cvsLotkaVolterra_ASA — todo
- [ ] cvodes_rs example cvsParticle_dns — todo
- [ ] cvodes_rs example cvsPendulum_dns — todo
- [ ] cvodes_rs example cvsRoberts_ASAi_dns — todo
- [ ] cvodes_rs example cvsRoberts_ASAi_dns_constraints — todo
- [ ] cvodes_rs example cvsRoberts_FSA_dns — todo
- [ ] cvodes_rs example cvsRoberts_FSA_dns_Switch — todo
- [ ] cvodes_rs example cvsRoberts_FSA_dns_constraints — todo
- [ ] cvodes_rs example cvsRoberts_dns — todo
- [ ] cvodes_rs example cvsRoberts_dnsL — todo
- [ ] cvodes_rs example cvsRoberts_dns_constraints — todo
- [ ] cvodes_rs example cvsRoberts_dns_uw — todo
- [ ] ida_rs example idaAnalytic_mels — todo
- [ ] ida_rs example idaFoodWeb_bnd — todo
- [ ] ida_rs example idaFoodWeb_kry — todo
- [ ] ida_rs example idaHeat2D_bnd — todo
- [ ] ida_rs example idaHeat2D_kry — todo
- [ ] ida_rs example idaKrylovDemo_ls — todo
- [ ] ida_rs example idaRoberts_dns — todo
- [ ] ida_rs example idaSlCrank_dns — todo
- [ ] idas_rs example idasAkzoNob_ASAi_dns — todo
- [ ] idas_rs example idasAkzoNob_dns — todo
- [ ] idas_rs example idasAnalytic_mels — todo
- [ ] idas_rs example idasFoodWeb_bnd — todo
- [ ] idas_rs example idasHeat2D_bnd — todo
- [ ] idas_rs example idasHeat2D_kry — todo
- [ ] idas_rs example idasHessian_ASA_FSA — todo
- [ ] idas_rs example idasKrylovDemo_ls — todo
- [ ] idas_rs example idasRoberts_ASAi_dns — todo
- [ ] idas_rs example idasRoberts_FSA_dns — todo
- [ ] idas_rs example idasRoberts_dns — todo
- [ ] idas_rs example idasSlCrank_FSA_dns — todo
- [ ] idas_rs example idasSlCrank_dns — todo
- [ ] kinsol_rs example kinAnalytic_fp — todo
- [ ] kinsol_rs example kinFerTron_dns — todo
- [ ] kinsol_rs example kinFoodWeb_kry — todo
- [ ] kinsol_rs example kinKrylovDemo_ls — todo
- [ ] kinsol_rs example kinLaplace_bnd — todo
- [ ] kinsol_rs example kinLaplace_picard_bnd — todo
- [ ] kinsol_rs example kinLaplace_picard_kry — todo
- [ ] kinsol_rs example kinRoberts_fp — todo
- [ ] kinsol_rs example kinRoboKin_dns — todo
