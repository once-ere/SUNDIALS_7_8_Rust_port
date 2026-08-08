# VERIFICATION — per-variant example matrix

One line per (example, args) reference variant parsed from the upstream
CMakeLists.txt files (199 total; tools/verify_examples.sh list regenerates
the tuple set). Status: todo | identical | last-digit(reason) |
excluded(reason) | ref-libm(reason) | OPEN(reason).

`OPEN`: the variant runs to completion but diverges from the reference in
solver-visible quantities (counters, converged values). Not yet diagnosed —
handed to the debug phase with the evidence recorded below the table. Never
"fixed" by tuning the example.

`ref-libm`: the shipped `.out` embeds the generating machine's glibc
transcendental rounding (sin/exp) inside the integration feedback loop and
cannot be byte-matched on this platform's libm; the port is verified
byte-identical to a pristine upstream-C build (reference config,
`-ffp-contract=off`) run locally. See the diurnal-family note below the
table.

| crate | example | args | reference .out | status |
|---|---|---|---|---|
| cvode_rs | cvAdvDiff_bnd | — | cvAdvDiff_bnd.out | IDENTICAL |
| cvode_rs | cvAnalytic_mels | — | cvAnalytic_mels.out | IDENTICAL |
| cvode_rs | cvDirectDemo_ls | — | cvDirectDemo_ls.out | IDENTICAL |
| cvode_rs | cvDisc_dns | — | cvDisc_dns.out | IDENTICAL |
| cvode_rs | cvDiurnal_kry_bp | — | cvDiurnal_kry_bp.out | ref-libm(glibc-2.27-era CR sin + modern exp) |
| cvode_rs | cvDiurnal_kry | — | cvDiurnal_kry.out | ref-libm(glibc>=2.28 sin/exp) |
| cvode_rs | cvKrylovDemo_ls | — | cvKrylovDemo_ls.out | ref-libm(pre-2.27 glibc CR sin/exp) |
| cvode_rs | cvKrylovDemo_ls | 1 | cvKrylovDemo_ls_1.out | ref-libm(same as no-arg variant) |
| cvode_rs | cvKrylovDemo_ls | 2 | cvKrylovDemo_ls_2.out | ref-libm(same as no-arg variant) |
| cvode_rs | cvKrylovDemo_prec | — | cvKrylovDemo_prec.out | IDENTICAL |
| cvode_rs | cvParticle_dns | — | cvParticle_dns.out | IDENTICAL |
| cvode_rs | cvPendulum_dns | — | cvPendulum_dns.out | exception: upstream .out anomaly |
| cvode_rs | cvRoberts_dns | — | cvRoberts_dns.out | IDENTICAL |
| cvode_rs | cvRoberts_dns_constraints | — | cvRoberts_dns_constraints.out | IDENTICAL |
| cvode_rs | cvRoberts_dns_negsol | — | cvRoberts_dns_negsol.out | exception: stale ref line 20 |
| cvode_rs | cvRoberts_dns_uw | — | cvRoberts_dns_uw.out | IDENTICAL |
| cvode_rs | cvRocket_dns | — | cvRocket_dns.out | IDENTICAL |
| cvode_rs | cvVdp_auto_nls | — | cvVdp_auto_nls.out | IDENTICAL |
| cvode_rs | cvKrylovDemo_ls | 0 1 | cvKrylovDemo_ls_0_1.out | ref-libm(same as no-arg variant) |
| cvode_rs | cvAdvDiff_bndL | — | cvAdvDiff_bndL.out | IDENTICAL (native band for LAPACK) |
| cvode_rs | cvRoberts_dnsL | — | cvRoberts_dnsL.out | last-digit (LAPACK->native dense) |
| cvode_rs | cvRoberts_block_klu | — | cvRoberts_block_klu.out | excluded(klu) |
| cvode_rs | cvRoberts_klu | — | cvRoberts_klu.out | excluded(klu) |
| cvode_rs | cvRoberts_sps | — | cvRoberts_sps.out | excluded(superlu) |
| cvodes_rs | cvsAdvDiff_ASAi_bnd | — | cvsAdvDiff_ASAi_bnd.out | todo |
| cvodes_rs | cvsAdvDiff_FSA_non | -sensi sim t | cvsAdvDiff_FSA_non_-sensi_sim_t.out | todo |
| cvodes_rs | cvsAdvDiff_FSA_non | -sensi stg t | cvsAdvDiff_FSA_non_-sensi_stg_t.out | todo |
| cvodes_rs | cvsAdvDiff_bnd | — | cvsAdvDiff_bnd.out | todo |
| cvodes_rs | cvsAnalytic_mels | — | cvsAnalytic_mels.out | todo |
| cvodes_rs | cvsAnalytic_mels | cvodes.max_order 3 | cvsAnalytic_mels_cvodes.max_order_3.out | todo |
| cvodes_rs | cvsDirectDemo_ls | — | cvsDirectDemo_ls.out | todo |
| cvodes_rs | cvsDiurnal_FSA_kry | -sensi sim t | cvsDiurnal_FSA_kry_-sensi_sim_t.out | todo |
| cvodes_rs | cvsDiurnal_FSA_kry | -sensi stg t | cvsDiurnal_FSA_kry_-sensi_stg_t.out | todo |
| cvodes_rs | cvsDiurnal_kry | — | cvsDiurnal_kry.out | todo |
| cvodes_rs | cvsDiurnal_kry_bp | — | cvsDiurnal_kry_bp.out | todo |
| cvodes_rs | cvsFoodWeb_ASAi_kry | — | cvsFoodWeb_ASAi_kry.out | todo |
| cvodes_rs | cvsFoodWeb_ASAp_kry | — | cvsFoodWeb_ASAp_kry.out | todo |
| cvodes_rs | cvsHessian_ASA_FSA | — | cvsHessian_ASA_FSA.out | todo |
| cvodes_rs | cvsKrylovDemo_ls | — | cvsKrylovDemo_ls.out | todo |
| cvodes_rs | cvsKrylovDemo_ls | 1 | cvsKrylovDemo_ls_1.out | todo |
| cvodes_rs | cvsKrylovDemo_ls | 2 | cvsKrylovDemo_ls_2.out | todo |
| cvodes_rs | cvsKrylovDemo_prec | — | cvsKrylovDemo_prec.out | todo |
| cvodes_rs | cvsLotkaVolterra_ASA | — | cvsLotkaVolterra_ASA.out | todo |
| cvodes_rs | cvsParticle_dns | — | cvsParticle_dns.out | todo |
| cvodes_rs | cvsPendulum_dns | — | cvsPendulum_dns.out | todo |
| cvodes_rs | cvsRoberts_ASAi_dns | — | cvsRoberts_ASAi_dns.out | todo |
| cvodes_rs | cvsRoberts_ASAi_dns_constraints | — | cvsRoberts_ASAi_dns_constraints.out | todo |
| cvodes_rs | cvsRoberts_FSA_dns | -sensi sim t | cvsRoberts_FSA_dns_-sensi_sim_t.out | todo |
| cvodes_rs | cvsRoberts_FSA_dns | -sensi stg1 t | cvsRoberts_FSA_dns_-sensi_stg1_t.out | todo |
| cvodes_rs | cvsRoberts_FSA_dns_Switch | — | cvsRoberts_FSA_dns_Switch.out | todo |
| cvodes_rs | cvsRoberts_FSA_dns_constraints | -sensi stg1 t | cvsRoberts_FSA_dns_constraints_-sensi_stg1_t.out | todo |
| cvodes_rs | cvsRoberts_dns | — | cvsRoberts_dns.out | todo |
| cvodes_rs | cvsRoberts_dns_constraints | — | cvsRoberts_dns_constraints.out | todo |
| cvodes_rs | cvsRoberts_dns_uw | — | cvsRoberts_dns_uw.out | todo |
| cvodes_rs | cvsKrylovDemo_ls | 0 1 | cvsKrylovDemo_ls_0_1.out | todo |
| cvodes_rs | cvsAdvDiff_bndL | — | cvsAdvDiff_bndL.out | todo |
| cvodes_rs | cvsRoberts_dnsL | — | cvsRoberts_dnsL.out | todo |
| cvodes_rs | cvsRoberts_ASAi_klu | — | cvsRoberts_ASAi_klu.out | excluded(klu) |
| cvodes_rs | cvsRoberts_FSA_klu | -sensi stg1 t | cvsRoberts_FSA_klu_-sensi_stg1_t.out | excluded(klu) |
| cvodes_rs | cvsRoberts_klu | — | cvsRoberts_klu.out | excluded(klu) |
| cvodes_rs | cvsRoberts_ASAi_sps | — | cvsRoberts_ASAi_sps.out | excluded(superlu) |
| cvodes_rs | cvsRoberts_FSA_sps | -sensi stg1 t | cvsRoberts_FSA_sps_-sensi_stg1_t.out | excluded(superlu) |
| cvodes_rs | cvsRoberts_sps | — | cvsRoberts_sps.out | excluded(superlu) |
| kinsol_rs | kinAnalytic_fp | — | kinAnalytic_fp.out | IDENTICAL |
| kinsol_rs | kinAnalytic_fp | --damping_fp 0.5 | kinAnalytic_fp_--damping_fp_0.5.out | IDENTICAL |
| kinsol_rs | kinAnalytic_fp | --damping_fn | kinAnalytic_fp_--damping_fn.out | IDENTICAL |
| kinsol_rs | kinAnalytic_fp | --m_aa 2 | kinAnalytic_fp_--m_aa_2.out | IDENTICAL |
| kinsol_rs | kinAnalytic_fp | --m_aa 2 --delay_aa 2 | kinAnalytic_fp_--m_aa_2_--delay_aa_2.out | IDENTICAL |
| kinsol_rs | kinAnalytic_fp | --m_aa 2 --damping_aa 0.5 | kinAnalytic_fp_--m_aa_2_--damping_aa_0.5.out | IDENTICAL |
| kinsol_rs | kinAnalytic_fp | --m_aa 2 --damping_fn | kinAnalytic_fp_--m_aa_2_--damping_fn.out | IDENTICAL |
| kinsol_rs | kinAnalytic_fp | --m_aa 3 --depth_fn | kinAnalytic_fp_--m_aa_3_--depth_fn.out | IDENTICAL |
| kinsol_rs | kinAnalytic_fp | --m_aa 2 --orth_aa 1 | kinAnalytic_fp_--m_aa_2_--orth_aa_1.out | IDENTICAL |
| kinsol_rs | kinAnalytic_fp | --m_aa 2 --orth_aa 2 | kinAnalytic_fp_--m_aa_2_--orth_aa_2.out | IDENTICAL |
| kinsol_rs | kinAnalytic_fp | --m_aa 2 --orth_aa 3 | kinAnalytic_fp_--m_aa_2_--orth_aa_3.out | IDENTICAL |
| kinsol_rs | kinFerTron_dns | — | kinFerTron_dns.out | IDENTICAL |
| kinsol_rs | kinFoodWeb_kry | — | kinFoodWeb_kry.out | IDENTICAL |
| kinsol_rs | kinKrylovDemo_ls | — | kinKrylovDemo_ls.out | IDENTICAL |
| kinsol_rs | kinLaplace_bnd | — | kinLaplace_bnd.out | IDENTICAL |
| kinsol_rs | kinLaplace_picard_bnd | — | kinLaplace_picard_bnd.out | IDENTICAL |
| kinsol_rs | kinLaplace_picard_kry | — | kinLaplace_picard_kry.out | IDENTICAL |
| kinsol_rs | kinRoberts_fp | — | kinRoberts_fp.out | IDENTICAL |
| kinsol_rs | kinRoberts_fp | kinsol.m_aa 1 | kinRoberts_fp_kinsol.m_aa_1.out | IDENTICAL |
| kinsol_rs | kinRoboKin_dns | — | kinRoboKin_dns.out | exception: stale ref (SUN_TABLE_WIDTH 28); values identical |
| kinsol_rs | kinFerTron_klu | — | kinFerTron_klu.out | excluded(klu) |
| kinsol_rs | kinRoboKin_slu | — | kinRoboKin_slu.out | excluded(superlu) |
| ida_rs | idaAnalytic_mels | — | idaAnalytic_mels.out | todo |
| ida_rs | idaAnalytic_mels | ida.scalar_tolerances 1e-3 1e-8 | idaAnalytic_mels_ida.scalar_tolerances_1e-3_1e-8.out | todo |
| ida_rs | idaFoodWeb_bnd | — | idaFoodWeb_bnd.out | todo |
| ida_rs | idaFoodWeb_kry | — | idaFoodWeb_kry.out | todo |
| ida_rs | idaHeat2D_bnd | — | idaHeat2D_bnd.out | todo |
| ida_rs | idaHeat2D_kry | — | idaHeat2D_kry.out | todo |
| ida_rs | idaKrylovDemo_ls | — | idaKrylovDemo_ls.out | todo |
| ida_rs | idaKrylovDemo_ls | 1 | idaKrylovDemo_ls_1.out | todo |
| ida_rs | idaKrylovDemo_ls | 2 | idaKrylovDemo_ls_2.out | todo |
| ida_rs | idaRoberts_dns | — | idaRoberts_dns.out | todo |
| ida_rs | idaSlCrank_dns | — | idaSlCrank_dns.out | todo |
| ida_rs | idaHeat2D_klu | — | idaHeat2D_klu.out | excluded(klu) |
| ida_rs | idaRoberts_klu | — | idaRoberts_klu.out | excluded(klu) |
| ida_rs | idaRoberts_sps | — | idaRoberts_sps.out | excluded(superlu) |
| idas_rs | idasAkzoNob_ASAi_dns | — | idasAkzoNob_ASAi_dns.out | exception: ref trailing whitespace stripped; values identical |
| idas_rs | idasAkzoNob_dns | — | idasAkzoNob_dns.out | IDENTICAL |
| idas_rs | idasAnalytic_mels | — | idasAnalytic_mels.out | IDENTICAL |
| idas_rs | idasAnalytic_mels | idas.init_step 1e-5 | idasAnalytic_mels_idas.init_step_1e-5.out | IDENTICAL |
| idas_rs | idasFoodWeb_bnd | — | idasFoodWeb_bnd.out | OPEN(last-digit: hused col, t=0.7/1.0) |
| idas_rs | idasHeat2D_bnd | — | idasHeat2D_bnd.out | IDENTICAL |
| idas_rs | idasHeat2D_kry | — | idasHeat2D_kry.out | IDENTICAL |
| idas_rs | idasHessian_ASA_FSA | — | idasHessian_ASA_FSA.out | IDENTICAL |
| idas_rs | idasKrylovDemo_ls | — | idasKrylovDemo_ls.out | IDENTICAL |
| idas_rs | idasKrylovDemo_ls | 1 | idasKrylovDemo_ls_1.out | IDENTICAL |
| idas_rs | idasKrylovDemo_ls | 2 | idasKrylovDemo_ls_2.out | IDENTICAL |
| idas_rs | idasRoberts_ASAi_dns | — | idasRoberts_ASAi_dns.out | IDENTICAL |
| idas_rs | idasRoberts_FSA_dns | -sensi stg t | idasRoberts_FSA_dns_-sensi_stg_t.out | IDENTICAL |
| idas_rs | idasRoberts_dns | — | idasRoberts_dns.out | IDENTICAL |
| idas_rs | idasSlCrank_dns | — | idasSlCrank_dns.out | OPEN(nre/nni off by 1, G @ digit 11) |
| idas_rs | idasSlCrank_FSA_dns | — | idasSlCrank_FSA_dns.out | OPEN(nst 263 vs 233; dG/dp digit 5) |
| idas_rs | idasRoberts_ASAi_klu | — | idasRoberts_ASAi_klu.out | excluded(klu) |
| idas_rs | idasRoberts_FSA_klu | -sensi stg t | idasRoberts_FSA_klu_-sensi_stg_t.out | excluded(klu) |
| idas_rs | idasRoberts_klu | — | idasRoberts_klu.out | excluded(klu) |
| idas_rs | idasRoberts_ASAi_sps | — | idasRoberts_ASAi_sps.out | excluded(superlu) |
| idas_rs | idasRoberts_FSA_sps | -sensi stg t | idasRoberts_FSA_sps_-sensi_stg_t.out | excluded(superlu) |
| idas_rs | idasRoberts_sps | — | idasRoberts_sps.out | excluded(superlu) |
| arkode_rs | ark_analytic | — | ark_analytic.out | todo |
| arkode_rs | ark_analytic | arkode.scalar_tolerances 1e-6 1e-8 arkode.table_names ARKODE_ESDIRK547L2SA_7_4_5 ARKODE_ERK_NONE | ark_analytic_arkode.scalar_tolerances_1e-6_1e-8_arkode.table_names_ARKODE_ESDIRK547L2SA_7_4_5_ARKODE_ERK_NONE.out | todo |
| arkode_rs | ark_advection_diffusion_reaction_splitting | — | ark_advection_diffusion_reaction_splitting.out | todo |
| arkode_rs | ark_analytic_lsrk | — | ark_analytic_lsrk.out | todo |
| arkode_rs | ark_analytic_lsrk_varjac | — | ark_analytic_lsrk_varjac.out | todo |
| arkode_rs | ark_analytic_lsrk_domeigest | — | ark_analytic_lsrk_domeigest.out | todo |
| arkode_rs | ark_analytic_lsrk_domeigest | arkid.dom_eig_est_init_preprocess_iters 1 sundomeigestimator.max_iters 1 | ark_analytic_lsrk_domeigest_arkid.dom_eig_est_init_preprocess_iters_1_sundomeigestimator.max_iters_1.out | todo |
| arkode_rs | ark_analytic_mels | — | ark_analytic_mels.out | todo |
| arkode_rs | ark_analytic_nonlin | — | ark_analytic_nonlin.out | todo |
| arkode_rs | ark_analytic_partitioned | forcing | ark_analytic_partitioned_forcing.out | todo |
| arkode_rs | ark_analytic_partitioned | splitting | ark_analytic_partitioned_splitting.out | todo |
| arkode_rs | ark_analytic_partitioned | splitting ARKODE_SPLITTING_BEST_2_2_2 | ark_analytic_partitioned_splitting_ARKODE_SPLITTING_BEST_2_2_2.out | todo |
| arkode_rs | ark_analytic_partitioned | splitting ARKODE_SPLITTING_RUTH_3_3_2 | ark_analytic_partitioned_splitting_ARKODE_SPLITTING_RUTH_3_3_2.out | todo |
| arkode_rs | ark_analytic_partitioned | splitting ARKODE_SPLITTING_YOSHIDA_8_6_2 | ark_analytic_partitioned_splitting_ARKODE_SPLITTING_YOSHIDA_8_6_2.out | todo |
| arkode_rs | ark_analytic_ssprk | — | ark_analytic_ssprk.out | todo |
| arkode_rs | ark_brusselator_1D_mri | — | ark_brusselator_1D_mri.out | todo |
| arkode_rs | ark_brusselator_fp | — | ark_brusselator_fp.out | todo |
| arkode_rs | ark_brusselator_lsrk_domeigest | — | ark_brusselator_lsrk_domeigest.out | todo |
| arkode_rs | ark_brusselator_lsrk_externaldomeigest | — | ark_brusselator_lsrk_externaldomeigest.out | todo |
| arkode_rs | ark_brusselator_mri | — | ark_brusselator_mri.out | todo |
| arkode_rs | ark_brusselator | — | ark_brusselator.out | todo |
| arkode_rs | ark_brusselator1D_imexmri | 0 0.001 | ark_brusselator1D_imexmri_0_0.001.out | todo |
| arkode_rs | ark_brusselator1D_imexmri | 2 0.001 | ark_brusselator1D_imexmri_2_0.001.out | todo |
| arkode_rs | ark_brusselator1D_imexmri | 3 0.001 | ark_brusselator1D_imexmri_3_0.001.out | todo |
| arkode_rs | ark_brusselator1D_imexmri | 4 0.001 | ark_brusselator1D_imexmri_4_0.001.out | todo |
| arkode_rs | ark_brusselator1D_imexmri | 5 0.001 | ark_brusselator1D_imexmri_5_0.001.out | todo |
| arkode_rs | ark_brusselator1D_imexmri | 6 0.001 | ark_brusselator1D_imexmri_6_0.001.out | todo |
| arkode_rs | ark_brusselator1D_imexmri | 7 0.001 | ark_brusselator1D_imexmri_7_0.001.out | todo |
| arkode_rs | ark_brusselator1D | — | ark_brusselator1D.out | todo |
| arkode_rs | ark_conserved_exp_entropy_ark | 1 0 | ark_conserved_exp_entropy_ark_1_0.out | todo |
| arkode_rs | ark_conserved_exp_entropy_ark | 1 1 | ark_conserved_exp_entropy_ark_1_1.out | todo |
| arkode_rs | ark_conserved_exp_entropy_erk | 1 | ark_conserved_exp_entropy_erk_1.out | todo |
| arkode_rs | ark_damped_harmonic_symplectic | — | ark_damped_harmonic_symplectic.out | todo |
| arkode_rs | ark_dissipated_exp_entropy | 1 0 | ark_dissipated_exp_entropy_1_0.out | todo |
| arkode_rs | ark_dissipated_exp_entropy | 1 1 | ark_dissipated_exp_entropy_1_1.out | todo |
| arkode_rs | ark_harmonic_symplectic | — | ark_harmonic_symplectic.out | todo |
| arkode_rs | ark_heat1D_adapt | — | ark_heat1D_adapt.out | todo |
| arkode_rs | ark_heat1D | — | ark_heat1D.out | todo |
| arkode_rs | ark_kepler | --stepper ERK --step-mode adapt | ark_kepler_--stepper_ERK_--step-mode_adapt.out | todo |
| arkode_rs | ark_kepler | --stepper ERK --step-mode fixed --count-orbits | ark_kepler_--stepper_ERK_--step-mode_fixed_--count-orbits.out | todo |
| arkode_rs | ark_kepler | --stepper SPRK --step-mode fixed --count-orbits --use-compensated-sums | ark_kepler_--stepper_SPRK_--step-mode_fixed_--count-orbits_--use-compensated-sums.out | todo |
| arkode_rs | ark_kepler | --stepper SPRK --step-mode fixed --method ARKODE_SPRK_EULER_1_1 --tf 50 --check-order --nout 1 | ark_kepler_--stepper_SPRK_--step-mode_fixed_--method_ARKODE_SPRK_EULER_1_1_--tf_50_--check-order_--nout_1.out | todo |
| arkode_rs | ark_kepler | --stepper SPRK --step-mode fixed --method ARKODE_SPRK_LEAPFROG_2_2 --tf 50 --check-order --nout 1 | ark_kepler_--stepper_SPRK_--step-mode_fixed_--method_ARKODE_SPRK_LEAPFROG_2_2_--tf_50_--check-order_--nout_1.out | todo |
| arkode_rs | ark_kepler | --stepper SPRK --step-mode fixed --method ARKODE_SPRK_MCLACHLAN_2_2 --tf 50 --check-order --nout 1 | ark_kepler_--stepper_SPRK_--step-mode_fixed_--method_ARKODE_SPRK_MCLACHLAN_2_2_--tf_50_--check-order_--nout_1.out | todo |
| arkode_rs | ark_kepler | --stepper SPRK --step-mode fixed --method ARKODE_SPRK_MCLACHLAN_3_3 --tf 50 --check-order --nout 1 | ark_kepler_--stepper_SPRK_--step-mode_fixed_--method_ARKODE_SPRK_MCLACHLAN_3_3_--tf_50_--check-order_--nout_1.out | todo |
| arkode_rs | ark_kepler | --stepper SPRK --step-mode fixed --method ARKODE_SPRK_MCLACHLAN_4_4 --tf 50 --check-order --nout 1 | ark_kepler_--stepper_SPRK_--step-mode_fixed_--method_ARKODE_SPRK_MCLACHLAN_4_4_--tf_50_--check-order_--nout_1.out | todo |
| arkode_rs | ark_kepler | --stepper SPRK --step-mode fixed --method ARKODE_SPRK_MCLACHLAN_5_6 --tf 50 --check-order --nout 1 | ark_kepler_--stepper_SPRK_--step-mode_fixed_--method_ARKODE_SPRK_MCLACHLAN_5_6_--tf_50_--check-order_--nout_1.out | todo |
| arkode_rs | ark_kepler | --stepper SPRK --step-mode fixed --method ARKODE_SPRK_PSEUDO_LEAPFROG_2_2 --tf 50 --check-order --nout 1 | ark_kepler_--stepper_SPRK_--step-mode_fixed_--method_ARKODE_SPRK_PSEUDO_LEAPFROG_2_2_--tf_50_--check-order_--nout_1.out | todo |
| arkode_rs | ark_kepler | --stepper SPRK --step-mode fixed --method ARKODE_SPRK_RUTH_3_3 --tf 50 --check-order --nout 1 | ark_kepler_--stepper_SPRK_--step-mode_fixed_--method_ARKODE_SPRK_RUTH_3_3_--tf_50_--check-order_--nout_1.out | todo |
| arkode_rs | ark_kepler | --stepper SPRK --step-mode fixed --method ARKODE_SPRK_YOSHIDA_6_8 --tf 50 --check-order --nout 1 | ark_kepler_--stepper_SPRK_--step-mode_fixed_--method_ARKODE_SPRK_YOSHIDA_6_8_--tf_50_--check-order_--nout_1.out | todo |
| arkode_rs | ark_kepler | — | ark_kepler.out | todo |
| arkode_rs | ark_kpr_mri | 0 1 0.005 | ark_kpr_mri_0_1_0.005.out | todo |
| arkode_rs | ark_kpr_mri | 1 0 0.01 | ark_kpr_mri_1_0_0.01.out | todo |
| arkode_rs | ark_kpr_mri | 1 1 0.002 | ark_kpr_mri_1_1_0.002.out | todo |
| arkode_rs | ark_kpr_mri | 2 4 0.002 | ark_kpr_mri_2_4_0.002.out | todo |
| arkode_rs | ark_kpr_mri | 3 2 0.001 | ark_kpr_mri_3_2_0.001.out | todo |
| arkode_rs | ark_kpr_mri | 4 3 0.001 | ark_kpr_mri_4_3_0.001.out | todo |
| arkode_rs | ark_kpr_mri | 5 4 0.001 | ark_kpr_mri_5_4_0.001.out | todo |
| arkode_rs | ark_kpr_mri | 6 5 0.001 | ark_kpr_mri_6_5_0.001.out | todo |
| arkode_rs | ark_kpr_mri | 7 2 0.002 | ark_kpr_mri_7_2_0.002.out | todo |
| arkode_rs | ark_kpr_mri | 8 3 0.001 -100 100 0.5 1 | ark_kpr_mri_8_3_0.001_-100_100_0.5_1.out | todo |
| arkode_rs | ark_kpr_mri | 9 3 0.001 -100 100 0.5 1 | ark_kpr_mri_9_3_0.001_-100_100_0.5_1.out | todo |
| arkode_rs | ark_kpr_mri | 10 4 0.001 -100 100 0.5 1 | ark_kpr_mri_10_4_0.001_-100_100_0.5_1.out | todo |
| arkode_rs | ark_kpr_mri | 11 2 0.001 | ark_kpr_mri_11_2_0.001.out | todo |
| arkode_rs | ark_kpr_mri | 12 3 0.005 | ark_kpr_mri_12_3_0.005.out | todo |
| arkode_rs | ark_kpr_mri | 13 4 0.01 | ark_kpr_mri_13_4_0.01.out | todo |
| arkode_rs | ark_KrylovDemo_prec | — | ark_KrylovDemo_prec.out | todo |
| arkode_rs | ark_KrylovDemo_prec | 1 | ark_KrylovDemo_prec_1.out | todo |
| arkode_rs | ark_KrylovDemo_prec | 2 | ark_KrylovDemo_prec_2.out | todo |
| arkode_rs | ark_lotka_volterra_ASA | --check-freq 1 | ark_lotka_volterra_ASA_--check-freq_1.out | todo |
| arkode_rs | ark_lotka_volterra_ASA | --check-freq 5 | ark_lotka_volterra_ASA_--check-freq_5.out | todo |
| arkode_rs | ark_onewaycouple_mri | — | ark_onewaycouple_mri.out | todo |
| arkode_rs | ark_reaction_diffusion_mri | — | ark_reaction_diffusion_mri.out | todo |
| arkode_rs | ark_robertson_constraints | — | ark_robertson_constraints.out | todo |
| arkode_rs | ark_robertson_root | — | ark_robertson_root.out | todo |
| arkode_rs | ark_robertson | — | ark_robertson.out | todo |
| arkode_rs | ark_twowaycouple_mri | — | ark_twowaycouple_mri.out | todo |
| arkode_rs | ark_brusselator_fp | 1 | ark_brusselator_fp_1.out | todo |

## Documented exceptions

- **cvPendulum_dns**: the upstream reference `cvPendulum_dns.out` prints
  `atol = 1.00e-5` (single-digit exponent) on its 5 header lines while the
  C source formats both tolerances with `%8.2e` — no conforming C `printf`
  produces a one-digit exponent, so the shipped reference cannot be
  reproduced by its own source. The port prints `1.00e-05` (what `%8.2e`
  yields); those 5 lines (10 diff lines) are the only divergence accepted
  for this variant once the remainder verifies.
- **cvRoberts_dnsL**: LAPACK dense solver replaced by the native dense
  solver per the port plan; different factorization arithmetic gives
  last-digit drift in printed `y` values (§1 documented exception class).
- **cvRoberts_dns_negsol**: reference line 20 (`netf = 59     ncfn`, 5-space
  gap) is unproducible by the example's single `%-6ld` format string, which
  yields the 8-space gap seen on the same run's line 41 — the shipped line
  predates the current PrintFinalStats format. Port output matches the
  current C source; the 2-diff-line exception is accepted.
- **Deterministic `pow`** (not an exception — a fix): `SUNRpowerR` ports the
  ARM optimized-routines `pow` (via musl, MIT; the glibc >= 2.28 algorithm
  that generated the references) instead of calling platform libm — Apple
  libm `pow` is 1 ulp off glibc on rare arguments inside the step-size
  heuristics, which forked `cvDirectDemo_ls`, `cvParticle_dns`, and
  `cvVdp_auto_nls` before the port. All three are byte-IDENTICAL with it.
- **kinRoboKin_dns**: the only kinsol variant that calls `KINPrintAllStats`
  in `SUN_OUTPUTFORMAT_TABLE`. All 16 stat lines differ by exactly one space
  before the `=`: the shipped `.out` puts `=` at column 29 (name field padded
  to 28), while `src/sundials/sundials_utils.h:31` defines
  `SUN_TABLE_WIDTH 29` and `sunfprintf_long` formats `"%-*s = %ld\n"`, giving
  column 30. Every printed value is byte-identical (verified: 0
  non-whitespace diffs over all 46 lines). The shipped 7.8.0 reference tree
  is self-inconsistent on this point — `ark_kepler.out` has `=` at column 29
  while `ark_kepler_--stepper_ERK_--step-mode_adapt.out` has it at column 30
  for the same `Current time` field, and every cvode/cvodes/ida/idas
  reference uses column 30 — so a subset of `.out` files predates the
  `SUN_TABLE_WIDTH` 28 -> 29 change. The port follows the shipped header;
  matching the stale reference would require contradicting it.
- **idasAkzoNob_ASAi_dns**: 3 diff lines, zero value differences. (a) The
  `G:` line: the C source is
  `printf("G:          %24.16f \n", Ith(q, 1));` — note the space before
  `\n` — so the port emits a trailing space and the reference does not.
  (b) The file ends with
  `printf("------...------\n\n")`, so the port emits a final blank line the
  reference lacks. The reference has been trailing-whitespace-normalized;
  the stripping is NOT systematic — the sibling `idasAkzoNob_dns.c` has the
  byte-identical `G:` printf and its shipped `idasAkzoNob_dns.out` line 37
  DOES retain the trailing space (that variant is IDENTICAL). Port output
  matches the C source character-for-character.

## OPEN divergences handed to the debug phase (2026-08-07)

Three idas_rs variants run to completion but do not match. None is a
formatting or example-setup issue (every one of them is byte-identical up
to the first solver-visible quantity), so per §6 they are recorded, not
guessed at.

- **idasFoodWeb_bnd** — 4 diff lines. The entire trajectory table matches:
  all `c_bl`/`c_tr` species values, `nst = 239` and order `k = 1` are
  byte-identical at every output time, and `hused` matches at
  t = 1e-8 … 4e-1. Only the last column (`hused`, `IDAGetLastStep`, `%12.4e`)
  differs at the final two output times: ref `6.2655e-01` vs port
  `6.2656e-01` at t = 7.0e-1 and t = 1.0e+0. Relative delta ~1.6e-5, i.e. a
  value sitting on the `%.4e` rounding boundary (~0.626555); a 1-ulp
  difference in the last accepted step size flips the printed digit. Not a
  formatting bug: both C `printf` and the port's `fmt_e` round the exact
  binary double correctly, so identical doubles print identically. Suspect
  the step-size heuristic arithmetic on the final step; `nst` never diverges.
- **idasSlCrank_dns** — 6 diff lines, no sensitivities involved. All 26
  trajectory rows (q, dq, lambda, nst, k, h) are byte-identical through
  t = 10.0, and `nst = 251`, `nje = 39`, `netf = 1`, `ncfn = 20`,
  `nsf = 0` all match. Only `nre` (1066 -> 1065) and `nni` (676 -> 675) are
  each low by exactly one, and the quadrature `G` differs at digit 11
  (3.3366160662909388 vs 3.3366160663381925, rel. 1.4e-11). Signature: one
  nonlinear iteration fewer somewhere that does not perturb the step
  sequence — most likely in `IDACalcIC` or the first step. Cross-check
  against `ida_rs`'s `idaSlCrank_dns` was not possible: the ida_rs examples
  do not currently compile (unrelated, in-flight phase).
- **idasSlCrank_FSA_dns** — 22 diff lines, the largest divergence.
  `nst` 233 -> 263, `nre` 1180 -> 1203, `nje` 46 -> 44, `nni` 720 -> 763,
  `ncfn` 26 -> 23, `nsf` 1 -> 2; `G` differs at digit 8 and the four
  `dG/dp` blocks differ in the 5th significant digit of the second
  component (-3.6375e-01 vs -3.6376e-01 / -3.6373e-01). This example calls
  `IDASensInit(..., fS = None, ...)`, i.e. the INTERNAL difference-quotient
  sensitivity residual, and reaches the live `ida_p` from its user data via
  a `Weak<RefCell<IDAMemRec>>` handle rather than the shared-handle
  `SensParams` pattern that ARCHITECTURE §8 fixes for CVODES. Prime suspect
  is therefore the IDAS parameter-aliasing path (`IDASetSensParams` copies
  into `ida_mem.ida_p` instead of sharing the caller's array); the debug
  phase should decide whether to give IDAS the same `SensParams` contract
  CVODES has. Note the no-sensitivity `idasSlCrank_dns` above also diverges
  slightly, so part of this may be the same underlying base-integrator
  issue rather than sensitivity-specific.

## Diurnal-family reference-libm exception (2026-08-06)

The six cvode_rs variants marked `ref-libm` (cvDiurnal_kry, cvDiurnal_kry_bp,
cvKrylovDemo_ls x4) all solve the 2-species diurnal problem, whose RHS/Jtimes
evaluate `sin`/`exp` inside the integration feedback loop. Evidence that the
mismatch is the reference environment, not the port:

1. **Port == upstream C.** Pristine upstream 7.8.0 C sources compiled locally
   (clang and gcc, `-O3 -DNDEBUG -ffp-contract=off`, logging 2, monitoring on,
   profiling off, error checks off) produce output byte-identical to the Rust
   port for all six variants — including the same divergence from the shipped
   `.out` (e.g. cvDiurnal_kry t=2.88e4: both give nst=311/order 3 vs shipped
   nst=307/order 4).
2. **Shipped `.out`s reproduced by libm substitution.** Linking the same
   pristine C build against the reference platform's libm implementations
   reproduces each shipped `.out` byte-for-byte:
   - cvDiurnal_kry.out (regenerated 2024-09-10, LLNL commit bb6cf3e7): glibc
     dbl-64 `sin`+`exp` (glibc >= 2.28 era, IBM s_sin.c + Nagy e_exp.c).
   - cvDiurnal_kry_bp.out (same commit, different CI node): correctly-rounded
     `sin` (pre-2.28 glibc IBM sin with mp fallback) + modern (>= 2.27) glibc
     `exp` — the glibc 2.27 signature.
   - cvKrylovDemo_ls*.out (regenerated 2020-05-19, commit 56289b71): fully
     correctly-rounded `sin`/`exp` (pre-2.27 glibc, e.g. RHEL7 2.17) — all
     four argv variants match byte-for-byte.
3. **Mutual inconsistency.** The three `.out`s require three different `sin`
   implementations (glibc-2.28+ for kry, <= 2.27 CR for bp/ls), so no single
   libm — and therefore no faithful port — can byte-match all of them
   simultaneously. First consequential deviation for cvDiurnal_kry: Apple
   libm and glibc >= 2.28 sin differ by 1 ulp at x = 0x1.27b7ca8e314fp-3
   (om*t, t ~= 1986 s); 47 of ~7600 sin calls differ over the run; the
   ulp-level trajectory drift first flips a step-size/order decision between
   nst 277 (t = 2.16e4, still byte-identical) and the t = 2.88e4 checkpoint.

Acceptance for these six variants is therefore byte-identity against the
locally-built pristine upstream C binary (satisfied), not the shipped `.out`.
