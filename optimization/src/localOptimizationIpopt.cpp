#include "localOptimizationIpopt.h"
#include "localOptimizationIpoptTNLP.h"
#include "logging.h"
#include "optimizationOptions.h"
#include "optimizationProblem.h"
#include "parpeException.h"
#include <alloca.h>
#include <cassert>
#include <cmath>
#include <csignal>
#include <cstdbool>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>

namespace parpe {


// https://www.coin-or.org/Ipopt/documentation/node40.html
// grep -A1 -r "roptions->Add" ../ThirdParty/Ipopt-3.12.7
const std::vector<std::string> strOpts = {
    "accept_every_trial_step",
    "adaptive_mu_globalization",
    "adaptive_mu_kkt_norm_type",
    "adaptive_mu_restore_previous_iterate",
    "alpha_for_y",
    "check_derivatives_for_naninf",
    "constraint_violation_norm_type",
    "corrector_type",
    "dependency_detection_with_rhs",
    "dependency_detector",
    "derivative_test",
    "derivative_test_print_all",
    "evaluate_orig_obj_at_resto_trial",
    "expect_infeasible_problem",
    "fast_step_computation",
    "fixed_mu_oracle",
    "fixed_variable_treatment",
    "flexible_penalty_function",
    "hessian_approximation",
    "hessian_approximation_space",
    "hessian_constant",
    "honor_original_bounds",
    "inexact_algorithm",
    "inexact_linear_system_scaling",
    "inexact_step_decomposition",
    "inf_pr_output",
    "jac_c_constant",
    "jac_d_constant",
    "jacobian_approximation",
    "limited_memory_aug_solver",
    "limited_memory_initialization",
    "limited_memory_special_for_resto",
    "limited_memory_update_type",
    "linear_scaling_on_demand",
    "linear_solver",
    "linear_system_scaling",
    "line_search_method",
    "ma27_ignore_singularity",
    "ma27_skip_inertia_check",
    "ma57_automatic_scaling",
    "ma77_order",
    "ma86_order",
    "ma86_scaling",
    "ma97_dump_matrix",
    "ma97_order",
    "ma97_scaling",
    "ma97_scaling1",
    "ma97_scaling2",
    "ma97_scaling3",
    "ma97_solve_blas3",
    "ma97_switch1",
    "ma97_switch2",
    "ma97_switch3",
    "magic_steps",
    "mehrotra_algorithm",
    "modify_hessian_with_slacks",
    "mu_allow_fast_monotone_decrease",
    "mu_oracle",
    "mu_strategy",
    "neg_curv_test_reg",
    "nlp_scaling_method",
    "option_file_name",
    "output_file",
    "pardiso_iterative",
    "pardiso_matching_strategy",
    "pardiso_order",
    "pardiso_redo_symbolic_fact_only_if_inertia_wrong",
    "pardiso_repeated_perturbation_means_singular",
    "pardiso_skip_inertia_check",
    "perturb_always_cd",
    "print_options_documentation",
    "print_options_latex_mode",
    "print_timing_statistics",
    "print_user_options",
    "quality_function_balancing_term",
    "quality_function_centrality",
    "quality_function_norm_type",
    "recalc_y",
    "replace_bounds",
    "sb","print_info_string",
    "skip_corr_if_neg_curv",
    "skip_corr_in_monotone_mode",
    "skip_finalize_solution_call",
    "start_with_resto",
    "suppress_all_output",
    "warm_start_entire_iterate",
    "warm_start_same_structure",
    "wsmp_iterative",
    "wsmp_no_pivoting",
    "wsmp_skip_inertia_check",
};

const std::vector<std::string> intOpts = {
    "acceptable_iter",
    "accept_after_max_steps",
    "accept_every_trial_step",
    "adaptive_mu_kkterror_red_iters",
    "debug_print_level",
    "derivative_test_first_index",
    "file_print_level",
    "filter_reset_trigger",
    "inexact_desired_pd_residual_iter",
    "inexact_normal_max_iter",
    "inexact_regularization_ls_count_trigger",
    "limited_memory_max_history",
    "limited_memory_max_skipping",
    "ma57_block_size",
    "ma57_node_amalgamation",
    "ma57_pivot_order",
    "ma57_small_pivot_flag",
    "ma77_buffer_lpage",
    "ma77_buffer_npage",
    "ma77_file_size",
    "ma77_maxstore",
    "ma77_nemin",
    "ma77_print_level",
    "ma86_nemin",
    "ma86_print_level",
    "ma97_nemin",
    "ma97_print_level",
    "max_filter_resets",
    "max_iter",
    "max_refinement_steps",
    "max_resto_iter",
    "max_soc",
    "max_soft_resto_iters",
    "min_refinement_steps",
    "mumps_mem_percent",
    "mumps_permuting_scaling",
    "mumps_pivot_order",
    "mumps_scaling",
    "num_linear_variables",
    "pardiso_iter_coarse_size",
    "pardiso_iter_max_levels",
    "pardiso_iter_max_row_fill",
    "pardiso_max_droptol_corrections",
    "pardiso_max_iter",
    "pardiso_max_iterative_refinement_steps",
    "pardiso_msglvl",
    "print_frequency_iter",
    "print_frequency_time",
    "print_level",
    "quality_function_max_section_steps",
    "soc_method",
    "watchdog_shortened_iter_trigger",
    "watchdog_trial_iter_max",
    "wsmp_max_iter",
    "wsmp_num_threads",
    "wsmp_ordering_option",
    "wsmp_ordering_option2",
    "wsmp_scaling",
    "wsmp_write_matrix_iteration",
};

const std::vector<std::string> dblOpts = {
    "acceptable_compl_inf_tol",
    "acceptable_constr_viol_tol",
    "acceptable_dual_inf_tol",
    "acceptable_obj_change_tol",
    "acceptable_tol",
    "adaptive_mu_kkterror_red_fact",
    "adaptive_mu_monotone_init_factor",
    "adaptive_mu_safeguard_factor",
    "alpha_for_y_tol",
    "alpha_min_frac",
    "alpha_red_factor",
    "barrier_tol_factor",
    "bound_mult_reset_threshold",
    "bound_relax_factor",
    "compl_inf_tol",
    "constr_mult_reset_threshold",
    "constr_viol_tol",
    "corrector_compl_avrg_red_fact",
    "delta", "Multiplier for constraint violation in the switching rule.",
    "derivative_test_perturbation",
    "derivative_test_tol",
    "diverging_iterates_tol",
    "dual_inf_tol",
    "eta_phi",
    "expect_infeasible_problem_ctol",
    "expect_infeasible_problem_ytol",
    "filter_margin_fact",
    "filter_max_margin",
    "findiff_perturbation",
    "first_hessian_perturbation",
    "gamma_phi",
    "gamma_theta",
    "inexact_decomposition_activate_tol",
    "inexact_decomposition_inactivate_tol",
    "inexact_desired_pd_residual",
    "inexact_normal_tol",
    "jacobian_regularization_exponent",
    "jacobian_regularization_value",
    "kappa_d",
    "kappa_sigma",
    "kappa_soc",
    "limited_memory_init_val",
    "limited_memory_init_val_max",
    "limited_memory_init_val_min",
    "local_inf_Ac_tol",
    "ma27_la_init_factor",
    "ma27_liw_init_factor",
    "ma27_meminc_factor",
    "ma27_pivtol",
    "ma27_pivtolmax",
    "ma28_pivtol",
    "ma57_pivtol",
    "ma57_pivtolmax",
    "ma57_pre_alloc",
    "ma77_small",
    "ma77_static",
    "ma77_u",
    "ma77_umax",
    "ma86_small",
    "ma86_static",
    "ma86_u",
    "ma86_umax",
    "ma97_small",
    "ma97_u",
    "ma97_umax",
    "max_cpu_time",
    "max_hessian_perturbation",
    "min_hessian_perturbation",
    "mu_init", "Initial value for the barrier parameter.",
    "mu_linear_decrease_factor",
    "mu_max",
    "mu_max_fact",
    "mu_min",
    "mumps_dep_tol",
    "mumps_pivtol",
    "mumps_pivtolmax",
    "mu_superlinear_decrease_power",
    "mu_target",
    "neg_curv_test_tol",
    "nlp_lower_bound_inf",
    "nlp_scaling_constr_target_gradient",
    "nlp_scaling_max_gradient", "Maximum gradient after NLP scaling.",
    "nlp_scaling_min_value",
    "nlp_scaling_obj_target_gradient",
    "nlp_upper_bound_inf",
    "nu_inc",
    "nu_init",
    "nu_low_fact",
    "nu_low_init",
    "nu_update_inf_skip_tol",
    "obj_max_inc",
    "obj_scaling_factor",
    "pardiso_iter_dropping_factor",
    "pardiso_iter_dropping_schur",
    "pardiso_iter_inverse_norm_factor",
    "pardiso_iter_relative_tol",
    "perturb_dec_fact",
    "perturb_inc_fact",
    "perturb_inc_fact_first",
    "point_perturbation_radius",
    "print_frequency_time",
    "quality_function_section_qf_tol",
    "quality_function_section_sigma_tol",
    "recalc_y_feas_tol",
    "required_infeasibility_reduction",
    "residual_improvement_factor",
    "residual_ratio_max",
    "residual_ratio_singular",
    "resto_failure_feasibility_threshold",
    "resto_penalty_parameter",
    "resto_proximity_weight",
    "rho",
    "sigma_max",
    "sigma_min",
    "slack_move",
    "slack_scale_max",
    "s_max",
    "soft_resto_pderror_reduction_factor",
    "s_phi",
    "s_theta",
    "tau_min",
    "tcc_psi",
    "tcc_theta",
    "tcc_theta_mu_exponent",
    "tcc_zeta",
    "theta_max_fact",
    "theta_min_fact",
    "tiny_step_tol",
    "tiny_step_y_tol",
    "tol",
    "tt_eps2",
    "tt_eps3",
    "tt_kappa1",
    "tt_kappa2",
    "warm_start_bound_frac",
    "warm_start_bound_push",
    "warm_start_mult_bound_push",
    "warm_start_mult_init_max",
    "warm_start_slack_bound_frac",
    "warm_start_slack_bound_push",
    "warm_start_target_mu",
    "wsmp_inexact_droptol",
    "wsmp_inexact_fillin_limit",
    "wsmp_pivtol",
    "wsmp_pivtolmax",
    "wsmp_singularity_threshold",
};

#ifdef INSTALL_SIGNAL_HANDLER
extern volatile sig_atomic_t caughtTerminationSignal;
#endif

static_assert(sizeof(double) == sizeof(Number),
              "Sizeof IpOpt::Number != sizeof double");


using namespace Ipopt;


void setIpOptOption(const std::pair<const std::string, const std::string> &pair, SmartPtr<OptionsList>* o) {
    // for iterating over OptimizationOptions

    auto options = *o;

    const std::string &key = pair.first;
    const std::string &val = pair.second;

    bool success = true;
    if(std::find(strOpts.begin(), strOpts.end(), key) != strOpts.end())
        success = options->SetStringValue(key, val);
    else if(std::find(intOpts.begin(), intOpts.end(), key) != intOpts.end())
        success = options->SetIntegerValue(key, std::stoi(val));
    else if(std::find(dblOpts.begin(), dblOpts.end(), key) != dblOpts.end())
        success = options->SetNumericValue(key, std::stod(val));
    else {
        logmessage(LOGLVL_WARNING, "Ignoring unknown optimization option %s.", key.c_str());
        return;
    }

    RELEASE_ASSERT(success, "Problem setting IpOpt option");

    logmessage(LOGLVL_DEBUG, "Set optimization option %s to %s.", key.c_str(), val.c_str());
}

void setIpOptOptions(SmartPtr<OptionsList> optionsIpOpt,
                     OptimizationProblem *problem) {

    if (problem->getOptimizationOptions().printToStdout) {
        optionsIpOpt->SetIntegerValue("print_level", 5);
        optionsIpOpt->SetStringValue("print_user_options", "yes");
    } else {
        optionsIpOpt->SetIntegerValue("print_level", 0);
        optionsIpOpt->SetStringValue("print_user_options", "no");
        optionsIpOpt->SetStringValue("sb",
                                     "yes"); // suppress copyright message
    }

    // TODO: move to HDF5 file
    optionsIpOpt->SetStringValue("hessian_approximation", "limited-memory");
    optionsIpOpt->SetStringValue("limited_memory_update_type", "bfgs");

    optionsIpOpt->SetIntegerValue(
                "max_iter", problem->getOptimizationOptions().maxOptimizerIterations);

    // set IpOpt options from OptimizationOptions
    problem->getOptimizationOptions().for_each<SmartPtr<OptionsList> *>(setIpOptOption, &optionsIpOpt);
}

OptimizerIpOpt::OptimizerIpOpt() {}

std::tuple<int, double, std::vector<double> > OptimizerIpOpt::optimize(OptimizationProblem *problem) {

    ApplicationReturnStatus status = Unrecoverable_Exception;

    std::vector<double> finalParameters;
    double finalCost = NAN;

    { // ensure all IpOpt objects are destroyed before mutex is unlocked

        // lock because we pass control to IpOpt
        auto lock = ipOptGetLock();

        auto optimizationController = problem->getReporter();

        try {
            SmartPtr<TNLP> mynlp =
                    new LocalOptimizationIpoptTNLP(*problem, *optimizationController);
            SmartPtr<IpoptApplication> app = IpoptApplicationFactory();
            app->RethrowNonIpoptException(true);

            setIpOptOptions(app->Options(), problem);

            // TODO: for output redirection
            // app->Jnlst()->AddJournal();

            status = app->Initialize();
            assert(status == Solve_Succeeded);
            status = app->OptimizeTNLP(mynlp);

            if(status == Invalid_Number_Detected) {
                // TODO: print where
            }
        } catch (IpoptException& e) {
            logmessage(LOGLVL_ERROR, "IpOpt exception: %s",  e.Message().c_str());
        } catch (std::exception& e) {
            logmessage(LOGLVL_ERROR, "Unknown exception occured during optimization: %s", e.what());
        } catch (...) {
            logmessage(LOGLVL_ERROR, "Unknown exception occured during optimization");
        }
        finalCost = optimizationController->getFinalCost();
        finalParameters = optimizationController->getFinalParameters();

    }

    // TODO: need smarter way to decide if should retry or not
    //    if((int)status < Not_Enough_Degrees_Of_Freedom) {
    //        // should exit, retrying probably makes no sense
    //        throw ParPEException(std::string("Unrecoverable IpOpt problem - see messages above. Code ") + std::to_string(status));
    //    }

    return std::tuple<int, double, std::vector<double> >((int)status < Maximum_Iterations_Exceeded, finalCost, finalParameters);
}

} // namespace parpe
