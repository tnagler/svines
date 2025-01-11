#include <vinecopulib-wrappers.hpp>
#include "Rcpp.h"
#include "svines.hpp"

using namespace vinecopulib;

SVinecop
svinecop_wrap(const Rcpp::List& svinecop_r)
{
  size_t p = svinecop_r["p"];

  std::vector<size_t> out_vertices = svinecop_r["out_vertices"];
  std::vector<size_t> in_vertices = svinecop_r["in_vertices"];

  // omit R-vine matrix check, already done in R
  auto cs_structure = rvine_structure_wrap(svinecop_r["cs_structure"], false);

  // extract pair-copulas
  auto pair_copulas = pair_copulas_wrap(svinecop_r["pair_copulas"],
                                        cs_structure.get_dim() * (p + 1));

  std::vector<std::string> var_types = svinecop_r["var_types"];

  return SVinecop(
    pair_copulas, cs_structure, p, out_vertices, in_vertices, var_types);
}

Rcpp::List
svinecop_wrap(const SVinecop& svinecop_cpp, bool is_fitted)
{
  auto vine_structure =
    rvine_structure_wrap(svinecop_cpp.get_rvine_structure());
  auto cs_structure = rvine_structure_wrap(svinecop_cpp.get_cs_structure());
  std::vector<size_t> out_vertices = svinecop_cpp.get_out_vertices();
  std::vector<size_t> in_vertices = svinecop_cpp.get_in_vertices();
  size_t p = svinecop_cpp.get_p();

  auto pair_copulas = pair_copulas_wrap(
    svinecop_cpp.get_all_pair_copulas(), svinecop_cpp.get_dim(), false);

  double npars = svinecop_cpp.get_npars();
  double threshold = svinecop_cpp.get_threshold();
  double loglik = NAN;
  auto var_types = svinecop_cpp.get_var_types();
  var_types.resize(svinecop_cpp.get_cs_dim());

  if (is_fitted)
    loglik = svinecop_cpp.get_loglik();
  return Rcpp::List::create(Rcpp::Named("pair_copulas") = pair_copulas,
                            Rcpp::Named("structure") = vine_structure,
                            Rcpp::Named("var_types") = var_types,
                            Rcpp::Named("npars") = npars,
                            Rcpp::Named("loglik") = loglik,
                            Rcpp::Named("threshold") = threshold,
                            Rcpp::Named("p") = p,
                            Rcpp::Named("out_vertices") = out_vertices,
                            Rcpp::Named("in_vertices") = in_vertices,
                            Rcpp::Named("cs_structure") = cs_structure);
}

// [[Rcpp::export()]]
Rcpp::List
svinecop_create_cpp(const Rcpp::List& svine_r)
{
  SVinecop tv_cpp = svinecop_wrap(svine_r);
  return svinecop_wrap(tv_cpp, false);
}

// [[Rcpp::export()]]
Rcpp::List
svinecop_select_cpp(const Eigen::MatrixXd& data,
                    size_t p,
                    const std::vector<std::string>& var_types,
                    const std::vector<size_t>& out_vertices,
                    const std::vector<size_t>& in_vertices,
                    bool is_structure_provided,
                    Rcpp::List& structure,
                    const std::vector<std::string>& family_set,
                    std::string par_method,
                    std::string nonpar_method,
                    double mult,
                    int trunc_lvl,
                    std::string tree_criterion,
                    double threshold,
                    std::string selection_criterion,
                    const Eigen::VectorXd& weights,
                    double psi0,
                    bool select_trunc_lvl,
                    bool select_threshold,
                    bool preselect_families,
                    bool show_trace,
                    size_t num_threads)
{
  std::vector<BicopFamily> fam_set(family_set.size());
  for (unsigned int fam = 0; fam < fam_set.size(); ++fam) {
    fam_set[fam] = to_cpp_family(family_set[fam]);
  }
  FitControlsVinecop fit_controls(fam_set);
  fit_controls.set_family_set(fam_set);
  fit_controls.set_parametric_method(par_method);
  fit_controls.set_nonparametric_method(nonpar_method);
  fit_controls.set_nonparametric_mult(mult);
  fit_controls.set_trunc_lvl(trunc_lvl);
  fit_controls.set_tree_criterion(tree_criterion);
  fit_controls.set_threshold(threshold);
  fit_controls.set_selection_criterion(selection_criterion);
  fit_controls.set_weights(weights);
  fit_controls.set_psi0(psi0);
  fit_controls.set_preselect_families(preselect_families);
  fit_controls.set_select_threshold(select_threshold);
  fit_controls.set_select_trunc_lvl(select_trunc_lvl);
  fit_controls.set_show_trace(show_trace);
  fit_controls.set_num_threads(num_threads);
  
  SVinecop svine(var_types.size(), p, var_types);
  // if (is_structure_provided) {
  //   svine = SVinecop(rvine_structure_wrap(structure, false),
  //                    p,
  //                    out_vertices,
  //                    in_vertices,
  //                    var_types);
  //   svine.select_families(data, fit_controls);
  // } else {
    svine.select_all(data); //, fit_controls);
  // }

  return svinecop_wrap(svine, true);
}

// [[Rcpp::export()]]
double
svinecop_loglik_cpp(const Eigen::MatrixXd& u,
                    const Rcpp::List& svinecop_r,
                    size_t cores)
{
  return svinecop_wrap(svinecop_r).loglik(u, cores);
}


// [[Rcpp::export()]]
Eigen::MatrixXd
svinecop_sim_cpp(const Rcpp::List& svinecop_r,
                 const size_t n,
                 const size_t rep,
                 const Eigen::MatrixXd& data,
                 const bool qrng,
                 const size_t cores,
                 const std::vector<int>& seeds)
{
  auto sv_cpp = svinecop_wrap(svinecop_r);

  // make sure everything is random, but reproducible
  std::vector<std::vector<int>> new_seeds(rep);
  {
    auto tmp_seeds = tools_stats::simulate_uniform(rep, 10, false, seeds);
    for (size_t i = 0; i < rep; i++) {
      new_seeds[i].resize(10);
      for (int k = 0; k < 10; k++)
        new_seeds[i][k] = 
          std::floor(tmp_seeds(i, k) * std::numeric_limits<int>::max());
    }
  }

  auto cs_dim = sv_cpp.get_cs_dim();
  Eigen::MatrixXd sim(n, cs_dim * rep);
  
  RcppThread::parallelFor(
    0,
    rep,
    [&](size_t r) {
      if (data.size() == 0) {
        sim.block(0, cs_dim * r, n, cs_dim) =
          sv_cpp.as_continuous().simulate(n, qrng, new_seeds[r]);
      } else {
        sim.block(0, cs_dim * r, n, cs_dim) =
        sv_cpp.as_continuous().simulate_ahead(n, data, qrng, new_seeds[r]);
      }
    },
    cores);

  return sim;
}

// [[Rcpp::export()]]
Eigen::MatrixXd
svinecop_pseudo_residuals_cpp(const Eigen::MatrixXd& u,
                              const Rcpp::List& svinecop_r,
                              const size_t num_threads)
{
  return svinecop_wrap(svinecop_r).pseudo_residuals(u, num_threads);
}


// [[Rcpp::export()]]
Eigen::MatrixXd
svinecop_scores_cpp(const Eigen::MatrixXd& u,
                    const Rcpp::List& svinecop_r,
                    const size_t num_threads)
{
  return svinecop_wrap(svinecop_r).scores(u, true, num_threads);
}

// [[Rcpp::export()]]
Eigen::MatrixXd
svinecop_hessian_cpp(const Eigen::MatrixXd& u,
                     const Rcpp::List& svinecop_r,
                     const size_t num_threads)
{
  return svinecop_wrap(svinecop_r).hessian_exp(u, true, num_threads);
}

// [[Rcpp::export()]]
Rcpp::List
with_parameters_cop_cpp(const Rcpp::List& svinecop_r,
                        const Eigen::VectorXd parameters)
{
  auto svc = svinecop_wrap(svinecop_r);
  auto d0 = svc.get_cs_dim();
  auto d = svc.get_dim();
  auto p = svc.get_p();
  auto pcs = svc.get_all_pair_copulas();

  size_t i = 0;
  for (size_t t = 0; t < d - 1; t++) {
    for (size_t e = 0; e < std::min(d0, d - t - 1); e++) {
      if (pcs[t][e].get_family() == BicopFamily::indep)
        continue;
      auto lb = pcs[t][e].get_parameters_lower_bounds();
      auto ub = pcs[t][e].get_parameters_upper_bounds();
      auto new_pars =
        parameters.segment(i, lb.size()).cwiseMax(lb).cwiseMin(ub);
      pcs[t][e].set_parameters(new_pars);
      for (size_t lag = 1; lag <= p; lag++) {
        if (e + d0 * lag < d - t - 1) {
          pcs[t][e + d0 * lag] = pcs[t][e];
        }
      }
      i += lb.size();
    }
  }

  svc = SVinecop(pcs,
                 svc.get_cs_structure(),
                 p,
                 svc.get_out_vertices(),
                 svc.get_in_vertices());
  return svinecop_wrap(svc, false);
}