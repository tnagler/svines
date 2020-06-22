#define BOOST_NO_AUTO_PTR 1

#include "Rcpp.h"
#include "svine.hpp"
#include <vinecopulib-wrappers.hpp>

using namespace vinecopulib;

Svine
svinecop_wrap(const Rcpp::List& svinecop_r)
{
  size_t p = svinecop_r["p"];

  std::vector<size_t> in_vertices = svinecop_r["in_vertices"];
  std::vector<size_t> out_vertices = svinecop_r["out_vertices"];

  // omit R-vine matrix check, already done in R
  auto cs_structure = rvine_structure_wrap(svinecop_r["cs_structure"], false);

  // extract pair-copulas
  auto pair_copulas = pair_copulas_wrap(svinecop_r["pair_copulas"],
                                        cs_structure.get_dim() * (p + 1));

  std::vector<std::string> var_types = svinecop_r["var_types"];

  return Svine(
    pair_copulas, cs_structure, p, in_vertices, out_vertices, var_types);
}

Rcpp::List
svinecop_wrap(const Svine& svinecop_cpp, bool is_fitted)
{
  auto vine_structure =
    rvine_structure_wrap(svinecop_cpp.get_rvine_structure());
  auto cs_structure = rvine_structure_wrap(svinecop_cpp.get_cs_structure());
  std::vector<size_t> in_vertices = svinecop_cpp.get_in_vertices();
  std::vector<size_t> out_vertices = svinecop_cpp.get_out_vertices();
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
                            Rcpp::Named("in_vertices") = in_vertices,
                            Rcpp::Named("p") = p,
                            Rcpp::Named("out_vertices") = out_vertices,
                            Rcpp::Named("cs_structure") = cs_structure);
}

// [[Rcpp::export()]]
Rcpp::List
svinecop_create_cpp(const Rcpp::List& svine_r)
{
  Svine tv_cpp = svinecop_wrap(svine_r);
  return svinecop_wrap(tv_cpp, false);
}

// [[Rcpp::export()]]
Rcpp::List
svinecop_select_cpp(const Eigen::MatrixXd& data,
                    size_t p,
                    const std::vector<std::string>& var_types,
                    std::vector<size_t> in_vertices,
                    std::vector<size_t> out_vertices,
                    bool is_structure_provided,
                    Rcpp::List& structure,
                    std::vector<std::string> family_set,
                    std::string par_method,
                    std::string nonpar_method,
                    double mult,
                    int truncation_level,
                    std::string tree_criterion,
                    double threshold,
                    std::string selection_criterion,
                    const Eigen::VectorXd& weights,
                    double psi0,
                    bool select_truncation_level,
                    bool select_threshold,
                    bool preselect_families,
                    bool show_trace,
                    size_t num_threads)
{
  std::vector<BicopFamily> fam_set(family_set.size());
  for (unsigned int fam = 0; fam < fam_set.size(); ++fam) {
    fam_set[fam] = to_cpp_family(family_set[fam]);
  }

  FitControlsVinecop fit_controls(fam_set,
                                  par_method,
                                  nonpar_method,
                                  mult,
                                  truncation_level,
                                  tree_criterion,
                                  threshold,
                                  selection_criterion,
                                  weights,
                                  psi0,
                                  preselect_families,
                                  select_truncation_level,
                                  select_threshold,
                                  show_trace,
                                  num_threads);

  Svine svine(data.cols(), p, var_types);
  if (is_structure_provided) {
    svine = Svine(rvine_structure_wrap(structure, false),
                  p,
                  in_vertices,
                  out_vertices,
                  var_types);
    svine.select_families(data, fit_controls);
  } else {
    svine.select_all(data, fit_controls);
  }

  return svinecop_wrap(svine, TRUE);
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
Eigen::VectorXd
svinecop_cond_cdf_cpp(const Eigen::MatrixXd& u,
                      size_t conditioned,
                      const Rcpp::List& svinecop_r,
                      size_t cores)
{
  return svinecop_wrap(svinecop_r).cond_cdf(u, conditioned, cores);
}

// [[Rcpp::export()]]
Eigen::MatrixXd
svinecop_sim_cpp(const Rcpp::List& svinecop_r,
                 const size_t n,
                 const bool qrng,
                 std::vector<int> seeds)
{
  return svinecop_wrap(svinecop_r).simulate(n, qrng, seeds);
}

// [[Rcpp::export()]]
Eigen::MatrixXd
svinecop_sim_conditional_cpp(const Rcpp::List& svinecop_r,
                             const size_t n,
                             const Eigen::MatrixXd& data,
                             const bool qrng,
                             size_t cores,
                             const std::vector<int>& seeds)
{
  return svinecop_wrap(svinecop_r)
    .simulate_conditional(n, data, qrng, cores, seeds);
}

// [[Rcpp::export()]]
Eigen::MatrixXd
svinecop_sim_ahead_cpp(const Rcpp::List& svinecop_r,
                       const size_t n_ahead,
                       const Eigen::MatrixXd& data,
                       const bool qrng,
                       const std::vector<int>& seeds)
{
  return svinecop_wrap(svinecop_r).simulate_ahead(n_ahead, data, qrng, seeds);
}
