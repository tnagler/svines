#pragma once

#include "svine_structure.hpp"
#include "svine_selector.hpp"
#include "tools_stl.hpp"
#include <vinecopulib/misc/tools_optimization.hpp>
#include <vinecopulib/misc/tools_stats.hpp>
#include <vinecopulib/misc/tools_stl.hpp>
#include <vinecopulib/vinecop/class.hpp>
#include <vinecopulib/vinecop/tools_select.hpp>
#include <wdm/eigen.hpp>

// TODO: weights, discrete ?

namespace vinecopulib {

// --------------------- S-vine ----------------------------------

class SVine : public Vinecop
{
public:
  SVine(size_t cs_dim, size_t p, const std::vector<std::string>& var_types = {})
    : SVine(RVineStructure(tools_stl::seq_int(1, cs_dim)),
            p,
            tools_stl::rev(tools_stl::seq_int(1, cs_dim)),
            tools_stl::rev(tools_stl::seq_int(1, cs_dim)),
            var_types)
  {}

  SVine(const RVineStructure& cs_struct,
        size_t p,
        std::vector<size_t> in_vertices,
        std::vector<size_t> out_vertices,
        const std::vector<std::string>& var_types = {})
    : cs_dim_(cs_struct.get_dim())
    , p_(p)
    , in_vertices_(in_vertices)
    , out_vertices_(out_vertices)
    , svine_struct_(SVineStructure(cs_struct, p, in_vertices, out_vertices))
  {
    if (var_types.size() == 0) {
      var_types_ = std::vector<std::string>(svine_struct_.get_dim(), "c");
    } else {
      var_types_ = tools_stl::rep(var_types, (p + 1));
    }
    d_ = svine_struct_.get_dim();
    this->set_var_types(var_types_);
    threshold_ = 0.0;
    loglik_ = NAN;
    rvine_structure_ = svine_struct_;
    pair_copulas_ = make_pair_copula_store(d_);
  }

  SVine(const std::vector<std::vector<Bicop>>& pair_copulas,
        const RVineStructure& cs_struct,
        size_t p,
        std::vector<size_t> in_vertices,
        std::vector<size_t> out_vertices,
        const std::vector<std::string>& var_types = {})
    : SVine(cs_struct, p, in_vertices, out_vertices, var_types)
  {
    pair_copulas_ = pair_copulas;
  }

  size_t get_p() const { return p_; }

  size_t get_cs_dim() const { return cs_dim_; }

  std::vector<size_t> get_in_vertices() const { return in_vertices_; }

  std::vector<size_t> get_out_vertices() const { return out_vertices_; }

  RVineStructure get_cs_structure() const
  {
    return svine_struct_.get_cs_structure();
  }

  SVineStructure get_svine_structure() const { return svine_struct_; }

  void select_families(
    const Eigen::MatrixXd& data,
    const FitControlsVinecop& controls = FitControlsVinecop())
  {
    tools_eigen::check_if_in_unit_cube(data);
    check_data_dim(data);

    if (rvine_structure_.get_trunc_lvl() > 0) {
      auto vt0 = tools_stl::span(var_types_, 0, cs_dim_);
      tools_select::SVineFamilySelector selector(
        data,
        svine_struct_.get_cs_structure(),
        controls,
        in_vertices_,
        out_vertices_,
        vt0);

      selector.select_all_trees(data);
      for (size_t lag = 1; lag <= p_; lag++) {
        selector.add_lag();
        selector.select_all_trees(selector.data());
      }

      finalize_fit(selector);
      loglik_ = this->loglik(data);
    }
  }

  void select_all(const Eigen::MatrixXd& data,
                  const FitControlsVinecop& controls = FitControlsVinecop())
  {
    tools_eigen::check_if_in_unit_cube(data);
    check_data_dim(data);

    auto vt0 = tools_stl::span(var_types_, 0, cs_dim_);
    tools_select::SVineStructureSelector selector(data, controls, vt0);
    selector.select_all_trees(data);
    for (size_t lag = 1; lag <= p_; lag++) {
      selector.add_lag();
    }
    finalize_fit(selector);
    loglik_ = this->loglik(data);
  }

  Eigen::MatrixXd simulate(const size_t n,
                           const bool qrng = false,
                           const std::vector<int>& seeds = std::vector<int>())
  {
    auto U = tools_stats::simulate_uniform(n, cs_dim_, qrng, seeds);

    // initialize first p + 1 lags
    Eigen::MatrixXd sim(n, cs_dim_);
    Eigen::MatrixXd Ui(1, d_);
    for (size_t i = 0; i <= p_; i++) {
      Ui.row(0).segment(i * cs_dim_, cs_dim_) = U.row(i);
    }
    Eigen::MatrixXd V = inverse_rosenblatt(Ui);
    for (size_t i = 0; i <= p_; i++) {
      sim.row(i) = V.block(0, i * cs_dim_, 1, cs_dim_);
    }

    // simulate conditional on previous observations
    for (size_t i = p_ + 1; i < n; i++) {
      Ui.leftCols(d_ - cs_dim_) = get_last_cpits(sim.topRows(i));
      Ui.rightCols(cs_dim_) = U.row(i);
      sim.row(i) = inverse_rosenblatt(Ui).rightCols(cs_dim_);
    }

    return sim;
  }

  Eigen::MatrixXd simulate_conditional(
    size_t n,
    const Eigen::MatrixXd& data,
    const bool qrng = false,
    const size_t num_threads = 1,
    const std::vector<int>& seeds = std::vector<int>())
  {
    check_cond_data(data);

    Eigen::MatrixXd U(n, d_);
    if (p_ > 0) {
      U.leftCols(d_ - cs_dim_) = get_last_cpits(data).replicate(n, 1);
    }
    U.rightCols(cs_dim_) =
      tools_stats::simulate_uniform(n, cs_dim_, qrng, seeds);

    return inverse_rosenblatt(U, num_threads).rightCols(cs_dim_);
  }

  Eigen::MatrixXd simulate_ahead(
    size_t n_ahead,
    const Eigen::MatrixXd& data,
    const bool qrng = false,
    const std::vector<int>& seeds = std::vector<int>())
  {
    check_cond_data(data);

    Eigen::MatrixXd U(n_ahead + p_, cs_dim_);
    U.bottomRows(n_ahead) =
      tools_stats::simulate_uniform(n_ahead, cs_dim_, qrng, seeds);

    Eigen::MatrixXd Ui(1, d_);
    Eigen::MatrixXd sim(n_ahead + p_, cs_dim_);
    sim.topRows(p_) = data.bottomRows(p_);
    for (size_t i = p_; i < p_ + n_ahead; i++) {
      Ui.leftCols(d_ - cs_dim_) = get_last_cpits(sim.topRows(i));
      Ui.rightCols(cs_dim_) = U.row(i);
      sim.row(i) = inverse_rosenblatt(Ui).rightCols(cs_dim_);
    }

    return sim.bottomRows(n_ahead);
  }

  Eigen::VectorXd pdf(const Eigen::MatrixXd&, const size_t = 1) const
  {
    throw std::runtime_error("pdf not meaningful for S-vines; use loglik().");
  }

  double loglik(const Eigen::MatrixXd& u, const size_t num_threads = 1)
  {
    if (static_cast<size_t>(u.cols()) != cs_dim_)
      throw std::runtime_error("dimension of data and model don't match.");
    size_t n = u.rows();

    // iid model, can return loglik directly
    if ((p_ == 0) || (n == 1)) {
      rvine_structure_ = svine_struct_.get_cs_structure();
      return Vinecop::loglik(u, num_threads);
    }

    // first compute substraction component (otherwise some contributions
    // are counted twice)
    size_t p_tmp = std::min(n, p_) - 1;
    rvine_structure_ = SVineStructure(
      svine_struct_.get_cs_structure(), p_tmp, in_vertices_, out_vertices_);
    d_ = cs_dim_ * (1 + p_tmp);
    auto u_spr = u;
    for (size_t lag = 0; lag < p_tmp; ++lag) {
      u_spr = spread_lag(u_spr, cs_dim_);
    }

    n = u_spr.rows();
    double ll = 0.0;
    if (n > 2) {
      ll -=
        Vinecop::loglik(u_spr.bottomRows(n - 1).topRows(n - 2), num_threads);
    } else {
      ll -= Vinecop::loglik(u_spr.bottomRows(n - 1), num_threads);
    }

    // add loglik *as if* it was iid with cs_dim * (1 + p) vars
    u_spr = spread_lag(u_spr, cs_dim_);
    rvine_structure_ = SVineStructure(
      svine_struct_.get_cs_structure(), p_, in_vertices_, out_vertices_);
    d_ = cs_dim_ * (1 + p_);
    ll += Vinecop::loglik(u_spr, num_threads);

    return ll;
  }

  Eigen::VectorXd cond_cdf(const Eigen::MatrixXd& u,
                           size_t conditioned,
                           const size_t num_threads = 1) const
  {
    if (static_cast<size_t>(u.cols()) != cs_dim_)
      throw std::runtime_error("dimension of data and model don't match.");
    if (static_cast<size_t>(u.rows()) <= p_)
      throw std::runtime_error("insufficient number of time points.");

    auto v = u;
    for (size_t lag = 0; lag < p_; ++lag) {
      v = spread_lag(v, cs_dim_);
      conditioned += cs_dim_;
    }

    size_t n = v.rows();
    Eigen::VectorXd seq = Eigen::VectorXd::LinSpaced(100, 1e-10, 1 - 1e-10);
    Eigen::MatrixXd vv;
    Eigen::VectorXd out(n);
    for (size_t i = 0; i < n; ++i) {
      vv = v.row(i).replicate(100, 1);
      vv.col(conditioned) = seq;
      auto pdf = Vinecop::pdf(vv, num_threads);
      out(i) = pdf.head(std::ceil(v(i, conditioned) * 100)).sum();
      out(i) /= pdf.sum();
    }
    return out;
  }

  Eigen::VectorXi get_num_pars()
  {
    Eigen::VectorXi nums(cs_dim_ * cs_dim_ * p_ + cs_dim_ * (cs_dim_ - 1) / 2);
    size_t i = 0;
    for (size_t t = 0; t < cs_dim_ * (1 + p_) - 1; ++t) {
      for (size_t e = 0; e < cs_dim_; ++e) {
        if (e < pair_copulas_[t].size()) {
          if (pair_copulas_[t][e].get_family() == BicopFamily::tll) {
            nums(i++) = 0;
          } else {
            nums(i++) = pair_copulas_[t][e].get_parameters().size();
          }
        }
      }
    }

    return nums;
  }

protected:
  Eigen::MatrixXd get_last_cpits(const Eigen::MatrixXd& data)
  {
    auto cpits = Eigen::MatrixXd();

    if (p_ > 0) {
      // only most recent observations are used
      Eigen::MatrixXd cond_vals = data.bottomRows(p_);

      // spread past observations into one row with d_ - cs_dim_ columns
      for (size_t lag = 1; lag < p_; lag++) {
        cond_vals = spread_lag(cond_vals, cs_dim_);
      }

      // construct sub-model for last p_ lags
      d_ -= cs_dim_;
      rvine_structure_ = SVineStructure(
        svine_struct_.get_cs_structure(), p_ - 1, in_vertices_, out_vertices_);

      // initialize Ui with rosenblatt of past observations
      cpits = rosenblatt(cond_vals);

      // restore original model
      rvine_structure_ = svine_struct_;
      d_ += cs_dim_;
    }

    return cpits;
  }

  void finalize_fit(const tools_select::SVineFamilySelector& selector)
  {
    in_vertices_ = selector.get_in_vertices();
    out_vertices_ = selector.get_out_vertices();
    svine_struct_ = SVineStructure(
      selector.get_cs_structure(), p_, in_vertices_, out_vertices_);
    rvine_structure_ = svine_struct_;
    Vinecop::finalize_fit(selector);
    var_types_ = selector.get_var_types();
  }

  void finalize_fit(tools_select::SVineStructureSelector& selector)
  {
    selector.finalize(std::numeric_limits<size_t>::max());

    in_vertices_ = selector.get_in_vertices();
    out_vertices_ = selector.get_out_vertices();
    Vinecop::finalize_fit(selector);
    svine_struct_ = SVineStructure(
      selector.get_cs_structure(), p_, in_vertices_, out_vertices_);
    rvine_structure_ = svine_struct_;
    var_types_ = selector.get_var_types();
  }

  void check_data_dim(const Eigen::MatrixXd& data) const
  {
    if (cs_dim_ != static_cast<size_t>(data.cols())) {
      std::stringstream msg;
      msg << "wrong number of columns." << std::endl
          << "expected: " << cs_dim_ << std::endl
          << "provided: " << data.cols() << std::endl;
      throw std::runtime_error(msg.str());
    }
  }

  void check_cond_data(const Eigen::MatrixXd& data) const
  {
    check_data_dim(data);
    if (static_cast<size_t>(data.rows()) < p_) {
      std::stringstream msg;
      msg << "need at least p observations to condition on;" << std::endl
          << "expected: >= " << p_ << std::endl
          << "actual: " << data.rows() << std::endl;
      throw std::runtime_error(msg.str());
    }
  }

  size_t cs_dim_;
  size_t p_;
  std::vector<size_t> in_vertices_;
  std::vector<size_t> out_vertices_;
  SVineStructure svine_struct_;
};
}
