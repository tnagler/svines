#pragma once

#include "svine_structure.hpp"
#include "svine_selector.hpp"
#include <vinecopulib/vinecop/class.hpp>

// TODO: weights, discrete ?

namespace vinecopulib {

class SVinecop : public Vinecop
{
public:
  SVinecop(size_t cs_dim,
           size_t p,
           const std::vector<std::string>& var_types = {});

  SVinecop(const RVineStructure& cs_struct,
           size_t p,
           std::vector<size_t> in_vertices,
           std::vector<size_t> out_vertices,
           const std::vector<std::string>& var_types = {});

  SVinecop(const std::vector<std::vector<Bicop>>& pair_copulas,
           const RVineStructure& cs_struct,
           size_t p,
           std::vector<size_t> in_vertices,
           std::vector<size_t> out_vertices,
           const std::vector<std::string>& var_types = {});

  size_t get_p() const;

  size_t get_cs_dim() const;

  std::vector<size_t> get_in_vertices() const;

  std::vector<size_t> get_out_vertices() const;

  RVineStructure get_cs_structure() const;

  SVineStructure get_svine_structure() const;

  void select_families(
    const Eigen::MatrixXd& data,
    const FitControlsVinecop& controls = FitControlsVinecop());

  void select_all(const Eigen::MatrixXd& data,
                  const FitControlsVinecop& controls = FitControlsVinecop());

  Eigen::MatrixXd simulate(const size_t n,
                           const bool qrng = false,
                           const std::vector<int>& seeds = std::vector<int>());

  Eigen::MatrixXd simulate_conditional(
    size_t n,
    const Eigen::MatrixXd& data,
    const bool qrng = false,
    const size_t num_threads = 1,
    const std::vector<int>& seeds = std::vector<int>());

  Eigen::MatrixXd simulate_ahead(
    size_t n_ahead,
    const Eigen::MatrixXd& data,
    const bool qrng = false,
    const std::vector<int>& seeds = std::vector<int>());

  Eigen::VectorXd pdf(const Eigen::MatrixXd&, const size_t = 1) const;

  double loglik(const Eigen::MatrixXd& u, const size_t num_threads = 1);

  Eigen::VectorXd cond_cdf(const Eigen::MatrixXd& u,
                           size_t conditioned,
                           const size_t num_threads = 1) const;

  Eigen::VectorXi get_num_pars();

protected:
  Eigen::MatrixXd get_last_cpits(const Eigen::MatrixXd& data);

  void finalize_fit(const tools_select::SVineFamilySelector& selector);

  void finalize_fit(tools_select::SVineStructureSelector& selector);

  void check_data_dim(const Eigen::MatrixXd& data) const;

  void check_cond_data(const Eigen::MatrixXd& data) const;

  size_t cs_dim_;
  size_t p_;
  std::vector<size_t> in_vertices_;
  std::vector<size_t> out_vertices_;
  SVineStructure svine_struct_;
};
}

#include "implementation/svinecop.ipp"