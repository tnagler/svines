#pragma once

#include <vinecopulib/vinecop/class.hpp>
#include <vinecopulib/vinecop/tools_select.hpp>
#include <wdm/eigen.hpp>

namespace vinecopulib {

Eigen::MatrixXd
spread_lag(const Eigen::MatrixXd& data, size_t cs_dim);

namespace tools_select {

//! @brief Abstract selector class.
class SVineSelector
{
public:
  SVineSelector(const Eigen::MatrixXd& data,
                std::vector<size_t> in_vertices,
                std::vector<size_t> out_vertices,
                const std::vector<std::string>& var_types);

  SVineSelector(const Eigen::MatrixXd& data,
                const std::vector<std::string>& var_types);

  std::vector<size_t> get_in_vertices() const;

  std::vector<size_t> get_out_vertices() const;

  RVineStructure get_cs_structure() const;

protected:
  void duplicate_vertex(size_t v, VineTree& tree);

  void duplicate_edge(EdgeIterator e, VineTree& tree);

  void check_in_out_vertices() const;
  void check_controls(const FitControlsVinecop& controls);

  size_t cs_dim_;
  size_t lag_;
  std::vector<size_t> in_vertices_;
  std::vector<size_t> out_vertices_;
  Eigen::MatrixXd data_;
  RVineStructure cs_struct_;
};

//! @brief Selector class for copulas + structure.
class SVineStructureSelector
  : public VinecopSelector
  , public SVineSelector
{
public:
  SVineStructureSelector(const Eigen::MatrixXd& data,
                         const FitControlsVinecop& controls,
                         const std::vector<std::string>& var_types);

  ~SVineStructureSelector() = default;

  std::vector<std::string> get_var_types() const;

  void select_all_trees(const Eigen::MatrixXd& data);

  void add_lag();

  VineTree make_base_tree(const Eigen::MatrixXd& data);

  void duplicate_edges(VineTree& old_tree, VineTree& new_tree, size_t t);

  void add_allowed_connections(VineTree& tree, size_t t);

  void finalize_svine(size_t trunc_lvl);

  void finalize(size_t trunc_lvl);

protected:
  double compute_fit_id(const EdgeProperties& e) override;
};

//! @brief Selector class for copulas only.
class SVineFamilySelector
  : public VinecopSelector
  , public SVineSelector
{
public:
  SVineFamilySelector(const Eigen::MatrixXd& data,
                      const RVineStructure& cs_struct,
                      const FitControlsVinecop& controls,
                      std::vector<size_t> in_vertices,
                      std::vector<size_t> out_vertices,
                      const std::vector<std::string>& var_types);

  ~SVineFamilySelector() = default;

  std::vector<std::string> get_var_types() const;

  void select_tree(size_t t) override;

  void add_lag();

  Eigen::MatrixXd data();

protected:
  void finalize(size_t);

  double compute_fit_id(const EdgeProperties& e) override;

  void flip_pcs(std::vector<VineTree>& trees);
};

}
}

#include "implementation/svine_selector.ipp"