#pragma once

#include "svines/svine_structure.hpp"
#include <vinecopulib/misc/tools_stl.hpp>
#include <vinecopulib/vinecop/tools_select.hpp>
#include <wdm/eigen.hpp>

namespace vinecopulib {

//! @param Removes one row and adds more column so data correspond to one time
//! lag more.
//!
//! @param data the data.
//! @param cs_dim cross-sectional dimension.
inline Eigen::MatrixXd
spread_lag(const Eigen::MatrixXd& data, size_t cs_dim)
{
  if (data.rows() < 2) {
    throw std::runtime_error("insufficient number of observations");
  }
  if (data.cols() % cs_dim != 0) {
    throw std::runtime_error("number of columns is not a multiple of cs_dim");
  }
  size_t n = data.rows() - 1;
  Eigen::MatrixXd newdata(n, data.cols() + cs_dim);
  newdata << data.topRows(n), data.rightCols(cs_dim).bottomRows(n);
  return newdata;
}

namespace tools_select {

// ------- SVineSelector

SVineSelector::SVineSelector(const Eigen::MatrixXd& data,
                             std::vector<size_t> out_vertices,
                             std::vector<size_t> in_vertices,
                             const std::vector<std::string>& var_types)
  : cs_dim_(var_types.size())
  , lag_(0)
  , out_vertices_(out_vertices)
  , in_vertices_(in_vertices)
  , data_(data)
{
  check_out_in_vertices();
}

SVineSelector::SVineSelector(const Eigen::MatrixXd& data,
                             const std::vector<std::string>& var_types)
  : cs_dim_(var_types.size())
  , lag_(0)
  , data_(data)
{}


inline std::vector<size_t>
SVineSelector::get_out_vertices() const
{
  return out_vertices_;
}

inline std::vector<size_t>
SVineSelector::get_in_vertices() const
{
  return in_vertices_;
}

inline RVineStructure
SVineSelector::get_cs_structure() const
{
  return cs_struct_;
}

inline void
SVineSelector::duplicate_vertex(size_t v, VineTree& tree)
{
  auto v_new = boost::add_vertex(tree);
  auto shift = [this](std::vector<size_t> index) {
    for (auto& i : index)
      i = i + cs_dim_ * lag_;
    return index;
  };

  // copy structure information
  tree[v_new].conditioned = shift(tree[v].conditioned);
  tree[v_new].conditioning = shift(tree[v].conditioning);
  tree[v_new].all_indices = shift(tree[v].all_indices);
  tree[v_new].prev_edge_indices = shift(tree[v].prev_edge_indices);
  tree[v_new].var_types = tree[v].var_types;

  // copy data and remove rows
  size_t n = tree[v].hfunc1.rows() - 1;
  tree[v_new].hfunc1 = tree[v].hfunc1.bottomRows(n);
  tree[v].hfunc1.conservativeResize(n);
  if (tree[v].hfunc1_sub.size()) {
    tree[v_new].hfunc1_sub = tree[v].hfunc1_sub.bottomRows(n);
    tree[v].hfunc1_sub.conservativeResize(n);
  }
  if (tree[v].hfunc2.size() > 1) {
    tree[v_new].hfunc2 = tree[v].hfunc2.bottomRows(n);
    tree[v].hfunc2.conservativeResize(n);
    if (tree[v].hfunc2_sub.size()) {
      tree[v_new].hfunc2_sub = tree[v].hfunc2_sub.bottomRows(n);
      tree[v].hfunc2_sub.conservativeResize(n);
    }
  }
}

inline void
SVineSelector::duplicate_edge(EdgeIterator e, VineTree& tree)
{
  size_t v1 = boost::source(e, tree);
  size_t v2 = boost::target(e, tree);
  auto e_new = boost::add_edge(v1 + cs_dim_, v2 + cs_dim_, tree);
  tree[e_new.first].pair_copula = tree[e].pair_copula;
  tree[e_new.first].fit_id = tree[e].fit_id;
  auto shift = [this](std::vector<size_t> index) {
    for (auto& i : index)
      i = i + cs_dim_ * lag_;
    return index;
  };
  tree[e_new.first].conditioned = shift(tree[e].conditioned);
  tree[e_new.first].conditioning = shift(tree[e].conditioning);
  tree[e_new.first].all_indices = shift(tree[e].all_indices);
  tree[e_new.first].var_types = tree[e].var_types;
}

inline void
SVineSelector::check_out_in_vertices() const
{
  auto d = cs_dim_;
  if (!tools_stl::is_same_set(in_vertices_, tools_stl::seq_int(1, d)))
    throw std::runtime_error(
      "in_vertices must contain numbers 1, ..., cs_dim.");
  if (!tools_stl::is_same_set(out_vertices_, tools_stl::seq_int(1, d)))
    throw std::runtime_error(
      "out_vertices must contain numbers 1, ..., cs_dim.");
}

inline void
SVineSelector::check_controls(const FitControlsVinecop& controls)
{
  if (controls.get_select_trunc_lvl()) {
    throw std::runtime_error("Cannot select truncation level for S-vines.");
  }
  if (controls.get_trunc_lvl() < std::numeric_limits<int>::max()) {
    throw std::runtime_error("S-vines cannot be truncated.");
  }
}

// ---------- SVineStructureSelector

SVineStructureSelector::SVineStructureSelector(
  const Eigen::MatrixXd& data,
  const FitControlsVinecop& controls,
  const std::vector<std::string>& var_types)
  : VinecopSelector(data, controls, var_types)
  , SVineSelector(data, var_types)
{
  check_controls(controls);
  out_vertices_.resize(cs_dim_);
  in_vertices_.resize(cs_dim_);
}

std::vector<std::string>
SVineStructureSelector::get_var_types() const
{
  return var_types_;
}

void
SVineStructureSelector::select_all_trees(const Eigen::MatrixXd& data)
{
  if (cs_dim_ > 1) {
    VinecopSelector::select_all_trees(data);
  }
  trees_ = trees_opt_;
}

void
SVineStructureSelector::add_lag()
{
  controls_.set_trunc_lvl(std::numeric_limits<size_t>::max());
  lag_++;
  d_ += cs_dim_;
  data_ = spread_lag(data_, cs_dim_);
  if (controls_.get_weights().size())
    controls_.set_weights(controls_.get_weights().head(data_.rows()));
  auto vt0 = var_types_;
  vt0.resize(cs_dim_);
  var_types_ = tools_stl::cat(var_types_, vt0);
  trees_.resize(d_);
  trees_opt_.resize(d_);

  vine_struct_ = RVineStructure(tools_stl::seq_int(1, d_), 1, false);

  trees_[0] = make_base_tree(data_);
  // add vertices and edges for lagged variable
  for (size_t t = 1; t < trees_.size(); t++) {
    auto old_tree = trees_[t];
    auto new_tree = edges_as_vertices(trees_[t - 1]);

    duplicate_edges(old_tree, new_tree, t);
    trees_opt_[t] = new_tree;
    min_spanning_tree(new_tree);
    add_edge_info(new_tree);
    select_pair_copulas(new_tree, trees_opt_[t]);

    trees_[t] = new_tree;
    if (controls_.get_show_trace()) {
      std::cout << "** Tree: " << t - 1 << std::endl;
      print_pair_copulas_of_tree(t - 1);
    }
  }
}

VineTree
SVineStructureSelector::make_base_tree(const Eigen::MatrixXd& data)
{
  VineTree base_tree(d_);
  auto order = vine_struct_.get_order();
  auto disc_cols = get_disc_cols(var_types_);

  // a star connects the root node (d) with all other nodes
  for (size_t target = 0; target < d_; ++target) {
    tools_interface::check_user_interrupt(target % 10000 == 0);
    // add edge and extract edge iterator
    auto e = add_edge(d_, target, base_tree).first;
    // inititialize hfunc1 with actual data for variable "target"
    // data need are reordered to correspond to natural order (neccessary
    // when structure is fixed)
    base_tree[e].hfunc1 = data.col(order[target] - 1);
    if (var_types_[order[target] - 1] == "d") {
      base_tree[e].hfunc1_sub = data.col(d_ + disc_cols[order[target] - 1]);
      base_tree[e].var_types = { "d", "d" };
    }

    // identify edge with variable "target" and initialize sets
    base_tree[e].conditioned.reserve(2);
    base_tree[e].conditioned.push_back(order[target] - 1);
    base_tree[e].conditioning.reserve(d_ - 2);
    base_tree[e].all_indices = base_tree[e].conditioned;
  }

  return base_tree;
}

void
SVineStructureSelector::duplicate_edges(VineTree& old_tree,
                                        VineTree& new_tree,
                                        size_t t)
{
  // copy existing edges
  for (auto e : boost::edges(old_tree)) {
    size_t v0 = boost::source(e, old_tree);
    size_t v1 = boost::target(e, old_tree);
    auto e_new = boost::add_edge(v0, v1, new_tree).first;
    new_tree[e_new].pair_copula = old_tree[e].pair_copula;
    new_tree[e_new].fit_id = old_tree[e].fit_id;
    new_tree[e_new].conditioned = old_tree[e].conditioned;
    new_tree[e_new].conditioning = old_tree[e].conditioning;
    new_tree[e_new].all_indices = old_tree[e].all_indices;
    new_tree[e_new].var_types = old_tree[e].var_types;
  }

  if (t >= (lag_ - 1) * cs_dim_) {
    if ((cs_dim_ > 1) || (t >= lag_)) {
      add_allowed_connections(new_tree, t);
    }
  }

  // add edges for new lag
  auto shift = [this](std::vector<size_t> index) {
    for (auto& i : index)
      i = i + cs_dim_;
    return index;
  };
  long int ei = 0;
  long int e_start = (lag_ - 1) * cs_dim_ - t;
  long int ei_max = lag_ * cs_dim_ - t;
  for (auto e : boost::edges(new_tree)) {
    if (ei++ < e_start)
      continue;
    if (ei > ei_max)
      continue;

    size_t v0 = boost::source(e, old_tree);
    size_t v1 = boost::target(e, old_tree);
    auto e_new = boost::add_edge(v0 + cs_dim_, v1 + cs_dim_, new_tree).first;
    new_tree[e_new].pair_copula = old_tree[e].pair_copula;
    new_tree[e_new].fit_id = old_tree[e].fit_id;
    new_tree[e_new].conditioned = shift(old_tree[e].conditioned);
    new_tree[e_new].conditioning = shift(old_tree[e].conditioning);
    new_tree[e_new].all_indices = shift(old_tree[e].all_indices);
    new_tree[e_new].var_types = old_tree[e].var_types;
  }
}

void
SVineStructureSelector::add_allowed_connections(VineTree& tree, size_t t)
{
  auto add_edge = [&tree, this](size_t v0, size_t v1) {
    auto crit = std::fabs(wdm::wdm(get_pc_data(v0, v1, tree), "kendall")(0, 1));
    auto w = 1 - crit;
    auto e = boost::add_edge(v0, v1, w, tree).first;
    tree[e].weight = w;
    tree[e].crit = 1 - w;
    tree[e].fit_id = 0.0;
  };
  auto nv = boost::num_vertices(tree);
  if (t == 1) {
    for (size_t v0 = 0; v0 < cs_dim_ - t + 1; ++v0) {
      for (size_t v1 = cs_dim_ - t + 1; v1 < nv; ++v1) {
        add_edge(v0, v1);
      }
    }
  } else {
    long int vstart = 0;
    if (t < cs_dim_) {
      vstart = cs_dim_ - t + 1;
    }
    for (size_t v0 = vstart; v0 < std::min(cs_dim_, nv); ++v0) {
      for (size_t v1 = 0; v1 < nv; ++v1) {
        if (v0 == v1)
          continue;
        find_common_neighbor(v0, v1, tree);
        if (find_common_neighbor(v0, v1, tree) > -1) {
          if (!boost::edge(v1, v0, tree).second) {
            add_edge(v0, v1);
          }
        }
      }
    }
  }
}

void
SVineStructureSelector::finalize_svine(size_t trunc_lvl)
{
  using namespace tools_stl;
  trees_opt_ = trees_;
  pair_copulas_ = make_pair_copula_store(d_, trunc_lvl);
  TriangularArray<size_t> mat(d_, trunc_lvl);
  std::vector<size_t> order(d_);

  if (d_ > cs_dim_) {
    size_t min_c, max_c, min_D, max_D;
    for (size_t t = 0; t < cs_dim_; t++) {
      for (auto e : boost::edges(trees_[t + 1])) {
        min_c = min_vec(trees_[t + 1][e].conditioned);
        max_c = max_vec(trees_[t + 1][e].conditioned);
        if (max_c / cs_dim_ != 1)
          continue;
        if (t > 0) {
          min_D = min_vec(trees_[t + 1][e].conditioning);
          max_D = max_vec(trees_[t + 1][e].conditioning);
          if (max_D / cs_dim_ == 0)
            out_vertices_[t] = min_c;
          if ((min_D / cs_dim_ == 1) && (min_c / cs_dim_ == 0))
            in_vertices_[t] = max_c % cs_dim_;
        } else {
          if (min_c / cs_dim_ == 0) {
            out_vertices_[t] = min_c;
            in_vertices_[t] = max_c % cs_dim_;
          }
        }
      }
    }
  }

  if (trunc_lvl > 0) {
    std::vector<size_t> ned_set;
    std::vector<size_t> ning_set;

    // fill matrix column by column
    for (size_t col = 0; col < d_ - 1; ++col) {
      tools_interface::check_user_interrupt();
      // matrix above trunc_lvl is left empty
      size_t t =
        std::max(std::min(trunc_lvl, d_ - 1 - col), static_cast<size_t>(1));
      // start with highest tree in this column
      for (auto e : boost::edges(trees_[t])) {
        // find an edge that contains a leaf
        size_t v0 = boost::source(e, trees_[t]);
        size_t v1 = boost::target(e, trees_[t]);
        size_t min_deg = std::min(boost::out_degree(v0, trees_[t]),
                                  boost::out_degree(v1, trees_[t]));
        if (min_deg > 1) {
          continue; // not a leaf
        }

        // check if edge contains in_vertex
        ned_set = trees_[t][e].conditioned;
        if ((in_vertices_[t % cs_dim_] != (ned_set[0] % cs_dim_)) &&
            (in_vertices_[t % cs_dim_] != (ned_set[1] % cs_dim_)))
          continue;

        // find position of in_vertex in the edge
        size_t pos = (ned_set[1] % cs_dim_ == in_vertices_[t % cs_dim_]);
        if (pos == 1) {
          trees_[t][e].pair_copula.flip();
        }

        // fill diagonal entry with in vertex
        order[col] = ned_set[pos];

        // entry in row t-1 is other index of the edge
        mat(t - 1, col) = ned_set[1 - pos];

        // assign fitted pair copula to appropriate entry, see
        // `Vinecop::get_pair_copula()`.
        if (trunc_lvl > 0) {
          pair_copulas_[t - 1][col] = trees_[t][e].pair_copula;
        }

        // initialize running set with full conditioning set of this edge
        ning_set = trees_[t][e].conditioning;

        // remove edge (must not be reused in another column!)
        boost::remove_edge(v0, v1, trees_[t]);
        break;
      }

      // fill column bottom to top
      for (size_t k = 1; k < t; ++k) {
        auto check_set = cat(order[col], ning_set);
        for (auto e : boost::edges(trees_[t - k])) {
          // search for an edge in lower tree that shares all
          // indices in the conditioning set + diagonal entry
          if (!is_same_set(trees_[t - k][e].all_indices, check_set)) {
            continue;
          }
          // found suitable edge ->
          // next matrix entry is conditioned variable of new edge
          // that's not equal to the diagonal entry of this column
          auto e_new = trees_[t - k][e];
          ptrdiff_t pos = (order[col] == e_new.conditioned[1]);
          if (pos == 1) {
            e_new.pair_copula.flip();
          }
          mat(t - k - 1, col) = e_new.conditioned[std::abs(1 - pos)];

          // assign fitted pair copula to appropriate entry, see
          // Vinecop::get_pair_copula().
          pair_copulas_[t - 1 - k][col] = e_new.pair_copula;

          // start over with conditioned set of next edge
          ning_set = e_new.conditioning;

          // remove edge (must not be reused in another column!)
          size_t v0 = boost::source(e, trees_[t - k]);
          size_t v1 = boost::target(e, trees_[t - k]);
          boost::remove_edge(v0, v1, trees_[t - k]);
          break;
        }
      }
    }

    // The last column contains a single element which must be different
    // from all other diagonal elements. Based on the properties of an
    // R-vine matrix, this must be the element next to it.

    order[d_ - 1] = mat(0, d_ - 2);

    // change to user-facing format
    // (variable index starting at 1 instead of 0)
    for (size_t i = 0; i < std::min(d_ - 1, trunc_lvl); ++i) {
      for (size_t j = 0; j < d_ - i - 1; ++j) {
        mat(i, j) += 1;
      }
    }
    for (size_t i = 0; i < d_; i++) {
      order[i] += 1;
    }
    for (size_t i = 0; i < cs_dim_; i++) {
      in_vertices_[i] += 1;
      out_vertices_[i] += 1;
    }

  } else {
    // order doesn't matter for truncated
    order = tools_stl::seq_int(1, d_);
  }

  // return as RVineStructure
  vine_struct_ = RVineStructure(order, mat);
}

void
SVineStructureSelector::finalize(size_t trunc_lvl)
{
  if (d_ == cs_dim_) {
    trees_ = trees_opt_;
    auto mat = vine_struct_.get_matrix();
    cs_struct_ = RVineStructure(mat.block(0, d_ - cs_dim_, cs_dim_, cs_dim_).eval());
    in_vertices_ = tools_stl::rev(cs_struct_.get_order());
    out_vertices_ = in_vertices_;
  } else {
    finalize_svine(trunc_lvl);
    auto mat = vine_struct_.get_matrix();
    cs_struct_ = RVineStructure(mat.block(0, d_ - cs_dim_, cs_dim_, cs_dim_).eval());
  }
}

double
SVineStructureSelector::compute_fit_id(const EdgeProperties& e)
{
  size_t min_c = std::min(e.conditioned[0], e.conditioned[1]);
  size_t max_c = std::max(e.conditioned[0], e.conditioned[1]);
  size_t lagdiff = (max_c - min_c) / cs_dim_;
  return (min_c % cs_dim_) + (max_c % cs_dim_) * std::pow(cs_dim_, 2) +
         lagdiff * std::pow(cs_dim_, 4);
}

// ------------- SVineFamilySelector

SVineFamilySelector::SVineFamilySelector(
  const Eigen::MatrixXd& data,
  const RVineStructure& cs_struct,
  const FitControlsVinecop& controls,
  std::vector<size_t> out_vertices,
  std::vector<size_t> in_vertices,
  const std::vector<std::string>& var_types)
  : VinecopSelector(data, cs_struct, controls, var_types)
  , SVineSelector(data, out_vertices, in_vertices, var_types)
{
  check_controls(controls);
  cs_struct_ = SVineStructure(cs_struct, 0, out_vertices, in_vertices);
}

inline std::vector<std::string>
SVineFamilySelector::get_var_types() const
{
  return var_types_;
}

inline void
SVineFamilySelector::select_tree(size_t t)
{
  auto new_tree = edges_as_vertices(trees_[t]);
  remove_edge_data(trees_[t]);
  add_allowed_edges(new_tree);
  if (boost::num_vertices(new_tree) > 0) {
    add_edge_info(new_tree); // for pc estimation and next tree
    if (controls_.get_selection_criterion() == "mbicv") {
      // adjust prior probability to tree level
      controls_.set_psi0(std::pow(psi0_, t + 1));
    }
    if (trees_opt_.size() > t + 1) {
      select_pair_copulas(new_tree, trees_opt_[t + 1]);
    } else {
      select_pair_copulas(new_tree);
    }
  }
  // make sure there is space for new tree
  trees_.resize(t + 2);
  trees_[t + 1] = new_tree;
}

inline void
SVineFamilySelector::add_lag()
{
  lag_++;
  d_ += cs_dim_;
  auto vt0 = var_types_;
  vt0.resize(cs_dim_);
  var_types_ = tools_stl::cat(var_types_, vt0);

  // add vertices and edges for lagged variable
  for (size_t t = 1; t < trees_.size(); t++) {
    auto old_tree = trees_[t];
    for (auto v : boost::vertices(old_tree))
      duplicate_vertex(v, trees_[t]);
    
    long int ei = 0;
    long int e_start = (lag_ - 1) * cs_dim_ - t;
    long int ei_max = lag_ * cs_dim_ - t;
    for (auto e : boost::edges(old_tree)) {
      if (ei++ < e_start)
        continue;
      if (ei > ei_max)
        continue;
      duplicate_edge(e, trees_[t]);
    }
  }

  // update trees and structure
  trees_opt_ = trees_;
  trees_ = std::vector<VineTree>(1);
  vine_struct_ = SVineStructure(cs_struct_, lag_, out_vertices_, in_vertices_);
  data_ = spread_lag(data_, cs_dim_);
  if (controls_.get_weights().size())
    controls_.set_weights(controls_.get_weights().head(data_.rows()));
  controls_.set_trunc_lvl(std::numeric_limits<size_t>::max());
}

inline Eigen::MatrixXd
SVineFamilySelector::data()
{
  return data_;
}

inline void SVineFamilySelector::finalize(size_t)
{
  pair_copulas_ = make_pair_copula_store(d_, d_);
  std::vector<size_t> order = vine_struct_.get_order();
  size_t d = order.size();
  auto struct_array = vine_struct_.get_struct_array();
  for (size_t t = 1; t < trees_.size(); t++) {
    for (auto e : boost::edges(trees_[t])) {
      std::vector<size_t> check_set(2);
      for (size_t k = 0; k < d; k++) {
        // conditioned starts at 0 -> substract -1
        check_set = { order[k] - 1,
                      static_cast<size_t>(struct_array(t - 1, k) - 1) };
        if (tools_stl::is_same_set(trees_[t][e].conditioned, check_set)) {
          pair_copulas_[t - 1][k] = trees_[t][e].pair_copula;
          break;
        }
      }
    }
  }
}

inline double
SVineFamilySelector::compute_fit_id(const EdgeProperties& e)
{
  size_t min_c = std::min(e.conditioned[0], e.conditioned[1]);
  size_t max_c = std::max(e.conditioned[0], e.conditioned[1]);
  size_t lagdiff = (max_c - min_c) / cs_dim_;
  return (min_c % cs_dim_) + (max_c % cs_dim_) * std::pow(cs_dim_, 2) +
         lagdiff * std::pow(cs_dim_, 4);
}

inline void
SVineFamilySelector::flip_pcs(std::vector<VineTree>& trees)
{
  std::vector<size_t> order = cs_struct_.get_order();
  size_t d = order.size();
  auto struct_array = cs_struct_.get_struct_array();
  for (size_t t = 1; t < trees.size(); t++) {
    for (auto e : boost::edges(trees[t])) {
      std::vector<size_t> check_set(2);
      for (size_t k = 0; k < d; k++) {
        // conditioned starts at 0 -> substract -1
        check_set = { order[k] - 1,
                      static_cast<size_t>(struct_array(t - 1, k) - 1) };
        if (tools_stl::is_same_set(trees_[t][e].conditioned, check_set)) {
          break;
        }
      }
      if (trees[t][e].conditioned[0] != check_set[0]) {
        trees[t][e].pair_copula.flip();
      }
    }
  }
}

}
}
