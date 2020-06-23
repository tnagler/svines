#pragma once

#include <svines/tools_stl.hpp>
#include <vinecopulib/misc/tools_stl.hpp>
#include <vinecopulib/vinecop/rvine_structure.hpp>

namespace vinecopulib {

//! Default constructor
SVineStructure::SVineStructure()
  : RVineStructure()
{}

//! @brief Creates an S-vine structure.
//!
//! @param cs_dim cross-sectional dimension.
//! @param p Markov order.
SVineStructure::SVineStructure(size_t cs_dim, size_t p)
  : SVineStructure(RVineStructure(tools_stl::seq_int(1, cs_dim)),
                   p,
                   tools_stl::seq_int(1, cs_dim),
                   tools_stl::seq_int(1, cs_dim))
{}

//! @brief Creates an S-vine structure.
//!
//! @param cs_struct cross-sectional structure.
//! @param p Markov order.
//! @param out_vertices out-vertices.
//! @param in_vertices in-vertices.
SVineStructure::SVineStructure(const RVineStructure& cs_struct,
                               size_t p,
                               std::vector<size_t> out_vertices,
                               std::vector<size_t> in_vertices)
  : p_(p)
  , out_vertices_(out_vertices)
  , in_vertices_(in_vertices)
{
  check_out_in_vertices(cs_struct, out_vertices, in_vertices);
  cs_struct_ = reorder_structure(cs_struct, in_vertices);
  order_ = expand_order(cs_struct_.get_order(), p);
  struct_array_ = build_s_vine_array(cs_struct_, p, out_vertices, in_vertices);

  RVineStructure new_struct;
  try {
    new_struct = RVineStructure(order_, struct_array_);
  } catch (const std::exception& e) {
    throw std::runtime_error(
      "out_vertices are not compatible with cross-sectional structure.");
  }
  d_ = new_struct.get_dim();
  trunc_lvl_ = new_struct.get_trunc_lvl();
  struct_array_ = new_struct.get_struct_array(true);
  min_array_ = new_struct.get_min_array();
  needed_hfunc1_ = new_struct.get_needed_hfunc1();
  needed_hfunc2_ = new_struct.get_needed_hfunc2();
}

//! @brief Gets the Markov order.
inline size_t
SVineStructure::get_p() const
{
  return p_;
}

//! @brief Gets the cross-sectional dimension.
inline size_t
SVineStructure::get_cs_dim() const
{
  return cs_struct_.get_dim();
}


//! @brief Gets the out-vertices.
inline std::vector<size_t>
SVineStructure::get_out_vertices() const
{
  return out_vertices_;
}

//! @brief Gets the in-vertices.
inline std::vector<size_t>
SVineStructure::get_in_vertices() const
{
  return in_vertices_;
}

//! @brief Gets the cross-sectional structure
inline RVineStructure
SVineStructure::get_cs_structure() const
{
  return cs_struct_;
}

//! @brief Checks whether in- and out-vertices are compatible.
//!
//! @param cs_struct cross-sectional structure.
//! @param in_vertices in-vertices.
//! @param out_vertices out-vertices.
inline void
SVineStructure::check_out_in_vertices(const RVineStructure& cs_struct,
                                      std::vector<size_t> out_vertices,
                                      std::vector<size_t> in_vertices) const
{
  auto d = cs_struct.get_dim();
  if (!tools_stl::is_same_set(in_vertices, tools_stl::seq_int(1, d)))
    throw std::runtime_error(
      "in_vertices must contain numbers 1, ..., cs_dim.");
  if (!tools_stl::is_same_set(out_vertices, tools_stl::seq_int(1, d)))
    throw std::runtime_error(
      "out_vertices must contain numbers 1, ..., cs_dim.");
}

//! @brief Expands the cross-sectional order vector to that of a p-Markov model.
//!
//! @param order cross-sectional order.
//! @param p Markov order.
inline std::vector<size_t>
SVineStructure::expand_order(const std::vector<size_t>& order, size_t p) const
{
  size_t cs_dim = order.size();
  size_t d = cs_dim * (p + 1);
  std::vector<size_t> new_order(d);
  for (size_t i = 0; i < d; i++) {
    new_order[i] = order[i % cs_dim] + ((d - 1 - i) / cs_dim) * cs_dim;
  }

  return new_order;
}

//! @brief Initial guess for new column in reordered cross-sectional structure 
//! (see Daniel's thesis).
//!
//! @param old_diag old diagonal elements.
//! @param old_struct old cross-sectional structure array.
//! @param new_el new diagonal element. 
//! @param column the column to fill.
inline std::vector<size_t>
SVineStructure::sup_diag(const std::vector<size_t>& old_diag,
                         const TriangularArray<size_t>& old_struct,
                         size_t new_el,
                         size_t column) const
{
  size_t d = old_diag.size();
  std::vector<size_t> x(d - 1);

  size_t i = 0;
  while (old_diag[i] != new_el) {
    x[i] = old_diag[i];
    i++;
  }

  size_t pivot_col = i;
  while (i < d - 1) {
    x[i] = old_struct(d - 2 - i, pivot_col);
    i++;
  }

  return tools_stl::span(tools_stl::rev(x), 0, d - 1 - column);
}

//! @brief Reorders a cross-sectional structure so it can be used with a 
//! compatible vector of in-vertices.
//!
//! @param structure a vine structure.
//! @param in_vertices a compatible vector of in-vertices.
inline RVineStructure
SVineStructure::reorder_structure(const RVineStructure& structure,
                                  std::vector<size_t> in_vertices) const
{
  using namespace tools_stl;
  size_t d = structure.get_dim();
  if (structure.get_trunc_lvl() < d - 1) {
    throw std::runtime_error("S-vines cannot be truncated.");
  }

  auto old_struct = structure.get_struct_array();
  auto new_struct = old_struct;

  // prepare objects
  auto old_order = structure.get_order();
  auto new_order = rev(in_vertices);

  // loop through all columns
  for (size_t i = 0; i < d - 1; i++) {
    auto new_column = sup_diag(old_order, old_struct, new_order[i], i);
    auto diag_until = span(new_order, 0, i);
    auto diag_after = span(new_order, i, d - i);
    auto max_col = d - 1; // find_position(new_order[i], old_order);
    for (size_t t = 0; t < new_column.size(); t++) {
      // Check whether an element in this column is already contained in
      // the diagonal to the left. If so, we need to find another node
      // that is connected to the diagonal entry of column i. We search
      // for such an edge in the old structure, but only to the left of
      // the column, where the element appeared on the diagonal.
      if (is_member(new_column[t], diag_until)) {
        bool found_node = false;
        for (size_t j = 0; j <= max_col; j++) {
          if (new_order[i] == old_struct(t, j)) {
            if (is_member(old_order[j], diag_after)) {
              new_column[t] = old_order[j];
              found_node = true;
            }
          } else if (new_order[i] == old_order[j]) {
            if (is_member(old_struct(t, j), diag_after)) {
              new_column[t] = old_struct(t, j);
              found_node = true;
            }
          }
          for (size_t k = 0; k < new_column.size(); ++k)
            new_struct(k, i) = new_column[k];
          if (found_node) {
            // The new entry may already be contained in this
            // column. We need to check the next rows for that
            // as well.
            diag_until = cat(new_column[t], diag_until);
            break;
          }
        }
      }
    }
    for (size_t k = 0; k < new_column.size(); ++k)
      new_struct(k, i) = new_column[k];
  }

  // this must always hold beacuse the first in vertex comes last on the
  // diagonal:
  if (d > 1)
    new_struct(0, d - 2) = new_order[d - 1];

  RVineStructure new_rvine;
  try {
    new_rvine = RVineStructure(new_order, new_struct);
  } catch (const std::exception& e) {
    throw std::runtime_error(
      "in_vertices are not compatible with cross-sectional structure.");
  }

  return new_rvine;
}


//! @brief Builds a S-vine array with given cross-section and in-/out-vertices.
//!
//! @param cs_struct cross-sectional structure.
//! @param p Markov order.
//! @param in_vertices in-vertices.
//! @param out_vertices out-vertices.
inline TriangularArray<size_t>
SVineStructure::build_s_vine_array(const RVineStructure& cs_struct,
                                   size_t p,
                                   std::vector<size_t> out_vertices,
                                   std::vector<size_t> in_vertices) const
{
  size_t cs_dim = cs_struct.get_dim();
  size_t d = cs_dim * (p + 1);
  auto diag = cs_struct.get_order();

  RVineStructure new_struct = cs_struct;
  if (diag[cs_dim - 1] != in_vertices[0]) {
    new_struct = reorder_structure(new_struct, in_vertices);
  }

  auto struct_array = cs_struct.get_struct_array();

  TriangularArray<size_t> strct(d);
  // copy cross-sectional structure
  for (size_t i = 0; i < cs_dim - 1; i++) {
    for (size_t j = 0; j < cs_dim - 1 - i; j++) {
      strct(i, j) = struct_array(i, j) + cs_dim * p;
    }
  }

  // fill parallelograms
  std::vector<size_t> par = out_vertices;
  for (size_t lag = 1; lag <= p; lag++) {
    for (size_t i = 0; i < cs_dim; i++) {
      for (size_t j = 0; j < cs_dim; j++) {
        strct(i + cs_dim * lag - j - 1, j) = par[i] + cs_dim * (p - lag);
      }
    }
  }

  // copy to other lags
  for (size_t lag = 1; lag <= p; lag++) {
    for (size_t j = 0; j < cs_dim; j++) {
      for (size_t i = 0; i < d - 1 - (j + cs_dim * lag); i++) {
        strct(i, j + cs_dim * lag) = strct(i, j) - cs_dim * lag;
      }
    }
  }

  return strct;
}

}