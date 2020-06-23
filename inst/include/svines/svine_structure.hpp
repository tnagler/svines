
#pragma once

#include <vinecopulib/vinecop/rvine_structure.hpp>

namespace vinecopulib {

class SVineStructure : public RVineStructure
{
public:
  SVineStructure();

  SVineStructure(size_t cs_dim, size_t p);

  SVineStructure(const RVineStructure& cs_struct,
                 size_t p,
                 std::vector<size_t> out_vertices,
                 std::vector<size_t> in_vertices);

  size_t get_p() const;

  size_t get_cs_dim() const;

  std::vector<size_t> get_out_vertices() const;

  std::vector<size_t> get_in_vertices() const;

  RVineStructure get_cs_structure() const;

private:
  void check_out_in_vertices(const RVineStructure& cs_struct,
                             std::vector<size_t> out_vertices,
                             std::vector<size_t> in_vertices) const;

  std::vector<size_t> expand_order(const std::vector<size_t>& order,
                                   size_t p) const;

  std::vector<size_t> sup_diag(const std::vector<size_t>& old_diag,
                               const TriangularArray<size_t>& old_struct,
                               size_t new_el,
                               size_t column) const;

    RVineStructure reorder_structure(const RVineStructure& structure,
                                     std::vector<size_t> in_vertices) const;

    TriangularArray<size_t> build_s_vine_array(
      const RVineStructure& cs_struct,
      size_t p,
      std::vector<size_t> out_vertices,
      std::vector<size_t> in_vertices) const;

  private:
    size_t p_;
    std::vector<size_t> out_vertices_;
    std::vector<size_t> in_vertices_;
    RVineStructure cs_struct_;
};

}

#include <svines/implementation/svine_structure.ipp>