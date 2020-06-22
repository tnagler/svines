#pragma once
#include <algorithm>
#include <vector>
#include <vinecopulib/misc/tools_stl.hpp>

namespace vinecopulib {

namespace tools_stl {

//! @brief Returns a segment of a vector.
template<class T>
std::vector<T>
span(std::vector<T> x, size_t start, size_t len)
{
  x.erase(x.begin(), x.begin() + std::min(x.size(), start));
  x.resize(std::min(x.size(), len));
  return x;
}

//! Reverses the order of elements in a vector.
template<class T>
std::vector<T>
rev(std::vector<T> x)
{
  tools_stl::reverse(x);
  return x;
}

//! @brief Repeats a vector multiple times.
template<class T>
std::vector<T>
rep(std::vector<T> x, size_t times)
{
  std::vector<T> y = x;
  for (size_t t = 1; t < times; ++t)
    y = cat(y, x);
  return y;
}

//! @brief Computes the minimum element in a vector.
template<typename T>
T
min_vec(std::vector<T>& x)
{
  auto i = std::min_element(x.begin(), x.end());
  return x[std::distance(x.begin(), i)];
}

//! @brief Computes the maximum element in a vector.
template<typename T>
T
max_vec(std::vector<T>& x)
{
  auto i = std::max_element(x.begin(), x.end());
  return x[std::distance(x.begin(), i)];
}

}
}