/* XTHELPERS - small convenience factories for allocating xtensor arrays
   with a runtime-determined shape, used throughout PolConvert to replace
   the old new[]/malloc()-based dynamic arrays.

             Copyright (C) 2013-2022  Ivan Marti-Vidal
             Nordic Node of EU ALMA Regional Center (Onsala, Sweden)
             Max-Planck-Institut fuer Radioastronomie (Bonn, Germany)
             University of Valencia (Spain)

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
*/

#ifndef __XTHELPERS_H__
#define __XTHELPERS_H__

#include <algorithm>
#include <cstddef>
#include <xtensor/xtensor.hpp>

// Allocates a 1-D xtensor of length n (uninitialized contents, same as new T[n]).
// Since xt::xtensor<T,1>::operator[](i) forwards to operator()(i), the result
// is a drop-in replacement for a `T*` obtained via `new T[n]` at all read/write
// call sites of the form arr[i].
template <typename T>
inline xt::xtensor<T,1> xt_alloc(std::size_t n) {
  return xt::xtensor<T,1>(typename xt::xtensor<T,1>::shape_type({n}));
}

// Allocates a 2-D xtensor of shape (n0, n1) (uninitialized contents).
template <typename T>
inline xt::xtensor<T,2> xt_alloc(std::size_t n0, std::size_t n1) {
  return xt::xtensor<T,2>(typename xt::xtensor<T,2>::shape_type({n0, n1}));
}

// Allocates a 3-D xtensor of shape (n0, n1, n2) (uninitialized contents).
template <typename T>
inline xt::xtensor<T,3> xt_alloc(std::size_t n0, std::size_t n1, std::size_t n2) {
  return xt::xtensor<T,3>(typename xt::xtensor<T,3>::shape_type({n0, n1, n2}));
}

// Allocates a 4-D xtensor of shape (n0, n1, n2, n3) (uninitialized contents).
template <typename T>
inline xt::xtensor<T,4> xt_alloc(std::size_t n0, std::size_t n1, std::size_t n2, std::size_t n3) {
  return xt::xtensor<T,4>(typename xt::xtensor<T,4>::shape_type({n0, n1, n2, n3}));
}

// Grows (or shrinks) a 1-D xtensor to newSize, preserving the first
// min(oldSize, newSize) elements. This is a drop-in replacement for the
// old `arr = (T*) realloc(arr, newSize*sizeof(T));` growable-buffer pattern.
template <typename T>
inline void xt_grow(xt::xtensor<T,1>& arr, std::size_t newSize) {
  xt::xtensor<T,1> tmp = xt_alloc<T>(newSize);
  std::size_t copyN = std::min(arr.size(), newSize);
  if (copyN > 0) {
    std::copy(arr.data(), arr.data() + copyN, tmp.data());
  }
  arr = std::move(tmp);
}

#endif
