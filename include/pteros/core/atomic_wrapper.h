/*
 *
 *                This source code is part of
 *                    ******************
 *                    ***   Pteros   ***
 *                    ******************
 *                 molecular modeling library
 *
 * Copyright (c) 2009-2013, Semen Yesylevskyy
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of Artistic License:
 *
 * Please note, that Artistic License is slightly more restrictive
 * then GPL license in terms of distributing the modified versions
 * of this software (they should be approved first).
 * Read http://www.opensource.org/licenses/artistic-license-2.0.php
 * for details. Such license fits scientific software better then
 * GPL because it prevents the distribution of bugged derivatives.
 *
*/

#ifndef ATOMIC_WRAPPER_H_INCLUDED
#define ATOMIC_WRAPPER_H_INCLUDED

#include <atomic>

// Wrapper over std::atomic to make it usable in std::vector
template <typename T>
struct atomic_wrapper {
  std::atomic<T> _a;

  atomic_wrapper():_a(){}

  atomic_wrapper(const std::atomic<T> &a):_a(a.load()){}

  atomic_wrapper(const atomic_wrapper &other):_a(other._a.load()){}

  atomic_wrapper &operator=(const atomic_wrapper &other){
    _a.store(other._a.load());
  }

  T load(){
      return _a.load();
  }

  void store(const T& v){
      _a.store(v);
  }
};

#endif
