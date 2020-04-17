/*
 * This file is a part of
 *
 * ============================================
 * ###   Pteros molecular modeling library  ###
 * ============================================
 *
 * https://github.com/yesint/pteros
 *
 * (C) 2009-2020, Semen Yesylevskyy
 *
 * All works, which use Pteros, should cite the following papers:
 *  
 *  1.  Semen O. Yesylevskyy, "Pteros 2.0: Evolution of the fast parallel
 *      molecular analysis library for C++ and python",
 *      Journal of Computational Chemistry, 2015, 36(19), 1480–1488.
 *      doi: 10.1002/jcc.23943.
 *
 *  2.  Semen O. Yesylevskyy, "Pteros: Fast and easy to use open-source C++
 *      library for molecular analysis",
 *      Journal of Computational Chemistry, 2012, 33(19), 1632–1636.
 *      doi: 10.1002/jcc.22989.
 *
 * This is free software distributed under Artistic License:
 * http://www.opensource.org/licenses/artistic-license-2.0.php
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


