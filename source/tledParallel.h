// =========================================================================
// File:       tledParallel.h
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    February 2013
// 
// Copyright (c) 2013, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================
#ifndef tledParallel_H
#define tledParallel_H

#ifdef BOOST_DISABLE_THREADS
#undef _USE_BOOST_
#endif

#ifdef _USE_BOOST_
#include <boost/thread.hpp>

/**
 * \brief An API for CPU-parallel algorithms
 */
class tledParallel {
  /**
   * \name Thread Number Control
   * @{
   */
private:
  int m_NumThreads;

public:
  /** @{ */
  /** Number of threads used. Defaults to value of boost::thread::hardware_concurrency() */
  int GetNumberOfThreads(void) const { return m_NumThreads; }
  void SetNumberOfThreads(const int numThreads) { m_NumThreads = numThreads; }
  /** @} */

  /** Number of available CPU cores */
  static int GetNumberOfHardwareThreads(void) { return boost::thread::hardware_concurrency(); }
  /** @} */

  /**
   * \name Computation
   * @{
   */
protected:  
  /** Joins all thread objects in the vector, and destroys them */
  void JoinThreads(std::vector<boost::thread*> &rvp_threads) const;

  /** Merges results obtained by threads */
  void MergeResults(float *p_merged, float *p_threadResults, const int numElements, const int numComps) const;

  /** Splits a set evenly between the threads of the object */
  int GetThreadBlockSize(const int totalSize) const { return totalSize/this->GetNumberOfThreads() + (totalSize%this->GetNumberOfThreads() > 0); }
  /** @} */

  /**
   * \name Construction, Destruction
   * @{
   */
protected:
  /** Sets number of threads to the number of available cores */
  tledParallel(void) : m_NumThreads(GetNumberOfHardwareThreads()) {}
  tledParallel(const int numThreads) : m_NumThreads(numThreads) {}

public:
  virtual ~tledParallel(void) {}
  /** @} */
};

#endif
#endif
