// =========================================================================
// File:       tledParallel.cpp
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

#include "tledParallel.h"

#ifdef _USE_BOOST_

void tledParallel::JoinThreads(std::vector<boost::thread*> &rvp_threads) const {
  for (std::vector<boost::thread*>::iterator ip_t = rvp_threads.begin(); ip_t < rvp_threads.end(); ip_t++) {
    (*ip_t)->join();
    delete *ip_t;
  }
  rvp_threads.clear();
}


static void _MergeWorker(float *p_out, const float *pc_in0Begin, const float *pc_in0End, const float *pc_in1Begin, const float *pc_in1End) {
  float *p_dst;
  float const *pc_src0, *pc_src1;

  assert(pc_in0End - pc_in0Begin == pc_in1End - pc_in1Begin);
  for (p_dst = p_out, pc_src0 = pc_in0Begin, pc_src1 = pc_in1Begin; pc_src0 < pc_in0End; pc_src0++, pc_src1++, p_dst++) {
    *p_dst = (*pc_src0) + (*pc_src1);
  }
}

void tledParallel::MergeResults(float *p_merged, float *p_threadResults, const int numElements, const int numComps) const {
  const int resultBlockSize = numComps*numElements;

  int oldBlockSize = this->GetNumberOfThreads();

  while (oldBlockSize > 2) {
    int mergeBlock = oldBlockSize/2 + (oldBlockSize%2 > 0);
    std::vector<boost::thread*> vp_threads;    
    int threadsPerBlock;

    threadsPerBlock = this->GetNumberOfThreads()/(mergeBlock - (2*mergeBlock - 1 >= oldBlockSize));
    assert(threadsPerBlock > 1 && threadsPerBlock < this->GetNumberOfThreads());
    for (int i = 0; i < mergeBlock; i++) {
      if (i + mergeBlock < oldBlockSize) {
	float *p_dst = p_threadResults + i*resultBlockSize;
	float const *pc_src = p_threadResults + (i + mergeBlock)*resultBlockSize;

	for (int t = 0; t < threadsPerBlock; t++) {
	  const int offStart =  (t*resultBlockSize)/threadsPerBlock;
	  const int offEnd = ((t + 1)*resultBlockSize)/threadsPerBlock;

	  assert(offStart >= 0 && offStart < offEnd && offEnd <= resultBlockSize);
	  assert(t + 1 < threadsPerBlock || offEnd == resultBlockSize);
	  assert(t > 0 || offStart == 0);
	  vp_threads.push_back(new boost::thread(_MergeWorker, p_dst + offStart, p_dst + offStart, p_dst + offEnd, pc_src + offStart, pc_src + offEnd));					       
	}
      }
    }
    assert((int)vp_threads.size() <= this->GetNumberOfThreads());

    this->JoinThreads(vp_threads);
    oldBlockSize = mergeBlock;
  }

  assert(oldBlockSize == 2);

  {
    std::vector<boost::thread*> vp_threads;    
    float const *pc_src0 = p_threadResults, *pc_src1 = p_threadResults + resultBlockSize;

    for (int t = 0; t < this->GetNumberOfThreads(); t++) {
      const int offStart = (t*resultBlockSize)/this->GetNumberOfThreads();
      const int offEnd = ((t + 1)*resultBlockSize)/this->GetNumberOfThreads();

      vp_threads.push_back(new boost::thread(_MergeWorker, p_merged + offStart, pc_src0 + offStart, pc_src0 + offEnd, pc_src1 + offStart, pc_src1 + offEnd));
    }

    this->JoinThreads(vp_threads);
  }
}



#endif
