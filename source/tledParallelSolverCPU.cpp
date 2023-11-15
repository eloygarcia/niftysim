// =========================================================================
// File:       tledParallelSolverCPU.cpp
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    January 2013
// 
// Copyright (c) 2013, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================
#include "tledParallelSolverCPU.h"

#ifdef _USE_BOOST_
#include "tledShellSolverCPU.h"

#include <algorithm>
#include <vector>
#include <cassert>

tledParallelSolverCPU::~tledParallelSolverCPU() {
  delete[] mp_ThreadResultBuffer;
}

void tledParallelSolverCPU::Init(tledModel* p_model) {
  tledSolverCPU::Init(p_model);
  mp_ThreadResultBuffer = new float[this->GetNumberOfDOFs()*this->GetNumberOfNodes()*this->GetNumberOfThreads()];
}


static void _ElementForceWorker(float *p_F, tledElement **pp_elBegin, tledElement **pp_elEnd, const float U[], const int szOut) {
  std::fill(p_F, p_F + szOut, 0.0f);
  for (tledElement **pp_el = pp_elBegin; pp_el < pp_elEnd; pp_el++) {
    (*pp_el)->ComputeElementForces(U, p_F);
  }  
}

void tledParallelSolverCPU::ComputeNewForces() {
  if (this->GetNumberOfElements() > 0) {
    const int elBlockSize = this->GetThreadBlockSize(this->GetNumberOfElements());

    float *p_threadForces = mp_ThreadResultBuffer;
    float *p_currThreadForces = p_threadForces;
    std::vector<boost::thread*> vp_threads;

    for (tledElement **pp_elBegin = this->GetElements(); pp_elBegin < this->GetElements() + this->GetNumberOfElements();) {
      tledElement **pp_elEnd = std::min(pp_elBegin + elBlockSize, this->GetElements() + this->GetNumberOfElements());

      vp_threads.push_back(new boost::thread(_ElementForceWorker, p_currThreadForces, pp_elBegin, pp_elEnd, this->GetTimeStepper().GetCurrentDisplacements(), this->GetNumberOfDOFs()*this->GetNumberOfNodes()));
      assert((int)vp_threads.size() <= this->GetNumberOfThreads());

      pp_elBegin = pp_elEnd;
      p_currThreadForces += this->GetNumberOfDOFs()*this->GetNumberOfNodes();
    }

    this->JoinThreads(vp_threads);
    this->MergeResults(this->GetInternalForceBuffer(), p_threadForces, this->GetNumberOfNodes(), this->GetNumberOfDOFs());
  } else {
    std::fill(this->GetInternalForceBuffer(), this->GetInternalForceBuffer() + this->GetNumberOfDOFs()*this->GetNumberOfNodes(), 0.0f);    
  }
}

static void _ANPElementPressureWorker(float *p_Pa, tledElement **pp_elBegin, tledElement **pp_elEnd, const float U[], const int szOut) {
  std::fill(p_Pa, p_Pa + szOut, 0.0f);
  for (tledElement **pp_el = pp_elBegin; pp_el < pp_elEnd; pp_el++) {
    (*pp_el)->ComputeElementPressure(U, p_Pa);
  }  
}

static void _ANPElementForceWorker(float *p_F, tledElement **pp_elBegin, tledElement **pp_elEnd, const float Pa[], const float Va[], const int szOut) {
  std::fill(p_F, p_F + szOut, 0.0f);
  for (tledElement **pp_el = pp_elBegin; pp_el < pp_elEnd; pp_el++) {
    (*pp_el)->ComputeModifiedElementForces(Pa, Va, p_F);
  }  
}

void tledParallelSolverCPU::ComputeNewForcesANP() {
  if (this->GetNumberOfElements() > 0) {
    const int elBlockSize = this->GetThreadBlockSize(this->GetNumberOfElements());

    std::vector<boost::thread*> vp_threads;

    {
      float *p_threadPressures = mp_ThreadResultBuffer;
      float *p_currThreadPressures = p_threadPressures;
    
      for (tledElement **pp_elBegin = this->GetElements(); pp_elBegin < this->GetElements() + this->GetNumberOfElements();) {
	tledElement **pp_elEnd = std::min(pp_elBegin + elBlockSize, this->GetElements() + this->GetNumberOfElements());

	vp_threads.push_back(new boost::thread(_ANPElementPressureWorker, p_currThreadPressures, pp_elBegin, pp_elEnd, this->GetTimeStepper().GetCurrentDisplacements(), this->GetNumberOfNodes()));
	assert((int)vp_threads.size() <= this->GetNumberOfThreads());

	pp_elBegin = pp_elEnd;
	p_currThreadPressures += this->GetNumberOfNodes();
      }

      this->JoinThreads(vp_threads);
      this->MergeResults(this->GetPressureBuffer(), p_threadPressures, this->GetNumberOfNodes(), 1);
    }

    {
      float *p_threadForces = mp_ThreadResultBuffer;
      float *p_currThreadForces = p_threadForces;
      std::vector<boost::thread*> vp_threads;

      for (tledElement **pp_elBegin = this->GetElements(); pp_elBegin < this->GetElements() + this->GetNumberOfElements();) {
	tledElement **pp_elEnd = std::min(pp_elBegin + elBlockSize, this->GetElements() + this->GetNumberOfElements());

	vp_threads.push_back(new boost::thread(_ANPElementForceWorker, p_currThreadForces, pp_elBegin, pp_elEnd, this->GetPressureBuffer(), this->GetNodeVolumes(), this->GetNumberOfDOFs()*this->GetNumberOfNodes()));
	assert((int)vp_threads.size() <= this->GetNumberOfThreads());

	pp_elBegin = pp_elEnd;
	p_currThreadForces += this->GetNumberOfDOFs()*this->GetNumberOfNodes();
      }

      this->JoinThreads(vp_threads);
      this->MergeResults(this->GetInternalForceBuffer(), p_threadForces, this->GetNumberOfNodes(), this->GetNumberOfDOFs());
    }
  } else {
    std::fill(this->GetPressureBuffer(), this->GetPressureBuffer() + this->GetNumberOfNodes(), 0.0f);
    std::fill(this->GetInternalForceBuffer(), this->GetInternalForceBuffer() + this->GetNumberOfDOFs()*this->GetNumberOfNodes(), 0.0f);
  }
}

#endif
