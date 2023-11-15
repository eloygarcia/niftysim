// =========================================================================
// File:       tledParallelShellSolverCPU.cpp
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
#include "tledParallelShellSolverCPU.h"
#include "tledElementShellBSTP1.h"
#include "tledElementMembraneSimpleLinear.h"
#include "tledElementMembraneNonLinear.h"

#include <cassert>
#include <vector>

#ifdef _USE_BOOST_
template <class TShellElement>
static void _ElementForceWorker(float *p_F, typename std::vector<TShellElement>::iterator i_elsBegin, typename std::vector<TShellElement>::iterator i_elsEnd, const float u[], const int numNodes) {
  std::fill(p_F, p_F + 3*numNodes, 0.0f);
  for (typename std::vector<TShellElement>::iterator i_el = i_elsBegin; i_el < i_elsEnd; i_el++) {
    i_el->ComputeElementForces(p_F, u);
  }
}

static void _AddWorker(float *p_outStart, float *p_outEnd, const float *pc_addStart) {
  float const *pc_src = pc_addStart;

  for (float *p_dst = p_outStart; p_dst < p_outEnd; p_dst++, pc_src++) *p_dst += *pc_src;
}

template <class TShellElement>
void tledParallelShellSolverCPU::ElementSetImpl<TShellElement>::ComputeForces(float *p_F, const float u[]) {
  std::vector<boost::thread*> vp_threads;

  {
    const int numThreadEls = this->GetThreadBlockSize(this->m_Elements.size());

    float *p_tf = this->mp_ThreadResultBuffer;

    assert(this->mp_ThreadResultBuffer != NULL);
    for (typename std::vector<TShellElement>::iterator i_el = this->m_Elements.begin(); i_el < this->m_Elements.end(); i_el += numThreadEls, p_tf += 3*this->m_NumNodes) {
      vp_threads.push_back(new boost::thread(_ElementForceWorker<TShellElement>, p_tf, i_el, std::min(this->m_Elements.end(), i_el + numThreadEls), u, this->m_NumNodes));
      assert((int)vp_threads.size() <= this->GetNumberOfThreads());
    }
    this->JoinThreads(vp_threads);
    this->MergeResults(this->mp_ThreadResultBuffer, this->mp_ThreadResultBuffer, this->m_NumNodes, 3);
  }
  
  {
    const int szThreadBlock = this->GetThreadBlockSize(3*this->m_NumNodes);

    for (int offStart = 0; offStart < 3*this->m_NumNodes; offStart += szThreadBlock) {
      vp_threads.push_back(new boost::thread(_AddWorker, p_F + offStart, p_F + std::min(offStart + szThreadBlock, 3*this->m_NumNodes), mp_ThreadResultBuffer + offStart));      
    } 
    this->JoinThreads(vp_threads);
  }
}

template <class TShellElement>
static void _ElementThicknessWorker(float *p_T, const typename std::vector<TShellElement>::const_iterator ic_elsBegin, const typename std::vector<TShellElement>::const_iterator ic_elsEnd, const float u[]) {
  float *p_dst = p_T;

  for (typename std::vector<TShellElement>::const_iterator ic_el = ic_elsBegin; ic_el < ic_elsEnd; ic_el++, p_dst++) {
    *p_dst = ic_el->ComputeCurrentThickness(u);
  }
}

template <class TShellElement>
void tledParallelShellSolverCPU::ElementSetImpl<TShellElement>::ComputeElementThicknesses(float *p_dst, const float U[]) const {
  const int numThreadEls = this->GetThreadBlockSize(this->m_Elements.size());

  float *p_tT = p_dst;
  std::vector<boost::thread*> vp_threads;

  for (typename std::vector<TShellElement>::const_iterator ic_el = this->m_Elements.begin(); ic_el < this->m_Elements.end(); ic_el += numThreadEls, p_tT += numThreadEls) {
    vp_threads.push_back(new boost::thread(_ElementThicknessWorker<TShellElement>, p_tT, ic_el, std::min(this->m_Elements.end(), ic_el + numThreadEls), U));
    assert((int)vp_threads.size() <= this->GetNumberOfThreads());
  }
  this->JoinThreads(vp_threads);
}

void tledParallelShellSolverCPU::Init(tledSolver &r_solver, const tledModel &model) {
  assert(mp_ThreadResultBuffer == NULL);
  this->m_NumNodes = model.GetMesh()->GetNumNodes();
  mp_ThreadResultBuffer = new float[this->m_NumNodes*3*this->GetNumberOfThreads()];

  Superclass::Init(r_solver, model);
}

tledParallelShellSolverCPU::~tledParallelShellSolverCPU() {
  if (mp_ThreadResultBuffer != NULL) delete[] mp_ThreadResultBuffer;
}

tledShellSolver::ElementSet* tledParallelShellSolverCPU::CreateElementSetFromIndices3(tledShellMaterial &r_mat, const tledSurface &surface, const std::vector<int> &elInds) {
  typedef tledShellMesh<3> __Surface;

  const tledShellMesh<3> &mesh = dynamic_cast<const __Surface&>(surface);

  if (r_mat.HasBendingStiffness()) {
    return new ElementSetImpl<tledElementShellBSTP1>(r_mat, this->mp_ThreadResultBuffer, mesh, elInds, this->m_NumNodes);
  } else {
    if (!r_mat.IsNonLinear()) return new ElementSetImpl<tledElementMembraneSimpleLinear<3> >(r_mat, this->mp_ThreadResultBuffer, mesh, elInds, this->m_NumNodes);
    else return new ElementSetImpl<tledElementMembraneNonLinear<3> >(r_mat, this->mp_ThreadResultBuffer, mesh, elInds, this->m_NumNodes);
  }
}
#endif
