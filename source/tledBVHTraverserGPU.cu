// =========================================================================
// File:       tledBVHTraverserGPU.cu
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    May 2013
// 
// Copyright (c) 2013, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================
#ifndef tledBVHTraverserGPU_CU
#define tledBVHTraverserGPU_CU

#include "tledBVHTraverserGPU.h"
#include "tledCUDAMemoryBlock.h"
#include "tledComplexCUDAHelpers.h"

#include <thrust/functional.h>

tledBVHTraverserGPU::tledBVHTraverserGPU() : m_NumEdgeEdge(0), m_NumNodeFacet(0) {
  mp_EdgeEdgeResultBuffer = mp_NodeFacetResulBuffer = NULL;
  m_NumNodeFacet = m_NumEdgeEdge = 0;
}

tledBVHTraverserGPU::~tledBVHTraverserGPU() {
}

void tledBVHTraverserGPU::ResetResults(void) {
#ifndef NDEBUG
  if (mp_EdgeEdgeResultBuffer != NULL) {
    assert(!mp_EdgeEdgeResultBuffer->IsActive());
    mp_EdgeEdgeResultBuffer = NULL;
  }

  if (mp_NodeFacetResulBuffer != NULL) {
    assert(!mp_NodeFacetResulBuffer->IsActive());
    mp_NodeFacetResulBuffer = NULL;
  }
#endif

  this->SetNumberOfEdgeEdgeContacts(0);
  this->SetNumberOfNodeFacetContacts(0);
}

#endif
