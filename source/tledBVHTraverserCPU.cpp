// =========================================================================
// File:       tledBVHTraverserCPU.cpp
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    November 2012
// 
// Copyright (c) 2012, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================

#include "tledBVHTraverserCPU.h"
#include "tledUnstructuredContactManager.h"

std::vector<int> tledBVHTraverserCPU::s_ContactNodeIndices, tledBVHTraverserCPU::s_ContactEdgeIndices;
std::vector<std::pair<int, int> > tledBVHTraverserCPU::svv_NodeFacetNarrowPhasePairs, tledBVHTraverserCPU::svv_EdgeEdgeNarrowPhasePairs;
std::vector<tledBVHTraverserCPU::MasterEdgeContactData> tledBVHTraverserCPU::s_MasterEdges;
std::vector<tledBVHTraverserCPU::MasterFacetContactData> tledBVHTraverserCPU::s_MasterFacets;

void tledBVHTraverserCPU::FinishBroadPhase() {
  this->GetNodeFacetNarrowPhasePairs() = tledHelper::MakeSortedUnique(this->GetNodeFacetNarrowPhasePairs(), std::equal_to<std::pair<int, int> >(), NarrowPhaseOrdering());
  this->GetEdgeEdgeNarrowPhasePairs() = tledHelper::MakeSortedUnique(this->GetEdgeEdgeNarrowPhasePairs(), std::equal_to<std::pair<int, int> >(), NarrowPhaseOrdering());
}
