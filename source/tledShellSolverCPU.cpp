// =========================================================================
// File:       tledShellSolver.cpp
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    March 2012
// 
// Copyright (c) 2012, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================

#include "tledHelper.h"
#include "tledElementMembraneSimpleLinear.h"
#include "tledElementMembraneNonLinear.h"
#include "tledShellSolverCPU.h"
#include "tledMembraneMaterialLinear.h"
#include "tledElementShellBSTP1.h"

#include <typeinfo>
#include <cstdlib>
#include <iostream>
#include <iterator>

template <class TShellElement>
void tledShellSolverCPU::ElementSetImpl<TShellElement>::ComputeForces(float *p_F, const float u[]) {
  for (typename std::vector<TShellElement>::iterator i_el = this->m_Elements.begin(); i_el < this->m_Elements.end(); i_el++) i_el->ComputeElementForces(p_F, u);
}

void tledShellSolverCPU::ComputeNewForces(float *p_F, const float u[]) {
  for (std::vector<Superclass::ElementSet*>::iterator ip_elSet = mvp_ShellElementSets.begin(); ip_elSet < mvp_ShellElementSets.end(); ip_elSet++) {
    dynamic_cast<ElementSet*>(*ip_elSet)->ComputeForces(p_F, u);
  }
}

template <>
tledShellSolver::ElementSetImpl<tledElementShellBSTP1, tledShellSolverCPU::ElementSet>::ElementSetImpl(tledShellMaterial &r_mat, const tledElementShellBSTP1::Surface &surface, const std::vector<int> &elInds) : tledShellSolverCPU::ElementSet(r_mat, elInds) {
  tledSurfaceTopology<tledShellMesh<3> > topology(surface);
  std::vector<int> reverseElSetIndexMap(surface.GetNumberOfFacets(), -1);

  topology.ComputeEdges();
  m_Elements.resize(elInds.size());

  for (std::vector<int>::const_iterator ic_e = elInds.begin(); ic_e < elInds.end(); ic_e++) reverseElSetIndexMap[*ic_e] = ic_e - elInds.begin();
  for (std::vector<tledElementShellBSTP1>::iterator i_el = m_Elements.begin(); i_el < m_Elements.end(); i_el++) {
    i_el->SetMaterial(r_mat);
    i_el->InitialiseElement(topology, elInds[i_el-m_Elements.begin()]);
  }

  for (std::vector<tledElementShellBSTP1>::iterator i_el = m_Elements.begin(); i_el < m_Elements.end(); i_el++) {
    i_el->InitialiseEdgeElements(m_Elements, topology, elInds[i_el-m_Elements.begin()], reverseElSetIndexMap);
  }

  for (std::vector<tledElementShellBSTP1>::iterator i_el = m_Elements.begin(); i_el < m_Elements.end(); i_el++) i_el->FinaliseInitialisation();
}

template <>
void tledShellSolver::ElementSetImpl<tledElementShellBSTP1, tledShellSolverCPU::ElementSet>::ClampNodes(const std::vector<bool> &clampedNodeMask) {
  for (std::vector<tledElementShellBSTP1>::iterator i_el = m_Elements.begin(); i_el < m_Elements.end(); i_el++) {
      i_el->ClampNodes(clampedNodeMask);
    }
}

tledShellSolver::ElementSet* tledShellSolverCPU::CreateElementSetFromIndices3(tledShellMaterial &r_mat, const tledSurface &surface, const std::vector<int> &elInds) {
  typedef tledShellMesh<3> __Surface;

  const tledShellMesh<3> &mesh = dynamic_cast<const __Surface&>(surface);

  if (r_mat.HasBendingStiffness()) {
    return new ElementSetImpl<tledElementShellBSTP1>(r_mat, mesh, elInds);
  } else {
    if (!r_mat.IsNonLinear()) return new ElementSetImpl<tledElementMembraneSimpleLinear<3> >(r_mat, mesh, elInds);
    else return new ElementSetImpl<tledElementMembraneNonLinear<3> >(r_mat, mesh, elInds);
  }
}

void tledShellSolverCPU::ComputeElementSetBendingMoments(float *p_M, const int elSetIndex, const float uCurr[]) const {
  static_cast<const ElementSet*>(mvp_ShellElementSets[elSetIndex])->ComputeElementBendingMoments(p_M, uCurr);
}
