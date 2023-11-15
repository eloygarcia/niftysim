// =========================================================================
// File:       tledShellSolver.cpp
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    February 2012
// 
// Copyright (c) 2012, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================
#include <cstdlib>

#include "tledShellSolver.h"

tledShellSolver::~tledShellSolver() {
  for (std::vector<ElementSet*>::iterator ip_elSet = mvp_ShellElementSets.begin(); ip_elSet < mvp_ShellElementSets.end(); ip_elSet++) delete *ip_elSet;
  if (mp_Surface) delete mp_Surface;  
}

tledShellSolver::tledShellSolver() {
  mp_Surface = NULL;
  mp_MainSolver = NULL;
}

void tledShellSolver::Init(tledSolver &solver, const tledModel &model) {
  const std::string surfaceType = model.GetShellMeshType();

  mp_MainSolver = &solver;
  m_Dt = (float)model.GetTimeStep();
  if (surfaceType == "T3") {
    mp_Surface = model.GetShellMesh<3>();
    GetSurface<tledShellMesh<3> >().SetNodeVector(mp_MainSolver->GetMesh()->GetAllNodeCds(), mp_MainSolver->GetMesh()->GetNumNodes());
  } else if (surfaceType == "SURFACE" && (std::string(model.GetElType()) == "T4" || std::string(model.GetElType()) == "T4ANP")) {
    mp_Surface = model.GetShellMesh<3>();
  } else std::abort();

  m_ElementElementSetIndices.clear();
  m_ElementElementSetIndices.insert(m_ElementElementSetIndices.end(), mp_Surface->GetNumberOfFacets(), std::pair<int, int>(-1, -1));

  m_MShell.insert(m_MShell.end(), model.GetMesh()->GetNumNodes(), 0.0f);
  for (int elSetInd = 0; elSetInd < model.GetNumberOfShellElementSets(); elSetInd++) {
    std::vector<int> clampedNodes;

    mvp_ShellElementSets.push_back(CreateElementSet(*mp_Surface, model, elSetInd));
    mvp_ShellElementSets.back()->ComputeMass(&m_MShell.front());
    clampedNodes = model.GetClampedShellNodes(elSetInd);
    if (clampedNodes.size() > 0) {
      std::vector<bool> clampedNodeMask;

      clampedNodeMask.insert(clampedNodeMask.end(), mp_Surface->GetNumberOfNodes(), false);
      for (std::vector<int>::const_iterator ic_n = clampedNodes.begin(); ic_n < clampedNodes.end(); ic_n++) clampedNodeMask[*ic_n] = true;
      mvp_ShellElementSets[elSetInd]->ClampNodes(clampedNodeMask);
    }
  }
}

tledShellSolver::ElementSet* tledShellSolver::CreateElementSet(const tledSurface &surface, const tledModel &model, const int elSetInd) {
  tledShellMaterial &r_mat = *model.GetShellMaterial(elSetInd);

  if (std::string(model.GetShellMeshType()) == "SURFACE") {
    for (int e = 0; e < surface.GetNumberOfFacets(); e++) m_ElementElementSetIndices[e] = std::pair<int, int>(0, e);

    if (std::string(model.GetElType()) == "T4" || std::string(model.GetElType()) == "T4ANP") return CreateElementSetFromIndices3(r_mat, surface, tledSequenceGenerator::MakeSequence(0, surface.GetNumberOfFacets()));
    else std::cerr << "Quads not yet supported.\n";
  } else {
    std::vector<int> elSetInds = model.GetShellElementSet(elSetInd);
    
    if (elSetInds.size() == 1 && elSetInds.front() == -1) elSetInds = tledSequenceGenerator::MakeSequence(0, surface.GetNumberOfFacets());
    for (std::vector<int>::const_iterator ic_elInd = elSetInds.begin(); ic_elInd < elSetInds.end(); ic_elInd++) {
      m_ElementElementSetIndices[*ic_elInd] = std::pair<int, int>(elSetInd, ic_elInd - elSetInds.begin());
    }

    if (std::string(model.GetShellMeshType()) == "T3") {
      return CreateElementSetFromIndices3(r_mat, surface, elSetInds);
    } else if (std::string(model.GetShellMeshType()) == "Q4") std::cerr << "Quads not yet supported.\n";
    std::abort();    
  }

  return NULL;
}

const tledElementMembrane& tledShellSolver::GetElement(const int elIndex) const {
  const std::pair<int, int> &locInd = GetElementSetIndexForElement(elIndex);

  if (locInd.first == -1) return *(const tledElementMembrane*)NULL;
  else return mvp_ShellElementSets[locInd.first]->GetElement(locInd.second);
} 

void tledShellSolver::ComputeAllElementThicknesses(float *p_dst, const float U[]) const {
  int szBuffer = 0;
  std::vector<float> tmpBuffer;

  for (std::vector<ElementSet*>::const_iterator ipc_es = mvp_ShellElementSets.begin(); ipc_es < mvp_ShellElementSets.end(); ipc_es++) szBuffer = std::max(szBuffer, (*ipc_es)->GetNumberOfElements());
  tmpBuffer.resize(szBuffer);
  
  for (std::vector<ElementSet*>::const_iterator ipc_es = mvp_ShellElementSets.begin(); ipc_es < mvp_ShellElementSets.end(); ipc_es++) {
    (*ipc_es)->ComputeElementThicknesses(&tmpBuffer.front(), U);
    for (int lei = 0; lei < (*ipc_es)->GetNumberOfElements(); lei++) p_dst[(*ipc_es)->GetElementIndices()[lei]] = tmpBuffer[lei];
  }
}

void tledShellSolver::GetAllInitialElementThicknesses(float *p_dst) const {
  for (std::vector<ElementSet*>::const_iterator ipc_es = mvp_ShellElementSets.begin(); ipc_es < mvp_ShellElementSets.end(); ipc_es++) {
    const float t = (*ipc_es)->GetMaterial().GetThickness();

    for (std::vector<int>::const_iterator ic_ei = (*ipc_es)->GetElementIndices().begin(); ic_ei < (*ipc_es)->GetElementIndices().end(); ic_ei++) p_dst[*ic_ei] = t;
  }
}
