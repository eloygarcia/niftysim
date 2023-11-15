// =========================================================================
// File:       tledShellSolver.tpp
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

template <class TShellElement, class TElementSetInterface>
void tledShellSolver::ElementSetImpl<TShellElement, TElementSetInterface>::ComputeMass(float *p_mass) {
  for (int elInd = 0; elInd < this->GetNumberOfElements(); elInd++) this->GetElement(elInd).ComputeElementMass(p_mass);
}

template <class TShellElement, class TElementSetInterface>
void tledShellSolver::ElementSetImpl<TShellElement, TElementSetInterface>::ComputeElementThicknesses(float *p_dst, const float U[]) const {
  float *p_t = p_dst;

  for (typename std::vector<TShellElement>::const_iterator ic_el = this->m_Elements.begin(); ic_el < this->m_Elements.end(); ic_el++, p_t++) *p_t = ic_el->ComputeCurrentThickness(U);
}
