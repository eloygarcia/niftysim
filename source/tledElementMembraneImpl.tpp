// =========================================================================
// File:       tledElementMembraneImpl.tpp
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

template <const int t_numFacetVertices>
void tledElementMembraneImpl<t_numFacetVertices>::InitialiseElement(const Surface &surface, const int facetIndex) {
  const int *facet = surface.GetFacet(facetIndex).NodeIndices;

  std::copy(facet, facet + t_numFacetVertices, m_FacetVertexIndices);
}

template <const int t_numFacetVertices>
void tledElementMembraneImpl<t_numFacetVertices>::ComputeElementMass(float *p_dst) const {
  const float mass = GetDensity()*GetArea()*GetThickness()/t_numFacetVertices;

#ifndef _GPU_
  assert(!std::isnan(mass));
#endif
  for (int vInd = 0; vInd < t_numFacetVertices; vInd++) {
    p_dst[m_FacetVertexIndices[vInd]] += mass;
  }
}

template <const int t_numFacetVertices>
template <const int t_numStrainComponents, const int t_numStressComponents>
void tledElementMembraneImpl<t_numFacetVertices>::ComputeElementForcesTemplate(float *p_dst, const float U[]) {
  float strain[t_numStrainComponents], stress[t_numStrainComponents];

  this->ComputeStrain(strain, U);
  this->GetMaterial().ComputeStress(stress, strain);
  this->ComputeForces(p_dst, stress);
}

template <const int t_numFacetVertices>
void tledElementMembraneImpl<t_numFacetVertices>::ComputeElementBasis(float (*p_basis)[3], float &r_area, const float X[]) const {
  std::abort();
}

template <>
void tledElementMembraneImpl<3>::ComputeElementBasis(float (*p_basis)[3], float &r_area, const float X[]) const;

#ifdef _GPU_
template <const int t_numFacetVertices>
void tledElementMembraneImpl<t_numFacetVertices>::InitGPU(tledElementMembrane::GPUElement &r_dst) {
  r_dst.Area = this->GetArea();
  r_dst.Thickness = this->GetMaterial().GetThickness();

  r_dst.ElementNodeIndices.x = this->GetFacetVertexIndices()[0];
  r_dst.ElementNodeIndices.y = this->GetFacetVertexIndices()[1];
  r_dst.ElementNodeIndices.z = this->GetFacetVertexIndices()[2];
  if (t_numFacetVertices > 3) r_dst.ElementNodeIndices.w = this->GetFacetVertexIndices()[3];
}
#endif
