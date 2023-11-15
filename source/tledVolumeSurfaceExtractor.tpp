// =========================================================================
// File:       tledVolumeSurfaceExtractor.tpp
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

template <class TSurface>
void tledVolumeSurfaceExtractor<TSurface>::Create() {
  this->GetMeshTopology().ComputeNodeAdjacency();
  this->GetMeshTopology().ComputeFacets();
  m_SurfaceFacetIndices = this->GetMeshTopology().GetSurfaceFacets();

  this->InitFacets();
  this->InitNodes();
  this->OrientFacets();
}

template <class TSurface>
void tledVolumeSurfaceExtractor<TSurface>::InitFacets() {
  this->GetOutput().SetNumberOfFacets(GetSurfaceFacetIndices().size());
  for (int f = 0; f < this->GetOutput().GetNumberOfFacets(); f++) {
    const int sFInd = this->GetSurfaceFacetIndices()[f];
    const typename tledMeshTopology<NumberOfElementVertices>::Facet &topoFacet = this->GetMeshTopology().GetFacet(sFInd);

    std::copy(topoFacet.begin(), topoFacet.end(), this->GetOutput().GetFacet(f).NodeIndices);
  }
}

template <class TSurface>
void tledVolumeSurfaceExtractor<TSurface>::OrientFacets() {
  using namespace tledVectorArithmetic;  
  typedef typename TSurface::Facet __Facet;
  
  for (int f = 0; f < this->GetOutput().GetNumberOfFacets(); f++) {
    const int fGlobInd = this->GetSurfaceFacetIndices()[f];

    float cent[3], normal[3];
    __Facet &r_facet = this->GetOutput().GetFacet(f);
    float const *pc_v0;

    pc_v0 = this->GetOutput().GetNodeCoordinates(r_facet.NodeIndices[0]);
    if (Dot(this->GetOutput().ComputeFacetNormal(normal, f), Sub(cent, this->GetMesh().ComputeCentroid(cent, this->GetMeshTopology().GetFacetElements(fGlobInd).first), pc_v0)) > 0) {
      if (__Facet::NumberOfVertices == 3) {
	std::iter_swap(r_facet.NodeIndices + 1, r_facet.NodeIndices + 2);
      } else {
	std::reverse(r_facet.NodeIndices, r_facet.NodeIndices + __Facet::NumberOfVertices);
      }
    }
    assert(Dot(this->GetOutput().ComputeFacetNormal(normal, f), Sub(cent, this->GetMesh().ComputeCentroid(cent, this->GetMeshTopology().GetFacetElements(fGlobInd).first), pc_v0)) < 0);
  }
}
