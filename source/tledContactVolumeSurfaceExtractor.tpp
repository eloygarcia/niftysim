// =========================================================================
// File:       tledContactVolumeSurfaceExtractor.tpp
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
void tledContactVolumeSurfaceExtractor<TSurface>::_BaseExtractor::InitNodes() {
  typedef typename Surface::Facet __Facet;

  std::vector<int> surface2VolumeMap, volume2SurfaceMap;
  
  /*
   * - Compile list of surface vertices
   * - Create map from global vertex indices to local (surface) indices
   * - Translate surface triangle definitions from global to local
   * - Create copies of surface nodes
   */
  surface2VolumeMap.clear();
  for (typename std::vector<__Facet>::const_iterator ic_facet = this->GetOutput().GetAllFacets().begin(); ic_facet < this->GetOutput().GetAllFacets().end(); ic_facet++) {
    surface2VolumeMap.insert(surface2VolumeMap.end(), ic_facet->NodeIndices, ic_facet->NodeIndices + __Facet::NumberOfVertices);
  }
  this->GetOutput().SetSurface2VolumeNodeMap(tledHelper::MakeSortedUnique(surface2VolumeMap));
  surface2VolumeMap.clear();

  {
    float *p_dst;
    std::vector<int>::const_iterator ic_nInd;

    this->GetOutput().SetNumberOfNodes(this->GetOutput().GetSurface2VolumeNodeMap().size());
    for (ic_nInd = this->GetOutput().GetSurface2VolumeNodeMap().begin(), p_dst = this->GetOutput().GetAllNodeCoordinates(); ic_nInd < this->GetOutput().GetSurface2VolumeNodeMap().end(); ic_nInd++, p_dst += 3) {
      std::copy(this->GetMesh().GetAllNodeCds() + 3*(*ic_nInd), this->GetMesh().GetAllNodeCds() + 3*(*ic_nInd) + 3, p_dst);
    }
  }

  volume2SurfaceMap.clear();
  volume2SurfaceMap.insert(volume2SurfaceMap.end(), this->GetMesh().GetNumNodes(), -1);
  for (std::vector<int>::const_iterator ic_nInd = this->GetOutput().GetSurface2VolumeNodeMap().begin(); ic_nInd < this->GetOutput().GetSurface2VolumeNodeMap().end(); ic_nInd++) {
    volume2SurfaceMap[*ic_nInd] = ic_nInd - this->GetOutput().GetSurface2VolumeNodeMap().begin();
  }
  this->GetOutput().SetVolume2SurfaceNodeMap(volume2SurfaceMap);
  volume2SurfaceMap.clear();

  for (int f = 0; f < this->GetOutput().GetNumberOfFacets(); f++) {
    __Facet &r_facet = this->GetOutput().GetFacet(f);

    for (int vInd = 0; vInd < __Facet::NumberOfVertices; vInd++) r_facet.NodeIndices[vInd] = this->GetOutput().MapVolume2SurfaceNode(r_facet.NodeIndices[vInd]);
  }
}  
  
template <class TSurface> 
void tledContactVolumeSurfaceExtractor<TSurface>::Create() {
  Superclass::Create();
}
