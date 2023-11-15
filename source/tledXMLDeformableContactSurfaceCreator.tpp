// =========================================================================
// File:       tledXMLDeformableContactSurfaceCreator.tpp
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

template <class TSurface>
void tledXMLDeformableContactSurfaceCreator<TSurface>::XMLImporter::ParseMapsData() {
  this->GetOutput().SetVolume2SurfaceNodeMap(GetXMLTextAsVector<int>(this->GetUniqueChild("Volume2SurfaceNodeMap", true)));
  this->GetOutput().SetSurface2VolumeNodeMap(GetXMLTextAsVector<int>(this->GetUniqueChild("Surface2VolumeNodeMap", true)));
}

template <class TSurface>
void tledXMLDeformableContactSurfaceCreator<TSurface>::XMLImporter::ParseNodeData() {
  Superclass::XMLImporter::ParseNodeData();
  
  {
    std::vector<int> numNodeFacets = GetXMLTextAsVector<int>(this->GetUniqueChild("NumberOfNodeFacetIndices", true));
    std::vector<int> nodeFacetInds = GetXMLTextAsVector<int>(this->GetUniqueChild("NodeFacetIndices", true));
    
    this->GetOutput().SetNodeFacetIndices(nodeFacetInds, numNodeFacets);
  }

  std::copy(this->GetOutput().GetAllNodeCoordinates(), this->GetOutput().GetAllNodeCoordinates() + 3*this->GetOutput().GetNumberOfNodes(), this->GetOutput().GetAllNodeCoordinates0().begin());
}

template <class TSurface>
void tledXMLDeformableContactSurfaceCreator<TSurface>::XMLImporter::Import() {
  this->GetOutput().ResetSaveCount();
  this->GetOutput().ResetUpdateCount();

  Superclass::XMLImporter::Import();
  this->ParseMapsData();
  this->GetOutput().LoadFromXMLPostloadHook();
}
