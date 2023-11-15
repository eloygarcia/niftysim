// =========================================================================
// File:       tledDeformableContactSurfaceXMLExporter.tpp
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    April 2013
// 
// Copyright (c) 2013, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================

template <class TContactSurface>
void tledDeformableContactSurfaceXMLExporter<TContactSurface>::WriteMaps() {
  this->CreateNumericListNode("Volume2SurfaceNodeMap", this->GetInput().GetVolume2SurfaceNodeMap());
  this->CreateNumericListNode("Surface2VolumeNodeMap", this->GetInput().GetSurface2VolumeNodeMap());
}

template <class TContactSurface>
void tledDeformableContactSurfaceXMLExporter<TContactSurface>::WriteNodes() {
  Superclass::WriteNodes();
  
  {
    std::vector<int> numNodeFacets;

    numNodeFacets.reserve(this->GetInput().GetNumberOfNodes());
    for (int n = 0; n < this->GetInput().GetNumberOfNodes(); n++) {
      numNodeFacets.push_back(this->GetInput().GetNumberOfNodeFacets(n));
    }
    this->CreateNumericListNode("NumberOfNodeFacetIndices", numNodeFacets);
  }

  {
    std::vector<int> nodeFacets;

    nodeFacets.reserve(8*this->GetInput().GetNumberOfNodes());
    for (int n = 0; n < this->GetInput().GetNumberOfNodes(); n++) {
      nodeFacets.insert(nodeFacets.end(), this->GetInput().GetNodeFacetIndices(n), this->GetInput().GetNodeFacetIndices(n) + this->GetInput().GetNumberOfNodeFacets(n));
    }
    this->CreateNumericListNode("NodeFacetIndices", nodeFacets);
  }
}

template <class TContactSurface>
void tledDeformableContactSurfaceXMLExporter<TContactSurface>::WriteBody() {
  Superclass::WriteBody();
  this->WriteMaps();
}
