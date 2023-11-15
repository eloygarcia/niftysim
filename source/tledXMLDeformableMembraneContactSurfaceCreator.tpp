// =========================================================================
// File:       tledXMLDeformableMembraneContactSurfaceCreator.tpp
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
void tledXMLDeformableMembraneContactSurfaceCreator<TSurface>::XMLImporter::ParseIndexBounds() {
  this->GetOutput().SetNumberOfMembraneNodes(this->template GetNumericElementValue<int>(this->GetUniqueChild("NumberOfMembraneNodes", true)));
  this->GetOutput().SetNumberOfMembraneElements(this->template GetNumericElementValue<int>(this->GetUniqueChild("NumberOfMembraneFacets", true)));
  this->GetOutput().SetNumberOfMembraneEdges(this->template GetNumericElementValue<int>(this->GetUniqueChild("NumberOfMembraneEdges", true)));

  this->GetOutput().SetMembraneFacetBaseIndex(this->template GetNumericElementValue<int>(this->GetUniqueChild("MembraneFacetBaseIndex", true)));
  this->GetOutput().SetMembraneNodeBaseIndex(this->template GetNumericElementValue<int>(this->GetUniqueChild("MembraneNodeBaseIndex", true)));
  this->GetOutput().SetMembraneEdgeBaseIndex(this->template GetNumericElementValue<int>(this->GetUniqueChild("MembraneEdgeBaseIndex", true)));
}

template <class TSurface>
void tledXMLDeformableMembraneContactSurfaceCreator<TSurface>::XMLImporter::ParseCoordinates0() {
  if (this->GetRootNode().nChildNode("Nodes0") > 0) {
    const XMLNode nodeNode = this->GetUniqueChild("Nodes0", true);
    const std::vector<float> nodes = GetXMLTextAsVector<float>(nodeNode);

    std::copy(nodes.begin(), nodes.end(), this->GetOutput().GetAllNodeCoordinates0().begin());
  }
}

template <class TSurface>
void tledXMLDeformableMembraneContactSurfaceCreator<TSurface>::XMLImporter::ParseNodeData() {
  Superclass::XMLImporter::ParseNodeData();
  this->ParseCoordinates0();
}

template <class TSurface>
void tledXMLDeformableMembraneContactSurfaceCreator<TSurface>::XMLImporter::Import() {
  Superclass::XMLImporter::Import();
  this->ParseIndexBounds();
}
