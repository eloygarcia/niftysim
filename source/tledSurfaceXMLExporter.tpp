// =========================================================================
// File:       tledSurfaceXMLExporter.tpp
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

template <class TSurface>
void tledSurfaceXMLExporter<TSurface>::WriteNodes() {
  XMLNode nodeNode = this->CreateNumericListNode("Nodes", this->GetInput().GetAllNodeCoordinates(), this->GetInput().GetNumberOfNodes()*3);
  
  this->AddNumericAttribute(nodeNode, "DOF", 3);
  this->AddNumericAttribute(nodeNode, "NumNodes", this->GetInput().GetNumberOfNodes());
}

template <class TSurface>
void tledSurfaceXMLExporter<TSurface>::WriteFacets() {
  std::vector<int> facetVtcs;
  XMLNode facetNode;

  facetVtcs.reserve(this->GetInput().GetNumberOfFacets()*Surface::Facet::NumberOfVertices);
  for (typename std::vector<typename Surface::Facet>::const_iterator ic_facet = this->GetInput().GetAllFacets().begin(); ic_facet < this->GetInput().GetAllFacets().end(); ic_facet++) {
    facetVtcs.insert(facetVtcs.end(), ic_facet->NodeIndices, ic_facet->NodeIndices + Surface::Facet::NumberOfVertices);
  }
  facetNode = this->CreateNumericListNode("Elements", facetVtcs);
  facetNode.addAttribute("Type", Surface::Facet::NumberOfVertices == 3? "T3" : "Q4");
  this->AddNumericAttribute(facetNode, "NumEls", this->GetInput().GetNumberOfFacets());
}

template <class TSurface>
void tledSurfaceXMLExporter<TSurface>::WriteBody() {
  this->WriteNodes();
  this->WriteFacets();
}
