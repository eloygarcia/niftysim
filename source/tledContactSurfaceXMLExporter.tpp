// =========================================================================
// File:       tledContactSurfaceXMLExporter.tpp
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
void tledContactSurfaceXMLExporter<TContactSurface>::WriteEdges() {
  std::vector<int> edges;

  edges.reserve(2*this->GetInput().GetNumberOfEdges());
  for (std::vector<std::pair<int, int> >::const_iterator ic_e = this->GetInput().GetAllEdges().begin(); ic_e < this->GetInput().GetAllEdges().end(); ic_e++) {
    edges.push_back(ic_e->first);
    edges.push_back(ic_e->second);
  }
  this->CreateNumericListNode("Edges", edges);  
}

template <class TContactSurface>
void tledContactSurfaceXMLExporter<TContactSurface>::WriteFacets() {
  Superclass::WriteFacets();

  {
    std::vector<int> facetEdges;

    facetEdges.reserve(Surface::Facet::NumberOfVertices*this->GetInput().GetNumberOfFacets());
    for (typename std::vector<typename Surface::Facet>::const_iterator ic_facet = this->GetInput().GetAllFacets().begin(); ic_facet < this->GetInput().GetAllFacets().end(); ic_facet++) {
      facetEdges.insert(facetEdges.end(), ic_facet->EdgeIndices, ic_facet->EdgeIndices + Surface::Facet::NumberOfVertices);
    }
    this->CreateNumericListNode("FacetEdgeIndices", facetEdges);
  }

  if (this->GetInput().GetMinDiameter() == this->GetInput().GetMinDiameter()) {
    this->CreateNumericNode("MinH", this->GetInput().GetMinDiameter());
  }

  if (this->GetInput().GetMaxDiameter() == this->GetInput().GetMaxDiameter()) {
    this->CreateNumericNode("MaxH", this->GetInput().GetMaxDiameter());
  }
}

template <class TContactSurface>
void tledContactSurfaceXMLExporter<TContactSurface>::WriteBody() {
  Superclass::WriteBody();
  this->WriteEdges();
}
