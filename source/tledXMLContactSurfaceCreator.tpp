// =========================================================================
// File:       tledXMLContactSurfaceCreator.tpp
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
void tledXMLContactSurfaceCreator<TSurface>::XMLImporter::ParseFacetData() {
  Superclass::XMLImporter::ParseFacetData();
  
  {
    typedef typename TSurface::Facet __Facet;

    std::vector<int> facetEdges = GetXMLTextAsVector<int>(this->GetUniqueChild("FacetEdgeIndices", true));

    if ((int)facetEdges.size() != __Facet::NumberOfVertices*this->GetOutput().GetNumberOfFacets()) {
      std::cerr << "Expected " << __Facet::NumberOfVertices*this->GetOutput().GetNumberOfFacets() << " facet edge-indices, found " << facetEdges.size() << std::endl;
      std::abort();
    }

    for (int f = 0; f < this->GetOutput().GetNumberOfFacets(); f++) {
      std::vector<int>::const_iterator ic_fe = facetEdges.begin() + __Facet::NumberOfVertices*f;

      std::copy(ic_fe, ic_fe + __Facet::NumberOfVertices, this->GetOutput().GetFacet(f).EdgeIndices);
    }
  }
}

template <class TSurface>
void tledXMLContactSurfaceCreator<TSurface>::XMLImporter::ParseEdgeData() {
  const XMLNode edgeNode = this->GetUniqueChild("Edges", true);
  
  std::vector<int> edges = GetXMLTextAsVector<int>(edgeNode);
  std::vector<int>::const_iterator ic_en = edges.begin();
  
  this->GetOutput().GetAllEdges().resize(edges.size()/2);
  for (std::vector<std::pair<int, int> >::iterator i_e = this->GetOutput().GetAllEdges().begin(); i_e < this->GetOutput().GetAllEdges().end(); i_e++) {
    i_e->first = *(ic_en++);
    i_e->second = *(ic_en++);
  }
}

template <class TSurface>
void tledXMLContactSurfaceCreator<TSurface>::XMLImporter::Import() {
  Superclass::XMLImporter::Import();
  this->ParseEdgeData();
}
