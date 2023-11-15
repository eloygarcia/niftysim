// =========================================================================
// File:       tledDeformableMembraneContactSurfaceXMLExporter.tpp
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

template <class TContactSurface>
void tledDeformableMembraneContactSurfaceXMLExporter<TContactSurface>::WriteIndexBounds() {
  this->CreateNumericNode("NumberOfMembraneFacets", this->GetInput().GetNumberOfMembraneElements());
  this->CreateNumericNode("NumberOfMembraneNodes", this->GetInput().GetNumberOfMembraneNodes());
  this->CreateNumericNode("NumberOfMembraneEdges", this->GetInput().GetNumberOfMembraneEdges());
  
  this->CreateNumericNode("MembraneFacetBaseIndex", this->GetInput().GetMembraneFacetBaseIndex());
  this->CreateNumericNode("MembraneNodeBaseIndex", this->GetInput().GetMembraneNodeBaseIndex());
  this->CreateNumericNode("MembraneEdgeBaseIndex", this->GetInput().GetMembraneEdgeBaseIndex());
}

template <class TContactSurface>
void tledDeformableMembraneContactSurfaceXMLExporter<TContactSurface>::WriteNodes0() {
  XMLNode nodeNode = this->CreateNumericListNode("Nodes0", this->GetInput().GetAllNodeCoordinates0(), this->GetInput().GetNumberOfNodes()*3);
  
  this->AddNumericAttribute(nodeNode, "DOF", 3);
  this->AddNumericAttribute(nodeNode, "NumNodes", this->GetInput().GetNumberOfNodes());
  nodeNode.addAttribute("DOF", "3");
}

template <class TContactSurface>
void tledDeformableMembraneContactSurfaceXMLExporter<TContactSurface>::WriteBody() {
  Superclass::WriteBody();
  this->WriteIndexBounds();
  this->WriteNodes0();
}
