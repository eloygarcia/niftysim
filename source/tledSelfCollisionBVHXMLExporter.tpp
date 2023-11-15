// =========================================================================
// File:       tledSelfCollisionBVHXMLExporter.tpp
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

template <class TBVH>
void tledSelfCollisionBVHXMLExporter<TBVH>::WriteAdjacencyData() {
  if (this->GetInput().GetNonAdjacentGeometryNodes().size() > 0) this->CreateNumericListNode("NonAdjacentGeometryNodes", this->GetInput().GetNonAdjacentGeometryNodes());
  if (this->GetInput().GetGeometryClusterSubtreeRootIndices().size() > 0) this->CreateNumericListNode("GeometryClusterSubtreeRootIndices", this->GetInput().GetGeometryClusterSubtreeRootIndices());
}

template <class TBVH>
void tledSelfCollisionBVHXMLExporter<TBVH>::WriteCollisionCandidates() {
  assert(this->GetInput().GetSelfCollisionCandidates().size() > 0);
  this->CreateNumericListNode("SelfCollisionCandidates", this->GetInput().GetSelfCollisionCandidates());
}

template <class TBVH>
void tledSelfCollisionBVHXMLExporter<TBVH>::WriteBody() {
  Superclass::WriteBody();
  this->WriteAdjacencyData();
  this->WriteCollisionCandidates();
}
