// =========================================================================
// File:       tledBVXMLExporter.tpp
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

template <class TBV>
void tledBVXMLExporter<TBV>::WriteBody() {
  this->CreateNumericNode("ParentIndex", this->GetInput().ParentIndex);  
  if (this->GetInput().PrimitiveIndex >= 0) {
    this->CreateNumericNode("PrimitiveIndex", this->GetInput().PrimitiveIndex);
  } else {
    this->CreateNumericListNode("ChildIndices", this->GetInput().ChildIndices, BoundingVolume::NumberOfChildBVs);
  }

  this->AddNumericAttribute(this->GetRootNode(), "NumberOfChildBVs", BoundingVolume::NumberOfChildBVs);
}
