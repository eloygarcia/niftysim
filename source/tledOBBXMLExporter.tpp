// =========================================================================
// File:       tledOBBXMLExporter.tpp
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    March 2014
// 
// Copyright (c) 2013, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================

template<class TBV>
void tledOBBXMLExporter<TBV>::WriteBody() {
  Superclass::WriteBody();
  this->CreateNumericListNode("Extents", this->GetInput().Extents, 3);
  this->CreateNumericListNode("Axes", &this->GetInput().Axes[0][0], 3*3);
  this->CreateNumericListNode("Centroid", this->GetInput().Centroid, 3);
}
