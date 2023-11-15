// =========================================================================
// File:       tledAABBXMLExporter.tpp
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

template<class TBV>
void tledAABBXMLExporter<TBV>::WriteBody() {
  Superclass::WriteBody();
  this->CreateNumericListNode("Bounds", &this->GetInput().Bounds[0][0], 2*3);
}
