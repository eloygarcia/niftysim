// =========================================================================
// File:       tledAABBXMLImporter.tpp
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

template <class TBV>
void tledAABBXMLImporter<TBV>::Import() {
  Superclass::Import();

  {
    std::vector<float> bnds = GetXMLTextAsVector<float>(this->GetUniqueChild("Bounds", true));

    std::copy(bnds.begin(), bnds.end(), &this->GetOutput().Bounds[0][0]);
  }
}
