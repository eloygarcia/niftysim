// =========================================================================
// File:       tledOBBXMLImporter.tpp
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    March 2014
// 
// Copyright (c) 2014, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================

template <class TBV>
void tledOBBXMLImporter<TBV>::Import() {
  Superclass::Import();

  {
    std::vector<float> axes = GetXMLTextAsVector<float>(this->GetUniqueChild("Axes", true));
    std::vector<float> cnt = GetXMLTextAsVector<float>(this->GetUniqueChild("Centroid", true));
    std::vector<float> exts = GetXMLTextAsVector<float>(this->GetUniqueChild("Extents", true));

    std::copy(axes.begin(), axes.end(), &this->GetOutput().Axes[0][0]);
    std::copy(cnt.begin(), cnt.end(), this->GetOutput().Centroid);
    std::copy(exts.begin(), exts.end(), this->GetOutput().Extents);
  }
}
