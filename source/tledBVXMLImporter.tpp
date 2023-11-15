// =========================================================================
// File:       tledBVXMLImporter.tpp
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
void tledBVXMLImporter<TBV>::Import() {
  this->GetOutput().ParentIndex = this->template GetNumericElementValue<int>(this->GetUniqueChild("ParentIndex", true));
  if (this->GetRootNode().nChildNode("PrimitiveIndex")) {
    this->GetOutput().PrimitiveIndex = this->template GetNumericElementValue<int>(this->GetUniqueChild("PrimitiveIndex", true));
    std::fill(this->GetOutput().ChildIndices, this->GetOutput().ChildIndices + BoundingVolume::NumberOfChildBVs, -1);
  } else if (this->GetRootNode().nChildNode("ChildIndices")) {
    std::vector<int> childIndices = GetXMLTextAsVector<int>(this->GetRootNode().getChildNode("ChildIndices"));

    std::copy(childIndices.begin(), childIndices.end(), this->GetOutput().ChildIndices);
    this->GetOutput().PrimitiveIndex = -1;
  } else {
    std::cerr << "Have neither a primitive index nor child indices.\n";
    std::abort();
  }
}
