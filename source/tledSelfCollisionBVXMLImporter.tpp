// =========================================================================
// File:       tledSelfCollisionBVXMLImporter.tpp
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

template <class TBaseImporter>
void tledSelfCollisionBVXMLImporter<TBaseImporter>::Import() {
  Superclass::Import();
  
  this->GetOutput().UpdateCounter = this->template GetNumericElementValue<int>(this->GetUniqueChild("UpdateCounter", true));
  this->GetOutput().VolinoAngle = this->template GetNumericElementValue<float>(this->GetUniqueChild("VolinoAngle", true));
  if (this->GetOutput().VolinoAngle >= 0 && this->GetOutput().VolinoAngle <= BoundingVolume::GetVolinoThresholdAngle()) {
    std::vector<float> a = GetXMLTextAsVector<float>(this->GetUniqueChild("VolinoAxis", true));

    std::copy(a.begin(), a.end(), this->GetOutput().VolinoAxis);
    this->GetOutput().SubtreeMinH = this->template GetNumericElementValue<float>(this->GetUniqueChild("SubtreeMinH", true));
  } else {
    this->GetOutput().SubtreeMinH = std::numeric_limits<float>::quiet_NaN();
    std::fill(this->GetOutput().VolinoAxis, this->GetOutput().VolinoAxis + 3, std::numeric_limits<float>::quiet_NaN());
  }
}
