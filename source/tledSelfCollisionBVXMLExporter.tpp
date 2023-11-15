// =========================================================================
// File:       tledSelfCollisionBVXMLExporter.tpp
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

template <class TBaseExporter>
void tledSelfCollisionBVXMLExporter<TBaseExporter>::WriteBody() {
  Superclass::WriteBody();
  this->CreateNumericNode("UpdateCounter", this->GetInput().UpdateCounter);

  if (this->GetInput().VolinoAngle >= 0 && this->GetInput().VolinoAngle <= BoundingVolume::GetVolinoThresholdAngle()) {
    this->CreateNumericListNode("VolinoAxis", this->GetInput().VolinoAxis, 3);
    this->CreateNumericNode("SubtreeMinH", this->GetInput().SubtreeMinH);
  } 
  this->CreateNumericNode("VolinoAngle", this->GetInput().VolinoAngle);
}
