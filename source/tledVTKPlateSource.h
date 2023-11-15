// =========================================================================
// File:       tledVTKPlateSource.h
// Purpose:    Model visualisation (VTK-based)
// Package:    NiftySim: Nonlinear finite element program
// Author:     Zeike Taylor
// Language:   C++
// Created:    November 2010
// 
// Copyright (c) 2010, University of Queensland. All rights reserved.
// MedTeQ Centre
// See the LICENSE.txt file in the root folder
// 
// ztaylor@itee.uq.edu.au
// =========================================================================

#ifdef _Visualisation_

#ifndef tledVTKPlateSource_H
#define tledVTKPlateSource_H

#include "tledHelper.h"

#include <vtkPolyData.h>
#include <vtkSmartPointer.h>

#include <vector>

class tledVTKPlateSource
{
public:
   tledVTKPlateSource() {;}
   tledVTKPlateSource(float* A, float* B, float* C);
   tledVTKPlateSource(std::vector<float> A, std::vector<float> B, std::vector<float> C);
   ~tledVTKPlateSource();
   
   vtkSmartPointer<vtkPolyData> GetOutput();
   
private:
   vtkSmartPointer<vtkPolyData> plt;
   float* a;
   float* b;
   float* c;
};

#endif // tledVTKPlateSource_H

#endif // _Visualisation_
