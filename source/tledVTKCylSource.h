// =========================================================================
// File:       tledVTKCylSource.h
// Purpose:    Model visualisation (VTK-based)
// Package:    NiftySim: Nonlinear finite element program
// Author:     Zeike Taylor
// Language:   C++
// Created:    August 2010
// 
// Copyright (c) 2010, University of Queensland. All rights reserved.
// MedTeQ Centre
// See the LICENSE.txt file in the root folder
// 
// ztaylor@itee.uq.edu.au
// =========================================================================

#ifdef _Visualisation_

#ifndef tledVTKCylSource_H
#define tledVTKCylSource_H

#include "tledHelper.h"

#include <vtkPolyData.h>
#include <vtkSmartPointer.h>

#include <vector>

class tledVTKCylSource
{
public:
   tledVTKCylSource() {;}
   tledVTKCylSource(float* a, float* A, float R, float L, int N);
   tledVTKCylSource(std::vector<float> a, std::vector<float> A, float R, float L, int N);
   ~tledVTKCylSource();
   
   vtkSmartPointer<vtkPolyData> GetOutput();
   
private:
   vtkSmartPointer<vtkPolyData> cyl;
   float* orig;
   float* axis;
   float rad;
   float len;
   int n;
};

#endif // tledVTKCylSource_H

#endif // _Visualisation_
