// =========================================================================
// File:       tledVTKUltrasoundProbeSource.h
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

#ifndef tledVTKUltrasoundProbeSource_H
#define tledVTKUltrasoundProbeSource_H

#include "tledHelper.h"

#include <vtkPolyData.h>
#include <vtkSmartPointer.h>

#include <vector>

class tledVTKUltrasoundProbeSource
{
public:
   tledVTKUltrasoundProbeSource() {;}
   tledVTKUltrasoundProbeSource(float* a, float* A, float R, float L, int Nth, int Nph);
   tledVTKUltrasoundProbeSource(std::vector<float> a, std::vector<float> A, float R, float L, int Nth, int Nph);
   ~tledVTKUltrasoundProbeSource();
   
   vtkSmartPointer<vtkPolyData> GetOutput();
   
private:
   vtkSmartPointer<vtkPolyData> prb;
   float* orig;
   float* axis;
   float rad;
   float len;
   int nth;
   int nph;
};

#endif // tledVTKUltrasoundProbeSource_H

#endif // _Visualisation_
