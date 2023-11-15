// =========================================================================
// File:       tledVTKPlateSource.cpp
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

#include "tledVTKPlateSource.h"

#include <vtkQuad.h>
#include <vtkPoints.h>
#include <vtkPlaneSource.h>

#include <cstring>

using namespace std;

tledVTKPlateSource::tledVTKPlateSource(float* A, float* B, float* C)
{
   a = new float[3];
   memcpy(a,A,sizeof(float)*3);
   b = new float[3];
   memcpy(b,B,sizeof(float)*3);
   c = new float[3];
   memcpy(c,C,sizeof(float)*3);
   plt = vtkSmartPointer<vtkPolyData>::New();
}

tledVTKPlateSource::tledVTKPlateSource(vector<float> A, vector<float> B, vector<float> C)
{
   a = new float[3];
   a[0] = A[0]; a[1] = A[1]; a[2] = A[2];
   b = new float[3];
   b[0] = B[0]; b[1] = B[1]; b[2] = B[2];
   c = new float[3];
   c[0] = C[0]; c[1] = C[1]; c[2] = C[2];
   plt = vtkSmartPointer<vtkPolyData>::New();
}

tledVTKPlateSource::~tledVTKPlateSource()
{
   delete[] a;
   delete[] b;
   delete[] c;
}

vtkSmartPointer<vtkPolyData> tledVTKPlateSource::GetOutput()
{
   // Compute 4th corner position
   float d[3];
   d[0] = (b[0] + c[0])/2; d[1] = (b[1] + c[1])/2; d[2] = (b[2] + c[2])/2; // midpoint between corners b & c
   // d-a is the vector from a to d --> add this to the midpoint to get the opposite corner
   d[0] += d[0] - a[0]; d[1] += d[1] - a[1]; d[2] += d[2] - a[2];
   
   // Assemble points
   vtkSmartPointer<vtkPoints> pnts = vtkSmartPointer<vtkPoints>::New();
   pnts->InsertNextPoint(a);
   pnts->InsertNextPoint(b);
   pnts->InsertNextPoint(d);
   pnts->InsertNextPoint(c);
   plt->SetPoints(pnts);
   
   // Construct the quad
   vtkSmartPointer<vtkQuad> face;
   plt->Allocate();
   face = vtkSmartPointer<vtkQuad>::New();
   face->GetPointIds()->SetId(0,0);
   face->GetPointIds()->SetId(1,1);
   face->GetPointIds()->SetId(2,2);
   face->GetPointIds()->SetId(3,3);
   
   plt->InsertNextCell(face->GetCellType(), face->GetPointIds());
      
   return plt;
}
