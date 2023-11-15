// =========================================================================
// File:       tledVTKCylSource.cpp
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

#include "tledVTKCylSource.h"

#include <vtkQuad.h>

#include <cstring>

using namespace std;

tledVTKCylSource::tledVTKCylSource(float* a, float* A, float R, float L, int N)
{
   orig = new float[3];
   memcpy(orig,a,sizeof(float)*3);
   axis = new float[3];
   memcpy(axis,A,sizeof(float)*3);
   rad = R;
   len = L;
   n = N;
   cyl = vtkSmartPointer<vtkPolyData>::New();
}

tledVTKCylSource::tledVTKCylSource(vector<float> a, vector<float> A, float R, float L, int N)
{
   orig = new float[3];
   orig[0] = a[0]; orig[1] = a[1]; orig[2] = a[2];
   axis = new float[3];
   axis[0] = A[0]; axis[1] = A[1]; axis[2] = A[2];
   rad = R;
   len = L;
   n = N;
   cyl = vtkSmartPointer<vtkPolyData>::New();
}

tledVTKCylSource::~tledVTKCylSource()
{
   delete[] orig;
   delete[] axis;
}

vtkSmartPointer<vtkPolyData> tledVTKCylSource::GetOutput()
{
   // Construct points in cylinder local coords ==========
   
   double pi = 3.141592653589793;
   double g = 2*pi/n;
   double th;
   double* px = new double[(n+1)*2];
   double* py = new double[(n+1)*2];
   double* pz = new double[(n+1)*2];
   for (int i = 0; i < n+1; i++)
   {
      th = i*g;
      px[i] = 0; px[i+n+1] = len;
      py[i] = rad*cos(th); py[i+n+1] = py[i];
      pz[i] = rad*sin(th); pz[i+n+1] = pz[i];
   }
   
   // Transform points to global coords ==========
   
   // Define the local system (a,b,c) (NB: non-unique)
   float* a = axis;
   int i = 0;
   float aMax = a[0];
   for (int j = 1; j < 3; j++)
   {
      if (a[j] > aMax)
      {
         aMax = a[j];
         i = j;
      }
   }
   int ind[5] = {0,1,2,0,1};
   float b[3] = {0,0,0};
   if ( (a[ind[i+1]] == 0) & (a[ind[i+2]] == 0) )
   {
      b[ind[i+1]] = 1;
   }
   else if ( a[ind[i+1]] == 0 )
   {
      b[ind[i+1]] = 1;
   }
   else if ( a[ind[i+2]] == 0 )
   {
      b[ind[i+2]] = 1;
   }
   else
   {
      b[ind[i+1]] = 1;
      b[ind[i+2]] = -a[ind[i+1]]/a[ind[i+2]];
      float nb = sqrt(b[0]*b[0] + b[1]*b[1] + b[2]*b[2]);
      b[0] = b[0]/nb; b[1] = b[1]/nb; b[2] = b[2]/nb;
   }
   float c[3] = {a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]}; // c = cross(a,b)
   // Transformation matrix
   float B[3*3] = {a[0],b[0],c[0],a[1],b[1],c[1],a[2],b[2],c[2]};
   // Transform the points
   for (int i = 0; i < (n+1)*2; i++)
   {
     float wkp[3] = {(float)px[i],(float)py[i],(float)pz[i]};
      float wkpD[3];
      for (int j = 0; j < 3; j++)
      {
         wkpD[j] = B[3*j]*wkp[0] + B[3*j+1]*wkp[1] + B[3*j+2]*wkp[2];
         wkpD[j] += orig[j];
      }
      px[i] = wkpD[0]; py[i] = wkpD[1]; pz[i] = wkpD[2];
   }
   
   // Assemble polydata ==========
   
   vtkSmartPointer<vtkPoints> pnts = vtkSmartPointer<vtkPoints>::New();
   for (int i = 0; i < (n+1)*2; i++)
   {
      pnts->InsertNextPoint(px[i],py[i],pz[i]);
   }
   cyl->SetPoints(pnts);
   vtkSmartPointer<vtkQuad> face;
   cyl->Allocate();
   for (int i = 0; i < n; i++)
   {
     face = vtkSmartPointer<vtkQuad>::New();
      face->GetPointIds()->SetId(0,i);
      face->GetPointIds()->SetId(1,i+1);
      face->GetPointIds()->SetId(2,i+n+2);
      face->GetPointIds()->SetId(3,i+n+1);
      
      cyl->InsertNextCell(face->GetCellType(), face->GetPointIds());
   }

   delete[] px;
   delete[] py;
   delete[] pz;
   
   return cyl;
}
