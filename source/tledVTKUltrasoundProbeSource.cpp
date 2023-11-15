// =========================================================================
// File:       tledVTKUltrasoundProbeSource.cpp
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

#include "tledVTKUltrasoundProbeSource.h"

#include <vtkQuad.h>

using namespace std;

tledVTKUltrasoundProbeSource::tledVTKUltrasoundProbeSource(float* a, float* A, float R, float L, int Nth, int Nph)
{
   orig = new float[3];
   memcpy(orig,a,sizeof(float)*3);
   axis = new float[3];
   memcpy(axis,A,sizeof(float)*3);
   rad = R;
   len = L;
   nth = Nth;
   nph = Nph;
   prb = vtkPolyData::New();
}

tledVTKUltrasoundProbeSource::tledVTKUltrasoundProbeSource(vector<float> a, vector<float> A, float R, float L, int Nth, int Nph)
{
   orig = new float[3];
   orig[0] = a[0]; orig[1] = a[1]; orig[2] = a[2];
   axis = new float[3];
   axis[0] = A[0]; axis[1] = A[1]; axis[2] = A[2];
   rad = R;
   len = L;
   nth = Nth;
   nph = Nph;
   prb = vtkSmartPointer<vtkPolyData>::New();
}

tledVTKUltrasoundProbeSource::~tledVTKUltrasoundProbeSource()
{
   delete[] orig;
   delete[] axis;
}

vtkSmartPointer<vtkPolyData> tledVTKUltrasoundProbeSource::GetOutput()
{
   // ******************** CONSTRUCT THE CYLINDRICAL SHAFT ********************
   
   // Construct points in cylinder local coords ==========
   
   double pi = 3.141592653589793;
   double g = 2*pi/nth;
   double th;
   double* px = new double[(nth+1)*2];
   double* py = new double[(nth+1)*2];
   double* pz = new double[(nth+1)*2];
   for (int i = 0; i < nth+1; i++)
   {
      th = i*g;
      px[i] = 0; px[i+nth+1] = len;
      py[i] = rad*cos(th); py[i+nth+1] = py[i];
      pz[i] = rad*sin(th); pz[i+nth+1] = pz[i];
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
   for (int i = 0; i < (nth+1)*2; i++)
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
   for (int i = 0; i < (nth+1)*2; i++)
      pnts->InsertNextPoint(px[i],py[i],pz[i]);
   prb->SetPoints(pnts);
   vtkSmartPointer<vtkQuad> face;
   prb->Allocate();
   for (int i = 0; i < nth; i++)
   {
      face = vtkSmartPointer<vtkQuad>::New();
      face->GetPointIds()->SetId(0,i);
      face->GetPointIds()->SetId(1,i+1);
      face->GetPointIds()->SetId(2,i+nth+2);
      face->GetPointIds()->SetId(3,i+nth+1);
      
      prb->InsertNextCell(face->GetCellType(), face->GetPointIds());
   }
   
   int numCylPnts = pnts->GetNumberOfPoints();
   
   // ******************** CONSTRUCT THE HEMISPHERE CAP ********************
   
   float dth = 2*pi/nth;
   float dph = pi/2/nph;
   float ph;
   delete[] px; delete[] py; delete[] pz;
   px = new double[(nth+1)*nph+1];
   py = new double[(nth+1)*nph+1];
   pz = new double[(nth+1)*nph+1];
   int cntr = 0;
   for (int i = 0; i < nph; i++)
   {
      ph = dph*i;
      for (int j = 0; j < nth+1; j++)
      {
         th = dth*j;
         px[cntr] = rad*sin(ph);
         py[cntr] = rad*cos(ph)*sin(th);
         pz[cntr] = rad*cos(ph)*cos(th);
         cntr++;
      }
   }
   px[cntr] = rad; py[cntr] = 0; pz[cntr] = 0;
   
   // Transform points to global coords ==========
   
   // Transform the points
   float offset[3] = {orig[0]+len*axis[0],orig[1]+len*axis[1],orig[2]+len*axis[2]};
   for (int i = 0; i < (nth+1)*nph+1; i++)
   {
      float wkp[3] = {(float)px[i],(float)py[i],(float)pz[i]};
      float wkpD[3];
      for (int j = 0; j < 3; j++)
      {
         wkpD[j] = B[3*j]*wkp[0] + B[3*j+1]*wkp[1] + B[3*j+2]*wkp[2];
         wkpD[j] += offset[j];
      }
      px[i] = wkpD[0]; py[i] = wkpD[1]; pz[i] = wkpD[2];
   }
   
   // Assemble polydata ==========
   
   for (int i = 0; i < (nth+1)*nph+1; i++)
      pnts->InsertNextPoint(px[i],py[i],pz[i]);
   prb->SetPoints(pnts);
   int baseNd = 0;
   for (int i = 0; i < nph-1; i++) // Faces, layers 1 to end-1
   {
      for (int j = 0; j < nth; j++)
      {
	 face = vtkSmartPointer<vtkQuad>::New();
         // NB: hemisphere point Ids should begin where the cyl Ids left off
         // --> hence the offset by numCylPnts
         face->GetPointIds()->SetId(0,numCylPnts+baseNd+j);
         face->GetPointIds()->SetId(1,numCylPnts+baseNd+j+1);
         face->GetPointIds()->SetId(2,numCylPnts+baseNd+j+nth+2);
         face->GetPointIds()->SetId(3,numCylPnts+baseNd+j+nth+1);
         
         prb->InsertNextCell(face->GetCellType(), face->GetPointIds());
      }
      baseNd += nth+1;
   }
   for (int j = 1; j < nth/2+1; j++) // Faces, last layer
   {
      face = vtkSmartPointer<vtkQuad>::New();
      face->GetPointIds()->SetId(0,numCylPnts+baseNd+2*j-2);
      face->GetPointIds()->SetId(1,numCylPnts+baseNd+2*j-1);
      face->GetPointIds()->SetId(2,numCylPnts+baseNd+2*j);
      face->GetPointIds()->SetId(3,numCylPnts+(nth+1)*nph);
      
      prb->InsertNextCell(face->GetCellType(), face->GetPointIds());
   }
   delete[] px; delete[] py; delete[] pz;
   
   return prb;
}
