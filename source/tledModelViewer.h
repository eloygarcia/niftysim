// =========================================================================
// File:       tledModelViewer.h
// Purpose:    Model visualisation (VTK-based)
// Package:    NiftySim: Nonlinear finite element program
// Author:     Zeike Taylor
// Language:   C++
// Created:    September 2010
// 
// Copyright (c) 2010, University of Queensland. All rights reserved.
// MedTeQ Centre
// See the LICENSE.txt file in the root folder
// 
// ztaylor@itee.uq.edu.au
// =========================================================================

#ifdef _Visualisation_

#ifndef tledModelViewer_H
#define tledModelViewer_H

#include "tledHelper.h"
#include "tledSolver.h"
#include "tledModel.h"
#include "tledSimulator.h"
#ifdef _GPU_
#include "tledContactManager.h"
#endif // _GPU_

#include <vtkUnstructuredGrid.h>
#include <vtkPolyData.h>
#include <vtkWarpVector.h>
#include <vtkSmartPointer.h>

class tledModelViewer
{
public:
   tledModelViewer(void);
   tledModelViewer(tledSimulator* simulator);
   ~tledModelViewer(void);

   void DisplayModel(void);
   void CreateVTKModels();
   
private:
   vtkSmartPointer<vtkUnstructuredGrid> mdlGrid;
   vtkSmartPointer<vtkPolyData> mp_Membrane;
   vtkSmartPointer<vtkWarpVector> mp_MembraneWarpVector;
   tledSimulator* Simulator;
   tledSolver* Solver;
   tledModel* Model;
   vtkSmartPointer<vtkWarpVector> mdlWarpVec;
#ifdef _GPU_
   tledContactManager* Contacts;
#endif // _GPU_
   int numContactCyls;
   int numContactPrbs;
   int numContactPlts;
   std::vector<vtkSmartPointer<vtkPolyData> > contactCylsPolyData;
   std::vector<vtkSmartPointer<vtkPolyData> > contactPrbsPolyData;
   std::vector<vtkSmartPointer<vtkPolyData> > contactPltsPolyData;
   std::vector<vtkSmartPointer<vtkWarpVector> > cylsWarpVec;
   std::vector<vtkSmartPointer<vtkWarpVector> > prbsWarpVec;
   std::vector<vtkSmartPointer<vtkWarpVector> > pltsWarpVec;
};

#endif // tledModelViewer_H

#endif // _Visualisation_
