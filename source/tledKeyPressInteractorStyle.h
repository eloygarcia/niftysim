// =========================================================================
// File:       tledKeyPressInteractorStyle.h
// Purpose:    Model visualisation (VTK-based)
// Package:    NiftySim: Nonlinear finite element program
// Author:     Zeike Taylor
// Language:   C++
// Created:    November 2009
// 
// Copyright (c) 2010, University of Queensland. All rights reserved.
// MedTeQ Centre
// See the LICENSE.txt file in the root folder
// 
// ztaylor@itee.uq.edu.au
// =========================================================================

#ifdef _Visualisation_

#ifndef tledKeyPressInteractorStyle_H
#define tledKeyPressInteractorStyle_H

#include <vtkActor.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkWarpVector.h>
#include <vtkSmartPointer.h>

#include <string>
#include <vector>

class tledKeyPressInteractorStyle : public vtkInteractorStyleTrackballCamera
{
protected:
  tledKeyPressInteractorStyle(void) : membraneWV(NULL) {}

public:
   static tledKeyPressInteractorStyle* New()
   {return new tledKeyPressInteractorStyle();}
   
   void SetDirec(std::string Direc) {direc = Direc;}
   void SetModelActor(vtkSmartPointer<vtkActor> mactor) {modelActor = mactor;}
   void SetMembraneActor(vtkSmartPointer<vtkActor> memActor) { membraneActor = memActor; }
   void SetCylActors(std::vector<vtkSmartPointer<vtkActor> > cactors) {cylActors = cactors;}
   void SetPrbActors(std::vector<vtkSmartPointer<vtkActor> > pactors) {prbActors = pactors;}
   void SetPltActors(std::vector<vtkSmartPointer<vtkActor> > pactors) {pltActors = pactors;}
   void SetRigidSurfaceActors(std::vector<vtkSmartPointer<vtkActor> > rsActors) {rigidSurfaceActors = rsActors;}
   void SetModelWarpVector(vtkSmartPointer<vtkWarpVector> warp) {modelWV = warp;}
   void SetMembraneWarpVector(vtkSmartPointer<vtkWarpVector> warp) {membraneWV = warp;}
   void SetCylWarpVectors(std::vector<vtkSmartPointer<vtkWarpVector> > warp) {cylsWV = warp;}
   void SetPrbWarpVectors(std::vector<vtkSmartPointer<vtkWarpVector> > warp) {prbsWV = warp;}
   void SetPltWarpVectors(std::vector<vtkSmartPointer<vtkWarpVector> > warp) {pltsWV = warp;}
   void SetRigidSurfaceWarpVectors(std::vector<vtkSmartPointer<vtkWarpVector> > warp) { rigidSurfWV = warp; }
   virtual void OnChar(void);
   
private:
   void CaptureImage(vtkSmartPointer<vtkRenderWindowInteractor> rwi);
   void PrintKeyPressOptions(void);
   
   vtkSmartPointer<vtkActor> modelActor;
   vtkSmartPointer<vtkActor> membraneActor;  
   std::vector<vtkSmartPointer<vtkActor> > cylActors;
   std::vector<vtkSmartPointer<vtkActor> > prbActors;
   std::vector<vtkSmartPointer<vtkActor> > pltActors;
   std::vector<vtkSmartPointer<vtkActor> > rigidSurfaceActors;
   vtkSmartPointer<vtkWarpVector> modelWV;
   vtkSmartPointer<vtkWarpVector> membraneWV;
   std::vector<vtkSmartPointer<vtkWarpVector> > cylsWV;
   std::vector<vtkSmartPointer<vtkWarpVector> > prbsWV;
   std::vector<vtkSmartPointer<vtkWarpVector> > pltsWV;
   std::vector<vtkSmartPointer<vtkWarpVector> > rigidSurfWV;
   std::string direc;
};

#endif // tledKeyPressInteractorStyle_H

#endif // _Visualisation_
