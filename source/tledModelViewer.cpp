// =========================================================================
// File:       tledModelViewer.cpp
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

#include "tledModelViewer.h"
#include "tledVTKUnstructuredContactSurfaceSource.h"
#include "tledMesh.h"
#include "tledKeyPressInteractorStyle.h"
#include "tledVTKMeshSource.h"
#include "tledVTKCylSource.h"
#include "tledVTKUltrasoundProbeSource.h"
#include "tledVTKPlateSource.h"
#include "tledShellSolver.h"

#if VTK_MAJOR_VERSION >= 6
#include <vtkOpenGLRenderer.h>
#endif
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkActor.h>
#include <vtkDataSetMapper.h>
#include <vtkProperty.h>
#include <vtkPoints.h>
#include <vtkTetra.h>
#include <vtkHexahedron.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>

using namespace std;

tledModelViewer::tledModelViewer(void)
{
   mdlGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
   mdlWarpVec = vtkSmartPointer<vtkWarpVector>::New();
   Solver = NULL;
   Model = NULL;
   mp_Membrane = NULL;
   mp_MembraneWarpVector = NULL;
   numContactCyls = 0;
   numContactPrbs = 0;
   numContactPlts = 0;
}

tledModelViewer::tledModelViewer(tledSimulator* simulator)
{
   Simulator = simulator;
   Solver = Simulator->GetSolver();
   Model = Simulator->GetModel();
   mdlGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
   mdlWarpVec = vtkSmartPointer<vtkWarpVector>::New();
   mp_Membrane = NULL;
   mp_MembraneWarpVector = NULL;
#ifdef _GPU_
   Contacts = Simulator->GetContactManager();
#endif // _GPU_
   numContactCyls = 0;
   numContactPrbs = 0;
   numContactPlts = 0;
}

tledModelViewer::~tledModelViewer(void)
{
}

void tledModelViewer::DisplayModel(void)
{
   vtkSmartPointer<tledKeyPressInteractorStyle> iStyle = vtkSmartPointer<tledKeyPressInteractorStyle>::New();

   // Model **********
   
   //create mapper
   vtkSmartPointer<vtkDataSetMapper> mdlMapper = vtkSmartPointer<vtkDataSetMapper>::New();
   mdlMapper->ImmediateModeRenderingOn();
   tledVTK6CompatSetInput(mdlMapper, mdlGrid);
      
   // create actor
   vtkSmartPointer<vtkActor> mdlActor = vtkSmartPointer<vtkActor>::New();
   mdlActor->SetMapper(mdlMapper);
   mdlActor->SetVisibility(1);
   mdlActor->GetProperty()->SetColor(1,0.8,0.8);
   mdlActor->GetProperty()->SetEdgeColor(0,0,0);
   mdlActor->GetProperty()->EdgeVisibilityOn();
   
   // Contact cylinders **********

   vector<vtkSmartPointer<vtkPolyDataMapper> > cylMappers;
   cylMappers.resize(numContactCyls);
   vector<vtkSmartPointer<vtkActor> > cylActors;
   cylActors.resize(numContactCyls);
   for (int cylNum = 0; cylNum < numContactCyls; cylNum++)
   {
      // create mappers
      cylMappers[cylNum] = vtkSmartPointer<vtkPolyDataMapper>::New();
      cylMappers[cylNum]->ImmediateModeRenderingOn();
      cylMappers[cylNum]->SetResolveCoincidentTopology(1);
      tledVTK6CompatSetInput(cylMappers[cylNum], contactCylsPolyData[cylNum]);
      
      // create actors
      cylActors[cylNum] = vtkSmartPointer<vtkActor>::New();
      cylActors[cylNum]->SetMapper(cylMappers[cylNum]);
      cylActors[cylNum]->SetVisibility(1);
      cylActors[cylNum]->GetProperty()->SetColor(0.8,0.8,1);
      cylActors[cylNum]->GetProperty()->SetEdgeColor(0,0,0);
      cylActors[cylNum]->GetProperty()->EdgeVisibilityOff();
      cylActors[cylNum]->GetProperty()->SetOpacity(0.8);
   }
   
   // Contact US probes **********
   
   vector<vtkSmartPointer<vtkPolyDataMapper> > prbMappers;
   prbMappers.resize(numContactPrbs);
   vector<vtkSmartPointer<vtkActor> > prbActors;
   prbActors.resize(numContactPrbs);
   for (int prbNum = 0; prbNum < numContactPrbs; prbNum++)
   {
      // create mappers
      prbMappers[prbNum] = vtkSmartPointer<vtkPolyDataMapper>::New();
      prbMappers[prbNum]->ImmediateModeRenderingOn();
      prbMappers[prbNum]->SetResolveCoincidentTopology(1);
      tledVTK6CompatSetInput(prbMappers[prbNum], contactPrbsPolyData[prbNum]);
      
      // create actors
      prbActors[prbNum] = vtkSmartPointer<vtkActor>::New();
      prbActors[prbNum]->SetMapper(prbMappers[prbNum]);
      prbActors[prbNum]->SetVisibility(1);
      prbActors[prbNum]->GetProperty()->SetColor(0.8,0.8,1);
      prbActors[prbNum]->GetProperty()->SetEdgeColor(0,0,0);
      prbActors[prbNum]->GetProperty()->EdgeVisibilityOff();
      prbActors[prbNum]->GetProperty()->SetOpacity(0.8);
   }
   
   // Contact plates **********
   
   vector<vtkSmartPointer<vtkPolyDataMapper> > pltMappers;
   pltMappers.resize(numContactPlts);
   vector<vtkSmartPointer<vtkActor> > pltActors;
   pltActors.resize(numContactPlts);
   for (int pltNum = 0; pltNum < numContactPlts; pltNum++)
   {
      // create mappers
      pltMappers[pltNum] = vtkSmartPointer<vtkPolyDataMapper>::New();
      pltMappers[pltNum]->ImmediateModeRenderingOn();
      pltMappers[pltNum]->SetResolveCoincidentTopology(1);
      tledVTK6CompatSetInput(pltMappers[pltNum], contactPltsPolyData[pltNum]);
      
      // create actors
      pltActors[pltNum] = vtkSmartPointer<vtkActor>::New();
      pltActors[pltNum]->SetMapper(pltMappers[pltNum]);
      pltActors[pltNum]->SetVisibility(1);
      pltActors[pltNum]->GetProperty()->SetColor(0.8,0.8,1);
      pltActors[pltNum]->GetProperty()->SetEdgeColor(0,0,0);
      pltActors[pltNum]->GetProperty()->EdgeVisibilityOff();
      pltActors[pltNum]->GetProperty()->SetOpacity(0.8);
   }
   
   // Render **********
   
   // a renderer and render window
#if VTK_MAJOR_VERSION >= 6
   vtkSmartPointer<vtkOpenGLRenderer> Renderer = vtkSmartPointer<vtkOpenGLRenderer>::New();
#else
   vtkSmartPointer<vtkRenderer> Renderer = vtkSmartPointer<vtkRenderer>::New();
#endif
   vtkSmartPointer<vtkRenderWindow> RenderWindow = vtkSmartPointer<vtkRenderWindow>::New();
   RenderWindow->AddRenderer(Renderer);
   
   // an interactor
   vtkSmartPointer<vtkRenderWindowInteractor> RenderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
   RenderWindowInteractor->SetRenderWindow(RenderWindow);

   /*
    * Create and add actors for rigid contact surfaces.
    */
   if (Simulator->GetContactManager()->DoUnstructuredContacts()) {
     const tledUnstructuredContactManager &unstructuredMeshContacts = Simulator->GetContactManager()->GetUnstructuredContactManager();

     std::vector<vtkSmartPointer<vtkWarpVector> > vsp_rigidSurfaceWarpVecs;
     std::vector<vtkSmartPointer<vtkActor> > vsp_rigidSurfaceActors;

     Simulator->GetContactManager()->GetUnstructuredContactManager().ResetRigidSurfaces();
     for (int sInd = 0; sInd < unstructuredMeshContacts.GetNumberOfRigidSurfaces(); sInd++) {
       tledVTKUnstructuredContactSurfaceSource converter;
       vtkSmartPointer<vtkActor> sp_actor;
       vtkSmartPointer<vtkDataSetMapper> sp_mapper;

       sp_actor = vtkSmartPointer<vtkActor>::New();
       sp_actor->GetProperty()->SetColor(0.7, 0.7, 0.7);
       sp_actor->GetProperty()->SetEdgeColor(0, 0, 0);
       sp_actor->GetProperty()->EdgeVisibilityOn();
       sp_actor->GetProperty()->SetOpacity(0.7);

       sp_mapper = vtkSmartPointer<vtkDataSetMapper>::New();
       sp_mapper->ImmediateModeRenderingOn();

       converter.SetInput(unstructuredMeshContacts.GetRigidSurface<tledRigidContactSurface>(sInd));
       converter.Update();
       if (unstructuredMeshContacts.GetRigidSurface<tledRigidContactSurface>(sInd).IsMoving()) {
	 std::vector<float> disps(3*unstructuredMeshContacts.GetRigidSurface<tledRigidContactSurface>(sInd).GetNumberOfNodes());

	 unstructuredMeshContacts.GetRigidSurface<tledRigidContactSurface>(sInd).GetFinalDisplacement(&disps.front());
	 converter.AddNodeVectorAttribute("Displacements", &disps.front());
	 converter.GetOutput()->GetPointData()->SetActiveVectors("Displacements");

	 vsp_rigidSurfaceWarpVecs.push_back(vtkSmartPointer<vtkWarpVector>::New());
	 tledVTK6CompatSetInput(vsp_rigidSurfaceWarpVecs.back(), converter.GetOutput());
	 vsp_rigidSurfaceWarpVecs.back()->SetScaleFactor(1);
	 vsp_rigidSurfaceWarpVecs.back()->Update();
	 tledVTK6CompatSetInput(sp_mapper, vsp_rigidSurfaceWarpVecs.back()->GetPolyDataOutput());
       } else {
	 tledVTK6CompatSetInput(sp_mapper, converter.GetOutput());
       }

       sp_actor->SetMapper(sp_mapper);
       sp_actor->SetVisibility(1);

       Renderer->AddActor(sp_actor);
       vsp_rigidSurfaceActors.push_back(sp_actor);
     }
     iStyle->SetRigidSurfaceActors(vsp_rigidSurfaceActors);
     iStyle->SetRigidSurfaceWarpVectors(vsp_rigidSurfaceWarpVecs);
   }

   // Membrane **********
   if (mp_Membrane != NULL) {
     vtkSmartPointer<vtkDataSetMapper> memMapper = vtkSmartPointer<vtkDataSetMapper>::New();
     vtkSmartPointer<vtkActor> memActor = vtkSmartPointer<vtkActor>::New();

     memMapper->ImmediateModeRenderingOn();
     tledVTK6CompatSetInput(memMapper, mp_Membrane);
     memActor->SetMapper(memMapper);
     memActor->SetVisibility(1);
     memActor->GetProperty()->SetColor(1,0.8,0.8);
     memActor->GetProperty()->SetEdgeColor(0,0,0);
     memActor->GetProperty()->EdgeVisibilityOn();
     
     Renderer->AddActor(memActor);
     iStyle->SetMembraneActor(memActor);
   }

   // add the actors to the scene
   Renderer->AddActor(mdlActor);
   for (int cylNum = 0; cylNum < numContactCyls; cylNum++)
      Renderer->AddActor(cylActors[cylNum]);
   for (int prbNum = 0; prbNum < numContactPrbs; prbNum++)
      Renderer->AddActor(prbActors[prbNum]);
   for (int pltNum = 0; pltNum < numContactPlts; pltNum++)
      Renderer->AddActor(pltActors[pltNum]);
   
//    Renderer->SetBackground(1,1,1); // Background color white
   Renderer->SetBackground(0,0,0); // Background color black
   
   // render an image (lights and cameras are created automatically)
   RenderWindow->Render();
   
   // interactor style
   iStyle->SetModelActor(mdlActor);
   iStyle->SetModelWarpVector(mdlWarpVec);

   if (mp_MembraneWarpVector != NULL) iStyle->SetMembraneWarpVector(mp_MembraneWarpVector);

   if (numContactCyls > 0)
   {
      iStyle->SetCylActors(cylActors);
      iStyle->SetCylWarpVectors(cylsWarpVec);
   }
   if (numContactPrbs > 0)
   {
      iStyle->SetPrbActors(prbActors);
      iStyle->SetPrbWarpVectors(prbsWarpVec);
   }
   if (numContactPlts > 0)
   {
      iStyle->SetPltActors(pltActors);
      iStyle->SetPltWarpVectors(pltsWarpVec);
   }
   iStyle->SetDirec(Model->GetDirectory());
 
   RenderWindowInteractor->SetInteractorStyle(iStyle);
 
   // begin mouse interaction
   cout << "Press \"h\" for keypress options" << endl;

   RenderWindowInteractor->Start();
}

void tledModelViewer::CreateVTKModels()
{
   if (!Simulator)
   {
     tledFatalError("Simulator not set");
     return;
   }
   
   // CREATE MODEL GRID ********************
   
   tledMesh* Mesh = Solver->GetMesh();
   
   // Create grid
   tledVTKMeshSource mtg;
   mtg.SetInput(*Mesh);
   mtg.Update();
   // Add displacement vectors
   mtg.AddNodeVectorAttribute("Displacements", Solver->GetAllDisps());
   vtkSmartPointer<vtkUnstructuredGrid> initGrid = mtg.GetOutput();
   initGrid->GetPointData()->SetActiveVectors("Displacements");
   
   tledVTK6CompatSetInput(mdlWarpVec, initGrid);
   mdlWarpVec->SetScaleFactor(1);
   mdlWarpVec->Update();
   mdlGrid = mdlWarpVec->GetUnstructuredGridOutput();

   /* 
    * Membranes
    */
   if (Model->GetNumberOfShellElementSets() > 0) {
     tledVTKSurfaceSource vtkAdapter;
     vtkSmartPointer<vtkPolyData> sp_initMem;

     /* Nodes vector is shared with main mesh, hence we can largely apply the same code as for the mesh */
     vtkAdapter.SetInput(Simulator->GetSolver()->GetShellSolver().GetSurface());
     vtkAdapter.Update();
     vtkAdapter.AddNodeVectorAttribute("Displacements", Solver->GetAllDisps());
     sp_initMem = vtkAdapter.GetOutput();
     sp_initMem->GetPointData()->SetActiveVectors("Displacements");

     mp_MembraneWarpVector = vtkSmartPointer<vtkWarpVector>::New();     
     tledVTK6CompatSetInput(mp_MembraneWarpVector, sp_initMem);
     mp_MembraneWarpVector->SetScaleFactor(1);
     mp_MembraneWarpVector->Update();
     mp_Membrane = mp_MembraneWarpVector->GetPolyDataOutput();     
   }
   
#ifdef _GPU_
   // CREATE CYLS DATA ********************
   
   numContactCyls = Contacts->GetNumContactCyls();
   contactCylsPolyData.resize(numContactCyls);
   cylsWarpVec.resize(numContactCyls);
   for (int cylNum = 0; cylNum < numContactCyls; cylNum++)
   {
      // Create initial polydata
      tledContactCylinder* Cyl = Contacts->GetContactCyl(cylNum);
      vector<float> orig = Cyl->GetStartOriginV();
      vector<float> axis = Cyl->GetStartAxisV();
      float R = Cyl->GetStartRadius();
      float L = Cyl->GetStartLength();
      tledVTKCylSource* cylSrc = new tledVTKCylSource(orig,axis,R,L,40);
      vtkSmartPointer<vtkPolyData> initData = cylSrc->GetOutput();
      int numPnts = initData->GetNumberOfPoints();
      
      // Add displacement vectors
      // (create a deformed cylinder then subtract the original point coords)
      vector<float> origDisps = Cyl->GetOriginDispV();
      orig[0] += origDisps[0]; orig[1] += origDisps[1]; orig[2] += origDisps[2];
      R += Cyl->GetRadiusChng();
      tledVTKCylSource* cylSrcDef = new tledVTKCylSource(orig,axis,R,L,40);
      vtkSmartPointer<vtkPolyData> defData = cylSrcDef->GetOutput();
      vtkSmartPointer<vtkDoubleArray> cylDispVec = vtkSmartPointer<vtkDoubleArray>::New();
      cylDispVec->SetNumberOfComponents(3);
      double p0[3];
      double p1[3];
      double d[3];
      for (int i = 0; i < numPnts; i++)
      {
         initData->GetPoint(i,&p0[0]);
         defData->GetPoint(i,&p1[0]);
         d[0] = p1[0] - p0[0]; d[1] = p1[1] - p0[1]; d[2] = p1[2] - p0[2];
         cylDispVec->InsertNextTuple(&d[0]);
      }
      initData->GetPointData()->SetVectors(cylDispVec);
      
      cylsWarpVec[cylNum] = vtkSmartPointer<vtkWarpVector>::New();
      tledVTK6CompatSetInput(cylsWarpVec[cylNum], initData);
      cylsWarpVec[cylNum]->SetScaleFactor(1);
      contactCylsPolyData[cylNum] = cylsWarpVec[cylNum]->GetPolyDataOutput();

      delete cylSrcDef;
      delete cylSrc;
   }
   
   // CREATE US PROBES DATA ********************
   
   numContactPrbs = Contacts->GetNumContactPrbs();
   contactPrbsPolyData.resize(numContactPrbs);
   prbsWarpVec.resize(numContactPrbs);
   for (int prbNum = 0; prbNum < numContactPrbs; prbNum++)
   {
      // Create initial polydata
      tledContactUSProbe* Prb = Contacts->GetContactPrb(prbNum);
      vector<float> orig = Prb->GetStartOriginV();
      vector<float> axis = Prb->GetStartAxisV();
      float R = Prb->GetStartRadius();
      float L = Prb->GetStartLength();
      tledVTKUltrasoundProbeSource* prbSrc = new tledVTKUltrasoundProbeSource(orig,axis,R,L,40,16);
      vtkSmartPointer<vtkPolyData> initData = prbSrc->GetOutput();
      int numPnts = initData->GetNumberOfPoints();
               
      // Add displacement vectors
      // (create a deformed probe then subtract the original point coords)
      vector<float> origDisps = Prb->GetOriginDispV();
      orig[0] += origDisps[0]; orig[1] += origDisps[1]; orig[2] += origDisps[2];
      R += Prb->GetRadiusChng();
      tledVTKUltrasoundProbeSource* prbSrcDef = new tledVTKUltrasoundProbeSource(orig,axis,R,L,40,16);
      vtkSmartPointer<vtkPolyData> defData = prbSrcDef->GetOutput();
      vtkSmartPointer<vtkDoubleArray> prbDispVec = vtkSmartPointer<vtkDoubleArray>::New();
      prbDispVec->SetNumberOfComponents(3);
      double p0[3];
      double p1[3];
      double d[3];
      for (int i = 0; i < numPnts; i++)
      {
         initData->GetPoint(i,&p0[0]);
         defData->GetPoint(i,&p1[0]);
         d[0] = p1[0] - p0[0]; d[1] = p1[1] - p0[1]; d[2] = p1[2] - p0[2];
         prbDispVec->InsertNextTuple(&d[0]);
      }
      initData->GetPointData()->SetVectors(prbDispVec);
      
      prbsWarpVec[prbNum] = vtkSmartPointer<vtkWarpVector>::New();
      tledVTK6CompatSetInput(prbsWarpVec[prbNum], initData);
      prbsWarpVec[prbNum]->SetScaleFactor(1);
      contactPrbsPolyData[prbNum] = prbsWarpVec[prbNum]->GetPolyDataOutput();

      delete prbSrcDef;
      delete prbSrc;
   }
   
   // CREATE PLATES DATA ********************
   
   numContactPlts = Contacts->GetNumContactPlts();
   contactPltsPolyData.resize(numContactPlts);
   pltsWarpVec.resize(numContactPlts);
   for (int pltNum = 0; pltNum < numContactPlts; pltNum++)
   {
      // Create initial polydata
      tledContactPlate* Plt = Contacts->GetContactPlt(pltNum);
      vector<float> a = Plt->GetStartCrnrAV();
      vector<float> b = Plt->GetStartCrnrBV();
      vector<float> c = Plt->GetStartCrnrCV();
      tledVTKPlateSource* pltSrc = new tledVTKPlateSource(a,b,c);
      vtkSmartPointer<vtkPolyData> initData = pltSrc->GetOutput();
      int numPnts = initData->GetNumberOfPoints();

      // Add displacement vectors
      // (create a deformed plate then subtract the original point coords)
      vector<float> Disps = Plt->GetDispV();
      a[0] += Disps[0]; a[1] += Disps[1]; a[2] += Disps[2];
      b[0] += Disps[0]; b[1] += Disps[1]; b[2] += Disps[2];
      c[0] += Disps[0]; c[1] += Disps[1]; c[2] += Disps[2];
      tledVTKPlateSource* pltSrcDef = new tledVTKPlateSource(a,b,c);
      vtkSmartPointer<vtkPolyData> defData = pltSrcDef->GetOutput();
      vtkSmartPointer<vtkDoubleArray> pltDispVec = vtkSmartPointer<vtkDoubleArray>::New();
      pltDispVec->SetNumberOfComponents(3);
      double p0[3];
      double p1[3];
      double d[3];
      for (int i = 0; i < numPnts; i++)
      {
         initData->GetPoint(i,&p0[0]);
         defData->GetPoint(i,&p1[0]);
         d[0] = p1[0] - p0[0]; d[1] = p1[1] - p0[1]; d[2] = p1[2] - p0[2];
         pltDispVec->InsertNextTuple(&d[0]);
      }
      initData->GetPointData()->SetVectors(pltDispVec);
      
      pltsWarpVec[pltNum] = vtkSmartPointer<vtkWarpVector>::New();
      tledVTK6CompatSetInput(pltsWarpVec[pltNum], initData);
      pltsWarpVec[pltNum]->SetScaleFactor(1);
      contactPltsPolyData[pltNum] = pltsWarpVec[pltNum]->GetPolyDataOutput();

      delete pltSrc;
      delete pltSrcDef;
   }   
#endif // _GPU_
}


#endif // _Visualisation_
