// =========================================================================
// File:       tledKeyPressInteractorStyle.cpp
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

#include "tledKeyPressInteractorStyle.h"
#include "tledHelper.h"

#include <vtkProperty.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkWindowToImageFilter.h>
#include <vtkPNGWriter.h>
#include <vtkRenderWindow.h>

using namespace std;

void tledKeyPressInteractorStyle::OnChar()
{
   vtkRenderWindowInteractor *rwi = this->Interactor;
   char ch = rwi->GetKeyCode() ;
   
   bool vis;
   double opacity;
   bool warp;
   switch (ch)
   {
      case 'u':
         warp = modelWV->GetScaleFactor() != 0.0;
         modelWV->SetScaleFactor((double)!warp);
	 modelWV->Update();
	 if (membraneWV != NULL) {
	   membraneWV->SetScaleFactor(!warp? 1.0 : 0.0);
	   membraneWV->Update();
	 }
         for (int i = 0; i < (int)cylActors.size(); i++) {
	   cylsWV[i]->SetScaleFactor(!warp? 1.0 : 0.0);
	   cylsWV[i]->Update();
	 }
         for (int i = 0; i < (int)prbActors.size(); i++) {
	   prbsWV[i]->SetScaleFactor(!warp? 1.0 : 0.0);
	   prbsWV[i]->Update();
	 }
         for (int i = 0; i < (int)pltActors.size(); i++) {
	   pltsWV[i]->SetScaleFactor(!warp? 1.0 : 0.0);
	   pltsWV[i]->Update();
	 }
	 for (int i = 0; i < int(rigidSurfWV.size()); i++) {
	   rigidSurfWV[i]->SetScaleFactor(!warp? 1.0 : 0.0);
	   rigidSurfWV[i]->Update();
	 }
         break;
      case 'x': // Flip about x-axis
         modelActor->RotateX(180);
	 if (membraneActor != NULL) membraneActor->RotateX(180);
         for (int i = 0; i < (int)cylActors.size(); i++)
            cylActors[i]->RotateX(180);
         for (int i = 0; i < (int)prbActors.size(); i++)
            prbActors[i]->RotateX(180);
         for (int i = 0; i < (int)pltActors.size(); i++)
            pltActors[i]->RotateX(180);
	 for (int i = 0; i < int(rigidSurfWV.size()); i++) {
	   rigidSurfaceActors[i]->RotateX(180);
	 }
         break;
      case 'y': // Flip about y-axis
         modelActor->RotateY(180);
	 if (membraneActor != NULL) membraneActor->RotateY(180);
         for (int i = 0; i < (int)cylActors.size(); i++)
            cylActors[i]->RotateY(180);
         for (int i = 0; i < (int)prbActors.size(); i++)
            prbActors[i]->RotateY(180);
         for (int i = 0; i < (int)pltActors.size(); i++)
            pltActors[i]->RotateY(180);
	 for (int i = 0; i < int(rigidSurfWV.size()); i++) {
	   rigidSurfaceActors[i]->RotateY(180);
	 }
         break;
      case 'z': // Flip about z-axis
         modelActor->RotateZ(180);
	 if (membraneActor != NULL) membraneActor->RotateZ(180);
         for (int i = 0; i < (int)cylActors.size(); i++)
            cylActors[i]->RotateZ(180);
         for (int i = 0; i < (int)prbActors.size(); i++)
            prbActors[i]->RotateZ(180);
         for (int i = 0; i < (int)pltActors.size(); i++)
            pltActors[i]->RotateZ(180);
	 for (int i = 0; i < int(rigidSurfWV.size()); i++) {
	   rigidSurfaceActors[i]->RotateZ(180);
	 }
         break;
      case 'e': // Toggle model edge visibility
         vis = modelActor->GetProperty()->GetEdgeVisibility() != 0;
         modelActor->GetProperty()->SetEdgeVisibility((int)!vis);
	 if (membraneActor != NULL) membraneActor->GetProperty()->SetEdgeVisibility(int(!vis));
         break;
      case 'v': // Toggle model visibility
         vis = modelActor->GetVisibility() != 0;
         modelActor->SetVisibility((int)!vis);
	 if (membraneActor != NULL) membraneActor->SetVisibility(int(!vis));
         break;
      case 'o': // Toggle model opacity
         opacity = modelActor->GetProperty()->GetOpacity();
         if (opacity < 1.0)
            opacity = 1.0;
         else
            opacity = 0.5;
         modelActor->GetProperty()->SetOpacity(opacity);
	 if (membraneActor != NULL) membraneActor->GetProperty()->SetOpacity(opacity);
         break;
      case 'w': // Make model wireframe
         modelActor->GetProperty()->SetRepresentationToWireframe();
	 if (membraneActor != NULL) membraneActor->GetProperty()->SetRepresentationToWireframe();
         break;
      case 's': // Make model solid
         modelActor->GetProperty()->SetRepresentationToSurface();
	 if (membraneActor != NULL) membraneActor->GetProperty()->SetRepresentationToSurface();
         break;
      case 'c': // Capture image
         CaptureImage(rwi);
         break;
      case 'q': // Quit
         vtkInteractorStyleTrackballCamera::OnChar();
         break;
      case 'h': // Show keypress options
         PrintKeyPressOptions();
         break;
      default:
         vtkInteractorStyleTrackballCamera::OnChar();
   }
   rwi->Render();
}

void tledKeyPressInteractorStyle::CaptureImage(vtkSmartPointer<vtkRenderWindowInteractor> rwi)
{
   // Capture
   vtkWindowToImageFilter* windowToImageFilter = vtkWindowToImageFilter::New();
   windowToImageFilter->SetInput(rwi->GetRenderWindow());
   windowToImageFilter->Update();
   // Write
   vtkPNGWriter* writer = vtkPNGWriter::New();
   stringstream filename;
   filename << direc;
   filename << "screenshot.png";
   writer->SetFileName(filename.str().c_str());
   tledVTK6CompatSetInput(writer, windowToImageFilter->GetOutput());
   writer->Write();
   cout << "Image saved to: " << filename.str().c_str() << endl;
   // Clean up
   windowToImageFilter->Delete();
   writer->Delete();
}

void tledKeyPressInteractorStyle::PrintKeyPressOptions(void)
{
   cout << "\n========== Keypress options ==========" << endl;
   cout << "h: Display these options" << endl;
   cout << "u: Toggle between deformed and undeformed shapes" << endl;
   cout << "x: Rotate about the x-axis" << endl;
   cout << "y: Rotate about the y-axis" << endl;
   cout << "z: Rotate about the z-axis" << endl;
   cout << "w: Make model representation wireframe" << endl;
   cout << "s: Make model representation surface" << endl;
   cout << "v: Toggle model visibility" << endl;
   cout << "o: Toggle model opacity" << endl;
   cout << "e: Toggle model edge visibility" << endl;
   cout << "c: Capture an image of the current view" << endl;
   cout << "q: Quit" << endl;
   cout << endl;
}

#endif // _Visualisation_

