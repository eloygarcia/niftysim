// =========================================================================
// File:       tledContactManager.cu
// Purpose:    Manage all contact objects in the model
// Package:    NiftySim: Nonlinear finite element program
// Author:     Zeike Taylor
// Language:   C++
// Created:    July 2010
// 
// Copyright (c) 2010, University of Queensland. All rights reserved.
// MedTeQ Centre
// See the LICENSE.txt file in the root folder
// 
// ztaylor@itee.uq.edu.au
// =========================================================================

#ifdef _GPU_

#include "tledContactManager.h"
#include "tledDeviceDeclarations.h"
#include "tledCUDAHelpers.h"

using namespace std;

template <class TGPUContact, class TCPUContact>
void _FreeGPUContactMemory(TGPUContact **dpp_gpuContacts, const int szArray) {
  if (szArray == 0) return;
  for (int c = 0; c < szArray; c++) {
    TGPUContact *dp_obj;

    tledCUDAHelpers::CopyFromDevice<TGPUContact*>(&dp_obj, dpp_gpuContacts + c);
    TCPUContact::ReleaseGPUMemory(dp_obj);
  }
  tledCheckCUDAErrors(cudaFree(dpp_gpuContacts));
}

tledContactManager::~tledContactManager()
{
  if (DoUnstructuredContacts()) delete mp_UnstructuredContacts;
  if (d_GPUContacts != NULL) {
    _FreeGPUContactMemory<tledGPUContactCylinder, tledContactCylinder>(h_GPUContacts.Cyls, h_GPUContacts.NumContactCyls);
    _FreeGPUContactMemory<tledGPUContactUSProbe, tledContactUSProbe>(h_GPUContacts.Prbs, h_GPUContacts.NumContactPrbs);
    _FreeGPUContactMemory<tledGPUContactPlate, tledContactPlate>(h_GPUContacts.Plts, h_GPUContacts.NumContactPlts);
    tledCheckCUDAErrors(cudaFree(d_GPUContacts));
  }
}

tledUnstructuredContactManager* tledContactManager::CreateUnstructuredContactManager() {
  return new tledUnstructuredContactManager();
}

void tledContactManager::InitialiseUnstructuredContacts(tledSolver &r_solver, const bool useGPU) {
  if (Model == NULL) {
    tledFatalError("Have no model; cannot initiliase unstructured-mesh contacts.");
  }

  if (this->GetModel().DoDeformableDeformableContactHandling() || this->GetModel().GetNumberOfRigidContactSurfaces() > 0) {
    mp_UnstructuredContacts = this->CreateUnstructuredContactManager();
    mp_UnstructuredContacts->Init(this->GetModel(), r_solver, useGPU);
  }
}

tledContactManager::tledContactManager(tledModel* model) : mp_UnstructuredContacts(NULL)
{
   Model = model;

   // Set up contact cyls
   d_GPUContacts = NULL;
   numContactCyls = Model->GetNumContactCyls();
   h_GPUContacts.NumContactCyls = numContactCyls;
   if (numContactCyls > 0) {
     tledCUDAHelpers::AllocateDeviceMemory(h_GPUContacts.Cyls, numContactCyls);
     Cyls.reserve(numContactCyls);
     for (int cylNum = 0; cylNum < numContactCyls; cylNum++) {
	 tledContactCylinder wkCyl(Model->GetContactCylOrigin(cylNum),
				   Model->GetContactCylAxis(cylNum),
				   Model->GetContactCylRadius(cylNum),
				   Model->GetContactCylLength(cylNum),
				   Model->GetContactCylSlvs(cylNum),
				   Model->GetContactCylOrigDisp(cylNum),
				   Model->GetContactCylRadChange(cylNum),
				   Model->GetNumNodes());
	 tledGPUContactCylinder* d_GPUCyls = wkCyl.GetGPUContactCylinder();


	 Cyls.push_back(wkCyl);
	 tledCUDAHelpers::CopyToDevice<tledGPUContactCylinder*>(h_GPUContacts.Cyls + cylNum, &d_GPUCyls);
     }
   }

   // Set up contact ultrasound probes
   numContactPrbs = Model->GetNumContactPrbs();
   h_GPUContacts.NumContactPrbs = numContactPrbs;
   if (numContactPrbs > 0) {
     tledCUDAHelpers::AllocateDeviceMemory<tledGPUContactUSProbe*>(h_GPUContacts.Prbs, numContactPrbs);
     Prbs.reserve(numContactPrbs);

     for (int prbNum = 0; prbNum < numContactPrbs; prbNum++) {
       tledContactUSProbe wkPrb(Model->GetContactPrbOrigin(prbNum),
				Model->GetContactPrbAxis(prbNum),
				Model->GetContactPrbRadius(prbNum),
				Model->GetContactPrbLength(prbNum),
				Model->GetContactPrbSlvs(prbNum),
				Model->GetContactPrbOrigDisp(prbNum),
				Model->GetContactPrbRadChange(prbNum),
				Model->GetNumNodes());
       tledGPUContactUSProbe* d_GPUPrbs = wkPrb.GetGPUContactUSProbe();

       Prbs.push_back(wkPrb);
       tledCUDAHelpers::CopyToDevice<tledGPUContactUSProbe*>(&h_GPUContacts.Prbs[prbNum], &d_GPUPrbs);
     }
   }

   // Set up contact plates
   numContactPlts = Model->GetNumContactPlts();
   h_GPUContacts.NumContactPlts = numContactPlts;
   if (numContactPlts > 0) {
     tledCUDAHelpers::AllocateDeviceMemory(h_GPUContacts.Plts, numContactPlts);
     Plts.reserve(numContactPlts);
     for (int pltNum = 0; pltNum < numContactPlts; pltNum++) {
       tledContactPlate wkPlt(Model->GetContactPltCrnrA(pltNum),
			      Model->GetContactPltCrnrB(pltNum),
			      Model->GetContactPltCrnrC(pltNum),
			      Model->GetContactPltSlvs(pltNum),
			      Model->GetContactPltDisp(pltNum),
			      Model->GetNumNodes());
       tledGPUContactPlate* d_GPUPlts = wkPlt.GetGPUContactPlate();

       Plts.push_back(wkPlt);
       tledCUDAHelpers::CopyToDevice<tledGPUContactPlate*>(h_GPUContacts.Plts + pltNum, &d_GPUPlts);
     }
   }

   // Construct device pointer
   tledCUDAHelpers::AllocateDeviceMemory(d_GPUContacts);
   tledCUDAHelpers::CopyToDevice(d_GPUContacts, &h_GPUContacts);
}

void tledContactManager::Update(double TR)
{
   // Update all contact cylinder positions
   for (int cylNum = 0; cylNum < numContactCyls; cylNum++)
      Cyls[cylNum].Update(TR);
   
   // Update all contact ultrasound probe positions
   for (int prbNum = 0; prbNum < numContactPrbs; prbNum++)
      Prbs[prbNum].Update(TR);
      
   // Update all contact plate positions
   for (int pltNum = 0; pltNum < numContactPlts; pltNum++)
      Plts[pltNum].Update(TR);
      
   // No need to update kinematic contacts, since they are fixed
}

tledGPUContacts* tledContactManager::GetContactsDevicePointer(void)
{
   return d_GPUContacts;
}

void tledContactManager::SetContactPltDisp(int num, vector<float> disp)
{
   if (num >= numContactPlts)
   {
      cerr << "!!! Warning: requested contact plate number does not exist" << endl;
      return;
   }
   Plts[num].SetDisp(disp);
}

void tledContactManager::SetContactPltStartCrnrA(int num, vector<float> A)
{
   if (num >= numContactPlts)
   {
      cerr << "!!! Warning: requested contact plate number does not exist" << endl;
      return;
   }
   Plts[num].SetStartCrnrA(A);
}

void tledContactManager::SetContactPltStartCrnrB(int num, vector<float> B)
{
   if (num >= numContactPlts)
   {
      cerr << "!!! Warning: requested contact plate number does not exist" << endl;
      return;
   }
   Plts[num].SetStartCrnrB(B);
}

void tledContactManager::SetContactPltStartCrnrC(int num, vector<float> C)
{
   if (num >= numContactPlts)
   {
      cerr << "!!! Warning: requested contact plate number does not exist" << endl;
      return;
   }
   Plts[num].SetStartCrnrC(C);
}

void tledContactManager::SetContactCylDisp(int num, vector<float> disp)
{
   if (num >= numContactCyls)
   {
      cerr << "!!! Warning: requested contact cylinder number does not exist" << endl;
      return;
   }
   Cyls[num].SetOriginDisp(disp);
}

void tledContactManager::SetContactCylRadChange(int num, float dr)
{
   if (num >= numContactCyls)
   {
      cerr << "!!! Warning: requested contact cylinder number does not exist" << endl;
      return;
   }
   Cyls[num].SetRadiusChng(dr);
}

void tledContactManager::SetContactCylOrigin(int num, vector<float> orig)
{
   if (num >= numContactCyls)
   {
      cerr << "!!! Warning: requested contact cylinder number does not exist" << endl;
      return;
   }
   Cyls[num].SetStartOrigin(orig);
}

void tledContactManager::SetContactCylAxis(int num, vector<float> axis)
{
   if (num >= numContactCyls)
   {
      cerr << "!!! Warning: requested contact cylinder number does not exist" << endl;
      return;
   }
   Cyls[num].SetStartAxis(axis);
}

void tledContactManager::SetContactCylRadius(int num, float r)
{
   if (num >= numContactCyls)
   {
      cerr << "!!! Warning: requested contact cylinder number does not exist" << endl;
      return;
   }
   Cyls[num].SetStartRadius(r);
}

void tledContactManager::SetContactCylLength(int num, float l)
{
   if (num >= numContactCyls)
   {
      cerr << "!!! Warning: requested contact cylinder number does not exist" << endl;
      return;
   }
   Cyls[num].SetStartLength(l);
}

void tledContactManager::SetContactPrbDisp(int num, vector<float> disp)
{
   if (num >= numContactPrbs)
   {
      cerr << "!!! Warning: requested contact US probe number does not exist" << endl;
      return;
   }
   Prbs[num].SetOriginDisp(disp);
}

void tledContactManager::SetContactPrbRadChange(int num, float dr)
{
   if (num >= numContactPrbs)
   {
      cerr << "!!! Warning: requested contact US probe number does not exist" << endl;
      return;
   }
   Prbs[num].SetRadiusChng(dr);
}

void tledContactManager::SetContactPrbOrigin(int num, vector<float> orig)
{
   if (num >= numContactPrbs)
   {
      cerr << "!!! Warning: requested contact US probe number does not exist" << endl;
      return;
   }
   Prbs[num].SetStartOrigin(orig);
}

void tledContactManager::SetContactPrbAxis(int num, vector<float> axis)
{
   if (num >= numContactPrbs)
   {
      cerr << "!!! Warning: requested contact US probe number does not exist" << endl;
      return;
   }
   Prbs[num].SetStartAxis(axis);
}

void tledContactManager::SetContactPrbRadius(int num, float r)
{
   if (num >= numContactPrbs)
   {
      cerr << "!!! Warning: requested contact US probe number does not exist" << endl;
      return;
   }
   Prbs[num].SetStartRadius(r);
}

void tledContactManager::SetContactPrbLength(int num, float l)
{
   if (num >= numContactPrbs)
   {
      cerr << "!!! Warning: requested contact US probe number does not exist" << endl;
      return;
   }
   Prbs[num].SetStartLength(l);
}

#endif // _GPU_
