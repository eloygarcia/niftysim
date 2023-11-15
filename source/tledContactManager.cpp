// =========================================================================
// File:       tledContactManager.cpp
// Purpose:    Manage all contact objects in the model (CPU version)
// Package:    NiftySim: Nonlinear finite element program
// Author:     Zeike Taylor
// Language:   C++
// Created:    April 2011
// 
// Copyright (c) 2011, University of Queensland. All rights reserved.
// MedTeQ Centre
// See the LICENSE.txt file in the root folder
// 
// ztaylor@itee.uq.edu.au
// =========================================================================

#ifndef _GPU_

#include "tledSolver.h"
#include "tledDeformableDeformableContactSolver.h"
#include "tledDeformableRigidContactSolver.h"
#include "tledContactManager.h"

using namespace std;

tledContactManager::~tledContactManager() {
  if (DoUnstructuredContacts()) delete mp_UnstructuredContacts;
}

tledContactManager::tledContactManager(tledModel *p_model) : mp_UnstructuredContacts(NULL), numContactCyls(0), numContactPrbs(0), numContactPlts(0) {
  Model = p_model;
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

void tledContactManager::Update(double TR)
{
   // CPU version not implemented
}

void tledContactManager::SetContactPltDisp(int num, vector<float> disp)
{
   // CPU version not implemented
}

void tledContactManager::SetContactPltStartCrnrA(int num, vector<float> A)
{
   // CPU version not implemented
}

void tledContactManager::SetContactPltStartCrnrB(int num, vector<float> B)
{
   // CPU version not implemented
}

void tledContactManager::SetContactPltStartCrnrC(int num, vector<float> C)
{
   // CPU version not implemented
}

void tledContactManager::SetContactCylDisp(int num, vector<float> disp)
{
   // CPU version not implemented
}

void tledContactManager::SetContactCylRadChange(int num, float dr)
{
   // CPU version not implemented
}

void tledContactManager::SetContactCylOrigin(int num, vector<float> orig)
{
   // CPU version not implemented
}

void tledContactManager::SetContactCylAxis(int num, vector<float> axis)
{
   // CPU version not implemented
}

void tledContactManager::SetContactCylRadius(int num, float r)
{
   // CPU version not implemented
}

void tledContactManager::SetContactCylLength(int num, float l)
{
   // CPU version not implemented
}

void tledContactManager::SetContactPrbDisp(int num, vector<float> disp)
{
   // CPU version not implemented
}

void tledContactManager::SetContactPrbRadChange(int num, float dr)
{
   // CPU version not implemented
}

void tledContactManager::SetContactPrbOrigin(int num, vector<float> orig)
{
   // CPU version not implemented
}

void tledContactManager::SetContactPrbAxis(int num, vector<float> axis)
{
   // CPU version not implemented
}

void tledContactManager::SetContactPrbRadius(int num, float r)
{
   // CPU version not implemented
}

void tledContactManager::SetContactPrbLength(int num, float l)
{
   // CPU version not implemented
}

#endif // _GPU_
