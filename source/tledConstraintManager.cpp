// =========================================================================
// File:       tledConstraintManager.h
// Purpose:    Class for managing all constraints in an analysis
// Package:    NiftySim: Nonlinear finite element program
// Author:     Zeike Taylor
// Language:   C++
// Created:    June 2008
// 
// Copyright (c) 2010, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// ztaylor@itee.uq.edu.au
// =========================================================================


#include "tledConstraintManager.h"
#include "tledMatrixFunctions.h"

#include <typeinfo>
#include <string>
#include <cstring>
#include <cassert>

using namespace std;

tledConstraintManager::tledConstraintManager()
{
   Constraints.resize(0);
   FixConstraints.resize(0);
}

tledConstraintManager::tledConstraintManager(tledModel* p_model, tledSolver* p_solver) {
  this->Init(*p_solver, *p_model);
}

template <class TConstraint>
static void _InstantiateConstraint(std::vector<tledConstraint*> &r_constraintVec, const std::vector<int> &inds, const std::vector<float> &mags, const loadShape ls, const int DOF) {
  if (DOF == tledModel::CONSTRAINT_ALL_DOFS) {
    for (int c = 0; c < 3; c++) r_constraintVec.push_back(new TConstraint(inds, mags, ls, c));
  } else r_constraintVec.push_back(new TConstraint(inds, mags, ls, DOF));
}

void tledConstraintManager::Init(tledSolver &r_solver, const tledModel &model) {
   // Loop through constraints in xml file
   for (int i = 0; i < model.GetNumConstraints(); i++)
   {
      const char* currConstraintType = model.GetConstraintType(i);
      int DOF;
      enum loadShape ls;
      vector<int> Ind;
      vector<float> Mag;

      if (!strcmp(currConstraintType,"Fix")) {
         DOF = model.GetConstraintDOF(i);
         Ind = model.GetConstraintInd(i);
         // Instantiate new constraint
	 if (DOF == tledModel::CONSTRAINT_ALL_DOFS) {
	   for (int c = 0; c < 3; c++) FixConstraints.push_back(new tledFixConstraint(Ind, c));
	 } else FixConstraints.push_back(new tledFixConstraint(Ind, DOF));
      } else if (!strcmp(currConstraintType,"Disp")) {
         DOF = model.GetConstraintDOF(i);
         ls = model.GetConstraintLoadShape(i);
         Ind = model.GetConstraintInd(i);
         Mag = model.GetConstraintMag(i);
	 _InstantiateConstraint<tledDispConstraint>(Constraints, Ind, Mag, ls, DOF);
      } else if (!strcmp(currConstraintType,"Force")) {
         DOF = model.GetConstraintDOF(i);
         ls = model.GetConstraintLoadShape(i);
         Ind = model.GetConstraintInd(i);
         Mag = model.GetConstraintMag(i);
	 _InstantiateConstraint<tledForceConstraint>(Constraints, Ind, Mag, ls, DOF);
      } else if (!strcmp(currConstraintType,"Gravity")) {
         // NB: THIS IS NOT REALLY A GOOD DESCRIPTION OF GRAVITATIONAL LOADING
         // - IT MAKES USE OF THE FACT THAT MASS IS LUMPED AT THE NODES AND APPLIES A FORCE
         // PROPORTIONAL TO THIS. BETTER TO INTEGRATE A BODY FORCE OVER THE ELEMENT VOLUME
         // --> USE WITH CAUTION!
         ls = model.GetConstraintLoadShape(i);
         Ind = model.GetConstraintInd(i);
         float accelMag = model.GetGravityMag(i);
         vector<float> accelDir = model.GetGravityDirec(i);
         float normAccelDir = Norm(accelDir);
         if (normAccelDir != 1)
         {
            // Normalise direction vector
            for (int j = 0; j < 3; j++)
               accelDir[j] /= normAccelDir;
         }
         float* allMass = new float[model.GetNumNodes()];
         r_solver.GetMassVector(allMass);
         vector<float> Mass(Ind.size());
         for (int i = 0; i < (int)Ind.size(); i++)
            Mass[i] = allMass[Ind[i]];
         tledGravityConstraint* newConstraint = new tledGravityConstraint(Ind,Mass,accelMag,accelDir,ls);
         Constraints.push_back(newConstraint);
         
         delete[] allMass;
      } else if (!strcmp(currConstraintType,"Pressure")) {
         ls = model.GetConstraintLoadShape(i);
         const char* fType = model.GetPressureFaceType(i);
         int FaceType;
         if (!strcmp(fType,"Quad"))
            FaceType = 0;
         else if (!strcmp(fType,"Tri"))
            FaceType = 1;
         else
           tledLogErrorStream(tledHelper::FatalError() << "Invalid face type: " << fType << " in XML file.");

         vector<int> Faces = model.GetPressureFaceNodeInds(i);
         float Magnitude = model.GetPressureMagnitude(i);
         tledPressureConstraint* newConstraint = new tledPressureConstraint(Faces, FaceType, &r_solver, Magnitude, ls);
         Constraints.push_back(newConstraint);
      } else if (!strcmp(currConstraintType, "Traction")) {
	const std::string fTypeStr = model.GetTractionFaceType(i);
	
	int faceType;
	std::vector<int> faces;
	std::vector<float> tractions;

	if (fTypeStr == "Quad") faceType = 0;
	else if (fTypeStr == "Tri") faceType = 1;
	else {
	  tledLogErrorStream(tledHelper::FatalError() << "Invalid face type: " << fTypeStr << " in XML file.");
	}

	faces = model.GetTractionFaceNodeInds(i);
	tractions = model.GetTractionFaceTractions(i);
	ls = model.GetConstraintLoadShape(i);

	Constraints.push_back(new tledTractionConstraint(r_solver, faceType, faces, tractions, ls));
      } else {
	Constraints.push_back(this->CreateConstraint(currConstraintType, i, model));
      }
   }
}

tledConstraint* tledConstraintManager::CreateConstraint(const std::string &cType, const int cInd, const tledModel &model) {
  tledLogErrorStream(tledHelper::FatalError() << "Invalid constraint type: " << cType << " in XML file.");

  return NULL;
}

tledConstraintManager::~tledConstraintManager() {
  for (std::vector<tledConstraint*>::iterator ip_c = Constraints.begin(); ip_c < Constraints.end(); ip_c++) {
    delete *ip_c;
  }
  
  for (std::vector<tledFixConstraint*>::iterator ip_c = FixConstraints.begin(); ip_c < FixConstraints.end(); ip_c++) {
    delete *ip_c;
  }
}

vector<int>* tledConstraintManager::GetDispInd(int dof) {
   vector<int>* workInds;
   int size = 0;

   // Loop over constraints and compile relevant indices
   IndDisp[dof].clear();
   for (int i = 0; i < (int)Constraints.size(); i++) size += Constraints[i]->GetDispInd(dof)->size();
   IndDisp[dof].reserve(size);
   for (int i = 0; i < (int)Constraints.size(); i++) {
      workInds = Constraints[i]->GetDispInd(dof);
      IndDisp[dof].insert(IndDisp[dof].end(), workInds->begin(), workInds->end());
   }
   
   return &IndDisp[dof];
}

vector<float>* tledConstraintManager::GetDispVal(int dof, int step, double dt, double T) {
   vector<float>* workVals;

   // Loop over dispConstraints and compile relevant displacements
   ValDisp[dof].clear();
   ValDisp[dof].reserve(IndDisp[dof].size());
   for (int i = 0; i < (int)Constraints.size(); i++) {
      workVals = Constraints[i]->GetDispVal(dof,step,dt,T);
      ValDisp[dof].insert(ValDisp[dof].end(), workVals->begin(), workVals->end());
   }
   assert(ValDisp[dof].size() == IndDisp[dof].size());
   
   return &ValDisp[dof];
}

vector<int>* tledConstraintManager::GetForceInd(int dof) {
  vector<int>* workInds;
  int size = 0;

  // Loop over constraints and compile relevant indices   
  IndForce[dof].clear();
  for (int i = 0; i < (int)Constraints.size(); i++) size += Constraints[i]->GetForceInd(dof)->size();
  IndForce[dof].reserve(size);
  for (int i = 0; i < (int)Constraints.size(); i++) {
    workInds = Constraints[i]->GetForceInd(dof);
    IndForce[dof].insert(IndForce[dof].end(), workInds->begin(), workInds->end());
  }
   
  return &IndForce[dof];
}

vector<float>* tledConstraintManager::GetForceVal(int dof, int step, double dt, double T) {
  vector<float>* workVals;

  // Loop over dispConstraints and compile relevant displacements
  ValForce[dof].clear();
  ValForce[dof].reserve(IndForce[dof].size());
  for (int i = 0; i < (int)Constraints.size(); i++) {
    workVals = Constraints[i]->GetForceVal(dof,step,dt,T);
    ValForce[dof].insert(ValForce[dof].end(), workVals->begin(), workVals->end());
  }
  assert(ValForce[dof].size() == IndForce[dof].size());
  
  return &ValForce[dof];
}

vector<int>* tledConstraintManager::GetFixInd(int dof) {
   vector<int>* workInds;
   int size = 0;

   // Loop over constraints and compile relevant indices
   IndFix[dof].clear();
   for (int i = 0; i < (int)FixConstraints.size(); i++) size += FixConstraints[i]->GetFixInd(dof)->size();
   IndFix[dof].reserve(size);
   for (int i = 0; i < (int)FixConstraints.size(); i++) {
      workInds = FixConstraints[i]->GetFixInd(dof);
      IndFix[dof].insert(IndFix[dof].end(), workInds->begin(), workInds->end());
   }
   
   return &IndFix[dof];
}

tledConstraint* tledConstraintManager::GetConstraint(int num)
{
  if (num >= GetNumConstraints())
    {
      cerr << "!!! Warning: requested constraint number does not exist" << endl;
      return NULL;
   }
   
  if (num < (int)Constraints.size()) {
   // NB: it is up to the calling function to dynamically cast the returned pointer,
   // since tledConstraint is abstract
   return Constraints[num];
  } else return FixConstraints[num-Constraints.size()];
}

const tledConstraint* tledConstraintManager::GetConstraint(int num) const {
  if (num >= GetNumConstraints())
    {
      cerr << "!!! Warning: requested constraint number does not exist" << endl;
      return NULL;
   }
   
  if (num < (int)Constraints.size()) {
   // NB: it is up to the calling function to dynamically cast the returned pointer,
   // since tledConstraint is abstract
   return Constraints[num];
  } else return FixConstraints[num-Constraints.size()];
}  

const tledFixConstraint* tledConstraintManager::GetFixConstraint(int num) const
{
  if (num >= (int)FixConstraints.size())
   {
      cerr << "!!! Warning: requested fix constraint number does not exist" << endl;
      return NULL;
   }
   
   return FixConstraints[num];
}

tledFixConstraint* tledConstraintManager::GetFixConstraint(int num) 
{
  if (num >= (int)FixConstraints.size())
   {
      cerr << "!!! Warning: requested fix constraint number does not exist" << endl;
      return NULL;
   }
   
   return FixConstraints[num];
}

int tledConstraintManager::GetNumDispConstraints() const
{
   tledDispConstraint dc;
   int sum = 0;
   for (int i = 0; i < (int)Constraints.size(); i++)
   {
      if ( typeid(*(Constraints[i])).name() == typeid(dc).name() )
         sum++;
   }
   return sum;
}

int tledConstraintManager::GetNumForceConstraints() const
{
   tledForceConstraint fc;
   int sum = 0;
   for (int i = 0; i < (int)Constraints.size(); i++)
   {
      if ( typeid(*(Constraints[i])).name() == typeid(fc).name() )
         sum++;
   }
   return sum;
}

void tledConstraintManager::AddConstraint(tledConstraint* constraint)
{
  if (typeid(*constraint) == typeid(tledFixConstraint)) AddFixConstraint(dynamic_cast<tledFixConstraint*>(constraint));
  else Constraints.push_back(constraint);
}

void tledConstraintManager::AddFixConstraint(tledFixConstraint* constraint)
{
   FixConstraints.push_back(constraint);
}
