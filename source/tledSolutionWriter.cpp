// =========================================================================
// File:       tledSolutionWriter.cpp
// Purpose:    Compile model solutions for saving
// Package:    NiftySim: Nonlinear finite element program
// Author:     Zeike Taylor
// Language:   C++
// Created:    Dec 2009
// 
// Copyright (c) 2010, University of Queensland. All rights reserved.
// MedTeQ Centre
// See the LICENSE.txt file in the root folder
// 
// ztaylor@itee.uq.edu.au
// =========================================================================

#include "tledSolutionWriter.h"

#include <fstream>
#include <iterator>
#include <algorithm>

tledSolutionWriter::tledSolutionWriter(void)
{
}

tledSolutionWriter::tledSolutionWriter(tledSolver* Solver, tledModel* Model, int numSteps)
{
   solver = Solver;
   model = Model;
   m_Freq = Model->GetOutputFreq();
   if (m_Freq > 0)
      nSvs = numSteps/m_Freq;
   else
      nSvs = 0;
   Scnt = 0;
   U = NULL;
   F = NULL;
   S = NULL;
   E = NULL;
   // Check output variable requests
   saveDisplacements = false;
   saveForces = false;
   saveStresses = false;
   saveStrains = false;
   saveEKin = false;
   saveEStrain = false;
   int numOutputVars = Model->GetNumOutputVars();
   for (int i = 0; i < numOutputVars; i++)
   {
      const char* var = Model->GetOutputVar(i);
      if (!strcmp(var,"U"))
         saveDisplacements = true;
      else if (!strcmp(var,"F"))
         saveForces = true;
      else if (!strcmp(var,"S"))
         saveStresses = true;
      else if (!strcmp(var,"E"))
         saveStrains = true;
      else if (!strcmp(var,"EKin"))
         saveEKin = true;
      else if (!strcmp(var,"EStrain"))
         saveEStrain = true;
      else
	tledLogErrorStream(tledHelper::FatalError() << "Invalid output request: " << var);
   }
   // Define output variable arrays
   int numNodes = solver->GetMesh()->GetNumNodes();
   int numEls = solver->GetMesh()->GetNumEls();
   if (saveDisplacements)
   {
      U = new std::vector<float>[nSvs];
      for (int i = 0; i < nSvs; i++)
         U[i].resize(numNodes*3,0);
   }
   if (saveForces)
   {
      F = new std::vector<float>[nSvs];
      for (int i = 0; i < nSvs; i++)
         F[i].resize(numNodes*3,0);
   }
   if (saveStresses)
   {
      S = new std::vector<float>[nSvs];
      for (int i = 0; i < nSvs; i++)
         S[i].resize(numEls*6,0);
   }
   if (saveStrains)
   {
      E = new std::vector<float>[nSvs];
      for (int i = 0; i < nSvs; i++)
         E[i].resize(numEls*6,0);
   }
   this->SetFilePrefix(Model->GetDirectory());
}

tledSolutionWriter::~tledSolutionWriter(void)
{
   if (saveDisplacements)
      delete[] U;
   if (saveForces)
      delete[] F;
   if (saveStresses)
      delete[] S;
   if (saveStrains)
      delete[] E;
}  

void tledSolutionWriter::SetFilePrefix(const std::string &pfx) {
  m_FilePrefix = pfx;
  m_FNameU = pfx + "U.txt";
  m_FNameF = pfx + "F.txt";
  m_FNameS = pfx + "S.txt";
  m_FNameE = pfx + "E.txt";
}

void tledSolutionWriter::SaveSolution(int step)
{
   if ( (nSvs > 0) && (step%m_Freq == 0) && (step > 0) )
   {
     solver->PrepareOutput();
      if (saveDisplacements)
      {
         float* u = solver->GetAllDisps();
         for (int i = 0; i < solver->GetMesh()->GetNumNodes()*3; i++)
            (U[Scnt])[i] = u[i];
      }
      if (saveForces)
      {
         float* f = solver->GetAllForces();
         for (int i = 0; i < solver->GetMesh()->GetNumNodes()*3; i++)
            (F[Scnt])[i] = f[i];
      }
      if (saveStresses)
      {
         float* s = new float[solver->GetMesh()->GetNumEls()*6];
         solver->GetGaussPtSPKStress(s);
         for (int i = 0; i < solver->GetMesh()->GetNumEls()*6; i++)
            (S[Scnt])[i] = s[i];
	 delete[] s;
      }
      if (saveStrains)
      {
         float* e = new float[solver->GetMesh()->GetNumEls()*6];
         solver->GetGaussPtGreenStrain(e);
         for (int i = 0; i < solver->GetMesh()->GetNumEls()*6; i++)
            (E[Scnt])[i] = e[i];
	 delete[] e;
      }

      if (saveEStrain) {
	m_EStrain.push_back(solver->GetStrainEnergy());
      }

      if (saveEKin) {
	m_EKin.push_back(solver->GetKineticEnergy());
      }

      Scnt++;
   }
}

static void _SaveEnergy(const std::string &fBaseName, const std::string &dirPath, const std::vector<float> &energy) {
  const std::string fname = dirPath + fBaseName;
     
  std::ofstream fout(fname.c_str());
     
  if (fout.is_open()) {
    std::copy(energy.begin(), energy.end(), std::ostream_iterator<float>(fout, "\n"));
  } else {
    tledLogErrorStream(tledHelper::Warning() << "Could not open " << fname << " for writing.");
  }
}

void tledSolutionWriter::WriteSolution(void)
{
   if (saveDisplacements)
   {
      std::ofstream fileU(m_FNameU.c_str());
      if (!fileU)
      {
         tledLogErrorStream(tledHelper::Warning() << "Cannot open file " << m_FNameU << " --> results not saved");
      }
      for (int i = 0; i < solver->GetMesh()->GetNumNodes()*3; i++)
      {
         for (int j = 0; j < nSvs; j++)
         {
            fileU << (U[j])[i];
            fileU << " ";
         }
         fileU << "\n";
      }
      fileU.close();
   }
   
   if (saveForces)
   {
      std::ofstream fileF(m_FNameF.c_str());
      if (!fileF)
      {
         tledLogErrorStream(tledHelper::Warning() << "Cannot open file " << m_FNameF << " --> results not saved");
      }
      for (int i = 0; i < solver->GetMesh()->GetNumNodes()*3; i++)
      {
         for (int j = 0; j < nSvs; j++)
         {
            fileF << (F[j])[i];
            fileF << " ";
         }
         fileF << "\n";
      }
      fileF.close();
   }
   if (saveStresses)
   {
      std::ofstream fileS(m_FNameS.c_str());
      if (!fileS)
      {
         tledLogErrorStream(tledHelper::Warning() << "Cannot open file " << m_FNameS << " --> results not saved");
      }
      for (int i = 0; i < solver->GetMesh()->GetNumEls()*6; i++)
      {
         for (int j = 0; j < nSvs; j++)
         {
            fileS << (S[j])[i];
            fileS << " ";
         }
         fileS << "\n";
      }
      fileS.close();
   }
   if (saveStrains)
   {
      std::ofstream fileE(m_FNameE.c_str());
      if (!fileE)
      {
         tledLogErrorStream(tledHelper::Warning() << "Cannot open file " << m_FNameE << " --> results not saved");
      }
      for (int i = 0; i < solver->GetMesh()->GetNumEls()*6; i++)
      {
         for (int j = 0; j < nSvs; j++)
         {
            fileE << (E[j])[i];
            fileE << " ";
         }
         fileE << "\n";
      }
      fileE.close();
   }

   if (saveEKin) {
     _SaveEnergy("EKin.txt", this->GetFilePrefix(), m_EKin);
   }

   if (saveEStrain) {
     _SaveEnergy("EStrain.txt", this->GetFilePrefix(), m_EStrain);
   }
}

void tledSolutionWriter::InitialiseHistories(void)
{
   Scnt = 0;
   int numNodes = solver->GetMesh()->GetNumNodes();
   int numEls = solver->GetMesh()->GetNumEls();
   if (saveDisplacements)
   {
      for (int i = 0; i < nSvs; i++)
         U[i].resize(numNodes*3,0);
   }
   if (saveForces)
   {
      for (int i = 0; i < nSvs; i++)
         F[i].resize(numNodes*3,0);
   }
   if (saveStresses)
   {
      for (int i = 0; i < nSvs; i++)
         S[i].resize(numEls*6,0);
   }
   if (saveStrains)
   {
      for (int i = 0; i < nSvs; i++)
         E[i].resize(numEls*6,0);
   }

   if (saveEKin) {
     m_EKin.clear();
     m_EKin.reserve(nSvs);
   }

   if (saveEStrain) {
     m_EStrain.clear();
     m_EStrain.reserve(nSvs);
   }
}

void tledSolutionWriter::SetNumSteps(const int numSteps)
{
   if (m_Freq > 0)
      nSvs = numSteps/m_Freq;
   else
      nSvs = 0;
   Scnt = 0;
   int numNodes = solver->GetMesh()->GetNumNodes();
   int numEls = solver->GetMesh()->GetNumEls();
   if (saveDisplacements)
   {
      if (U)
         delete[] U;
      U = new std::vector<float>[nSvs];
      for (int i = 0; i < nSvs; i++)
         U[i].resize(numNodes*3,0);
   }
   if (saveForces)
   {
      if (F)
         delete[] F;
      F = new std::vector<float>[nSvs];
      for (int i = 0; i < nSvs; i++)
         F[i].resize(numNodes*3,0);
   }
   if (saveStresses)
   {
      if (S)
         delete[] S;
      S = new std::vector<float>[nSvs];
      for (int i = 0; i < nSvs; i++)
         S[i].resize(numEls*6,0);
   }
   if (saveStrains)
   {
      if (E)
         delete[] E;
      E = new std::vector<float>[nSvs];
      for (int i = 0; i < nSvs; i++)
         E[i].resize(numEls*6,0);
   }

   if (saveEKin) {
     m_EKin.clear();
     m_EKin.reserve(nSvs);
   }

   if (saveEStrain) {
     m_EStrain.clear();
     m_EStrain.reserve(nSvs);
   }
}
