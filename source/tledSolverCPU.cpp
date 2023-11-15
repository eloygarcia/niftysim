// =========================================================================
// File:       tledSolverCPU.cpp
// Purpose:    Main finite element object - CPU version
// Package:    NiftySim: Nonlinear finite element program
// Author:     Zeike Taylor
// Language:   C++
// Created:    March 2007
// 
// Copyright (c) 2010, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// ztaylor@itee.uq.edu.au
// =========================================================================

#include "tledSolverCPU.h"
#include "tledMatrixFunctions.h"
#include "tledShellSolverCPU.h"
#include "tledParallelShellSolverCPU.h"
#include "tledHelper.h"
#include "tledCentralDifferenceTimeStepperCPU.h"
#include "tledBVHTraverserCPU.h"
#include "tledNewmarkTimeStepperCPU.h"

#include <iostream>
#include <fstream>
#include <algorithm>
#include <cmath>

using namespace std;

tledSolverCPU::tledSolverCPU()
{
   Divergence=false;
   NumEls = 0;
   NumNodes = 0;
   mp_TimeStepper = NULL;
}

tledSolverCPU::~tledSolverCPU()
{
  if (mp_TimeStepper != NULL) delete mp_TimeStepper;
  for (int i = 0; i < NumEls; i++)
    {
      delete Elements[i];
    }
  if (NumEls > 0) delete[] Elements;

  if (NumNodes > 0 && NDOF > 0) {
    delete[] F;
    delete[] mp_EffectiveF;
    delete[] R;
    delete[] M;
    delete[] A;
    delete[] B;
    delete[] C;
    delete[] Va;
    delete[] Pa;
  }
}

void tledSolverCPU::Init(tledModel* Model)
{
   // Initialise variables
   Mesh = Model->GetMesh();
   NDOF = Model->GetNodeDOF();
   NumNodes = Model->GetNumNodes();
   NumEls = Model->GetNumEls();
   mp_EffectiveF = new float[NumNodes*NDOF]; 
   F = new float[NumNodes*NDOF]; memset(F,0,sizeof(float)*NumNodes*NDOF);
   R = new float[NumNodes*NDOF]; memset(R,0,sizeof(float)*NumNodes*NDOF);
   M = new float[NumNodes*NDOF];
   A = new float[NumNodes*NDOF]; memset(A,0,sizeof(float)*NumNodes*NDOF);
   B = new float[NumNodes*NDOF]; memset(B,0,sizeof(float)*NumNodes*NDOF);
   C = new float[NumNodes*NDOF]; memset(C,0,sizeof(float)*NumNodes*NDOF);
   Va = new float[NumNodes]; memset(Va,0,sizeof(float)*NumNodes);
   Pa = new float[NumNodes]; memset(Pa,0,sizeof(float)*NumNodes);
   this->InitialiseConstraints();

   // Initialise array of elements
   tledElement* pEl = NULL;
   if (NumEls > 0) Elements = new tledElement*[NumEls];
   int* ElChecker = new int[NumEls];
   memset(ElChecker,0,sizeof(int)*NumEls);
   const char* EType = Model->GetElType();
   if (!strcmp(EType,"T4ANP"))
      ANP = true;
   else
      ANP = false;
   for (int ElSetNum = 0; ElSetNum < Model->GetNumElSets(); ElSetNum++)	// Loop over element sets
   {
      const char* MatType = Model->GetMatType(ElSetNum);
      const float rho = Model->GetElSetMassDensity(ElSetNum);
      vector<int> ElSet = Model->GetElSet(ElSetNum);
      int NumElasticMatParams = Model->GetNumElasticParams(ElSetNum);
      int NumViscMatParams = 2*(Model->GetNumViscIsoTerms(ElSetNum) + Model->GetNumViscVolTerms(ElSetNum));
      int NumMatParams;

      assert(!std::isnan(rho));
      if (NumViscMatParams > 0)
         NumMatParams = NumElasticMatParams + NumViscMatParams + 3;	// 3 extra for dt, Ni, and Nv
      else
         NumMatParams = NumElasticMatParams;
      float* MatParams = new float[NumMatParams];
      Model->GetCombinedMatParams(ElSetNum,MatParams);
      for (int i = 0; i < (int)ElSet.size(); i++)
      {
         if (ElChecker[ElSet[i]] != 0)
            cerr << "\n!!! Element " << ElSet[i] << " listed in more than 1 element set" << endl;
         else
         {
            if (!strcmp(EType,"T4"))
               pEl = new tledElementT4(Mesh, ElSet[i], MatType, MatParams, rho);
            else if (!strcmp(EType,"H8"))
               pEl = new tledElementH8(Mesh, ElSet[i], MatType, MatParams, rho,Model->GetHGKappa());
            else if (!strcmp(EType,"T4ANP"))
            {
               pEl = new tledElementT4ANP(Mesh, ElSet[i], MatType, MatParams, rho);
               float Ve = pEl->GetVolume();
               vector<int> EInd = Mesh->GetElNodeInds(ElSet[i]);
               for (int a = 0; a < 4; a++)
               {
                  Va[EInd[a]] += Ve/4;
               }
            }
            else
               cerr << "\n!!! Invalid element type specified: " << EType << endl;
            Elements[ElSet[i]] = pEl;
            ElChecker[ElSet[i]] = 1;
         }
      }
      delete[] MatParams;
   }
   // Check that all elements instantiated
   vector<int> BadEls;
   bool SomeBadEls = false;
   for (int i = 0; i < NumEls; i++)
   {
      if (ElChecker[i] == 0)
      {
         BadEls.push_back(i);
         SomeBadEls = true;
      }
   }
   if (SomeBadEls == true)
   {
      cerr << "\n!!! Some elements not assigned to an element list:";
      for (int i = 0; i < (int)BadEls.size(); i++)
         cout << " " << BadEls[i];
      cout << endl;

      if ((int)BadEls.size() == NumEls) return;
   }

   delete[] ElChecker;

   if (Model->GetNumberOfShellElementSets() > 0) {
#ifdef _USE_BOOST_
     this->SetShellSolver(new tledParallelShellSolverCPU());
#else
     this->SetShellSolver(new tledShellSolverCPU());
#endif
     this->GetShellSolver().Init(*this, *Model);
   } 

   CompileMass();
   Dt = Model->GetTimeStep();
   m_Alpha = Model->GetDampingCoeff();

   InstantiateTimeStepper();
}

void tledSolverCPU::SetContactManager(tledContactManager* contacts) { 
  tledSolver::SetContactManager(contacts);
}

void tledSolverCPU::SetFixed(vector<int>* IndX, vector<int>* IndY, vector<int>* IndZ)
{
   IndFX = IndX;
   IndFY = IndY;
   IndFZ = IndZ;
}

void tledSolverCPU::SetDisps(vector<int>* IndX, vector<int>* IndY, vector<int>* IndZ)
{
   IndDX = IndX;
   IndDY = IndY;
   IndDZ = IndZ;
}

void tledSolverCPU::SetDisps(vector<int>* IndX, vector<float>* UX, vector<int>* IndY, vector<float>* UY, vector<int>* IndZ, vector<float>* UZ)
{
   // Check sizes of index and displacement vectors
   if ( (IndX->size()!=UX->size()) || (IndY->size()!=UY->size()) || (IndZ->size()!=UZ->size()) )
   {
      tledLogErrorStream(tledHelper::FatalError() << "Input vector sizes are incompatible with index arrays.");
   }
   else
   {
      IndDX = IndX;
      IndDY = IndY;
      IndDZ = IndZ;
      UDX = UX;
      UDY = UY;
      UDZ = UZ;
   }
}

void tledSolverCPU::SetDisps(vector<float>* UX, vector<float>* UY, vector<float>* UZ)
{
   // Check sizes of displacement vectors

  if (IndDX->size() != UX->size() || IndDY->size() != UY->size() || IndDZ->size() != UZ->size()) { 
    tledLogErrorStream(tledHelper::FatalError() << "Input vector sizes are incompatible with index arrays.");
  } else {
    UDX = UX;
    UDY = UY;
    UDZ = UZ;
  }
}

void tledSolverCPU::SetExtForces(vector<int>* IndX, vector<int>* IndY, vector<int>* IndZ)
{
   IndRX = IndX;
   IndRY = IndY;
   IndRZ = IndZ;
}

void tledSolverCPU::SetExtForces(vector<int>* IndX, vector<float>* FX, vector<int>* IndY, vector<float>* FY, vector<int>* IndZ, vector<float>* FZ)
{
   // Check sizes of index and displacement vectors
   if (IndX->size() != FX->size() || IndY->size() != FY->size() || IndZ->size() != FZ->size()) {
     tledLogErrorStream(tledHelper::FatalError() << "Input vector sizes are incompatible with index arrays.");
   } else {
     IndRX = IndX;
     IndRY = IndY;
     IndRZ = IndZ;
     RX = FX;
     RY = FY;
     RZ = FZ;
   }
}

void tledSolverCPU::SetExtForces(vector<float>* FX, vector<float>* FY, vector<float>* FZ)
{
   // Check sizes of displacement vectors
   if (RX->size() != FX->size() || RY->size() != FY->size() || RZ->size() != FZ->size()) {
     tledLogErrorStream(tledHelper::FatalError() << "Input vector sizes are incompatible with constraint settings.");
   } else {
      RX = FX;
      RY = FY;
      RZ = FZ;
   }
}

static void _SetComponentToZero(float *p_disps, const std::vector<int> &inds, const int componentInd) {
  for (std::vector<int>::const_iterator ic_i = inds.begin(); ic_i < inds.end(); ic_i++) {
    p_disps[(*ic_i)*3+componentInd] = 0;
  }
}

void tledSolverCPU::ApplyFixed()
{
   _SetComponentToZero(mp_TimeStepper->GetNextDisplacements(), *IndFX, 0);
   _SetComponentToZero(mp_TimeStepper->GetNextDisplacements(), *IndFY, 1);
   _SetComponentToZero(mp_TimeStepper->GetNextDisplacements(), *IndFZ, 2);
}

void tledSolverCPU::ApplyDisps()
{
   int i;

   _SetComponentToZero(mp_TimeStepper->GetNextDisplacements(), *IndDX, 0);
   _SetComponentToZero(mp_TimeStepper->GetNextDisplacements(), *IndDY, 1);
   _SetComponentToZero(mp_TimeStepper->GetNextDisplacements(), *IndDZ, 2);

   for (i = 0; i < (int)(IndDX->size()); i++)
   {
     mp_TimeStepper->GetNextDisplacements()[(*IndDX)[i]*3] += (*UDX)[i];
   }
   for (i = 0; i < (int)(IndDY->size()); i++)
   {
     mp_TimeStepper->GetNextDisplacements()[(*IndDY)[i]*3+1] += (*UDY)[i];
   }
   for (i = 0; i < (int)(IndDZ->size()); i++)
   {
     mp_TimeStepper->GetNextDisplacements()[(*IndDZ)[i]*3+2] += (*UDZ)[i];
   }
}

void tledSolverCPU::ApplyExtForces()
{
   int i;
   memset(R,0,sizeof(float)*NumNodes*NDOF);
   for (i = 0; i < (int)(IndRX->size()); i++)
   {
      R[(*IndRX)[i]*3] += (*RX)[i];
   }
   for (i = 0; i < (int)(IndRY->size()); i++)
   {
      R[(*IndRY)[i]*3+1] += (*RY)[i];
   }
   for (i = 0; i < (int)(IndRZ->size()); i++)
   {
      R[(*IndRZ)[i]*3+2] += (*RZ)[i];
   }
}

float* tledSolverCPU::GetAllExtForces(float *p_dst) const {
  std::copy(R, R + NDOF*NumNodes, p_dst);

  return p_dst;
}

void tledSolverCPU::PerformStep()
{
   if (ANP == true)
      ComputeNewForcesANP();
   else
      ComputeNewForces();

   if (this->HasMembrane()) static_cast<tledShellSolverCPU&>(this->GetShellSolver()).ComputeNewForces(F, mp_TimeStepper->GetCurrentDisplacements());

   if (this->GetContactManager().DoUnstructuredContacts()) {
     tledUnstructuredContactManager &r_unstructuredMeshContacts = this->GetContactManager().GetUnstructuredContactManager();
     int it = 0;
     bool hadContacts;

     do {
       ComputeNewDisps();
       r_unstructuredMeshContacts.Update();

       hadContacts = false;

       if (r_unstructuredMeshContacts.GetNumberOfRigidSurfaces() > 0) {
	 hadContacts |= r_unstructuredMeshContacts.ComputeDeformableRigidContactResponses(F, mp_TimeStepper->GetNextDisplacements(), mp_TimeStepper->GetCurrentDisplacements());
	 hadContacts |= (it == 0 && tledBVHTraverserCPU::HaveActiveContacts());
       }

       if (r_unstructuredMeshContacts.DoDeformableDeformableContactHandling()) {
	 hadContacts |= r_unstructuredMeshContacts.ComputeDeformableDeformableContactResponses(F, mp_TimeStepper->GetNextDisplacements(), mp_TimeStepper->GetCurrentDisplacements());
	 hadContacts |= (it == 0 && tledBVHTraverserCPU::HaveActiveContacts());
       }
     } while (r_unstructuredMeshContacts.DoMultiPassContactHandling() && hadContacts && (it += 1) < 4);    
   }
   
   ComputeNewDisps();
   ApplyFixed();	// Apply BCs for the next step
   ApplyDisps();	// Apply displacements for the next step
   ApplyExtForces();	// Apply external forces for the next step

   if (this->GetContactManager().DoUnstructuredContacts()) this->GetContactManager().GetUnstructuredContactManager().FinishContactHandling();
   mp_TimeStepper->FinishTimeStep();
}

void tledSolverCPU::UnsetDivergence() {
   Divergence = false;
}

void tledSolverCPU::GetDivergence(bool* Div) {  
  Divergence = false;
  for (float const *pc_u = mp_TimeStepper->GetCurrentDisplacements(); pc_u < mp_TimeStepper->GetCurrentDisplacements() + 3*NumNodes && !Divergence; pc_u++) {
    Divergence |= std::isnan(*pc_u);
  }

  *Div = Divergence;
}

void tledSolverCPU::GetForces(vector<int>* NodeInd, vector<float>* Forces)
{
   // Populate returned forces vector
   for (int i = 0; i < (int)(NodeInd->size()); i++)
   {
      for (int j = 0; j < NDOF; j++)
      {
         Forces->push_back( F[(*NodeInd)[i]*3+j] );
      }
   }
}

void tledSolverCPU::GetDisps(vector<int>* NodeInd, vector<float>* Disps)
{
   for (int i = 0; i < (int)(NodeInd->size()); i++)
   {
      for (int j = 0; j < NDOF; j++)
      {
	Disps->push_back(mp_TimeStepper->GetCurrentDisplacements()[3*(*NodeInd)[i]+j] );
      }
   }
}

void tledSolverCPU::CompileMass(void)
{
   tledElement* pEl = NULL;
   memset(M,0,sizeof(float)*NumNodes*NDOF);
   for (int i = 0; i < NumEls; i++)
   {
      pEl = Elements[i];
      pEl->ComputeElementMass(M);
   }

   if (this->HasMembrane()) {
     float *p_m;
     std::vector<float>::const_iterator ic_mShell;

     for (ic_mShell = this->GetShellSolver().GetMass().begin(), p_m = M; ic_mShell < this->GetShellSolver().GetMass().end(); ic_mShell++) {
       *(p_m++) += *ic_mShell;
       *(p_m++) += *ic_mShell;
       *(p_m++) += *ic_mShell;
     }
   }
}

void tledSolverCPU::ComputeNewForces()
{
   memset(F,0,sizeof(float)*NumNodes*3);
   tledElement* pEl = NULL;
   for (int i = 0; i < NumEls; i++)
   {
      pEl = Elements[i];
      pEl->ComputeElementForces(mp_TimeStepper->GetCurrentDisplacements(), F);
   }
}

void tledSolverCPU::ComputeNewForcesANP()
{
   memset(F,0,sizeof(float)*NumNodes*3);
   memset(Pa,0,sizeof(float)*NumNodes);
   tledElement* pEl = NULL;
   for (int i = 0; i < NumEls; i++)	// 1st element loop
   {
      pEl = Elements[i];
      pEl->ComputeElementPressure(mp_TimeStepper->GetCurrentDisplacements(), Pa);
   }
   for (int i = 0; i < NumEls; i++)	// 2nd element loop
   {
      pEl = Elements[i];
      pEl->ComputeModifiedElementForces(Pa,Va,F);
   }
}

void tledSolverCPU::ComputeNewDisps() {
  float const *pc_r = R, *pc_f = F;
  float *p_ft = mp_EffectiveF;

  while (p_ft < mp_EffectiveF + 3*NumNodes) *(p_ft++) = *(pc_r++) - *(pc_f++);  

  mp_TimeStepper->EvolveDisplacements(mp_EffectiveF);
}

void tledSolverCPU::GetNodeVMStress(float* NodeSVM)
{
   memset(NodeSVM,0,sizeof(float)*NumNodes);
   int* NodeValence = new int[NumNodes];
   memset(NodeValence,0,sizeof(int)*NumNodes);
   vector<int> NodeList;
   float temp;
   float ElSVM;
   for (int i = 0; i < NumEls; i++)
   {
      NodeList = Mesh->GetElNodeInds(i);
      Elements[i]->GetVMStress(&ElSVM);
      for (int j = 0; j < (int)NodeList.size(); j++)
      {
         temp = NodeSVM[NodeList[j]]*NodeValence[NodeList[j]];
         temp += ElSVM;
         NodeValence[NodeList[j]] += 1;
         temp /= NodeValence[NodeList[j]];
         NodeSVM[NodeList[j]] = temp;
      }
   }

   delete NodeValence;
}

void tledSolverCPU::GetNodeVMStrain(float* NodeEVM)
{
   memset(NodeEVM,0,sizeof(float)*NumNodes);
   int* NodeValence = new int[NumNodes];
   memset(NodeValence,0,sizeof(int)*NumNodes);
   vector<int> NodeList;
   float temp;
   float ElEVM;
   for (int i = 0; i < NumEls; i++)
   {
      NodeList = Mesh->GetElNodeInds(i);
      Elements[i]->GetVMStrain(&ElEVM);
      for (int j = 0; j < (int)NodeList.size(); j++)
      {
         temp = NodeEVM[NodeList[j]]*NodeValence[NodeList[j]];
         temp += ElEVM;
         NodeValence[NodeList[j]] += 1;
         temp /= NodeValence[NodeList[j]];
         NodeEVM[NodeList[j]] = temp;
      }
   }

   delete NodeValence;
}

void tledSolverCPU::GetGaussPtSPKStress(float* S)
{
   memset(S,0,sizeof(float)*NumEls*6);
   for (int i = 0; i < NumEls; i++)
      Elements[i]->GetStress(&S[i*6]);
}

void tledSolverCPU::GetGaussPtGreenStrain(float* E)
{
   memset(E,0,sizeof(float)*NumEls*6);
   for (int i = 0; i < NumEls; i++)
      Elements[i]->GetStrain(&E[i*6]);
}

void tledSolverCPU::GetMassVector(float* mass) const 
{
   for (int i = 0; i < NumNodes; i++)
      mass[i] = M[3*i];
}

void tledSolverCPU::PrintNodalForces()
{
   float FX=0,FY=0,FZ=0;
   printf("-----------------------------------------------------\n");
   printf("Nodal Forces (ALL)\n");
   printf("Node\tFX\t\tFY\t\tFZ\n");
   for (int i = 0; i < NumNodes; i++)
   {
      printf("%i\t%+10.5e\t%+10.5e\t%+10.5e\n",i+1,F[i*3],F[i*3+1],F[i*3+2]);
      FX += F[i*3];
      FY += F[i*3+1];
      FZ += F[i*3+2];
   }
   printf("\nSums\t%+10.5e\t%+10.5e\t%+10.5e\n",FX,FY,FZ);
   printf("-----------------------------------------------------\n");
}

void tledSolverCPU::PrintDispNodalForces()
{
   float FX=0,FY=0,FZ=0;
   int* ind = new int[NumNodes]; memset(ind,0,sizeof(int)*NumNodes);
   for (int i = 0; i < (int)(IndDX->size()); i++)
   {
      ind[(*IndDX)[i]] = 1;
   }
   for (int i = 0; i < (int)(IndDY->size()); i++)
   {
      ind[(*IndDY)[i]] = 1;
   }
   for (int i = 0; i < (int)(IndDZ->size()); i++)
   {
      ind[(*IndDZ)[i]] = 1;
   }
   printf("-----------------------------------------------------\n");
   printf("Nodal Forces (LOADED NODES)\n");
   printf("Node\tFX\t\tFY\t\tFZ\n");
   for (int i = 0; i < NumNodes; i++)
   {
      if (ind[i] > 0)
      {
         printf("%i\t%+10.5e\t%+10.5e\t%+10.5e\n",
         i+1,F[3*i],F[3*i+1],F[3*i+2]);
         FX += F[3*i];
         FY += F[3*i+1];
         FZ += F[3*i+2];
      }
   }
   printf("\nSums\t%+10.5e\t%+10.5e\t%+10.5e\n",FX,FY,FZ);
   printf("-----------------------------------------------------\n");

   delete ind;
}

void tledSolverCPU::PrintDispNodalForceSums()
{
   float FX=0,FY=0,FZ=0;
   int* ind = new int[NumNodes]; memset(ind,0,sizeof(int)*NumNodes);
   for (int i = 0; i < (int)(IndDX->size()); i++)
   {
      ind[(*IndDX)[i]] = 1;
   }
   for (int i = 0; i < (int)(IndDY->size()); i++)
   {
      ind[(*IndDY)[i]] = 1;
   }
   for (int i = 0; i < (int)(IndDZ->size()); i++)
   {
      ind[(*IndDZ)[i]] = 1;
   }
   for (int i = 0; i < NumNodes; i++)
   {
      if (ind[i] > 0)
      {
         FX += F[3*i];
         FY += F[3*i+1];
         FZ += F[3*i+2];
      }
   }
   printf("-----------------------------------------------------\n");
   printf("Nodal Force Sums (LOADED NODES)\n");
   printf("\tFX\t\tFY\t\tFZ\n");
   printf("\nSums\t%+10.5e\t%+10.5e\t%+10.5e\n",FX,FY,FZ);
   printf("-----------------------------------------------------\n");

   delete ind;
}

void tledSolverCPU::PrintNodalDisps()
{
   printf("-----------------------------------------------------\n");
   printf("Nodal Displacements (ALL)\n");
   printf("Node\tUX\t\tUY\t\tUZ\n");
   for (int i = 0; i < NumNodes; i++)
   {
     printf("%i\t%+10.5e\t%+10.5e\t%+10.5e\n",i+1, mp_TimeStepper->GetCurrentDisplacements()[i*3], mp_TimeStepper->GetCurrentDisplacements()[i*3+1], mp_TimeStepper->GetCurrentDisplacements()[i*3+2]);
   }
   printf("-----------------------------------------------------\n");
}

void tledSolverCPU::PrintNodalForceSums()
{
   float FX=0,FY=0,FZ=0;
   printf("-----------------------------------------------------\n");
   printf("Nodal Force Sums (ALL)\n");
   printf("\tFX\t\tFY\t\tFZ\n");
   for (int i = 0; i < NumNodes; i++)
   {
      FX += F[i*3];
      FY += F[i*3+1];
      FZ += F[i*3+2];
   }
   printf("\nSums\t%+10.5e\t%+10.5e\t%+10.5e\n",FX,FY,FZ);
   printf("-----------------------------------------------------\n");
}

void tledSolverCPU::InitialiseSolutionVariables(void)
{
   memset(F,0,sizeof(float)*NumNodes*3);
   Divergence = false;
}

void tledSolverCPU::InitialiseConstraints(void)
{
   memset(R,0,sizeof(float)*NumNodes*3);
   IndFX = NULL;
   IndFY = NULL;
   IndFZ = NULL;
   IndDX = NULL;
   IndDY = NULL;
   IndDZ = NULL;
   IndRX = NULL;
   IndRY = NULL;
   IndRZ = NULL;
   UDX = NULL;
   UDY = NULL;
   UDZ = NULL;
   RX = NULL;
   RY = NULL;
   RZ = NULL;
}

void tledSolverCPU::InstantiateTimeStepper() {
  if (mp_TimeStepper != NULL) delete mp_TimeStepper;

  if (ANP) mp_TimeStepper = new tledNewmarkTimeStepperCPU(this->GetNumberOfNodes(), float(this->GetDt()), float(this->GetAlpha()), this->GetMassInternal());
  else mp_TimeStepper = new tledCentralDifferenceTimeStepperCPU(this->GetNumberOfNodes(), float(this->GetDt()), float(this->GetAlpha()), this->GetMassInternal());  
}

void tledSolverCPU::SetTimeStep(double dt)
{
   Dt = dt;
   // Update CD coeffs
   this->InstantiateTimeStepper();
   // Update visc integration params
   for (int i = 0; i < NumEls; i++)
     Elements[i]->GetMaterial()->SetTimeStep((float)Dt);
}

float tledSolverCPU::GetStrainEnergy(void)
{
   float e = 0;
   for (int i = 0; i < NumEls; i++)
      e += Elements[i]->GetStrainEnergy();
   return e;
}

float tledSolverCPU::GetKineticEnergy() {
  std::vector<float> deltaUs(3*NumNodes);
  float e = 0;

  mp_TimeStepper->GetCurrentDeltaDisplacements(&deltaUs.front());
  for (const float *pc_du = &deltaUs.front(), *pc_m = M; pc_m < M + 3*NumNodes; pc_m++, pc_du++) {
    e += *pc_du**pc_du**pc_m;
  }

  return (float)(e/(2*Dt*Dt));
}

void tledSolverCPU::SetAllDisps(float* U) {
  mp_TimeStepper->SetCurrentDisplacements(U);
}

void tledSolverCPU::ComputeStresses(void)
{
   if (ANP)
   {
      memset(Pa,0,sizeof(float)*NumNodes);
      for (int i = 0; i < NumEls; i++)	// 1st element loop
	Elements[i]->ComputeElementPressure(mp_TimeStepper->GetCurrentDisplacements(), Pa);
      for (int i = 0; i < NumEls; i++)	// 2nd element loop
         Elements[i]->ComputeModifiedStress(Pa,Va);
   }
   else
   {
      for (int i = 0; i < NumEls; i++)
	Elements[i]->ComputeStress(mp_TimeStepper->GetCurrentDisplacements());
   }
}

void tledSolverCPU::SetElementMatParams(int el, vector<float> params)
{
   Elements[el]->SetMaterialParams(params);
}

void tledSolverCPU::SetMultipleElementMatParams(vector<int> el, vector<float>* params)
{
   for (int i = 0; i < (int)el.size(); i++)
      Elements[el[i]]->SetMaterialParams(params[i]);
}

void tledSolverCPU::SetGeometry(vector<float> NodeCoords)
{
   Mesh->SetNodeCoordinates(NodeCoords);
   memset(Va,0,sizeof(float)*NumNodes);
   for (int el = 0; el < NumEls; el++)
   {
      Elements[el]->SetGeometry(Mesh);
      if (ANP)
      {
         float Ve = Elements[el]->GetVolume();
         vector<int> EInd = Mesh->GetElNodeInds(el);
         for (int a = 0; a < 4; a++)
            Va[EInd[a]] += Ve/4;
      }
   }
   CompileMass();

   InstantiateTimeStepper();
}
