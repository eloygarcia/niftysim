// =========================================================================
// File:       tledSolverGPU_ROM.h
// Purpose:    Main finite element object - GPU version using reduced order modelling
// Package:    NiftySim: Nonlinear finite element program
// Author:     Zeike Taylor
// Language:   C++
// Created:    December 2009
// 
// Copyright (c) 2010, University of Queensland. All rights reserved.
// MedTeQ Centre
// See the LICENSE.txt file in the root folder
// 
// ztaylor@itee.uq.edu.au
// =========================================================================
 
#ifdef _GPU_

#ifndef tledSolverGPU_ROM_H
#define tledSolverGPU_ROM_H

#include "tledSolver.h"
#include "tledCUDAHelpers.h"

/**
 * \brief Reduced-order modelling solver with GPU execution
 * \ingroup solver
 */
class tledSolverGPU_ROM : public tledSolver
{
public:
   tledSolverGPU_ROM();
   virtual ~tledSolverGPU_ROM();
   
   virtual void Init(tledModel* Model);
   virtual void SetFixed(std::vector<int>* IndX, std::vector<int>* IndY, std::vector<int>* IndZ);
   virtual void SetDisps(std::vector<int>* IndX, std::vector<int>* IndY, std::vector<int>* IndZ) {;}
   virtual void SetDisps(std::vector<int>* IndX, std::vector<float>* UX, std::vector<int>* IndY, std::vector<float>* UY,
                         std::vector<int>* IndZ, std::vector<float>* UZ);
   virtual void SetDisps(std::vector<float>* UX, std::vector<float>* UY, std::vector<float>* UZ) {;}
   virtual void SetExtForces(std::vector<int>* IndX, std::vector<int>* IndY, std::vector<int>* IndZ) {;}
   virtual void SetExtForces(std::vector<int>* IndX, std::vector<float>* FX, std::vector<int>* IndY, std::vector<float>* FY,
                             std::vector<int>* IndZ, std::vector<float>* FZ);
   virtual void SetExtForces(std::vector<float>* FX, std::vector<float>* FY, std::vector<float>* FZ) {;}
   virtual float* GetAllExtForces(float *p_dst) const;
   virtual void PerformStep();
   virtual void GetDivergence(bool* Div);
   virtual void UnsetDivergence(void);
   virtual void GetForces(std::vector<int>* NodeInd, std::vector<float>* Forces);
   virtual float* GetAllForces();
   virtual void GetDisps(std::vector<int>* NodeInd, std::vector<float>* Disps);
   virtual float* GetAllDisps();
   virtual float* GetAllNextDisps();
   virtual void GetNodeVMStress(float* NodeSVM);
   virtual void GetNodeVMStrain(float* NodeEVM);
   virtual void GetGaussPtSPKStress(float* S) {;}
   virtual void GetGaussPtGreenStrain(float* E) {;}
   virtual tledMesh* GetMesh() {return Mesh;}
   virtual const tledMesh* GetMesh() const {return Mesh;}
   virtual void GetMassVector(float* mass) const;
   virtual void PrepareOutput() {}
   virtual void PrintNodalForces();
   virtual void PrintDispNodalForces();
   virtual void PrintDispNodalForceSums();
   virtual void PrintNodalDisps();
   virtual void PrintNodalForceSums();
   
   virtual void SetContactManager(tledContactManager* contacts);
   
   virtual void InitialiseSolutionVariables(void);
   virtual void InitialiseConstraints(void);
   virtual void SetTimeStep(double dt);
   
   virtual float GetStrainEnergy(void);
   virtual float GetKineticEnergy(void);
   virtual void SetAllDisps(float* U);
   virtual void ComputeStresses(void); // Note: invalid with viscoelastic materials, since the deformation history is not taken into account
   
   virtual void SetElementMatParams(int el, std::vector<float> params){;}
   virtual void SetMultipleElementMatParams(std::vector<int> el, std::vector<float>* params){;}
   
   virtual void SetGeometry(std::vector<float> NodeCoords) {;}
   
protected:
   // FUNCTIONS
   void ComputeElementT4Variables(tledModel* Model, std::vector<std::vector<int2> >* NInd);
   void ComputeElementT4ANPVariables(tledModel* Model, std::vector<std::vector<int2> >* NInd);
   void ComputeElementH8Variables(tledModel* Model, std::vector<std::vector<int2> >* NInd);
   void ComputeLEparams(float* inMatParams, float4* outMatParams, float2* wkHGLame, float* ANPKappa);
   void ComputeNHparams(float* inMatParams, float4* outMatParams, float2* wkHGLame, float* ANPKappa);
   void ComputeTIparams(float* inMatParams, float4* outMatParams, float2* wkHGLame, float* ANPKappa);
   void ComputeViscparams(float* inMatViscParams, int Ni, int Nv, int2* N, float2* Ai, float2* Av);
   double ComputeHexVol(double x[8][3]);
   void ComputeH8DhDx(double x[8][3], double V, double fDhDx[8][3]);
   void ComputeBmat(double x[8][3], double B[8][3]);
   void ComputeCIJK(double C[8][8][8]);
   void ComputeNodalVolumes();
   void ComputeCDCoeffs();
   
   // HOST VARIABLES
   float4* h_DhDx;      // Shape function derivatives (global)
   int4* h_EInd;	// Element node indices
   float3* h_Ucurr;
   float3* h_Fint;	// Global nodal forces
   int2* h_NodeMap;	// Lookup array for FCds
   int2* h_FCds;	// Lookup array for FEl
   float4* h_HG;	// Hourglass parameters
   int4* h_FixMask;	// Binary mask identifying fixed DOFs
   int4* h_DispMask;
   float4* h_Uload;	// Displaced DOFs (defined as float4 from the outset since it will point to pinned memory)
   float4* h_Fext;		// External nodal forces
   float4* h_Vol_MType_K;	// Element volumes (.x), material types (.y), and bulk modulus (.z)
   float4* h_ElasticParams;	// Element material params
   int2* h_NumProny;		// Number of isochoric and volumetric Prony series terms
   float2* h_bkwdEulerIso;	// Backward Euler integration coefficients - Iso
   float2* h_bkwdEulerVol;	// Backward Euler integration coefficients - Vol
   float3* h_SPKa;         // Second Piola-Kirchhoff stress tensor [S11,S22,S33]
   float3* h_SPKb;         // Second Piola-Kirchhoff stress tensor [S12,S23,S13]
   float3* h_Ca;           // Right Cauchy-Green deformation tensor [C11,C22,C33]
   float3* h_Cb;           // Right Cauchy-Green deformation tensor [C12,C23,C13]
   float* h_Pa;		// Nodal pressures (ANP)
   float* h_Va;		// Nodal volumes (ANP)
   float* h_Mf;		// Mass std::vector (full)
   float* h_Phi;
   float* h_Mr;
   float* h_f;
   float* h_fr;
   float* h_tmp;
   float* h_UnextV;
   float4* l_U;          // Local copy of device displacements for output
   float* UOutput;
   
   // DEVICE VARIABLES
   float4* d_DhDx;
   int4* d_EInd;
   float4* d_HG;
   float4* d_Fint;
   float4* d_FEl;
   float4* d_Uprev;
   float4* d_Ucurr;
   float4* d_Unext;
   float* d_UnextV; // Vector form of Unext
   float4* d_Uload;
   float4* d_Fext;
   int2* d_FCds;
   int2* d_NodeMap;
   int4* d_FixMask;
   int4* d_DispMask;
   float4* d_Vol_MType_K;
   float4* d_ElasticParams;
   int2* d_NumProny;
   float2* d_bkwdEulerIso;
   float2* d_bkwdEulerVol;
   float4* d_StressStateIso;	// State variables for viscoelastic materials
   float4* d_StressStateVol;
   float4* d_SPKa;
   float4* d_SPKb;
   float4* d_Ca;
   float4* d_Cb;
   float* d_Pa;
   float* d_Va;
   float* d_Mfc;  // Compact version of mass vector
   bool* d_Divergence;
   float* d_f;   // Effective load
   float4* d_NodeCds;
   tledGPUContacts* d_Contacts;
   
   // UTILITY VARIABLES
   tledMesh* Mesh;         // Geometric information
   tledContactManager* Contacts;
   const char* EType;	// Element type
   int NumNodes;           // Number of nodes
   int NumEls;		// Number of elements
   int NDOF;		// Number of DOFs per node
   double Dt;		// Time step
   double Ttot;            // Total analysis time
   float alpha;            // Damping coefficient
   float D;		// Mass density
   int FCdsLength;         // Length of FCds array
   int numElasticParamsF4; // Number of float4 variables required to store each el's elastic mat params
   int2 maxNumViscTerms;	// Max number of isochoric and volumetric visco terms for any element in the mesh
   bool ANP;		// Flag to indicate use of ANP formulation
   float2* HGLame;      // Element Lame params, [Lambda,Mu], required for hourglass control
   float* SVM;    // List of element Von Mises effective stresses
   float* EVM;    // List of element Von Mises effective strains
   int numBasisVecs; // Number of vectors in the reduced basis
   float3 CDgamma;   // Central difference integration coeffs
};



#endif // _tledSolverGPU_ROM_H

#endif // _GPU_
