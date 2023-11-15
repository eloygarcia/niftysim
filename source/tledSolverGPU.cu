// =========================================================================
// File:       tledSolverGPU.cu
// Purpose:    Main finite element object - GPU version
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


#ifdef _GPU_

#include "tledSolverGPU.h"
#include "tledMatrixFunctions.h"
#include "tledTimer.h"
#include "tledDeviceDeclarations.h"
#include "tledShellSolverGPU.h"
#include "tledCUDAMemoryBlock.h"

#include <iostream>
#include <math.h>

#include "tledTimeStepperGPU.cu"
#include "tledShellSolverGPU.cu"
#include "tledSolverGPU_kernels.cu"

using namespace std;
using namespace tledSolverGPU_kernels;

tledSolverGPU::tledSolverGPU()
{
  h_HG = NULL;
  mp_TimeStepper = NULL;
}

tledSolverGPU::~tledSolverGPU()
{
   if (NumNodes == 0) return;

   // Host variables
   delete mp_TimeStepper;
   delete[] h_F;
   delete[] M;
   delete[] h_CD;
   delete[] h_NodeMap;
   delete[] h_FCds;

   if (NumEls > 0) {
     delete[] h_DhDx;
     delete[] h_EInd;
     delete[] h_Vol_MType_K;
     delete[] h_ElasticParams;
     delete[] h_NumProny;
     delete[] h_bkwdEulerIso;
     delete[] h_bkwdEulerVol;
     delete[] h_SPKa;
     delete[] h_SPKb;
     delete[] h_Ca;
     delete[] h_Cb;
     delete[] HGLame;
     delete[] ElementVols;
     if(h_HG) delete[] h_HG;
   }

   delete[] h_Pa;
   delete[] h_Va;

   // Pinned host variables
   tledCheckCUDAErrors(cudaFreeHost(h_Uload));
   tledCheckCUDAErrors(cudaFreeHost(h_R));
   
   // Device variables
   if (NumEls > 0) {
     tledCheckCUDAErrors(cudaFree(d_DhDx));
     tledCheckCUDAErrors(cudaFree(d_EInd));
     tledCheckCUDAErrors(cudaFree(d_FEl));
     tledCheckCUDAErrors(cudaFree(d_FCds));
     tledCheckCUDAErrors(cudaFree(d_ElasticParams));
     tledCheckCUDAErrors(cudaFree(d_Vol_MType_K));
     if (d_HG != NULL) {
       tledCheckCUDAErrors(cudaFree(d_HG));
     }
     if (d_NumProny != NULL) {
       tledCheckCUDAErrors(cudaFree(d_NumProny));
     }
     if (d_bkwdEulerIso != NULL) {
       tledCheckCUDAErrors(cudaFree(d_bkwdEulerIso));
     }
     if (d_bkwdEulerVol != NULL) {
       tledCheckCUDAErrors(cudaFree(d_bkwdEulerVol));
     }
     if (d_StressStateIso != NULL) {
       tledCheckCUDAErrors(cudaFree(d_StressStateIso));
     }
     if (d_StressStateVol != NULL) {
       tledCheckCUDAErrors(cudaFree(d_StressStateVol));
     }
     if (d_Pa != NULL) {
       tledCheckCUDAErrors(cudaFree(d_Pa));
     }
     if (d_Va != NULL) {
       tledCheckCUDAErrors(cudaFree(d_Va));
     }

     tledCheckCUDAErrors(cudaFree(d_SPKa));
     tledCheckCUDAErrors(cudaFree(d_SPKb));
     tledCheckCUDAErrors(cudaFree(d_Ca));
     tledCheckCUDAErrors(cudaFree(d_Cb));
   }
   tledCheckCUDAErrors(cudaFree(d_CD));
   tledCheckCUDAErrors(cudaFree(d_F));
   tledCheckCUDAErrors(cudaFree(d_FixMask));
   tledCheckCUDAErrors(cudaFree(d_NodeMap));
   tledCheckCUDAErrors(cudaFree(d_Uload));
   tledCheckCUDAErrors(cudaFree(d_DispMask));
   tledCheckCUDAErrors(cudaFree(d_NodeCds));
   tledCheckCUDAErrors(cudaFree(d_R));
   tledCheckCUDAErrors(cudaFree(d_Divergence));
}

const float4* tledSolverGPU::GetAllOnDeviceNextDisplacements(void) const {
  return const_cast<const tledTimeStepperGPU*>(mp_TimeStepper)->GetOnDeviceNextDisplacements();
}

const float4* tledSolverGPU::GetAllOnDeviceCurrentDisplacements(void) const {
  return const_cast<const tledTimeStepperGPU*>(mp_TimeStepper)->GetOnDeviceCurrentDisplacements();
}

void tledSolverGPU::Init(tledModel* Model)
{
  //   printf("\n");
  //   cudaDeviceProp deviceProp;
  //   cudaGetDeviceProperties(&deviceProp,0);
  //   printf("totalGlobalMem = %u\n",deviceProp.totalGlobalMem);
  //   printf("sharedMemPerBlock = %i\n",deviceProp.sharedMemPerBlock);
  //   printf("regsPerBlock = %i\n",deviceProp.regsPerBlock);
  //   printf("warpSize = %i\n",deviceProp.warpSize);
  //   printf("memPitch = %i\n",deviceProp.memPitch);
  //   printf("maxThreadsPerBlock = %i\n",deviceProp.maxThreadsPerBlock);
  //   printf("totalConstMem = %i\n",deviceProp.totalConstMem);
  //   printf("major = %i\n",deviceProp.major);
  //   printf("minor = %i\n",deviceProp.minor);
  //   printf("clockRate = %i\n",deviceProp.clockRate);
  //   printf("deviceOverlap = %i\n",deviceProp.deviceOverlap);
  //   printf("multiProcessorCount = %i\n",deviceProp.multiProcessorCount);
   
  // Initialise variables
  Mesh = Model->GetMesh();
  Dt = Model->GetTimeStep();
  Ttot = Model->GetTotalTime();
  alpha = Model->GetDampingCoeff();
  NDOF = Model->GetNodeDOF();
  EType = Model->GetElType();
  NumNodes = Model->GetNumNodes();
  NumEls = Model->GetNumEls();
  int NPE = Model->GetNodesPerEl();
  ANP = false;
  numElasticParamsF4 = (int)ceil((double)(Model->GetMaxNumElasticParams())/4);
  maxNumViscTerms.x = Model->GetMaxNumViscIsoTerms();
  maxNumViscTerms.y = Model->GetMaxNumViscVolTerms();

  if (NumNodes == 0) return;
   
  // Allocate host variables
  h_F = new float3[NumNodes]; memset(h_F,0,sizeof(float)*NumNodes*3);
  M = new float[NumNodes]; memset(M,0,sizeof(float)*NumNodes);
  h_CD = new float4[NumNodes]; memset(h_CD,0,sizeof(float)*NumNodes*4);
  h_NodeMap = new int2[NumNodes]; memset(h_NodeMap,0,sizeof(int)*NumNodes*2);
   
  m_ElementElSetIndex = std::vector<int>(NumEls, -1);
  m_ElSetDensities.clear();
  m_ElSetDensities.reserve(Model->GetNumElSets());
  for (int esInd = 0; esInd < Model->GetNumElSets(); esInd++) m_ElSetDensities.push_back(Model->GetElSetMassDensity(esInd));

  if (NumEls > 0) {
    h_EInd = new int4[NumEls*NPE/4]; memset(h_EInd,0,sizeof(int)*NumEls*NPE); // NPE/4 = 1 for T4, and 2 for H8
    h_DhDx = new float4[NumEls*NPE*3/4]; memset(h_DhDx,0,sizeof(float)*NumEls*NPE*3);
    h_Vol_MType_K = new float4[NumEls]; memset(h_Vol_MType_K,0,sizeof(float)*NumEls*4);
    h_ElasticParams = new float4[NumEls*numElasticParamsF4]; memset(h_ElasticParams,0,sizeof(float)*NumEls*4*numElasticParamsF4);
    h_NumProny = new int2[NumEls]; memset(h_NumProny,0,sizeof(int)*NumEls*2);
    h_bkwdEulerIso = new float2[NumEls*maxNumViscTerms.x]; memset(h_bkwdEulerIso,0,sizeof(float)*NumEls*2*maxNumViscTerms.x);
    h_bkwdEulerVol = new float2[NumEls*maxNumViscTerms.y]; memset(h_bkwdEulerVol,0,sizeof(float)*NumEls*2*maxNumViscTerms.y);
    h_SPKa = new float3[NumEls]; memset(h_SPKa,0,sizeof(float)*3*NumEls);
    h_SPKb = new float3[NumEls]; memset(h_SPKb,0,sizeof(float)*3*NumEls);
    h_Ca = new float3[NumEls]; memset(h_Ca,0,sizeof(float)*3*NumEls);
    h_Cb = new float3[NumEls]; memset(h_Cb,0,sizeof(float)*3*NumEls);
    HGLame = new float2[NumEls]; memset(HGLame,0,sizeof(float)*2*NumEls);
    ElementVols = new float[NumEls];
  }
  h_Va = NULL;
  h_Pa = NULL;   

  // Compute element variables
  vector< vector<int2> >* NInd = new vector< vector<int2> >;
  NInd->resize(NumNodes);
  if (!strcmp(EType,"T4"))  // 4-node tetrahedra
    {
      ComputeElementT4Variables(Model,NInd);
    }
  else if (!strcmp(EType,"T4ANP"))
    {
      ComputeElementT4ANPVariables(Model,NInd);
    }
  else // 8-node hexahedra
    {
      ComputeElementH8Variables(Model,NInd);
    }

  if (Model->GetNumberOfShellElementSets() > 0) {
    this->SetShellSolver(new tledShellSolverGPU());
    this->GetShellSolver().Init(*this, *Model);
    MatAdd(M, &this->GetShellSolver().GetMass().front(), this->GetShellSolver().GetMass().size(), 1, M);
  }

  // Compute node variables
  InstantiateTimeStepper();
  int* EPN = new int[NumNodes]; // Elements per node (different for each node)
  FCdsLength = 0;
  for (int node = 0; node < NumNodes; node++)
    {
      // Arrays for contructing NodeMap and FCds
      EPN[node] = (int)(*NInd)[node].size();
      FCdsLength += EPN[node];
    }
  assert(NumEls > 0 || FCdsLength == 0);
  if (FCdsLength > 0) h_FCds = new int2[FCdsLength];
  else h_FCds = NULL;
  int2* pFCds = h_FCds;
  int ind = 0;	// Keeps track of position in FCds
  // Second node loop to construct NodeMap and FCds
  for (int node = 0; node < NumNodes; node++)
    {
      h_NodeMap[node].x = ind;
      h_NodeMap[node].y = EPN[node];
      ind += EPN[node];
      for (int el = 0; el < EPN[node]; el++)
	{
	  *pFCds = ((*NInd)[node])[el];
	  pFCds++;
	}
    }

  float4* h_NodeCds = new float4[NumNodes];
  float* NodeCds = Mesh->GetAllNodeCds();
  for (int i = 0; i < NumNodes; i++)
    {
      h_NodeCds[i] = make_float4(NodeCds[3*i],NodeCds[3*i+1],NodeCds[3*i+2],0);
    }
   
  // Memory sizes
  int factorANP = 0;
  if (ANP)
    factorANP = 1;
  int factorHG = 0;
  if (!strcmp(EType,"H8"))
    factorHG = 1;

  const unsigned int MemSzDhDx = sizeof(float4)*NumEls*NPE*3/4;	// NPE*3/4 = 3 for T4, and 6 for H8
  const unsigned int MemSzFEl = sizeof(float4)*NumEls*NPE;
  const unsigned int MemSzEInd = sizeof(int4)*NumEls*NPE/4;	// NPE/4 = 1 for T4, and 2 for H8
  const unsigned int MemSzFCds = sizeof(int2)*FCdsLength;
  const unsigned int MemSzHG = sizeof(float4)*NumEls*16*factorHG;
  const unsigned int MemSzVol = sizeof(float4)*NumEls;
  const unsigned int MemSzElast = sizeof(float4)*NumEls*numElasticParamsF4;
  const unsigned int MemSzProny = sizeof(int2)*NumEls;
  const unsigned int MemSzEulerI = sizeof(float2)*NumEls*maxNumViscTerms.x;
  const unsigned int MemSzEulerV = sizeof(float2)*NumEls*maxNumViscTerms.y;
  const unsigned int MemSzStateI = sizeof(float4)*NumEls*2*maxNumViscTerms.x;	// Will produce 8 stresses for each Prony term (6 req'd)
  const unsigned int MemSzStateV = sizeof(float4)*NumEls*2*maxNumViscTerms.y;
  const unsigned int MemSzSPK = sizeof(float4)*NumEls;
  const unsigned int MemSzU = sizeof(float4)*NumNodes;
  const unsigned int MemSzNodeMap = sizeof(int2)*NumNodes;
  const unsigned int MemSzMask = sizeof(int4)*NumNodes;
  const unsigned int MemSzDiv = sizeof(bool);
  const unsigned int MemSzPa = sizeof(float)*NumNodes*factorANP;
   
  // Allocate pinned host memory
  tledCUDAHelpers::AllocateHostMemory(h_Uload, NumNodes);
  memset(h_Uload,0,sizeof(float)*NumNodes*4);
  cudaMallocHost((void**)&h_R,MemSzU);			// Imposed forces
  memset(h_R,0,sizeof(float)*NumNodes*4);
   
  // Allocate device memory
  if (NumEls > 0) {
    cudaMalloc((void**)&d_DhDx,MemSzDhDx);
    cudaMalloc((void**)&d_EInd,MemSzEInd);
    cudaMalloc((void**)&d_FEl,MemSzFEl);
    cudaMalloc((void**)&d_FCds,MemSzFCds);
    cudaMalloc((void**)&d_ElasticParams,MemSzElast);
    cudaMalloc((void**)&d_Vol_MType_K,MemSzVol);

    if (!strcmp(EType,"H8")) cudaMalloc((void**)&d_HG,MemSzHG);
    else d_HG = NULL;
    if (maxNumViscTerms.x + maxNumViscTerms.y > 0) cudaMalloc((void**)&d_NumProny,MemSzProny);
    else d_NumProny = NULL;

    if (maxNumViscTerms.x > 0) {
      cudaMalloc((void**)&d_bkwdEulerIso,MemSzEulerI);
      cudaMalloc((void**)&d_StressStateIso,MemSzStateI);
    } else {
      d_bkwdEulerIso = NULL;
      d_StressStateIso = NULL;
    }

    if (maxNumViscTerms.y > 0) {
      cudaMalloc((void**)&d_bkwdEulerVol,MemSzEulerV);
      cudaMalloc((void**)&d_StressStateVol,MemSzStateV);
    } else {
      d_bkwdEulerVol = NULL;
      d_StressStateVol = NULL;
    }       

    if (ANP) {
      cudaMalloc((void**)&d_Pa,MemSzPa);
      cudaMalloc((void**)&d_Va,MemSzPa);
    } else {
      d_Pa = NULL;
      d_Va = NULL;
    }

    cudaMalloc((void**)&d_SPKa,MemSzSPK);
    cudaMalloc((void**)&d_SPKb,MemSzSPK);
    cudaMalloc((void**)&d_Ca,MemSzSPK);
    cudaMalloc((void**)&d_Cb,MemSzSPK);
  }
  cudaMalloc((void**)&d_CD,MemSzU);
  cudaMalloc((void**)&d_F,MemSzU);
  cudaMalloc((void**)&d_NodeMap,MemSzNodeMap);
  tledCUDAHelpers::AllocateDeviceMemory(d_FixMask, NumNodes);
  tledCUDAHelpers::SetZero(d_FixMask, NumNodes);
  cudaMalloc((void**)&d_Uload,MemSzU);
  tledCUDAHelpers::AllocateDeviceMemory(d_DispMask, NumNodes);
  tledCUDAHelpers::SetZero(d_DispMask, NumNodes);
  cudaMalloc((void**)&d_R,MemSzU);
  cudaMalloc((void**)&d_NodeCds,MemSzU);
  cudaMalloc((void**)&d_Divergence,MemSzDiv);
  d_Contacts = NULL;
   
  // Copy host data to device
  if (NumEls > 0) {
    cudaMemcpy(d_DhDx,h_DhDx,MemSzDhDx,cudaMemcpyHostToDevice);
    cudaMemcpy(d_EInd,h_EInd,MemSzEInd,cudaMemcpyHostToDevice);
    cudaMemcpy(d_FCds,h_FCds,MemSzFCds,cudaMemcpyHostToDevice);
    cudaMemcpy(d_Vol_MType_K,h_Vol_MType_K,MemSzVol,cudaMemcpyHostToDevice);
    cudaMemcpy(d_ElasticParams,h_ElasticParams,MemSzElast,cudaMemcpyHostToDevice);

    if (!strcmp(EType,"H8"))
      cudaMemcpy(d_HG,h_HG,MemSzHG,cudaMemcpyHostToDevice);
    if (maxNumViscTerms.x + maxNumViscTerms.y > 0) cudaMemcpy(d_NumProny,h_NumProny,MemSzProny,cudaMemcpyHostToDevice);
    if (maxNumViscTerms.x > 0) cudaMemcpy(d_bkwdEulerIso,h_bkwdEulerIso,MemSzEulerI,cudaMemcpyHostToDevice);
    if (maxNumViscTerms.y > 0) cudaMemcpy(d_bkwdEulerVol,h_bkwdEulerVol,MemSzEulerV,cudaMemcpyHostToDevice);
    else d_bkwdEulerVol = NULL;
    if (ANP) {
      cudaMemcpy(d_Pa,h_Pa,MemSzPa,cudaMemcpyHostToDevice);
      cudaMemcpy(d_Va,h_Va,MemSzPa,cudaMemcpyHostToDevice);
    } 
  }

  cudaMemcpy(d_CD,h_CD,MemSzU,cudaMemcpyHostToDevice);
  cudaMemcpy(d_NodeMap,h_NodeMap,MemSzNodeMap,cudaMemcpyHostToDevice);
  cudaMemcpy(d_NodeCds,h_NodeCds,MemSzU,cudaMemcpyHostToDevice);

  cudaMemcpyToSymbol(c_NumNodes,&NumNodes,sizeof(int));
  cudaMemcpyToSymbol(c_NumEls,&NumEls,sizeof(int));
  cudaMemcpyToSymbol(c_NPE,&NPE,sizeof(int));
  cudaMemcpyToSymbol(c_maxNumViscTerms,&maxNumViscTerms,sizeof(int2));
      
  // Initialise some variables
  cudaMemset(d_Divergence,0,sizeof(bool));
  if (maxNumViscTerms.x + maxNumViscTerms.y > 0) {
    cudaMemset(d_StressStateIso,0,MemSzStateI);
    cudaMemset(d_StressStateVol,0,MemSzStateV);
  }
      
  // Bind textures
  if (NumEls > 0) {
    tledCheckCUDAErrors(cudaBindTexture(0,txDhDx,d_DhDx,MemSzDhDx));
    tledCheckCUDAErrors(cudaBindTexture(0,txEInd,d_EInd,MemSzEInd));
    if (!strcmp(EType,"H8")) {
      tledCheckCUDAErrors(cudaBindTexture(0,txHG,d_HG,MemSzHG));
    }
    tledCheckCUDAErrors(cudaBindTexture(0,txFCds,d_FCds,MemSzFCds));
    tledCheckCUDAErrors(cudaBindTexture(0,txFEl,d_FEl,MemSzFEl));
    tledCheckCUDAErrors(cudaBindTexture(0,txVol_MType_K,d_Vol_MType_K,MemSzVol));
    tledCheckCUDAErrors(cudaBindTexture(0,txElasticParams,d_ElasticParams,MemSzElast));
    if (maxNumViscTerms.x + maxNumViscTerms.y > 0)
      tledCheckCUDAErrors(cudaBindTexture(0,txNumProny,d_NumProny,MemSzProny));
    if (maxNumViscTerms.x > 0) {
      tledCheckCUDAErrors(cudaBindTexture(0,txEulerIso,d_bkwdEulerIso,MemSzEulerI));
      tledCheckCUDAErrors(cudaBindTexture(0,txStateIso,d_StressStateIso,MemSzStateI));
    }
    if (maxNumViscTerms.y > 0) {
      tledCheckCUDAErrors(cudaBindTexture(0,txEulerVol,d_bkwdEulerVol,MemSzEulerV));
      tledCheckCUDAErrors(cudaBindTexture(0,txStateVol,d_StressStateVol,MemSzStateV));
    }
    if (ANP) {
      tledCheckCUDAErrors(cudaBindTexture(0,txPa,d_Pa,MemSzPa));
      tledCheckCUDAErrors(cudaBindTexture(0,txVa,d_Va,MemSzPa));
    }
  }

  tledCheckCUDAErrors(cudaBindTexture(0,txNodeMap,d_NodeMap,MemSzNodeMap));   
  tledCheckCUDAErrors(cudaBindTexture(0,txCD,d_CD,MemSzU));
  tledCheckCUDAErrors(cudaBindTexture(0,txFixMask,d_FixMask,MemSzMask));
  tledCheckCUDAErrors(cudaBindTexture(0,txUload,d_Uload,MemSzU));
  tledCheckCUDAErrors(cudaBindTexture(0,txDispMask,d_DispMask,MemSzMask));
  tledCheckCUDAErrors(cudaBindTexture(0,txF4Fext,d_R,MemSzU));
  tledCheckCUDAErrors(cudaBindTexture(0,txNodeCds,d_NodeCds,MemSzU));
   
  delete NInd;
  delete[] EPN;
  delete[] h_NodeCds;
}

void tledSolverGPU::ComputeElementT4Variables(tledModel* Model, vector<vector<int2> >* NInd)
{
   vector<int> vEInd;
   double x[4][3];
   vector<float> NCds;
   double DhDr[4][3] = {-1, -1, -1,	// Shape function natural derivatives
                        1, 0, 0,
                        0, 1, 0,
                        0, 0, 1};
   double fDhDx[4][3];
   double J[3][3];	// Jacobian
   double detJ = 0;
   double invJ[3][3];
   double fVol = 0;
   double Mass = 0;
   int2 workInd;

   // Loop over element sets
   int* ElChecker = new int[NumEls];
   memset(ElChecker,0,sizeof(int)*NumEls);
   for (int ElSetNum = 0; ElSetNum < Model->GetNumElSets(); ElSetNum++)	// Loop over element sets
   {
      const char* MatType = Model->GetMatType(ElSetNum);
      const float D = this->GetElementSetDensity(ElSetNum);
      int currMType;
      // Get elastic params
      int currNumElasticParams = Model->GetNumElasticParams(ElSetNum);
      float* inElasticParams = new float[currNumElasticParams];
      Model->GetElasticParams(ElSetNum,inElasticParams);
      float4* outElasticParams = new float4[numElasticParamsF4];
      // Get visco params
      int currNumViscIsoTerms = Model->GetNumViscIsoTerms(ElSetNum);
      int currNumViscVolTerms = Model->GetNumViscVolTerms(ElSetNum);
      float* inViscParams = new float[2*(currNumViscIsoTerms+currNumViscVolTerms)];
      Model->GetViscoParams(ElSetNum,inViscParams);
      int2* N = new int2;
      float2* Ai = new float2[currNumViscIsoTerms];
      float2* Av = new float2[currNumViscVolTerms];

      assert(D == D);

      float2 wkHGLame;   // Lambda, Mu
      float ANPKappa; // Bulk modulus for use in ANP calcs
      if (!strcmp(MatType,"LE"))
      {
         currMType = 1;
         ComputeLEparams(inElasticParams,outElasticParams,&wkHGLame,&ANPKappa);
      }
      else if (!strcmp(MatType,"NH"))
      {
         currMType = 2;
         ComputeNHparams(inElasticParams,outElasticParams,&wkHGLame,&ANPKappa);
      }
      else if (!strcmp(MatType,"TI"))
      {
         currMType = 3;
         ComputeTIparams(inElasticParams,outElasticParams,&wkHGLame,&ANPKappa);
      }
      else if (!strcmp(MatType,"NHV"))
      {
         currMType = 4;
         ComputeNHparams(inElasticParams,outElasticParams,&wkHGLame,&ANPKappa);
         ComputeViscparams(inViscParams,currNumViscIsoTerms,currNumViscVolTerms,N,Ai,Av);
      }
      else if (!strcmp(MatType,"TIV"))
      {
         currMType = 5;
         ComputeTIparams(inElasticParams,outElasticParams,&wkHGLame,&ANPKappa);
         ComputeViscparams(inViscParams,currNumViscIsoTerms,currNumViscVolTerms,N,Ai,Av);
      }
      else if (!strcmp(MatType,"AB"))
      {
         currMType = 6;
         ComputeABparams(inElasticParams,outElasticParams,&wkHGLame,&ANPKappa);
      }
      else if (!strcmp(MatType,"PY"))
      {
         currMType = 7;
         ComputePYparams(inElasticParams,outElasticParams,&wkHGLame,&ANPKappa);
      }
      delete[] inElasticParams;
      delete[] inViscParams;

      vector<int> ElSet = Model->GetElSet(ElSetNum);
      // Loop over els in current set
      for (int ElSetPos = 0; ElSetPos < (int)ElSet.size(); ElSetPos++)
      {
         int el = ElSet[ElSetPos];

	 m_ElementElSetIndex[el] = ElSetNum;
         if (ElChecker[el] != 0)
            cerr << "\n!!! Element " << el << " listed in more than 1 element set" << endl;
         else
         {
            // Assemble material params
            for (int i = 0; i < numElasticParamsF4; i++)
            {
               h_ElasticParams[el*numElasticParamsF4 + i] = outElasticParams[i];
            }
            h_NumProny[el] = *N;
            for (int i = 0; i < currNumViscIsoTerms; i++)
            {
               h_bkwdEulerIso[el] = Ai[i];
            }
            for (int i = 0; i < currNumViscVolTerms; i++)
            {
               h_bkwdEulerVol[el] = Av[i];
            }

            // Assemble node coords
            vEInd = Model->GetElNodeInds(el);
            for (int i = 0; i < 4; i++)
            {
               NCds = Model->GetNodeCds(vEInd[i]);
               for (int j = 0; j < 3; j++)
               {
                  x[i][j] = (double)NCds[j];
               }
            }
            MatMult43T43(DhDr,x,J);
            // Compute element volume - note alternative, but equivalent formula
            MatDet33(J,&detJ);
            fVol = fabs(detJ/6);
            // Compute shape function global derivatives
            MatInv33(J,invJ);
            MatMult4333T(DhDr,invJ,fDhDx);
            // Compute element mass and add to nodal totals
            Mass = fVol*D/4;
            for (int i = 0; i < 4; i++)
            {
               M[vEInd[i]] += (float)Mass;
            }
            // Assemble NInd
            for (int node = 0; node < 4; node++)
            {
               workInd.x = el;
               workInd.y = node;
               (*NInd)[vEInd[node]].push_back(workInd);
            }

            // Collect final variables
            h_EInd[el].x = vEInd[0];
            h_EInd[el].y = vEInd[1];
            h_EInd[el].z = vEInd[2];
            h_EInd[el].w = vEInd[3];
            for (int i = 0; i < 3; i++)
            {
               h_DhDx[3*el+i].x = (float)fDhDx[0][i];
               h_DhDx[3*el+i].y = (float)fDhDx[1][i];
               h_DhDx[3*el+i].z = (float)fDhDx[2][i];
               h_DhDx[3*el+i].w = (float)fDhDx[3][i];
            }
            h_Vol_MType_K[el] = make_float4((float)fVol,(float)currMType,ANPKappa,0);
            ElementVols[el] = (float)fVol;
            HGLame[el] = wkHGLame;

            ElChecker[el] = 1;
         }
      }
      delete[] outElasticParams;
      delete N;
      delete[] Ai;
      delete[] Av;
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
   }

   delete[] ElChecker;
}

void tledSolverGPU::ComputeElementT4ANPVariables(tledModel* Model, vector<vector<int2> >* NInd)
{
   ComputeElementT4Variables(Model,NInd);
   ANP = true;
   h_Pa = new float[NumNodes]; memset(h_Pa,0,sizeof(float)*NumNodes);
   ComputeNodalVolumes();
}

void tledSolverGPU::ComputeElementH8Variables(tledModel* Model, vector<vector<int2> >* NInd)
{
   vector<int> vEInd;
   double x[8][3];
   vector<float> NCds;
   double J[3][3];
   double detJ = 0;
   double fVol = 0;
   const double a = 1./8;
   double DhDr[8][3] = {-a, -a, -a,	// Shape function natural derivatives
                        a, -a, -a,
                        a, a, -a,
                        -a, a, -a,
                        -a, -a, a,
                        a, -a, a,
                        a, a, a,
                        -a, a, a};
   double fDhDx[8][3];
   double Mass = 0;
   int2 workInd;

   h_HG = new float4[NumEls*16];
   double fHG[8][8];

   // Loop over element sets
   int* ElChecker = new int[NumEls];
   memset(ElChecker,0,sizeof(int)*NumEls);
   for (int ElSetNum = 0; ElSetNum < Model->GetNumElSets(); ElSetNum++)	// Loop over element sets
   {
      const char* MatType = Model->GetMatType(ElSetNum);
      const float D = this->GetElementSetDensity(ElSetNum);
      int currMType;
      // Get elastic params
      int currNumElasticParams = Model->GetNumElasticParams(ElSetNum);
      float* inElasticParams = new float[currNumElasticParams];
      Model->GetElasticParams(ElSetNum,inElasticParams);
      float4* outElasticParams = new float4[numElasticParamsF4];
      // Get visco params
      int currNumViscIsoTerms = Model->GetNumViscIsoTerms(ElSetNum);
      int currNumViscVolTerms = Model->GetNumViscVolTerms(ElSetNum);
      float* inViscParams = new float[2*(currNumViscIsoTerms+currNumViscVolTerms)];
      Model->GetViscoParams(ElSetNum,inViscParams);
      int2* N = new int2;
      float2* Ai = new float2[currNumViscIsoTerms];
      float2* Av = new float2[currNumViscVolTerms];

      assert(D == D);

      float2 wkHGLame;	// Lambda, Mu
      float ANPKappa;	// Bulk modulus for use in ANP calcs
      if (!strcmp(MatType,"LE"))
      {
         currMType = 1;
         ComputeLEparams(inElasticParams,outElasticParams,&wkHGLame,&ANPKappa);
      }
      else if (!strcmp(MatType,"NH"))
      {
         currMType = 2;
         ComputeNHparams(inElasticParams,outElasticParams,&wkHGLame,&ANPKappa);
      }
      else if (!strcmp(MatType,"TI"))
      {
         currMType = 3;
         ComputeTIparams(inElasticParams,outElasticParams,&wkHGLame,&ANPKappa);
      }
      else if (!strcmp(MatType,"NHV"))
      {
         currMType = 4;
         ComputeNHparams(inElasticParams,outElasticParams,&wkHGLame,&ANPKappa);
         ComputeViscparams(inViscParams,currNumViscIsoTerms,currNumViscVolTerms,N,Ai,Av);
      }
      else if (!strcmp(MatType,"TIV"))
      {
         currMType = 5;
         ComputeTIparams(inElasticParams,outElasticParams,&wkHGLame,&ANPKappa);
         ComputeViscparams(inViscParams,currNumViscIsoTerms,currNumViscVolTerms,N,Ai,Av);
      }
      else if (!strcmp(MatType,"AB"))
      {
         currMType = 6;
         ComputeABparams(inElasticParams,outElasticParams,&wkHGLame,&ANPKappa);
      }
      else if (!strcmp(MatType,"PY"))
      {
         currMType = 7;
         ComputePYparams(inElasticParams,outElasticParams,&wkHGLame,&ANPKappa);
      }
      delete[] inElasticParams;
      delete[] inViscParams;

      vector<int> ElSet = Model->GetElSet(ElSetNum);
      // Loop over els in current set
      for (int ElSetPos = 0; ElSetPos < (int)ElSet.size(); ElSetPos++)
      {
         int el = ElSet[ElSetPos];
         if (ElChecker[el] != 0)
            cerr << "\n!!! Element " << el << " listed in more than 1 element set" << endl;
         else
         {
            // Assemble material params
            for (int i = 0; i < numElasticParamsF4; i++)
            {
               h_ElasticParams[el*numElasticParamsF4 + i] = outElasticParams[i];
            }
            h_NumProny[el] = *N;
            for (int i = 0; i < currNumViscIsoTerms; i++)
            {
               h_bkwdEulerIso[el] = Ai[i];
            }
            for (int i = 0; i < currNumViscVolTerms; i++)
            {
               h_bkwdEulerVol[el] = Av[i];
            }

            // Assemble node coords
            vEInd = Model->GetElNodeInds(el);
            for (int i = 0; i < 8; i++)
            {
               NCds = Model->GetNodeCds(vEInd[i]);
               for (int j = 0; j < 3; j++)
               {
                  x[i][j] = (double)NCds[j];
               }
            }

            // Jacobian
            MatMult83T83(DhDr,x,J);
            MatDet33(J,&detJ);
            // Element volume
            fVol = ComputeHexVol(x);
            // Compute shape function global derivatives
            ComputeH8DhDx(x,fVol,fDhDx);
            // Compute element mass and add to nodal totals
            Mass = ((double)D)*fVol/8;
            for (int i = 0; i < 8; i++)
            {
               M[vEInd[i]] += (float)Mass;
            }
            // Assemble NInd
            for (int node = 0; node < 8; node++)
            {
               workInd.x = el;
               workInd.y = node;
               (*NInd)[vEInd[node]].push_back(workInd);
            }
            // Hourglass control parameters
            double a = 0;
            for (int i = 0; i < 8; i++)
            {
               for (int j = 0; j < 3; j++)
               {
                     a += fDhDx[i][j]*fDhDx[i][j];
               }
            }
            HGKappa = Model->GetHGKappa();
            double k = HGKappa*fVol*(wkHGLame.x+2*wkHGLame.y)*a/8;
            double Gamma[8][4] = {1,1,1,-1,
                                 -1,1,-1,1,
                                 1,-1,-1,-1,
                                 -1,-1,1,1,
                                 1,-1,-1,1,
                                 -1,-1,1,-1,
                                 1,1,1,1,
                                 -1,1,-1,-1};
            double A[8][8];
            MatMult8383T(fDhDx,x,A);
            double gamma[8][4];
            MatMult8884(A,Gamma,gamma);
            MatSubtract(&Gamma[0][0],&gamma[0][0],8,4,&gamma[0][0]);
            MatMult8484T(gamma,gamma,fHG);
            MatMultScalar(&fHG[0][0],8,8,k,&fHG[0][0]);

            // Collect final variables
            for (int i = 0; i < 2; i++)
            {
               h_EInd[2*el+i].x = vEInd[4*i];
               h_EInd[2*el+i].y = vEInd[4*i+1];
               h_EInd[2*el+i].z = vEInd[4*i+2];
               h_EInd[2*el+i].w = vEInd[4*i+3];
            }
            for (int i = 0; i < 3; i++)
            {
               h_DhDx[6*el+i].x = (float)fDhDx[0][i];
               h_DhDx[6*el+i].y = (float)fDhDx[1][i];
               h_DhDx[6*el+i].z = (float)fDhDx[2][i];
               h_DhDx[6*el+i].w = (float)fDhDx[3][i];

               h_DhDx[6*el+i+3].x = (float)fDhDx[4][i];
               h_DhDx[6*el+i+3].y = (float)fDhDx[5][i];
               h_DhDx[6*el+i+3].z = (float)fDhDx[6][i];
               h_DhDx[6*el+i+3].w = (float)fDhDx[7][i];
            }
            h_Vol_MType_K[el] = make_float4((float)(8*detJ),(float)currMType,ANPKappa,0);
            ElementVols[el] = (float)fVol;
            float4 b;
            for (int i = 0; i < 8; i++)
            {
               for (int j = 0; j < 2; j++)
               {
                  b.x = (float)fHG[i][4*j];
                  b.y = (float)fHG[i][4*j+1];
                  b.z = (float)fHG[i][4*j+2];
                  b.w = (float)fHG[i][4*j+3];
                  h_HG[16*el + 2*i+j] = b;
               }
            }
            HGLame[el] = wkHGLame;

            ElChecker[el] = 1;
         }
      }
      delete[] outElasticParams;
      delete N;
      delete[] Ai;
      delete[] Av;
   }
   // Check that all elements instantiated
   {
     bool someBadEls = false;

     for (int i = 0; i < NumEls; i++) if (ElChecker[i] == 0) {
	 if (!someBadEls) std::cerr << "\n!!! Some elements not assigned to an element list:\n"; 
	 std::cerr << i << std::endl;
	 someBadEls = true;
       }
   }

   delete[] ElChecker;
}

void tledSolverGPU::ComputeNodalVolumes()
{
   if (h_Va == NULL)
      h_Va = new float[NumNodes];
   memset(h_Va,0,sizeof(float)*NumNodes);
   for (int i = 0; i < NumEls; i++)
   {
      float Ve = h_Vol_MType_K[i].x;
      int4 EInd = h_EInd[i];
      h_Va[EInd.x] += Ve/4;
      h_Va[EInd.y] += Ve/4;
      h_Va[EInd.z] += Ve/4;
      h_Va[EInd.w] += Ve/4;
   }
}

double tledSolverGPU::ComputeHexVol(double x[8][3])
{
   // Calc CIJK first
   double C[8][8][8];
   ComputeCIJK(C);
   // Calc volume
   double V = 0;
   for (int I = 0; I < 8; I++)
   {
      for (int J = 0; J < 8; J++)
      {
         for (int K = 0; K < 8; K++)
         {
            V += x[I][0]*x[J][1]*x[K][2]*C[I][J][K];
         }
      }
   }
   return V;
}

void tledSolverGPU::ComputeH8DhDx(double x[8][3], double V, double fDhDx[8][3])
{
   // Calc B matrix
   double B[8][3];
   ComputeBmat(x,B);
   for (int i = 0; i < 8; i++)
   {
      for (int j = 0; j < 3; j++)
      {
         fDhDx[i][j] = B[i][j]/V;
      }
   }
}

void tledSolverGPU::ComputeBmat(double x[8][3], double B[8][3])
{
   // Calc CIJK first
   double C[8][8][8];
   ComputeCIJK(C);
   // Calc B
   memset(B,0,sizeof(double)*8*3);
   for (int I = 0; I < 8; I++)
   {
      for (int J = 0; J < 8; J++)
      {
         for (int K = 0; K < 8; K++)
         {
            B[I][0] += x[J][1]*x[K][2]*C[I][J][K];
            B[I][1] += x[J][2]*x[K][0]*C[I][J][K];
            B[I][2] += x[J][0]*x[K][1]*C[I][J][K];
         }
      }
   }
}

void tledSolverGPU::ComputeCIJK(double C[8][8][8])
{
   double a = (double)(1./12);
   double Ctemp[8*8*8] = 
   {0,0,0,0,0,0,0,0,
   0,0,-a,-a,a,a,0,0,
   0,a,0,-a,0,0,0,0,
   0,a,a,0,-a,0,0,-a,
   0,-a,0,a,0,-a,0,a,
   0,-a,0,0,a,0,0,0,
   0,0,0,0,0,0,0,0,
   0,0,0,a,-a,0,0,0,

   0,0,a,a,-a,-a,0,0,
   0,0,0,0,0,0,0,0,
   -a,0,0,-a,0,a,a,0,
   -a,0,a,0,0,0,0,0,
   a,0,0,0,0,-a,0,0,
   a,0,-a,0,a,0,-a,0,
   0,0,-a,0,0,a,0,0,
   0,0,0,0,0,0,0,0,

   0,-a,0,a,0,0,0,0,
   a,0,0,a,0,-a,-a,0,
   0,0,0,0,0,0,0,0,
   -a,-a,0,0,0,0,a,a,
   0,0,0,0,0,0,0,0,
   0,a,0,0,0,0,-a,0,
   0,a,0,-a,0,a,0,-a,
   0,0,0,-a,0,0,a,0,

   0,-a,-a,0,a,0,0,a,
   a,0,-a,0,0,0,0,0,
   a,a,0,0,0,0,-a,-a,
   0,0,0,0,0,0,0,0,
   -a,0,0,0,0,0,0,a,
   0,0,0,0,0,0,0,0,
   0,0,a,0,0,0,0,-a,
   -a,0,a,0,-a,0,a,0,

   0,a,0,-a,0,a,0,-a,
   -a,0,0,0,0,a,0,0,
   0,0,0,0,0,0,0,0,
   a,0,0,0,0,0,0,-a,
   0,0,0,0,0,0,0,0,
   -a,-a,0,0,0,0,a,a,
   0,0,0,0,0,-a,0,a,
   a,0,0,a,0,-a,-a,0,

   0,a,0,0,-a,0,0,0,
   -a,0,a,0,-a,0,a,0,
   0,-a,0,0,0,0,a,0,
   0,0,0,0,0,0,0,0,
   a,a,0,0,0,0,-a,-a,
   0,0,0,0,0,0,0,0,
   0,-a,-a,0,a,0,0,a,
   0,0,0,0,a,0,-a,0,

   0,0,0,0,0,0,0,0,
   0,0,a,0,0,-a,0,0,
   0,-a,0,a,0,-a,0,a,
   0,0,-a,0,0,0,0,a,
   0,0,0,0,0,a,0,-a,
   0,a,a,0,-a,0,0,-a,
   0,0,0,0,0,0,0,0,
   0,0,-a,-a,a,a,0,0,

   0,0,0,-a,a,0,0,0,
   0,0,0,0,0,0,0,0,
   0,0,0,a,0,0,-a,0,
   a,0,-a,0,a,0,-a,0,
   -a,0,0,-a,0,a,a,0,
   0,0,0,0,-a,0,a,0,
   0,0,a,a,-a,-a,0,0,
   0,0,0,0,0,0,0,0};

   memcpy(C,Ctemp,sizeof(double)*8*8*8);
}

void tledSolverGPU::ComputeLEparams(float* inMatParams, float4* outMatParams, float2* wkHGLame, float* ANPKappa)
{
   float E = inMatParams[0];
   float P = inMatParams[1];
   (*outMatParams).x = E*(1 - P)/(2*(1 + P)*(1 - 2*P));  // c1
   (*outMatParams).y = (*outMatParams).x*P/(1 - P);      // c2
   (*outMatParams).z = (*outMatParams).x + 2*(*outMatParams).y;   // c3
   (*outMatParams).w = (*outMatParams).x*(1 - 2*P)/(1 - P);    // c4

   // Compute HG Lame params
   (*wkHGLame).x = E*P/((1+P)*(1-2*P));
   (*wkHGLame).y = E/(2*(1+P));
   // Compute ANPKappa
   *ANPKappa = E/(3*(1-2*P));
}

void tledSolverGPU::ComputeNHparams(float* inMatParams, float4* outMatParams, float2* wkHGLame, float* ANPKappa)
{
   (*outMatParams).x = inMatParams[0]; // Mu
   (*outMatParams).y = inMatParams[1]; // K
   (*outMatParams).z = 0;     // unused
   (*outMatParams).w = 0;     // unused

   // Compute HG Lame params
   (*wkHGLame).x = (*outMatParams).y - 2*(*outMatParams).x/3;
   (*wkHGLame).y = (*outMatParams).x;
   // Compute ANPKappa
   *ANPKappa = inMatParams[1];
}

void tledSolverGPU::ComputeTIparams(float* inMatParams, float4* outMatParams, float2* wkHGLame, float* ANPKappa)
{
   outMatParams[0].x = inMatParams[0]; // Mu
   outMatParams[0].y = inMatParams[1]; // K
   outMatParams[0].z = inMatParams[2]; // Eta
   // Check that a is normalised
   float a[3] = {inMatParams[3],inMatParams[4],inMatParams[5]};
   float norma = sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2]);
   if (norma != 1.0f)
   {
      // Normalise the direction vector
      for (int i = 0; i < 3; i++)
         a[i] /= norma;
   }
   outMatParams[0].w = a[0];
   outMatParams[1].x = a[1];
   outMatParams[1].y = a[2];
   outMatParams[1].z = 0;  // unused
   outMatParams[1].w = 0;  // unused

   // Compute HG Lame params
   (*wkHGLame).x = outMatParams[0].y - 2*(outMatParams[0].x)/3;
   (*wkHGLame).y = outMatParams[0].x;
   // Compute ANPKappa
   *ANPKappa = inMatParams[1];
}

void tledSolverGPU::ComputeABparams(float* inMatParams, float4* outMatParams, float2* wkHGLame, float* ANPKappa)
{
   (*outMatParams).x = inMatParams[0]; // Mu
   (*outMatParams).y = inMatParams[1]; // Lm
   (*outMatParams).z = inMatParams[2]; // K
   (*outMatParams).w = 0;     // unused
   
   // Compute HG Lame params
   (*wkHGLame).x = (*outMatParams).z - 2*(*outMatParams).x/3;
   (*wkHGLame).y = (*outMatParams).x;
   // Compute ANPKappa
   *ANPKappa = inMatParams[2];
}

void tledSolverGPU::ComputePYparams(float* inMatParams, float4* outMatParams, float2* wkHGLame, float* ANPKappa)
{
   outMatParams[0].x = inMatParams[0]; // c10
   outMatParams[0].y = inMatParams[1]; // c01
   outMatParams[0].z = inMatParams[2]; // c20
   outMatParams[0].w = inMatParams[3]; // c02
   outMatParams[1].x = inMatParams[4]; // c11
   outMatParams[1].y = inMatParams[5]; // K
   outMatParams[1].z = 0;  // unused
   outMatParams[1].w = 0;  // unused
   
   // Compute HG Lame params
   float Mu = 2*(outMatParams[0].x+outMatParams[0].y); // mu = 2*(c10+c01)
   (*wkHGLame).x = outMatParams[1].y - 2*Mu/3;
   (*wkHGLame).y = Mu;
   // Compute ANPKappa
   *ANPKappa = inMatParams[5];
}

void tledSolverGPU::ComputeViscparams(float* inMatViscParams, int Ni, int Nv, int2* N, float2* Ai, float2* Av)
{
   *N = make_int2(Ni,Nv);

   // Set up isochoric terms
   for (int i = 0; i < Ni; i++)
   {
      Ai[i].x = (float)( Dt*inMatViscParams[2*i]/(Dt + inMatViscParams[2*i+1]) );
      Ai[i].y = (float)( inMatViscParams[2*i+1]/(Dt + inMatViscParams[2*i+1]) );
   }
   // Set up volumetric terms
   for (int i = 0; i < Nv; i++)
   {
      Av[i].x = (float)( Dt*inMatViscParams[2*Ni+2*i]/(Dt + inMatViscParams[2*Ni+2*i+1]) );
      Av[i].y = (float)( inMatViscParams[2*Ni+2*i+1]/(Dt + inMatViscParams[2*Ni+2*i+1]) );
   }
}

void tledSolverGPU::PerformStep()
{
   // Memory sizes
   const unsigned int MemSzU = sizeof(float4)*NumNodes;
   
   // Load new displacements
   cudaMemcpy(d_Uload,h_Uload,MemSzU,cudaMemcpyHostToDevice);
         
   // Load external forces
   cudaMemcpy(d_R,h_R,MemSzU,cudaMemcpyHostToDevice);
   
   // Set up execution parameters
   const unsigned int BlockSzP1 = 128;
   const unsigned int BlockSzP2 = 64;
   const unsigned int GridSzP1 = (int)ceil((double)NumEls/(double)BlockSzP1);	// Number of threads for each pass
   const unsigned int GridSzP2 = (int)ceil((double)NumNodes/(double)BlockSzP2);
   dim3 ThreadsP1(BlockSzP1,1,1);
   dim3 GridP1(GridSzP1,1,1);
   dim3 ThreadsP2(BlockSzP2,1,1);
   dim3 GridP2(GridSzP2,1,1);
   
   if (NumEls > 0) {
     if (!strcmp(EType,"T4"))
       {
	 ComputeNewForcesT4_kernel <<< GridP1,ThreadsP1 >>> (d_FEl,d_StressStateIso,d_StressStateVol,d_SPKa,d_SPKb,d_Ca,d_Cb);
       }
     else if (!strcmp(EType,"H8"))
       {
	 ComputeNewForcesH8_kernel <<< GridP1,ThreadsP1 >>> (d_FEl,d_StressStateIso,d_StressStateVol,d_SPKa,d_SPKb,d_Ca,d_Cb);
       }
     else // T4ANP
       {
	 ComputeElementPressuresT4ANP_kernel <<< GridP1,ThreadsP1 >>> (d_FEl);
	 ComputeNodalPressuresT4ANP_kernel <<< GridP2,ThreadsP2 >>> (d_Pa);
	 ComputeModifiedElementForcesT4ANP_kernel <<< GridP1,ThreadsP1 >>> (d_FEl,d_StressStateIso,d_StressStateVol,d_SPKa,d_SPKb,d_Ca,d_Cb);
       }
   }

   if (this->HasMembrane()) {
     dynamic_cast<tledShellSolverGPU&>(this->GetShellSolver()).ComputeNewForces();
     dynamic_cast<tledShellSolverGPU&>(this->GetShellSolver()).ComputeNodalForces(d_R);
   }     

   if (typeid(*mp_TimeStepper) == typeid(tledCentralDifferenceTimeStepperGPU)) {
     ComputeNewDisps_kernel<tledCentralDifferenceTimeStepperGPU::SolutionVariables> <<<GridP2, ThreadsP2>>> (d_F, static_cast<tledCentralDifferenceTimeStepperGPU::SolutionVariables*>(mp_TimeStepper->GetSolutionDeviceBuffer()), d_Contacts, d_Divergence);
   } else if (typeid(*mp_TimeStepper) == typeid(tledNewmarkTimeStepperGPU)) {
     ComputeNewDisps_kernel<tledNewmarkTimeStepperGPU::SolutionVariables> <<<GridP2, ThreadsP2>>> (d_F, static_cast<tledNewmarkTimeStepperGPU::SolutionVariables*>(mp_TimeStepper->GetSolutionDeviceBuffer()), d_Contacts, d_Divergence);
   } else {
     tledFatalError("Unsupported time ODE solver class.");
   }

#ifdef GPU_GP_CONTACT
   if (this->GetContactManager().DoUnstructuredContacts()) {
     tledUnstructuredContactManager &r_unstructuredMeshContacts = this->GetContactManager().GetUnstructuredContactManager();

     r_unstructuredMeshContacts.Update();

     if (r_unstructuredMeshContacts.GetNumberOfRigidSurfaces() > 0) {
       r_unstructuredMeshContacts.ComputeDeformableRigidContactResponses(d_R, this->GetAllOnDeviceNextDisplacements(), this->GetAllOnDeviceCurrentDisplacements());
     }

     if (r_unstructuredMeshContacts.DoDeformableDeformableContactHandling()) {
       r_unstructuredMeshContacts.ComputeDeformableDeformableContactResponses(d_R, this->GetAllOnDeviceNextDisplacements(), this->GetAllOnDeviceCurrentDisplacements());
     }

     if (typeid(*mp_TimeStepper) == typeid(tledCentralDifferenceTimeStepperGPU)) {
       ComputeNewDisps_kernel<tledCentralDifferenceTimeStepperGPU::SolutionVariables> <<<GridP2, ThreadsP2>>> (d_F, static_cast<tledCentralDifferenceTimeStepperGPU::SolutionVariables*>(mp_TimeStepper->GetSolutionDeviceBuffer()), d_Contacts, d_Divergence);
     } else if (typeid(*mp_TimeStepper) == typeid(tledNewmarkTimeStepperGPU)) {
       ComputeNewDisps_kernel<tledNewmarkTimeStepperGPU::SolutionVariables> <<<GridP2, ThreadsP2>>> (d_F, static_cast<tledNewmarkTimeStepperGPU::SolutionVariables*>(mp_TimeStepper->GetSolutionDeviceBuffer()), d_Contacts, d_Divergence);
     } else {
       std::abort();
     }

     r_unstructuredMeshContacts.FinishContactHandling();
   }
#endif

   mp_TimeStepper->FinishTimeStep();
}

void tledSolverGPU::GetDivergence(bool* Div)
{
   cudaMemcpy(Div,d_Divergence,sizeof(bool),cudaMemcpyDeviceToHost);
}

void tledSolverGPU::UnsetDivergence() {
  bool div;

  div = false;
  tledCUDAHelpers::CopyToDevice(d_Divergence, &div);
}

void tledSolverGPU::SetDisps(vector<int>* IndX, vector<float>* UX, vector<int>* IndY, vector<float>* UY,
                             vector<int>* IndZ, vector<float>* UZ)
{
  tledCUDAHostMemoryBlock &r_dispMask = tledCUDAHostMemoryBlock::GetNextFreeBufferWithSize<int4>(NumNodes);

   // Reset all vals to 0
   memset(h_Uload,0,sizeof(float)*NumNodes*4);

   if (IndX->size() != UX->size() || IndY->size() != UY->size() || IndZ->size() != UZ->size()) {
     tledLogErrorStream(tledHelper::FatalError() << "Input vector sizes are incompatible with index arrays.");
   }

   // Set new disps
   std::fill(r_dispMask.GetBuffer<int4>(), r_dispMask.GetBuffer<int4>() + NumNodes, make_int4(0, 0, 0, 0));
   for (int i = 0; i < (int)IndX->size(); i++)
   {
      h_Uload[(*IndX)[i]].x += (*UX)[i];
      r_dispMask.GetBuffer<int4>()[(*IndX)[i]].x = 1;
   }
   for (int i = 0; i < (int)IndY->size(); i++)
   {
      h_Uload[(*IndY)[i]].y += (*UY)[i];
      r_dispMask.GetBuffer<int4>()[(*IndY)[i]].y = 1;
   }
   for (int i = 0; i < (int)IndZ->size(); i++)
   {
      h_Uload[(*IndZ)[i]].z += (*UZ)[i];
      r_dispMask.GetBuffer<int4>()[(*IndZ)[i]].z = 1;
   }

   tledCUDAHelpers::CopyToDevice(d_DispMask, r_dispMask.GetBuffer<int4>(), NumNodes);
   r_dispMask.ToggleActive();
}

void tledSolverGPU::SetExtForces(vector<int>* IndX, vector<float>* FX, vector<int>* IndY, vector<float>* FY,
                                 vector<int>* IndZ, vector<float>* FZ)
{
   // Reset all vals to 0
   memset(h_R,0,sizeof(float)*NumNodes*4);

   if (IndX->size() != FX->size() || IndY->size() != FY->size() || IndZ->size() != FZ->size()) {
     tledLogErrorStream(tledHelper::FatalError() << "Input vector sizes are incompatible with index arrays.");
   } 

   // Set new disps
   for (int i = 0; i < (int)IndX->size(); i++)
      h_R[(*IndX)[i]].x += (*FX)[i];
   for (int i = 0; i < (int)IndY->size(); i++)
      h_R[(*IndY)[i]].y += (*FY)[i];
   for (int i = 0; i < (int)IndZ->size(); i++)
      h_R[(*IndZ)[i]].z += (*FZ)[i];
}

float* tledSolverGPU::GetAllExtForces(float *p_dst) const {
  for (int i = 0; i < NumNodes; i++) tledCUDAHelpers::ConvertFromFloatN(p_dst + 3*i, h_R[i]);

  return p_dst;
}

void tledSolverGPU::SetFixed(vector<int>* IndX, vector<int>* IndY, vector<int>* IndZ)
{
  tledCUDAHostMemoryBlock &r_fixMask = tledCUDAHostMemoryBlock::GetNextFreeBufferWithSize<int4>(NumNodes);

   // Reset all vals to 1
   for (int i = 0; i < NumNodes; i++)
   {
      r_fixMask.GetBuffer<int4>()[i].x = 1;
      r_fixMask.GetBuffer<int4>()[i].y = 1;
      r_fixMask.GetBuffer<int4>()[i].z = 1;
      r_fixMask.GetBuffer<int4>()[i].w = 1;	// NB: last value not used
   }
   // Set new disp indices to 0
   for (int i = 0; i < (int)IndX->size(); i++)
      r_fixMask.GetBuffer<int4>()[(*IndX)[i]].x = 0;
   for (int i = 0; i < (int)IndY->size(); i++)
      r_fixMask.GetBuffer<int4>()[(*IndY)[i]].y = 0;
   for (int i = 0; i < (int)IndZ->size(); i++)
      r_fixMask.GetBuffer<int4>()[(*IndZ)[i]].z = 0;

   // Send new fix mask to device
   tledCUDAHelpers::CopyToDevice(d_FixMask, r_fixMask.GetBuffer<int4>(), NumNodes);
   r_fixMask.ToggleActive();
}

void tledSolverGPU::SetContactManager(tledContactManager* contacts)
{   
   tledSolver::SetContactManager(contacts);

   int numContactCyls = this->GetContactManager().GetNumContactCyls();
   int numContactPrbs = this->GetContactManager().GetNumContactPrbs();
   int numContactPlts = this->GetContactManager().GetNumContactPlts();
   
   cudaMemcpyToSymbol(c_numContactCyls,&numContactCyls,sizeof(int));
   cudaMemcpyToSymbol(c_numContactPrbs,&numContactPrbs,sizeof(int));
   cudaMemcpyToSymbol(c_numContactPlts,&numContactPlts,sizeof(int));

   if (this->GetContactManager().GetNumContactObjects() > 0)
      d_Contacts = this->GetContactManager().GetContactsDevicePointer();
}

void tledSolverGPU::GetForces(vector<int>* NodeInd, vector<float>* Forces)
{
   // Get forces from device
   const unsigned int MemSzU = sizeof(float4)*NumNodes;
   float4* l_F = (float4*)malloc(MemSzU);
   cudaMemcpy(l_F,d_F,MemSzU,cudaMemcpyDeviceToHost);
   gpuV4ToV3(l_F, NumNodes, h_F);
   
   free(l_F);

   for (int i = 0; i < (int)(NodeInd->size()); i++)
   {
      Forces->push_back( h_F[3*(*NodeInd)[i]].x );
      Forces->push_back( h_F[3*(*NodeInd)[i]].y );
      Forces->push_back( h_F[3*(*NodeInd)[i]].z );
   }
}

float* tledSolverGPU::GetAllForces()
{
   // Get forces from device
   const unsigned int MemSzU = sizeof(float4)*NumNodes;
   float4* l_F = (float4*)malloc(MemSzU);
   cudaMemcpy(l_F,d_F,MemSzU,cudaMemcpyDeviceToHost);
   gpuV4ToV3(l_F, NumNodes, h_F);
   
   free(l_F);
   
   float* F = new float[NumNodes*3];
   for (int i = 0; i < NumNodes; i++)
   {
      F[3*i] = h_F[i].x;
      F[3*i+1] = h_F[i].y;
      F[3*i+2] = h_F[i].z;
   }
   return F;
}

void tledSolverGPU::GetDisps(vector<int>* NodeInd, vector<float>* Disps)
{
   // Get disps from device
  Disps->reserve(3*NodeInd->size());
  for (vector<int>::const_iterator ic_n = NodeInd->begin(); ic_n < NodeInd->end(); ic_n++) {
    Disps->push_back(mp_TimeStepper->GetCurrentDisplacements()[*ic_n*3]);
  }
}

float* tledSolverGPU::GetAllNextDisps() {
  return mp_TimeStepper->GetNextDisplacements();
}

float* tledSolverGPU::GetAllDisps() {
  return mp_TimeStepper->GetCurrentDisplacements();
}

void tledSolverGPU::GetNodeVMStress(float* NodeSVM)
{
   const unsigned int MemSzSPK = sizeof(float4)*NumEls;
   float4* l_SPKa = (float4*)malloc(MemSzSPK);
   float4* l_SPKb = (float4*)malloc(MemSzSPK);
   cudaMemcpy(l_SPKa,d_SPKa,MemSzSPK,cudaMemcpyDeviceToHost);
   cudaMemcpy(l_SPKb,d_SPKb,MemSzSPK,cudaMemcpyDeviceToHost);
   gpuV4ToV3(l_SPKa,NumEls,h_SPKa);
   gpuV4ToV3(l_SPKb,NumEls,h_SPKb);
   
   free(l_SPKa);
   free(l_SPKb);

   // Compute element VM stresses
   SVM = new float[NumEls];
   float Sv[6];
   for (int i = 0; i < NumEls; i++)
   {
      Sv[0] = h_SPKa[i].x; Sv[1] = h_SPKa[i].y; Sv[2] = h_SPKa[i].z;
      Sv[3] = h_SPKb[i].x; Sv[4] = h_SPKb[i].y; Sv[5] = h_SPKb[i].z;
      SVM[i] = sqrt( Sv[0]*Sv[0] + Sv[1]*Sv[1] + Sv[2]*Sv[2]
                     - Sv[0]*Sv[1] - Sv[1]*Sv[2] - Sv[0]*Sv[2]
                     + 3*(Sv[3]*Sv[3] + Sv[4]*Sv[4] + Sv[5]*Sv[5]) );
   }

   // Average over nodes
   memset(NodeSVM,0,sizeof(float)*NumNodes);
   int* NodeValence = new int[NumNodes];
   memset(NodeValence,0,sizeof(int)*NumNodes);
   vector<int> NodeList;
   float temp;
   for (int i = 0; i < NumEls; i++)
   {
      NodeList = Mesh->GetElNodeInds(i);
      for (int j = 0; j < (int)NodeList.size(); j++)
      {
         temp = NodeSVM[NodeList[j]]*NodeValence[NodeList[j]];
         temp += SVM[i];
         NodeValence[NodeList[j]] += 1;
         temp /= NodeValence[NodeList[j]];
         NodeSVM[NodeList[j]] = temp;
      }
   }

   delete NodeValence;
   delete SVM;
}

void tledSolverGPU::GetNodeVMStrain(float* NodeEVM)
{
   const unsigned int MemSzC = sizeof(float4)*NumEls;
   float4* l_Ca = (float4*)malloc(MemSzC);
   float4* l_Cb = (float4*)malloc(MemSzC);
   cudaMemcpy(l_Ca,d_Ca,MemSzC,cudaMemcpyDeviceToHost);
   cudaMemcpy(l_Cb,d_Cb,MemSzC,cudaMemcpyDeviceToHost);
   gpuV4ToV3(l_Ca,NumEls,h_Ca);
   gpuV4ToV3(l_Cb,NumEls,h_Cb);
   
   free(l_Ca);
   free(l_Cb);

   // Compute element VM strains
   EVM = new float[NumEls];
   float Ev[6];
   float L, M;
   for (int i = 0; i < NumEls; i++)
   {
      // E = (C-I)/2
      Ev[0] = (h_Ca[i].x - 1)/2; Ev[1] = (h_Ca[i].y - 1)/2; Ev[2] = (h_Ca[i].z - 1)/2;
      Ev[3] = h_Cb[i].x/2; Ev[4] = h_Cb[i].y/2; Ev[5] = h_Cb[i].z/2;
      // Get Poisson ratio for current element...
      L = HGLame[i].x;
      M = HGLame[i].y;
      float P = L/(2*(L+M));
      // Von Mises strain
      EVM[i] = sqrt( Ev[0]*Ev[0] + Ev[1]*Ev[1] + Ev[2]*Ev[2]
                     - Ev[0]*Ev[1] - Ev[1]*Ev[2] - Ev[0]*Ev[2]
                     + 3*(Ev[3]*Ev[3] + Ev[4]*Ev[4] + Ev[5]*Ev[5]) )/(1+P);
   }

   // Average over nodes
   memset(NodeEVM,0,sizeof(float)*NumNodes);
   int* NodeValence = new int[NumNodes];
   memset(NodeValence,0,sizeof(int)*NumNodes);
   vector<int> NodeList;
   float temp;
   for (int i = 0; i < NumEls; i++)
   {
      NodeList = Mesh->GetElNodeInds(i);
      for (int j = 0; j < (int)NodeList.size(); j++)
      {
         temp = NodeEVM[NodeList[j]]*NodeValence[NodeList[j]];
         temp += EVM[i];
         NodeValence[NodeList[j]] += 1;
         temp /= NodeValence[NodeList[j]];
         NodeEVM[NodeList[j]] = temp;
      }
   }

   delete NodeValence;
   delete EVM;
}

void tledSolverGPU::GetMassVector(float* mass) const
{
   memcpy(mass,M,sizeof(float)*NumNodes);
}

void tledSolverGPU::PrintNodalForces()
{
   // Get forces from device
   const unsigned int MemSzU = sizeof(float4)*NumNodes;
   float4* l_F = (float4*)malloc(MemSzU);
   cudaMemcpy(l_F,d_F,MemSzU,cudaMemcpyDeviceToHost);
   gpuV4ToV3(l_F, NumNodes, h_F);
   
   free(l_F);

   float FX=0,FY=0,FZ=0;
   printf("-----------------------------------------------------\n");
   printf("Nodal Forces (ALL)\n");
   printf("Node\tFX\t\tFY\t\tFZ\n");
   for (int i = 0; i < NumNodes; i++)
   {
      printf("%i\t%+10.5e\t%+10.5e\t%+10.5e\n",i+1,h_F[i].x,h_F[i].y,h_F[i].z);
      FX += h_F[i].x;
      FY += h_F[i].y;
      FZ += h_F[i].z;
   }
   printf("\nSums\t%+10.5e\t%+10.5e\t%+10.5e\n",FX,FY,FZ);
   printf("-----------------------------------------------------\n");
}

void tledSolverGPU::PrintDispNodalForces()
{
   // Get forces from device
   const unsigned int MemSzU = sizeof(float4)*NumNodes;
   float4* l_F = (float4*)malloc(MemSzU);
   cudaMemcpy(l_F,d_F,MemSzU,cudaMemcpyDeviceToHost);
   gpuV4ToV3(l_F, NumNodes, h_F);
   
   free(l_F);

   float FX=0,FY=0,FZ=0;
   printf("-----------------------------------------------------\n");
   printf("Nodal Forces (LOADED NODES)\n");
   printf("Node\tFX\t\tFY\t\tFZ\n");
   for (int i = 0; i < NumNodes; i++)
   {
      if ((h_Uload[i].x != 0)||(h_Uload[i].y != 0)||(h_Uload[i].z != 0))
      {
         printf("%i\t%+10.5e\t%+10.5e\t%+10.5e\n",i+1,h_F[i].x,h_F[i].y,h_F[i].z);
         FX += h_F[i].x;
         FY += h_F[i].y;
         FZ += h_F[i].z;
      }
   }
   printf("\nSums\t%+10.5e\t%+10.5e\t%+10.5e\n",FX,FY,FZ);
   printf("-----------------------------------------------------\n");
}

void tledSolverGPU::PrintDispNodalForceSums()
{
   // Get forces from device
   const unsigned int MemSzU = sizeof(float4)*NumNodes;
   float4* l_F = (float4*)malloc(MemSzU);
   cudaMemcpy(l_F,d_F,MemSzU,cudaMemcpyDeviceToHost);
   gpuV4ToV3(l_F, NumNodes, h_F);
   
   free(l_F);

   float FX=0,FY=0,FZ=0;
   printf("-----------------------------------------------------\n");
   printf("Nodal Force Sums (LOADED NODES)\n");
   printf("Node\tFX\t\tFY\t\tFZ\n");
   for (int i = 0; i < NumNodes; i++)
   {
      if ((h_Uload[i].x != 0)||(h_Uload[i].y != 0)||(h_Uload[i].z != 0))
      {
         FX += h_F[i].x;
         FY += h_F[i].y;
         FZ += h_F[i].z;
      }
   }
   printf("\nSums\t%+10.5e\t%+10.5e\t%+10.5e\n",FX,FY,FZ);
   printf("-----------------------------------------------------\n");
}

void tledSolverGPU::PrintNodalDisps()
{
   // Get disps from device
  const float *uCurr = GetAllDisps();

   printf("-----------------------------------------------------\n");
   printf("Nodal Displacements (ALL)\n");
   printf("Node\tUX\t\tUY\t\tUZ\n");
   for (int i = 0; i < NumNodes; i++) {
     printf("%i\t%+10.5e\t%+10.5e\t%+10.5e\n",i+1, uCurr[3*i], uCurr[3*i+1], uCurr[3*i+2]);
   }
   printf("-----------------------------------------------------\n");
}

void tledSolverGPU::PrintNodalForceSums()
{
}

void tledSolverGPU::InitialiseSolutionVariables(void)
{
   memset(h_F,0,sizeof(float)*NumNodes*3);
      
   cudaMemset(d_F,0,sizeof(float4)*NumNodes);

   cudaMemset(d_Divergence,0,sizeof(bool));
}

void tledSolverGPU::InitialiseConstraints(void)
{
  tledCUDAHostMemoryBlock &r_fixMask = tledCUDAHostMemoryBlock::GetNextFreeBufferWithSize<int4>(NumNodes);

   memset(h_R,0,sizeof(float)*NumNodes*4);
   cudaMemset(d_R,0,sizeof(float4)*NumNodes);
   // Reset all FixMask vals to 1
   for (int i = 0; i < NumNodes; i++)
   {
      r_fixMask.GetBuffer<int4>()[i].x = 1;
      r_fixMask.GetBuffer<int4>()[i].y = 1;
      r_fixMask.GetBuffer<int4>()[i].z = 1;
      r_fixMask.GetBuffer<int4>()[i].w = 1;	// NB: last value not used
   }
   tledCUDAHelpers::CopyToDevice(d_FixMask, r_fixMask.GetBuffer<int4>(), NumNodes);
   //Reset DispMask
   memset(h_Uload,0,sizeof(float)*NumNodes*4);
   tledCUDAHelpers::SetZero(d_DispMask, NumNodes);
}

void tledSolverGPU::SetTimeStep(double dt)
{
   double Dt_old = Dt;
   Dt = dt;
   // Update CD coeffs
   InstantiateTimeStepper();
   const unsigned int MemSz = sizeof(float4)*NumNodes;
   cudaMemcpy(d_CD,h_CD,MemSz,cudaMemcpyHostToDevice);
   // Update visc integration params
   for (int el = 0; el < NumEls; el++)
   {
      int N = h_NumProny[el].x;
      for (int i = 0; i < N; i++) // Iso terms
      {
         float A = h_bkwdEulerIso[maxNumViscTerms.x*el + i].x;
         float B = h_bkwdEulerIso[maxNumViscTerms.x*el + i].y;
         float t = float(Dt_old/(1/B-1));
         float a = float(A*(Dt_old+t)/Dt_old);
         h_bkwdEulerIso[maxNumViscTerms.x*el + i].x = float(Dt*a/(Dt+t)); // New vals
         h_bkwdEulerIso[maxNumViscTerms.x*el + i].y = float(t/(Dt+t));
      }
      N = h_NumProny[el].y;
      for (int i = 0; i < N; i++) // Vol terms
      {
         float A = h_bkwdEulerVol[maxNumViscTerms.y*el + i].x;
         float B = h_bkwdEulerVol[maxNumViscTerms.y*el + i].y;
         float t = float(Dt_old/(1/B-1));
         float a = float(A*(Dt_old+t)/Dt_old);
         h_bkwdEulerVol[maxNumViscTerms.y*el + i].x = float(Dt*a/(Dt+t)); // New vals
         h_bkwdEulerVol[maxNumViscTerms.y*el + i].y = float(t/(Dt+t));
      }
   }
   // NB: the coefficient values are changed, but not the number of them, i.e. the number of Prony terms is unchanged
   const unsigned int MemSzEulerI = sizeof(float2)*NumEls*maxNumViscTerms.x;
   const unsigned int MemSzEulerV = sizeof(float2)*NumEls*maxNumViscTerms.y;
   cudaMemcpy(d_bkwdEulerIso,h_bkwdEulerIso,MemSzEulerI,cudaMemcpyHostToDevice);
   cudaMemcpy(d_bkwdEulerVol,h_bkwdEulerVol,MemSzEulerV,cudaMemcpyHostToDevice);
}

float tledSolverGPU::GetStrainEnergy(void)
{   
   const unsigned int MemSzSE = sizeof(float)*NumEls;
   
   float* d_ElStrainEnergy;
   cudaMalloc((void**)&d_ElStrainEnergy,MemSzSE);
   
   // Set up execution parameters
   const unsigned int BlockSz = 128;
   const unsigned int GridSz = (int)ceil((double)NumEls/(double)BlockSz);	// Number of threads for each pass
   dim3 Threads(BlockSz,1,1);
   dim3 Grid(GridSz,1,1);
   
   // Execute
   ComputeElementStrainEnergy_kernel <<< Grid,Threads >>> (d_ElStrainEnergy,d_SPKa,d_SPKb,d_Ca,d_Cb);
   
   float* l_ElStrainEnergy = (float*)malloc(MemSzSE);
   cudaMemcpy(l_ElStrainEnergy,d_ElStrainEnergy,MemSzSE,cudaMemcpyDeviceToHost);
   float e = 0.0f;
   for (int i = 0; i < NumEls; i++)
      e += l_ElStrainEnergy[i];
   
   cudaFree(d_ElStrainEnergy);
   free(l_ElStrainEnergy);
   
   return e;
}

float tledSolverGPU::GetKineticEnergy() {  
  std::vector<float> deltaUs(3*NumNodes);
  std::vector<float>::const_iterator ic_du;
  float const *pc_m;
  float e = 0;
      
  mp_TimeStepper->GetCurrentDeltaDisplacements(&deltaUs.front());
  for (ic_du = deltaUs.begin(), pc_m = M; ic_du < deltaUs.end(); pc_m++) {
    e += *ic_du**ic_du**pc_m, ic_du++;
    e += *ic_du**ic_du**pc_m, ic_du++;
    e += *ic_du**ic_du**pc_m, ic_du++;
  }

  return e/float(2*Dt*Dt);
}

void tledSolverGPU::PrepareOutput() {
  mp_TimeStepper->RetrieveSolutionFromDevice();
}

void tledSolverGPU::SetAllDisps(float* U) {
  mp_TimeStepper->SetCurrentDisplacements(U);
}

void tledSolverGPU::ComputeStresses(void)
{
   // Set up execution parameters
   const unsigned int BlockSzP1 = 128;
   const unsigned int BlockSzP2 = 64;
   const unsigned int GridSzP1 = (int)ceil((double)NumEls/(double)BlockSzP1);	// Number of threads for each pass
   const unsigned int GridSzP2 = (int)ceil((double)NumNodes/(double)BlockSzP2);
   dim3 ThreadsP1(BlockSzP1,1,1);
   dim3 GridP1(GridSzP1,1,1);
   dim3 ThreadsP2(BlockSzP2,1,1);
   dim3 GridP2(GridSzP2,1,1);
   
   if (!strcmp(EType,"T4"))
   {
      ComputeElementStressT4_kernel <<< GridP1,ThreadsP1 >>> (d_StressStateIso,d_StressStateVol,d_SPKa,d_SPKb,d_Ca,d_Cb);
   }
   else if (!strcmp(EType,"H8"))
   {
      ComputeElementStressH8_kernel <<< GridP1,ThreadsP1 >>> (d_StressStateIso,d_StressStateVol,d_SPKa,d_SPKb,d_Ca,d_Cb);
   }
   else // T4ANP
   {
      ComputeElementPressuresT4ANP_kernel <<< GridP1,ThreadsP1 >>> (d_FEl);
      ComputeNodalPressuresT4ANP_kernel <<< GridP2,ThreadsP2 >>> (d_Pa);
      ComputeModifiedElementStressT4ANP_kernel <<< GridP1,ThreadsP1 >>> (d_StressStateIso,d_StressStateVol,d_SPKa,d_SPKb,d_Ca,d_Cb);
   }
}

void tledSolverGPU::SetElementMatParams(int el, vector<float> params)
{
   UpdateMaterialParams(el,params);

   cudaMemcpy(&d_Vol_MType_K[el],&h_Vol_MType_K[el],sizeof(float4),cudaMemcpyHostToDevice);
   cudaMemcpy(&d_ElasticParams[numElasticParamsF4*el],&h_ElasticParams[numElasticParamsF4*el],numElasticParamsF4*sizeof(float4),cudaMemcpyHostToDevice);
   cudaMemcpy(&d_bkwdEulerIso[maxNumViscTerms.x*el],&h_bkwdEulerIso[maxNumViscTerms.x*el],maxNumViscTerms.x*sizeof(float2),cudaMemcpyHostToDevice);
   cudaMemcpy(&d_bkwdEulerVol[maxNumViscTerms.y*el],&h_bkwdEulerVol[maxNumViscTerms.y*el],maxNumViscTerms.y*sizeof(float2),cudaMemcpyHostToDevice);
   
   if (!strcmp(EType,"H8"))
   {
      UpdateHGParams(el);
      cudaMemcpy(&d_HG[16*el],&h_HG[16*el],sizeof(float4)*16,cudaMemcpyHostToDevice);
   }
}

void tledSolverGPU::SetMultipleElementMatParams(vector<int> el, vector<float>* params)
{
   for (int i = 0; i < (int)el.size(); i++)
      UpdateMaterialParams(el[i],params[i]);

   cudaMemcpy(d_Vol_MType_K,h_Vol_MType_K,sizeof(float4)*NumEls,cudaMemcpyHostToDevice);
   cudaMemcpy(d_ElasticParams,h_ElasticParams,sizeof(float4)*NumEls*numElasticParamsF4,cudaMemcpyHostToDevice);
   cudaMemcpy(d_bkwdEulerIso,h_bkwdEulerIso,sizeof(float2)*NumEls*maxNumViscTerms.x,cudaMemcpyHostToDevice);
   cudaMemcpy(d_bkwdEulerVol,h_bkwdEulerVol,sizeof(float2)*NumEls*maxNumViscTerms.y,cudaMemcpyHostToDevice);
   
   if (!strcmp(EType,"H8"))
   {
      for (int i = 0; i < (int)el.size(); i++)
         UpdateHGParams(el[i]);
      cudaMemcpy(d_HG,h_HG,sizeof(float4)*16*NumEls,cudaMemcpyHostToDevice);
   }
}

void tledSolverGPU::UpdateMaterialParams(int el, vector<float> params)
{
   int MType = (int)h_Vol_MType_K[el].y;
   float* inElasticParams = NULL;
   float4* outElasticParams = new float4[numElasticParamsF4];
   float* inViscParams = NULL;
   int Ni=0;
   int Nv=0;
   int2* N = NULL;
   float2* Ai = NULL;
   float2* Av = NULL;
   float2 wkHGLame;
   float ANPKappa;
   int paramcounter;
   switch (MType)
   {
      case 1: // LE
         inElasticParams = new float[2];
         for (int i=0; i<2; i++) {inElasticParams[i] = params[i];}
         ComputeLEparams(inElasticParams,outElasticParams,&wkHGLame,&ANPKappa);
         break;
      case 2: // NH
         inElasticParams = new float[2];
         for (int i=0; i<2; i++) {inElasticParams[i] = params[i];}
         ComputeNHparams(inElasticParams,outElasticParams,&wkHGLame,&ANPKappa);
         break;
      case 3: // TI
         inElasticParams = new float[6];
         for (int i=0; i<6; i++) {inElasticParams[i] = params[i];}
         ComputeTIparams(inElasticParams,outElasticParams,&wkHGLame,&ANPKappa);
         break;
      case 4: // NHV
         inElasticParams = new float[2];
         paramcounter = 0;
         for (int i=0; i<2; i++) {inElasticParams[i] = params[paramcounter++];}
         ComputeNHparams(inElasticParams,outElasticParams,&wkHGLame,&ANPKappa);
         paramcounter++; // Next value is Dt -> not needed
         Ni = (int)params[paramcounter++];
         Nv = (int)params[paramcounter++];
         N = new int2;
         Ai = new float2[Ni];
         Av = new float2[Nv];
         inViscParams = new float[2*(Ni+Nv)];
         for (int i=0; i<2*Ni; i++) {inViscParams[i] = params[paramcounter++];}
         for (int i=0; i<2*Nv; i++) {inViscParams[i] = params[paramcounter++];}
         ComputeViscparams(inViscParams,Ni,Nv,N,Ai,Av);
         break;
      case 5: // TIV
         inElasticParams = new float[6];
         paramcounter = 0;
         for (int i=0; i<6; i++) {inElasticParams[i] = params[paramcounter++];}
         ComputeTIparams(inElasticParams,outElasticParams,&wkHGLame,&ANPKappa);
         paramcounter++; // Next value is Dt -> not needed
         Ni = (int)params[paramcounter++];
         Nv = (int)params[paramcounter++];
         N = new int2;
         Ai = new float2[Ni];
         Av = new float2[Nv];
         inViscParams = new float[2*(Ni+Nv)];
         for (int i=0; i<2*Ni; i++) {inViscParams[i] = params[paramcounter++];}
         for (int i=0; i<2*Nv; i++) {inViscParams[i] = params[paramcounter++];}
         ComputeViscparams(inViscParams,Ni,Nv,N,Ai,Av);
         break;
      case 6: // AB
         inElasticParams = new float[3];
         for (int i=0; i<3; i++) {inElasticParams[i] = params[i];}
         ComputeABparams(inElasticParams,outElasticParams,&wkHGLame,&ANPKappa);
         break;
      case 7: // PY
         inElasticParams = new float[6];
         for (int i=0; i<6; i++) {inElasticParams[i] = params[i];}
         ComputePYparams(inElasticParams,outElasticParams,&wkHGLame,&ANPKappa);
         break;
      default:
         ;
   }
   // Update global material params lists
   for (int i = 0; i < numElasticParamsF4; i++)
      h_ElasticParams[el*numElasticParamsF4 + i] = outElasticParams[i];
   if ((Ni>0)|(Nv>0))
      h_NumProny[el] = *N;
   for (int i = 0; i < Ni; i++)
      h_bkwdEulerIso[el] = Ai[i];
   for (int i = 0; i < Nv; i++)
      h_bkwdEulerVol[el] = Av[i];
   h_Vol_MType_K[el].z = ANPKappa;
   HGLame[el] = wkHGLame;
   
   // Clean up
   delete inElasticParams;
   delete outElasticParams;
   if (inViscParams) delete inViscParams;
   if (N) delete N;
   if (Ai) delete Ai;
   if (Av) delete Av;
}

void tledSolverGPU::UpdateHGParams(int el)
{
   vector<int> EInd = Mesh->GetElNodeInds(el);
   vector<float> NCds;
   double x[8][3];
   for (int i = 0; i < 8; i++)
   {
      NCds = Mesh ->GetNodeCds(EInd[i]);
      for (int j = 0; j < 3; j++)
      {
         x[i][j] = (double)NCds[j];
      }
   }

   double fDhDx[8][3];
   for (int i = 0; i < 3; i++)
   {
      fDhDx[0][i] = (double)h_DhDx[6*el+i].x;
      fDhDx[1][i] = (double)h_DhDx[6*el+i].y;
      fDhDx[2][i] = (double)h_DhDx[6*el+i].z;
      fDhDx[3][i] = (double)h_DhDx[6*el+i].w;
      
      fDhDx[4][i] = (double)h_DhDx[6*el+i+3].x;
      fDhDx[5][i] = (double)h_DhDx[6*el+i+3].y;
      fDhDx[6][i] = (double)h_DhDx[6*el+i+3].z;
      fDhDx[7][i] = (double)h_DhDx[6*el+i+3].w;
   }
   double a = 0;
   for (int i = 0; i < 8; i++)
   {
      for (int j = 0; j < 3; j++)
      {
         a += fDhDx[i][j]*fDhDx[i][j];
      }
   }
   double k = HGKappa*ElementVols[el]*(HGLame[el].x+2*HGLame[el].y)*a/8;
   double Gamma[8][4] = {1,1,1,-1,
      -1,1,-1,1,
      1,-1,-1,-1,
      -1,-1,1,1,
      1,-1,-1,1,
      -1,-1,1,-1,
      1,1,1,1,
      -1,1,-1,-1};
   double A[8][8];
   MatMult8383T(fDhDx,x,A);
   double gamma[8][4];
   MatMult8884(A,Gamma,gamma);
   MatSubtract(&Gamma[0][0],&gamma[0][0],8,4,&gamma[0][0]);
   double fHG[8][8];
   MatMult8484T(gamma,gamma,fHG);
   MatMultScalar(&fHG[0][0],8,8,k,&fHG[0][0]);

   float4 b;
   for (int i = 0; i < 8; i++)
   {
      for (int j = 0; j < 2; j++)
      {
         b.x = (float)fHG[i][4*j];
         b.y = (float)fHG[i][4*j+1];
         b.z = (float)fHG[i][4*j+2];
         b.w = (float)fHG[i][4*j+3];
         h_HG[16*el + 2*i+j] = b;
      }
   }
}

void tledSolverGPU::InstantiateTimeStepper() {
  if (mp_TimeStepper != NULL) delete mp_TimeStepper;
  if (ANP) mp_TimeStepper = new tledNewmarkTimeStepperGPU(NumNodes, float(Dt), alpha, M);
  else mp_TimeStepper = new tledCentralDifferenceTimeStepperGPU(NumNodes, float(Dt), alpha, M);
}

void tledSolverGPU::SetGeometry(vector<float> NodeCoords)
{
   Mesh->SetNodeCoordinates(NodeCoords);
   float4* h_NodeCds = new float4[NumNodes];
   for (int i = 0; i < NumNodes; i++)
      h_NodeCds[i] = make_float4(NodeCoords[3*i],NodeCoords[3*i+1],NodeCoords[3*i+2],0);
   memset(M,0,sizeof(float)*NumNodes);
   if (!strcmp(EType,"T4"))
      UpdateElementT4Geometry();
   else if (!strcmp(EType,"T4ANP"))
   {
      UpdateElementT4Geometry();
      ComputeNodalVolumes();
   }
   else // 8-node hexahedra
      UpdateElementH8Geometry();
   
   InstantiateTimeStepper();
   
   // Send new data to GPU
   int NPE = Mesh->GetNodesPerEl();
   int factorANP = 0;
   if (ANP)
      factorANP = 1;
   int factorHG = 0;
   if (!strcmp(EType,"H8"))
      factorHG = 1;
   const unsigned int MemSzDhDx = sizeof(float4)*NumEls*NPE*3/4;	// NPE*3/4 = 3 for T4, and 6 for H8
   const unsigned int MemSzVol = sizeof(float4)*NumEls;
   const unsigned int MemSzCD = sizeof(float4)*NumNodes;
   const unsigned int MemSzHG = sizeof(float4)*NumEls*16*factorHG;
   const unsigned int MemSzVa = sizeof(float)*NumNodes*factorANP;
   
   cudaMemcpy(d_DhDx,h_DhDx,MemSzDhDx,cudaMemcpyHostToDevice);
   cudaMemcpy(d_CD,h_CD,MemSzCD,cudaMemcpyHostToDevice);
   cudaMemcpy(d_HG,h_HG,MemSzHG,cudaMemcpyHostToDevice);
   cudaMemcpy(d_Vol_MType_K,h_Vol_MType_K,MemSzVol,cudaMemcpyHostToDevice);
   cudaMemcpy(d_Va,h_Va,MemSzVa,cudaMemcpyHostToDevice);
   cudaMemcpy(d_NodeCds,h_NodeCds,MemSzCD,cudaMemcpyHostToDevice);

   delete h_NodeCds;
}

void tledSolverGPU::UpdateElementT4Geometry()
{
   vector<int> vEInd;
   double x[4][3];
   vector<float> NCds;
   double DhDr[4][3] = {-1, -1, -1,	// Shape function natural derivatives
                        1, 0, 0,
                        0, 1, 0,
                        0, 0, 1};
   double fDhDx[4][3];
   double J[3][3];	// Jacobian
   double detJ = 0;
   double invJ[3][3];
   double fVol = 0;
   double Mass = 0;
   
   for (int el = 0; el < NumEls; el++)
   {
      vEInd = Mesh->GetElNodeInds(el);
      for (int i = 0; i < 4; i++)
      {
         NCds = Mesh->GetNodeCds(vEInd[i]);
         for (int j = 0; j < 3; j++)
            x[i][j] = (double)NCds[j];
      }
      MatMult43T43(DhDr,x,J);
      // Compute element volume - note alternative, but equivalent formula
      MatDet33(J,&detJ);
      fVol = fabs(detJ/6);
      // Compute shape function global derivatives
      MatInv33(J,invJ);
      MatMult4333T(DhDr,invJ,fDhDx);
      // Compute element mass and add to nodal totals
      Mass = fVol*this->GetElementDensity(el)/4;
      for (int i = 0; i < 4; i++)
         M[vEInd[i]] += (float)Mass;
      
      // Collect final variables
      for (int i = 0; i < 3; i++)
      {
         h_DhDx[3*el+i].x = (float)fDhDx[0][i];
         h_DhDx[3*el+i].y = (float)fDhDx[1][i];
         h_DhDx[3*el+i].z = (float)fDhDx[2][i];
         h_DhDx[3*el+i].w = (float)fDhDx[3][i];
      }
      h_Vol_MType_K[el].x = (float)fVol;
      ElementVols[el] = (float)fVol;
   }
}

void tledSolverGPU::UpdateElementH8Geometry()
{
   vector<int> vEInd;
   double x[8][3];
   vector<float> NCds;
   double J[3][3];
   double detJ = 0;
   double fVol = 0;
   const double a = 1./8;
   double DhDr[8][3] = {-a, -a, -a,	// Shape function natural derivatives
                        a, -a, -a,
                        a, a, -a,
                        -a, a, -a,
                        -a, -a, a,
                        a, -a, a,
                        a, a, a,
                        -a, a, a};
   double fDhDx[8][3];
   double Mass = 0;
   float2 wkHGLame;
   
   h_HG = new float4[NumEls*16];
   double fHG[8][8];
   
   for (int el = 0; el < NumEls; el++)
   {
      vEInd = Mesh->GetElNodeInds(el);
      for (int i = 0; i < 8; i++)
      {
         NCds = Mesh->GetNodeCds(vEInd[i]);
         for (int j = 0; j < 3; j++)
            x[i][j] = (double)NCds[j];
      }

      // Jacobian
      MatMult83T83(DhDr,x,J);
      MatDet33(J,&detJ);
      // Element volume
      fVol = ComputeHexVol(x);
      // Compute shape function global derivatives
      ComputeH8DhDx(x,fVol,fDhDx);
      // Compute element mass and add to nodal totals
      Mass = double(this->GetElementDensity(el))*fVol/8;
      for (int i = 0; i < 8; i++)
         M[vEInd[i]] += (float)Mass;
      
      // Hourglass control parameters
      double a = 0;
      for (int i = 0; i < 8; i++)
      {
         for (int j = 0; j < 3; j++)
            a += fDhDx[i][j]*fDhDx[i][j];
      }
      wkHGLame = HGLame[el];
      double k = HGKappa*fVol*(wkHGLame.x+2*wkHGLame.y)*a/8;
      double Gamma[8][4] = {1,1,1,-1,
                           -1,1,-1,1,
                           1,-1,-1,-1,
                           -1,-1,1,1,
                           1,-1,-1,1,
                           -1,-1,1,-1,
                           1,1,1,1,
                           -1,1,-1,-1};
      double A[8][8];
      MatMult8383T(fDhDx,x,A);
      double gamma[8][4];
      MatMult8884(A,Gamma,gamma);
      MatSubtract(&Gamma[0][0],&gamma[0][0],8,4,&gamma[0][0]);
      MatMult8484T(gamma,gamma,fHG);
      MatMultScalar(&fHG[0][0],8,8,k,&fHG[0][0]);
      
      // Collect final variables
      for (int i = 0; i < 3; i++)
      {
         h_DhDx[6*el+i].x = (float)fDhDx[0][i];
         h_DhDx[6*el+i].y = (float)fDhDx[1][i];
         h_DhDx[6*el+i].z = (float)fDhDx[2][i];
         h_DhDx[6*el+i].w = (float)fDhDx[3][i];

         h_DhDx[6*el+i+3].x = (float)fDhDx[4][i];
         h_DhDx[6*el+i+3].y = (float)fDhDx[5][i];
         h_DhDx[6*el+i+3].z = (float)fDhDx[6][i];
         h_DhDx[6*el+i+3].w = (float)fDhDx[7][i];
      }
      h_Vol_MType_K[el].x = (float)(8*detJ);
      ElementVols[el] = (float)fVol;
      float4 b;
      for (int i = 0; i < 8; i++)
      {
         for (int j = 0; j < 2; j++)
         {
            b.x = (float)fHG[i][4*j];
            b.y = (float)fHG[i][4*j+1];
            b.z = (float)fHG[i][4*j+2];
            b.w = (float)fHG[i][4*j+3];
            h_HG[16*el + 2*i+j] = b;
         }
      }
   }
}



#endif // _GPU_
