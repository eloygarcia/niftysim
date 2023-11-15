// =========================================================================
// File:       tledSolverGPU_ROM.cu
// Purpose:    Main finite element object - GPU version using reduced order modelling
// Package:    NiftySim: Nonlinear finite element program
// Author:     Zeike Taylor
// Language:   C++
// Created:    December 2009 (converted to cu file, April 2011)
// 
// Copyright (c) 2010, University of Queensland. All rights reserved.
// MedTeQ Centre
// See the LICENSE.txt file in the root folder
// 
// ztaylor@itee.uq.edu.au
// =========================================================================

#ifdef _GPU_

#include "tledSolverGPU_ROM.h"
#include "tledTimer.h"
#include "tledSolverGPU_ROM_kernels.cu"
#include "tledDeviceDeclarations.h"

#include <limits>
#include <math.h>

using namespace std;
using namespace tledSolverGPU_ROM_kernels;

tledSolverGPU_ROM::tledSolverGPU_ROM()
{
}

tledSolverGPU_ROM::~tledSolverGPU_ROM()
{
   // Host variables
   delete[] h_DhDx;
   delete[] h_EInd;
   delete[] h_Ucurr;
   delete[] h_Fint;
   delete[] h_Mf;
   delete[] h_NodeMap;
   delete[] h_FixMask;
   delete[] h_FCds;
   if(!h_HG) delete[] h_HG;
   delete[] h_Vol_MType_K;
   delete[] h_ElasticParams;
   delete[] h_NumProny;
   delete[] h_bkwdEulerIso;
   delete[] h_bkwdEulerVol;
   delete[] h_SPKa;
   delete[] h_SPKb;
   delete[] h_Ca;
   delete[] h_Cb;
   delete[] h_Pa;
   delete[] h_Va;
   delete[] HGLame;
   delete[] h_Mr;
   delete[] h_Phi;
   delete[] h_fr;
   delete[] h_tmp;
   delete[] l_U;
   delete[] UOutput;
   
   // Pinned host variables
   cudaFreeHost(h_Uload);
   cudaFreeHost(h_Fext);
   cudaFreeHost(h_DispMask);
   cudaFreeHost(h_f);
   cudaFreeHost(h_UnextV);
   
   // Device variables
   cudaFree(d_DhDx);
   cudaFree(d_EInd);
   cudaFree(d_HG);
   cudaFree(d_Fint);
   cudaFree(d_FEl);
   cudaFree(d_Uprev);
   cudaFree(d_Ucurr);
   cudaFree(d_Unext);
   cudaFree(d_UnextV);
   cudaFree(d_Uload);
   cudaFree(d_Fext);
   cudaFree(d_FCds);
   cudaFree(d_NodeCds);
   cudaFree(d_NodeMap);
   cudaFree(d_FixMask);
   cudaFree(d_DispMask);
   cudaFree(d_Vol_MType_K);
   cudaFree(d_ElasticParams);
   cudaFree(d_NumProny);
   cudaFree(d_bkwdEulerIso);
   cudaFree(d_bkwdEulerVol);
   cudaFree(d_StressStateIso);
   cudaFree(d_StressStateVol);
   cudaFree(d_SPKa);
   cudaFree(d_SPKb);
   cudaFree(d_Ca);
   cudaFree(d_Cb);
   cudaFree(d_Divergence);
   cudaFree(d_Pa);
   cudaFree(d_Va);
   cudaFree(d_f);
}

void tledSolverGPU_ROM::Init(tledModel* Model)
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
   D = Model->GetDensity();
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
   numBasisVecs = Model->GetNumBasisVectors();
   
   // Allocate host variables
   h_DhDx = new float4[NumEls*NPE*3/4]; memset(h_DhDx,0,sizeof(float)*NumEls*NPE*3);
   h_EInd = new int4[NumEls*NPE/4]; memset(h_EInd,0,sizeof(int)*NumEls*NPE); // NPE/4 = 1 for T4, and 2 for H8
   h_Ucurr = new float3[NumNodes]; memset(h_Ucurr,0,sizeof(float)*NumNodes*3);
   h_Fint = new float3[NumNodes]; memset(h_Fint,0,sizeof(float)*NumNodes*3);
   h_Mf = new float[NumNodes*3]; memset(h_Mf,0,sizeof(float)*NumNodes*3);
   h_NodeMap = new int2[NumNodes]; memset(h_NodeMap,0,sizeof(int)*NumNodes*2);
   h_FixMask = new int4[NumNodes]; memset(h_FixMask,0,sizeof(int)*NumNodes*4);
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
   h_Va = NULL;
   h_Pa = NULL;
   h_Mr = new float[numBasisVecs*numBasisVecs];
   h_fr = new float[numBasisVecs];
   h_tmp = new float[numBasisVecs];

   l_U = new float4[NumNodes];
   UOutput = new float[3*NumNodes];
   
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
   
   // Compute node variables
   int* EPN = new int[NumNodes]; // Elements per node (different for each node)
   FCdsLength = 0;
   for (int node = 0; node < NumNodes; node++)
   {
      // Arrays for contructing NodeMap and FCds
      EPN[node] = (int)(*NInd)[node].size();
      FCdsLength += EPN[node];
   }
   h_FCds = new int2[FCdsLength];
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
   
   // Compute central difference coeffs
   ComputeCDCoeffs();
   
   // Get reduced basis
   h_Phi = Model->GetReducedBasis(); // Get row-major array of basis vectors
   // Compute Mr = Phi'*M*Phi
   float* tmp = new float[NumNodes*3*numBasisVecs];
   MatMultdAB(h_Mf,3*NumNodes,h_Phi,NumNodes*3,numBasisVecs,tmp); // tmp = M*Phi
   MatMultAtB(h_Phi,3*NumNodes,numBasisVecs,tmp,NumNodes*3,numBasisVecs,h_Mr); //Mr = Phi'*tmp (= Phi'*M*Phi)
   delete tmp;
   // Invert Mr
   MatInverse(h_Mr,numBasisVecs);
   
   // Make compacted full mass vector
   float* Mfc = new float[NumNodes];
   for (int i = 0; i < NumNodes; i++)
      Mfc[i] = h_Mf[3*i];
   
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
   const unsigned int MemSzU = sizeof(float4)*NumNodes;
   const unsigned int MemSzNodeMap = sizeof(int2)*NumNodes;
   const unsigned int MemSzFCds = sizeof(int2)*FCdsLength;
   const unsigned int MemSzHG = sizeof(float4)*NumEls*16*factorHG;
   const unsigned int MemSzMask = sizeof(int4)*NumNodes;
   const unsigned int MemSzVol = sizeof(float4)*NumEls;
   const unsigned int MemSzElast = sizeof(float4)*NumEls*numElasticParamsF4;
   const unsigned int MemSzProny = sizeof(int2)*NumEls;
   const unsigned int MemSzEulerI = sizeof(float2)*NumEls*maxNumViscTerms.x;
   const unsigned int MemSzEulerV = sizeof(float2)*NumEls*maxNumViscTerms.y;
   const unsigned int MemSzStateI = sizeof(float4)*NumEls*2*maxNumViscTerms.x;	// Will produce 8 stresses for each Prony term (6 req'd)
   const unsigned int MemSzStateV = sizeof(float4)*NumEls*2*maxNumViscTerms.y;
   const unsigned int MemSzSPK = sizeof(float4)*NumEls;
   const unsigned int MemSzDiv = sizeof(bool);
   const unsigned int MemSzPa = sizeof(float)*NumNodes*factorANP;
   const unsigned int MemSzf = sizeof(float)*NumNodes*3;
   const unsigned int MemSzM = sizeof(float)*NumNodes;
   
   // Allocate pinned host memory
   cudaMallocHost((void**)&h_Uload,MemSzU);  // Imposed displacements
   memset(h_Uload,0,MemSzU);
   cudaMallocHost((void**)&h_DispMask,MemSzMask);  // Mask of displaced nodes
   memset(h_DispMask,0,MemSzMask);
   cudaMallocHost((void**)&h_Fext,MemSzU);   // Imposed forces
   memset(h_Fext,0,MemSzU);
   cudaMallocHost((void**)&h_f,MemSzf);   // Effective forces
   memset(h_f,0,MemSzf);
   cudaMallocHost((void**)&h_UnextV,MemSzf);   // Vector form of Unext
   memset(h_UnextV,0,MemSzf);
   
   // Allocate device memory
   cudaMalloc((void**)&d_DhDx,MemSzDhDx);
   cudaMalloc((void**)&d_EInd,MemSzEInd);
   cudaMalloc((void**)&d_Fint,MemSzU);
   cudaMalloc((void**)&d_Ucurr,MemSzU);
   cudaMalloc((void**)&d_Uprev,MemSzU);
   cudaMalloc((void**)&d_Unext,MemSzU);
   cudaMalloc((void**)&d_UnextV,MemSzf);
   cudaMalloc((void**)&d_FEl,MemSzFEl);
   cudaMalloc((void**)&d_NodeMap,MemSzNodeMap);
   cudaMalloc((void**)&d_FCds,MemSzFCds);
   cudaMalloc((void**)&d_FixMask,MemSzMask);
   cudaMalloc((void**)&d_HG,MemSzHG);
   cudaMalloc((void**)&d_Uload,MemSzU);
   cudaMalloc((void**)&d_DispMask,MemSzMask);
   cudaMalloc((void**)&d_Fext,MemSzU);
   cudaMalloc((void**)&d_Vol_MType_K,MemSzVol);
   cudaMalloc((void**)&d_ElasticParams,MemSzElast);
   cudaMalloc((void**)&d_NumProny,MemSzProny);
   cudaMalloc((void**)&d_bkwdEulerIso,MemSzEulerI);
   cudaMalloc((void**)&d_bkwdEulerVol,MemSzEulerV);
   cudaMalloc((void**)&d_StressStateIso,MemSzStateI);
   cudaMalloc((void**)&d_StressStateVol,MemSzStateV);
   cudaMalloc((void**)&d_SPKa,MemSzSPK);
   cudaMalloc((void**)&d_SPKb,MemSzSPK);
   cudaMalloc((void**)&d_Ca,MemSzSPK);
   cudaMalloc((void**)&d_Cb,MemSzSPK);
   cudaMalloc((void**)&d_Divergence,MemSzDiv);
   cudaMalloc((void**)&d_Pa,MemSzPa);
   cudaMalloc((void**)&d_Va,MemSzPa);
   cudaMalloc((void**)&d_f,MemSzf);
   cudaMalloc((void**)&d_Mfc,MemSzM);
   cudaMalloc((void**)&d_NodeCds,MemSzU);
   
   // Copy host data to device
   cudaMemcpy(d_DhDx,h_DhDx,MemSzDhDx,cudaMemcpyHostToDevice);
   cudaMemcpy(d_EInd,h_EInd,MemSzEInd,cudaMemcpyHostToDevice);
   cudaMemcpy(d_NodeMap,h_NodeMap,MemSzNodeMap,cudaMemcpyHostToDevice);
   cudaMemcpy(d_FCds,h_FCds,MemSzFCds,cudaMemcpyHostToDevice);
   cudaMemcpy(d_HG,h_HG,MemSzHG,cudaMemcpyHostToDevice);
   cudaMemcpy(d_Vol_MType_K,h_Vol_MType_K,MemSzVol,cudaMemcpyHostToDevice);
   cudaMemcpy(d_ElasticParams,h_ElasticParams,MemSzElast,cudaMemcpyHostToDevice);
   cudaMemcpy(d_NumProny,h_NumProny,MemSzProny,cudaMemcpyHostToDevice);
   cudaMemcpy(d_bkwdEulerIso,h_bkwdEulerIso,MemSzEulerI,cudaMemcpyHostToDevice);
   cudaMemcpy(d_bkwdEulerVol,h_bkwdEulerVol,MemSzEulerV,cudaMemcpyHostToDevice);
   cudaMemcpy(d_Pa,h_Pa,MemSzPa,cudaMemcpyHostToDevice);
   cudaMemcpy(d_Va,h_Va,MemSzPa,cudaMemcpyHostToDevice);
   cudaMemcpy(d_Mfc,Mfc,MemSzM,cudaMemcpyHostToDevice);
   cudaMemcpy(d_NodeCds,h_NodeCds,MemSzU,cudaMemcpyHostToDevice);
   cudaMemcpyToSymbol(c_NumNodes,&NumNodes,sizeof(int));
   cudaMemcpyToSymbol(c_NumEls,&NumEls,sizeof(int));
   cudaMemcpyToSymbol(c_NPE,&NPE,sizeof(int));
   cudaMemcpyToSymbol(c_maxNumViscTerms,&maxNumViscTerms,sizeof(int2));
   cudaMemcpyToSymbol(c_gamma,&CDgamma,sizeof(float3));
   
   // Initialise some variables
   cudaMemset(d_Ucurr,0,MemSzU);
   cudaMemset(d_Uprev,0,MemSzU);
   cudaMemset(d_Unext,0,MemSzU);
   cudaMemset(d_StressStateIso,0,MemSzStateI);
   cudaMemset(d_StressStateVol,0,MemSzStateV);
   cudaMemset(d_Divergence,0,MemSzDiv);
   
   // Bind textures
   cudaBindTexture(0,txDhDx,d_DhDx,MemSzDhDx);
   cudaBindTexture(0,txEInd,d_EInd,MemSzEInd);
   cudaBindTexture(0,txHG,d_HG,MemSzHG);
   cudaBindTexture(0,txNodeMap,d_NodeMap,MemSzNodeMap);
   cudaBindTexture(0,txFCds,d_FCds,MemSzFCds);
   cudaBindTexture(0,txFEl,d_FEl,MemSzFEl);
   cudaBindTexture(0,txFixMask,d_FixMask,MemSzMask);
   cudaBindTexture(0,txUload,d_Uload,MemSzU);
   cudaBindTexture(0,txDispMask,d_DispMask,MemSzMask);
   cudaBindTexture(0,txF4Fext,d_Fext,MemSzU);
   cudaBindTexture(0,txVol_MType_K,d_Vol_MType_K,MemSzVol);
   cudaBindTexture(0,txElasticParams,d_ElasticParams,MemSzElast);
   cudaBindTexture(0,txNumProny,d_NumProny,MemSzProny);
   cudaBindTexture(0,txEulerIso,d_bkwdEulerIso,MemSzEulerI);
   cudaBindTexture(0,txEulerVol,d_bkwdEulerVol,MemSzEulerV);
   cudaBindTexture(0,txStateIso,d_StressStateIso,MemSzStateI);
   cudaBindTexture(0,txStateVol,d_StressStateVol,MemSzStateV);
   cudaBindTexture(0,txPa,d_Pa,MemSzPa);
   cudaBindTexture(0,txVa,d_Va,MemSzPa);
   cudaBindTexture(0,txUnextV,d_UnextV,MemSzf);
   cudaBindTexture(0,txMfc,d_Mfc,MemSzM);
   cudaBindTexture(0,txNodeCds,d_NodeCds,MemSzU);
   
   delete NInd;
   delete EPN;
   delete Mfc;
   delete h_NodeCds;
}

void tledSolverGPU_ROM::ComputeElementT4Variables(tledModel* Model, vector<vector<int2> >* NInd)
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
               h_ElasticParams[el*numElasticParamsF4 + i] = outElasticParams[i];
            h_NumProny[el] = *N;
            for (int i = 0; i < currNumViscIsoTerms; i++)
               h_bkwdEulerIso[el] = Ai[i];
            for (int i = 0; i < currNumViscVolTerms; i++)
               h_bkwdEulerVol[el] = Av[i];

            // Assemble node coords
            vEInd = Model->GetElNodeInds(el);
            for (int i = 0; i < 4; i++)
            {
               NCds = Model->GetNodeCds(vEInd[i]);
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
            Mass = fVol*D/4;
            for (int i = 0; i < 4; i++)
            {
               h_Mf[vEInd[i]*3] += (float)Mass;
               h_Mf[vEInd[i]*3+1] += (float)Mass;
               h_Mf[vEInd[i]*3+2] += (float)Mass;
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
            HGLame[el] = wkHGLame;

            ElChecker[el] = 1;
         }
      }
      delete outElasticParams;
      delete N;
      delete Ai;
      delete Av;
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

   delete ElChecker;
}

void tledSolverGPU_ROM::ComputeElementT4ANPVariables(tledModel* Model, vector<vector<int2> >* NInd)
{
   ComputeElementT4Variables(Model,NInd);
   ANP = true;
   h_Pa = new float[NumNodes]; memset(h_Pa,0,sizeof(float)*NumNodes);
   ComputeNodalVolumes();
}

void tledSolverGPU_ROM::ComputeElementH8Variables(tledModel* Model, vector<vector<int2> >* NInd)
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
               h_ElasticParams[el*numElasticParamsF4 + i] = outElasticParams[i];
            h_NumProny[el] = *N;
            for (int i = 0; i < currNumViscIsoTerms; i++)
               h_bkwdEulerIso[el] = Ai[i];
            for (int i = 0; i < currNumViscVolTerms; i++)
               h_bkwdEulerVol[el] = Av[i];

            // Assemble node coords
            vEInd = Model->GetElNodeInds(el);
            for (int i = 0; i < 8; i++)
            {
               NCds = Model->GetNodeCds(vEInd[i]);
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
            Mass = ((double)D)*fVol/8;
            for (int i = 0; i < 8; i++)
            {
               h_Mf[vEInd[i]*3] += (float)Mass;
               h_Mf[vEInd[i]*3+1] += (float)Mass;
               h_Mf[vEInd[i]*3+2] += (float)Mass;
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
                  a += fDhDx[i][j]*fDhDx[i][j];
            }
            double k = Model->GetHGKappa()*fVol*(wkHGLame.x+2*wkHGLame.y)*a/8;
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
      delete outElasticParams;
      delete N;
      delete Ai;
      delete Av;
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

   delete ElChecker;
}

void tledSolverGPU_ROM::ComputeNodalVolumes()
{
   h_Va = new float[NumNodes]; memset(h_Va,0,sizeof(float)*NumNodes);
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

double tledSolverGPU_ROM::ComputeHexVol(double x[8][3])
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
            V += x[I][0]*x[J][1]*x[K][2]*C[I][J][K];
      }
   }
   return V;
}

void tledSolverGPU_ROM::ComputeH8DhDx(double x[8][3], double V, double fDhDx[8][3])
{
   // Calc B matrix
   double B[8][3];
   ComputeBmat(x,B);
   for (int i = 0; i < 8; i++)
   {
      for (int j = 0; j < 3; j++)
         fDhDx[i][j] = B[i][j]/V;
   }
}

void tledSolverGPU_ROM::ComputeBmat(double x[8][3], double B[8][3])
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

void tledSolverGPU_ROM::ComputeCIJK(double C[8][8][8])
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

void tledSolverGPU_ROM::ComputeLEparams(float* inMatParams, float4* outMatParams, float2* wkHGLame, float* ANPKappa)
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

void tledSolverGPU_ROM::ComputeNHparams(float* inMatParams, float4* outMatParams, float2* wkHGLame, float* ANPKappa)
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

void tledSolverGPU_ROM::ComputeTIparams(float* inMatParams, float4* outMatParams, float2* wkHGLame, float* ANPKappa)
{
   outMatParams[0].x = inMatParams[0]; // Mu
   outMatParams[0].y = inMatParams[1]; // K
   outMatParams[0].z = inMatParams[2]; // Eta
   outMatParams[0].w = inMatParams[3]; // a0
   outMatParams[1].x = inMatParams[4]; // a1
   outMatParams[1].y = inMatParams[5]; // a2
   outMatParams[1].z = 0;  // unused
   outMatParams[1].w = 0;  // unused

   // Compute HG Lame params
   (*wkHGLame).x = outMatParams[0].y - 2*(outMatParams[0].x)/3;
   (*wkHGLame).y = outMatParams[0].x;
   // Compute ANPKappa
   *ANPKappa = inMatParams[1];
}

void tledSolverGPU_ROM::ComputeViscparams(float* inMatViscParams, int Ni, int Nv, int2* N, float2* Ai, float2* Av)
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

void tledSolverGPU_ROM::PerformStep()
{
   const unsigned int MemSzU = sizeof(float4)*NumNodes;
   const unsigned int MemSzMask = sizeof(int4)*NumNodes;
   const unsigned int MemSzf = sizeof(float)*NumNodes*3;
   
   // Load new displacements
   cudaMemcpy(d_Uload,h_Uload,MemSzU,cudaMemcpyHostToDevice);
   cudaMemcpy(d_DispMask,h_DispMask,MemSzMask,cudaMemcpyHostToDevice);
   
   // Load external forces
   cudaMemcpy(d_Fext,h_Fext,MemSzU,cudaMemcpyHostToDevice);
   
   // Update displacements
   float4* d_Utemp;
   d_Utemp = d_Uprev;
   d_Uprev = d_Ucurr;
   d_Ucurr = d_Unext;
   d_Unext = d_Utemp;
   cudaBindTexture(0,txUprev,d_Uprev,MemSzU);
   cudaBindTexture(0,txUcurr,d_Ucurr,MemSzU);
   
   // Set up execution parameters
   const unsigned int BlockSzP1 = 128;
   const unsigned int BlockSzP2 = 16;
   const unsigned int BlockSzP3 = 16;
   const unsigned int GridSzP1 = (int)ceil((double)NumEls/(double)BlockSzP1);	// Number of threads for each pass
   const unsigned int GridSzP2 = (int)ceil((double)NumNodes/(double)BlockSzP2);
   const unsigned int GridSzP3 = (int)ceil((double)NumNodes/(double)BlockSzP3);
   dim3 ThreadsP1(BlockSzP1,1,1);
   dim3 GridP1(GridSzP1,1,1);
   dim3 ThreadsP2(BlockSzP2,1,1);
   dim3 GridP2(GridSzP2,1,1);
   dim3 ThreadsP3(BlockSzP3,1,1);
   dim3 GridP3(GridSzP3,1,1);

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

   ComputeEffectiveLoads_kernel <<< GridP2,ThreadsP2 >>> (d_Fint,d_f,d_Contacts);
   
   // fr = Phi'*f
   cudaMemcpy(h_f,d_f,MemSzf,cudaMemcpyDeviceToHost);
   MatMultAtB(h_Phi,NumNodes*3,numBasisVecs,h_f,NumNodes*3,1,h_tmp); // tmp = Phi'*f
   
   // fr = Mr*fr
   MatMultAB(h_Mr,numBasisVecs,numBasisVecs,h_tmp,numBasisVecs,1,h_fr); //fr = Mr*tmp
   
   // Unext = Phi*fr
   MatMultAB(h_Phi,NumNodes*3,numBasisVecs,h_fr,numBasisVecs,1,h_UnextV);
   cudaMemcpy(d_UnextV,h_UnextV,MemSzf,cudaMemcpyHostToDevice);
   
   // Unext = g1*Unext + g2*Ucurr + g3*Uprev (and set BCs)
   ComputeNewDisps_kernel <<< GridP3,ThreadsP3 >>> (d_UnextV,d_Unext,d_Divergence);
}

void tledSolverGPU_ROM::GetDivergence(bool* Div)
{
   cudaMemcpy(Div,d_Divergence,sizeof(bool),cudaMemcpyDeviceToHost);
}

void tledSolverGPU_ROM::UnsetDivergence() {
  bool div;

  div = false;
  tledCUDAHelpers::CopyToDevice(d_Divergence, &div);
}

void tledSolverGPU_ROM::SetDisps(vector<int>* IndX, vector<float>* UX, vector<int>* IndY, vector<float>* UY,
                             vector<int>* IndZ, vector<float>* UZ)
{
   // Reset all vals to 0
   memset(h_Uload,0,sizeof(float)*NumNodes*4);
   memset(h_DispMask,0,sizeof(int)*NumNodes*4);
   // Set new disps
   for (int i = 0; i < (int)IndX->size(); i++)
   {
      h_Uload[(*IndX)[i]].x += (*UX)[i];
      h_DispMask[(*IndX)[i]].x = 1;
   }
   for (int i = 0; i < (int)IndY->size(); i++)
   {
      h_Uload[(*IndY)[i]].y += (*UY)[i];
      h_DispMask[(*IndY)[i]].y = 1;
   }
   for (int i = 0; i < (int)IndZ->size(); i++)
   {
      h_Uload[(*IndZ)[i]].z += (*UZ)[i];
      h_DispMask[(*IndZ)[i]].z = 1;
   }
}

void tledSolverGPU_ROM::SetExtForces(vector<int>* IndX, vector<float>* FX, vector<int>* IndY, vector<float>* FY,
                                 vector<int>* IndZ, vector<float>* FZ)
{
   // Reset all vals to 0
   memset(h_Fext,0,sizeof(float)*NumNodes*4);
   // Set new disps
   for (int i = 0; i < (int)IndX->size(); i++)
      h_Fext[(*IndX)[i]].x += (*FX)[i];
   for (int i = 0; i < (int)IndY->size(); i++)
      h_Fext[(*IndY)[i]].y += (*FY)[i];
   for (int i = 0; i < (int)IndZ->size(); i++)
      h_Fext[(*IndZ)[i]].z += (*FZ)[i];
}

float* tledSolverGPU_ROM::GetAllExtForces(float *p_dst) const {
  for (int i = 0; i < NumNodes; i++) tledCUDAHelpers::ConvertFromFloatN(p_dst + 3*i, h_Fext[i]);

  return p_dst;
}

void tledSolverGPU_ROM::SetFixed(vector<int>* IndX, vector<int>* IndY, vector<int>* IndZ)
{
   // Reset all vals to 1
   for (int i = 0; i < NumNodes; i++)
   {
      h_FixMask[i].x = 1;
      h_FixMask[i].y = 1;
      h_FixMask[i].z = 1;
      h_FixMask[i].w = 1;	// NB: last value not used
   }
   // Set new disp indices to 0
   for (int i = 0; i < (int)IndX->size(); i++)
      h_FixMask[(*IndX)[i]].x = 0;
   for (int i = 0; i < (int)IndY->size(); i++)
      h_FixMask[(*IndY)[i]].y = 0;
   for (int i = 0; i < (int)IndZ->size(); i++)
      h_FixMask[(*IndZ)[i]].z = 0;

   // Send new fix mask to device
   cudaMemcpy(d_FixMask,h_FixMask,sizeof(int4)*NumNodes,cudaMemcpyHostToDevice);
}

void tledSolverGPU_ROM::ComputeCDCoeffs()
{
  CDgamma.x = float(2*Dt*Dt/(alpha*Dt+2));
  CDgamma.y = float(4/(alpha*Dt+2));
  CDgamma.z = 1 - CDgamma.y;
}

void tledSolverGPU_ROM::SetContactManager(tledContactManager* contacts)
{   
   Contacts = contacts;
   int numContactCyls = Contacts->GetNumContactCyls();
   int numContactPrbs = Contacts->GetNumContactPrbs();
   int numContactPlts = Contacts->GetNumContactPlts();
   
   cudaMemcpyToSymbol(c_numContactCyls,&numContactCyls,sizeof(int));
   cudaMemcpyToSymbol(c_numContactPrbs,&numContactPrbs,sizeof(int));
   cudaMemcpyToSymbol(c_numContactPlts,&numContactPlts,sizeof(int));

   if (Contacts->GetNumContactObjects() > 0)
      d_Contacts = Contacts->GetContactsDevicePointer();
}

void tledSolverGPU_ROM::GetForces(vector<int>* NodeInd, vector<float>* Forces)
{
   // Get forces from device
   const unsigned int MemSzU = sizeof(float4)*NumNodes;
   float4* l_F = (float4*)malloc(MemSzU);
   cudaMemcpy(l_F,d_Fint,MemSzU,cudaMemcpyDeviceToHost);
   gpuV4ToV3(l_F, NumNodes, h_Fint);
   
   free(l_F);

   for (int i = 0; i < (int)(NodeInd->size()); i++)
   {
      Forces->push_back( h_Fint[3*(*NodeInd)[i]].x );
      Forces->push_back( h_Fint[3*(*NodeInd)[i]].y );
      Forces->push_back( h_Fint[3*(*NodeInd)[i]].z );
   }
}

float* tledSolverGPU_ROM::GetAllForces()
{
   // Get forces from device
   const unsigned int MemSzU = sizeof(float4)*NumNodes;
   float4* l_F = (float4*)malloc(MemSzU);
   cudaMemcpy(l_F,d_Fint,MemSzU,cudaMemcpyDeviceToHost);
   gpuV4ToV3(l_F, NumNodes, h_Fint);
   
   free(l_F);
   
   float* F = new float[NumNodes*3];
   for (int i = 0; i < NumNodes; i++)
   {
      F[3*i] = h_Fint[i].x;
      F[3*i+1] = h_Fint[i].y;
      F[3*i+2] = h_Fint[i].z;
   }
   return F;
}

void tledSolverGPU_ROM::GetDisps(vector<int>* NodeInd, vector<float>* Disps)
{
   // Get disps from device
   const unsigned int MemSzU = sizeof(float4)*NumNodes;
   float4* l_Ucurr = (float4*)malloc(MemSzU);
   cudaMemcpy(l_Ucurr,d_Ucurr,MemSzU,cudaMemcpyDeviceToHost);
   gpuV4ToV3(l_Ucurr, NumNodes, h_Ucurr);
   
   free(l_Ucurr);
   
   for (int i = 0; i < (int)(NodeInd->size()); i++)
   {
      Disps->push_back( h_Ucurr[3*(*NodeInd)[i]].x );
      Disps->push_back( h_Ucurr[3*(*NodeInd)[i]].y );
      Disps->push_back( h_Ucurr[3*(*NodeInd)[i]].z );
   }
}

float* tledSolverGPU_ROM::GetAllNextDisps() {
  cudaMemcpy(l_U, d_Unext, sizeof(float4)*NumNodes, cudaMemcpyDeviceToHost);

  for (int i = 0; i < NumNodes; i++) {
    UOutput[3*i] = l_U[i].x;
    UOutput[3*i+1] = l_U[i].y;
    UOutput[3*i+2] = l_U[i].z;
  }

  return UOutput;
}

float* tledSolverGPU_ROM::GetAllDisps() {
  cudaMemcpy(l_U, d_Ucurr, sizeof(float4)*NumNodes, cudaMemcpyDeviceToHost);

  for (int i = 0; i < NumNodes; i++) {
    UOutput[3*i] = l_U[i].x;
    UOutput[3*i+1] = l_U[i].y;
    UOutput[3*i+2] = l_U[i].z;
  }

  return UOutput;
}

float tledSolverGPU_ROM::GetKineticEnergy() {
  std::cerr << "Kinetic energy not available.\n";

  return std::numeric_limits<float>::quiet_NaN();
}

void tledSolverGPU_ROM::GetNodeVMStress(float* NodeSVM)
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

void tledSolverGPU_ROM::GetNodeVMStrain(float* NodeEVM)
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

void tledSolverGPU_ROM::GetMassVector(float* mass) const
{
   for (int i = 0; i < NumNodes; i++)
      mass[i] = h_Mf[3*i];
}

void tledSolverGPU_ROM::PrintNodalForces()
{
   // Get forces from device
   const unsigned int MemSzU = sizeof(float4)*NumNodes;
   float4* l_F = (float4*)malloc(MemSzU);
   cudaMemcpy(l_F,d_Fint,MemSzU,cudaMemcpyDeviceToHost);
   gpuV4ToV3(l_F, NumNodes, h_Fint);
   
   free(l_F);

   float FX=0,FY=0,FZ=0;
   printf("-----------------------------------------------------\n");
   printf("Nodal Forces (ALL)\n");
   printf("Node\tFX\t\tFY\t\tFZ\n");
   for (int i = 0; i < NumNodes; i++)
   {
      printf("%i\t%+10.5e\t%+10.5e\t%+10.5e\n",i+1,h_Fint[i].x,h_Fint[i].y,h_Fint[i].z);
      FX += h_Fint[i].x;
      FY += h_Fint[i].y;
      FZ += h_Fint[i].z;
   }
   printf("\nSums\t%+10.5e\t%+10.5e\t%+10.5e\n",FX,FY,FZ);
   printf("-----------------------------------------------------\n");
}

void tledSolverGPU_ROM::PrintDispNodalForces()
{
   // Get forces from device
   const unsigned int MemSzU = sizeof(float4)*NumNodes;
   float4* l_F = (float4*)malloc(MemSzU);
   cudaMemcpy(l_F,d_Fint,MemSzU,cudaMemcpyDeviceToHost);
   gpuV4ToV3(l_F, NumNodes, h_Fint);
   
   free(l_F);

   float FX=0,FY=0,FZ=0;
   printf("-----------------------------------------------------\n");
   printf("Nodal Forces (LOADED NODES)\n");
   printf("Node\tFX\t\tFY\t\tFZ\n");
   for (int i = 0; i < NumNodes; i++)
   {
      if ((h_Uload[i].x != 0)||(h_Uload[i].y != 0)||(h_Uload[i].z != 0))
      {
         printf("%i\t%+10.5e\t%+10.5e\t%+10.5e\n",i+1,h_Fint[i].x,h_Fint[i].y,h_Fint[i].z);
         FX += h_Fint[i].x;
         FY += h_Fint[i].y;
         FZ += h_Fint[i].z;
      }
   }
   printf("\nSums\t%+10.5e\t%+10.5e\t%+10.5e\n",FX,FY,FZ);
   printf("-----------------------------------------------------\n");
}

void tledSolverGPU_ROM::PrintDispNodalForceSums()
{
   // Get forces from device
   const unsigned int MemSzU = sizeof(float4)*NumNodes;
   float4* l_F = (float4*)malloc(MemSzU);
   cudaMemcpy(l_F,d_Fint,MemSzU,cudaMemcpyDeviceToHost);
   gpuV4ToV3(l_F, NumNodes, h_Fint);
   
   free(l_F);

   float FX=0,FY=0,FZ=0;
   printf("-----------------------------------------------------\n");
   printf("Nodal Force Sums (LOADED NODES)\n");
   printf("Node\tFX\t\tFY\t\tFZ\n");
   for (int i = 0; i < NumNodes; i++)
   {
      if ((h_Uload[i].x != 0)||(h_Uload[i].y != 0)||(h_Uload[i].z != 0))
      {
         FX += h_Fint[i].x;
         FY += h_Fint[i].y;
         FZ += h_Fint[i].z;
      }
   }
   printf("\nSums\t%+10.5e\t%+10.5e\t%+10.5e\n",FX,FY,FZ);
   printf("-----------------------------------------------------\n");
}

void tledSolverGPU_ROM::PrintNodalDisps()
{
   // Get disps from device
   const unsigned int MemSzU = sizeof(float4)*NumNodes;
   float4* l_Ucurr = (float4*)malloc(MemSzU);
   cudaMemcpy(l_Ucurr,d_Ucurr,MemSzU,cudaMemcpyDeviceToHost);
   gpuV4ToV3(l_Ucurr, NumNodes, h_Ucurr);
   
   free(l_Ucurr);

   printf("-----------------------------------------------------\n");
   printf("Nodal Displacements (ALL)\n");
   printf("Node\tUX\t\tUY\t\tUZ\n");
   for (int i = 0; i < NumNodes; i++)
      printf("%i\t%+10.5e\t%+10.5e\t%+10.5e\n",i+1,h_Ucurr[i].x,h_Ucurr[i].y,h_Ucurr[i].z);
   printf("-----------------------------------------------------\n");
}

void tledSolverGPU_ROM::PrintNodalForceSums()
{
}

void tledSolverGPU_ROM::InitialiseSolutionVariables(void)
{
   memset(h_Ucurr,0,sizeof(float)*NumNodes*3);
   memset(h_Fint,0,sizeof(float)*NumNodes*3);

   cudaMemset(d_Uprev,0,sizeof(float4)*NumNodes);
   cudaMemset(d_Ucurr,0,sizeof(float4)*NumNodes);
   cudaMemset(d_Unext,0,sizeof(float4)*NumNodes);
   
   cudaMemset(d_Fint,0,sizeof(float4)*NumNodes);

   cudaMemset(d_Divergence,0,sizeof(bool));
}

void tledSolverGPU_ROM::InitialiseConstraints(void)
{
   memset(h_Fext,0,sizeof(float)*NumNodes*4);
   cudaMemset(d_Fext,0,sizeof(float4)*NumNodes);
   // Reset all FixMask vals to 1
   for (int i = 0; i < NumNodes; i++)
   {
      h_FixMask[i].x = 1;
      h_FixMask[i].y = 1;
      h_FixMask[i].z = 1;
      h_FixMask[i].w = 1;	// NB: last value not used
   }
   cudaMemcpy(d_FixMask,h_FixMask,sizeof(int4)*NumNodes,cudaMemcpyHostToDevice);
   //Reset DispMask
   memset(h_Uload,0,sizeof(float)*NumNodes*4);
   memset(h_DispMask,0,sizeof(int)*NumNodes*4);
   memset(h_f,0,sizeof(float)*NumNodes*3);
   memset(h_UnextV,0,sizeof(float)*NumNodes*3);
   // NB: not necessary to reset d_Uload, d_f, d_UnextV, h_fr, nor h_tmp since this is done whenever gpuPerformStep is called
}

void tledSolverGPU_ROM::SetTimeStep(double dt)
{
   double Dt_old = Dt;
   Dt = dt;
   // Update CD coeffs
   ComputeCDCoeffs();
   cudaMemcpyToSymbol(c_gamma,&CDgamma,sizeof(float3));
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

float tledSolverGPU_ROM::GetStrainEnergy(void)
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

void tledSolverGPU_ROM::SetAllDisps(float* U)
{
   for (int i = 0; i < NumNodes; i++)
   {
      h_Ucurr[i].x = U[3*i];
      h_Ucurr[i].y = U[3*i+1];
      h_Ucurr[i].z = U[3*i+2];
   }
   const unsigned int MemSz = sizeof(float4)*NumNodes;
   float4* l_Ucurr = (float4*)malloc(MemSz);
   gpuV3ToV4(h_Ucurr,NumNodes,l_Ucurr);
   cudaMemcpy(d_Ucurr,l_Ucurr,MemSz,cudaMemcpyHostToDevice);
   
   free(l_Ucurr);
}

void tledSolverGPU_ROM::ComputeStresses(void)
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

#endif // _GPU_
