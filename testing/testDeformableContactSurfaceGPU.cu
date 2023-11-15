// =========================================================================
// File:       testDeformableContactSurfaceGPU.cu
// Purpose:    tledDeformableContactSurfaceGPU unit test
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   CUDA
// Created:    August 2013
// 
// Copyright (c) 2013, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================
#ifdef GPU_GP_CONTACT
#include "tledUnitTest.h"
#include "tledCUDAUnitTest.h"
#include "tledCUDAHelpers.h"
#include "tledCUDAUnitTest.h"
#include "tledCUDAMemoryBlock.h"
#include "tledDeformableContactSurfaceGPU.h"
#include "tledCUDAMemoryBlock.h"

#include <limits>

#include "tledCUDAHelpers.cu"

static void _TestNormalsOnSphereGPU() {
  const float maxDev = 5.f*tledPi/180.f;
  const int indexStartOffset = 10;
  const int indexEndOffset = 7;

  tledDeformableContactSurfaceT3GPU surface(tledUnitTest::LoadMSHMesh(tledUnitTest::GetMeshPath("sphere.msh"), "T4"));
  
  surface.Init();

  {
    tledCUDADeviceMemoryBlock &r_nodeInds = tledCUDADeviceMemoryBlock::GetNextFreeBufferWithSize<int>(surface.GetNumberOfNodes() - indexStartOffset - indexEndOffset);
    std::vector<int> nodeInds = tledSequenceGenerator::MakeSequence(indexStartOffset, surface.GetNumberOfNodes() - indexEndOffset);
    
    tledCUDAHelpers::CopyToDevice(r_nodeInds.GetBuffer<int>(), nodeInds);    
    surface.ComputeNodeNormals(r_nodeInds.GetBuffer<int>(), nodeInds.size());
    tledDeviceSyncDebug;
    r_nodeInds.ToggleActive();
  }

  {
    std::vector<float3> nodeNormals(surface.GetNumberOfNodes());
    float3 const *dp_normals;
    float err = 0;
   
    dp_normals = surface.GetHostGPUSurface().NodeNormals;
    tledCUDAHelpers::CopyFromDevice(&nodeNormals.front(), dp_normals, surface.GetNumberOfNodes());
    for (int nInd = indexStartOffset; nInd < surface.GetNumberOfNodes() - indexEndOffset; nInd++) {
      using namespace tledVectorArithmetic;

      const float *x = surface.GetNodeCoordinates(nInd);

      float n[3], ref[3];

      tledCUDAHelpers::ConvertFromFloatN(n, nodeNormals[nInd]);
      ScalarDiv(ref, x, Norm(x));
      tledUnitTestAssert(std::fabs(ComputeAngleNormalised(n, ref)) < maxDev);
      err += std::fabs(ComputeAngleNormalised(n, ref));      
    }

    err /= surface.GetNumberOfNodes() - indexStartOffset - indexEndOffset;
    std::cout << "Avg. normal error: " << err << std::endl;
  }
}

static void _TestUpdateSave(const std::string &meshPath) {
  const int numTimeSteps = 100;
  const int histSizeMin = 1;
  const int histSizeMax = 30;
  const tledMesh solidMesh = tledUnitTest::LoadMSHMesh(meshPath, "T4");

  tledCUDADeviceMemoryBlock &r_devDisps = tledCUDADeviceMemoryBlock::GetNextFreeBufferWithSize<float4>(solidMesh.GetNumNodes());
  tledCUDAHostMemoryBlock &r_hostDisps = tledCUDAHostMemoryBlock::GetNextFreeBufferWithSize<float4>(solidMesh.GetNumNodes());
    
  for (int h = histSizeMin; h < histSizeMax; h++) {
    tledDeformableContactSurfaceT3GPU surface(solidMesh);
  
    surface.SetCoordinateHistorySize(h);
    surface.Init();
    for (int t = 0; t < numTimeSteps; t++) {
      for (float4 *p_dst = r_hostDisps.GetBuffer<float4>(); p_dst < r_hostDisps.GetBuffer<float4>() + solidMesh.GetNumNodes(); p_dst++) {
	p_dst->x = 1.f*t;
	p_dst->y = 2.f*t;
	p_dst->z = 4.f*t;
	p_dst->w = std::numeric_limits<float>::quiet_NaN();
      }   
      r_devDisps.CopyFromHost<float4>(r_hostDisps, solidMesh.GetNumNodes());

      surface.Update(r_devDisps.GetBuffer<float4>());
      tledDeviceSyncDebug;
      tledUnitTestAssert(surface.GetSaveCount() == t);
      
      tledCUDAHelpers::CopyFromDevice(r_hostDisps.GetBuffer<float3>(), static_cast<const tledDeformableContactSurfaceT3GPU::GPUSurface&>(surface.GetHostGPUSurface()).OldNodeCoordinates, surface.GetNumberOfNodes());
      for (int n = 0; n < surface.GetNumberOfNodes(); n++) {
	const float3 &xOld = r_hostDisps.GetBuffer<float3>()[n];

	float3 x0 = make_float3(surface.GetNodeCoordinates(n)[0], surface.GetNodeCoordinates(n)[1], surface.GetNodeCoordinates(n)[2]);
	
	if (t >= h) {
	  float3 off = make_float3((t - h), (t - h)*2, (t - h)*4);

	  tledUnitTestAssert(std::fabs(xOld.x - (x0.x + off.x)) < 1e-4f);
	  tledUnitTestAssert(std::fabs(xOld.y - (x0.y + off.y)) < 1e-4f);
	  tledUnitTestAssert(std::fabs(xOld.z - (x0.z + off.z)) < 1e-4f);
	} else {
	  tledUnitTestAssert(std::fabs(xOld.x - x0.x) < 1e-4f);
	  tledUnitTestAssert(std::fabs(xOld.y - x0.y) < 1e-4f);
	  tledUnitTestAssert(std::fabs(xOld.z - x0.z) < 1e-4f);
	}
      }
      surface.Save();
      tledDeviceSyncDebug;
    }
  }

  r_devDisps.ToggleActive();
  r_hostDisps.ToggleActive();  
}

int main(int argc, char *argv[]) {
  tledUnitTest::InitUnitTest();
  tledCUDAUnitTest::InitCUDATests();

  _TestNormalsOnSphereGPU();

  _TestUpdateSave(tledUnitTest::GetMeshPath("sphere.msh"));
  _TestUpdateSave(tledUnitTest::GetMeshPath("organic_shape.msh"));

  tledCUDAUnitTest::FinishCUDATests();
  tledUnitTestPrintSuccess;

  return EXIT_SUCCESS;
} /* main */

#else 
#include "tledUnitTest.h"
#include "tledHelper.h"

int main(int argc, char *argv[]) {
  tledLogErrorStream(tledHelper::Warning() << "Test disabled");

  return EXIT_SUCCESS;
} /* main */
#endif
