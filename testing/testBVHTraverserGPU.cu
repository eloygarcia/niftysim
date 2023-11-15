// =========================================================================
// File:       testBVHTraverserGPU.cu
// Purpose:    tledBVHTraverserGPU unit test
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   CUDA
// Created:    May 2013
// 
// Copyright (c) 2013, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================
#ifdef GPU_GP_CONTACT
#define __GPU_TEST_LINK_CONFLICT_NO_INCLUDE

#include "tledUnitTest.h"
#include "tledCUDAUnitTest.h"
#include "tledBVHTraverserGPU.h"
#include "tledBV.h"
#include "tledAABB.h"
#include "tledBV_kernels.h"
#include "tledNarrowConeSelfCollisionBVHGPU.h"
#include "tledDeformableContactSurfaceGPU.h"
#include "tledSelfCollisionBVHCreator.h"
#include "tledDeformableDeformableBVHTraverserImplGPU.h"
#include "tledGreedySelfCollisionBVHUpdaterGPU.h"
#include "tledUnstructuredContactManager.h"

#include <vector>
#include <algorithm>

typedef tledAABB<2> TestBV;

class TestTraverserSurfaceGPU : public tledDeformableContactSurfaceT3GPU {
public:
  typedef tledDeformableContactSurfaceT3GPU Superclass;
  typedef Superclass::Facet Facet;
  typedef Superclass::GPUSurface GPUSurface;

private:
  const int mc_NumFacetsPerSide;

public:
  int GetNumberOfFacetsPerSide(void) const { return mc_NumFacetsPerSide; }

  float GetDX(void) const { return 1.25f; }
  float GetX(void) const { return this->GetDX()*this->GetNumberOfFacetsPerSide(); }
  float GetY(void) const { return 1.f; }
  float GetZ(void) const { return 0.01f; }

  void FakeHistory(void);

private:
  void _InsertEdges(const int fInd);

public:
  TestTraverserSurfaceGPU(const int numFacetsPerSide);
};

typedef tledNarrowConeSelfCollisionBVHImplGPU<TestTraverserSurfaceGPU::Superclass, TestBV > TestBVHGPU;

class TestBVHCreator : public tledSelfCollisionBVHCreator<TestBVHGPU> {
public:
  typedef tledSelfCollisionBVHCreator<TestBVHGPU> Superclass;
  typedef TestBVHGPU BVH;
  typedef TestBV BoundingVolume;

private:
  const TestTraverserSurfaceGPU *mpc_Surface;

private:
  int _GenerateRecursive(const int startInd, const int endInd, const int parentInd);

protected:
  virtual void SplitBV(std::vector<int>::iterator p_childIndexBounds[][2], const int BVIndex, const std::vector<int>::iterator &pIndicesBegin, const std::vector<int>::iterator &pIndicesEnd) { tledFatalError("Not supported"); }
  virtual void GenerateMain(void);

public:
  TestBVHCreator(const TestTraverserSurfaceGPU &surface) : mpc_Surface(&surface) {}
};

class TestBroadPhaseBVHTraverserGPU : public tledDeformableDeformableBVHTraverserImplGPU<TestBVHGPU> {
public:
  typedef tledDeformableDeformableBVHTraverserImplGPU<TestBVHGPU> Superclass;

private:
  TestTraverserSurfaceGPU *mp_Surface;
  bool m_WasChecked;

public:
  void RunDetectionTest(void);
  void RunInitTest(void);
  virtual void PrepareNarrowPhase(const int2 *dpc_inBPList, const int *dpc_actualNumResultItems, const int maxBPResultItems);

public:
  TestBroadPhaseBVHTraverserGPU(TestBVHGPU &r_bvh, TestTraverserSurfaceGPU &r_surface) : Superclass(r_bvh), mp_Surface(&r_surface), m_WasChecked(false) {}
};

class TestNarrowPhaseBVHTraverserGPU : public tledDeformableDeformableBVHTraverserImplGPU<TestBVHGPU> {
public:
  typedef tledDeformableDeformableBVHTraverserImplGPU<TestBVHGPU> Superclass;

private:
  TestTraverserSurfaceGPU *mp_Surface;

protected:
  virtual void RunBroadPhase(void);
  
public:
  void CheckNarrowPhaseResults(void);

public:
  TestNarrowPhaseBVHTraverserGPU(TestBVHGPU &r_bvh, TestTraverserSurfaceGPU &r_surface) : Superclass(r_bvh), mp_Surface(&r_surface) {}
};

void TestTraverserSurfaceGPU::_InsertEdges(const int fInd) {
  for (int e = 0; e < 3; e++) {
    std::pair<int, int> edge;      

    edge.first = 3*fInd + std::min(e, (e + 1)%3);
    edge.second = 3*fInd + std::max(e, (e + 1)%3);

    this->GetFacet(fInd).EdgeIndices[e] = this->GetAllEdges().size();
    this->GetAllEdges().push_back(edge);
  }
}

TestTraverserSurfaceGPU::TestTraverserSurfaceGPU(const int numFacetsPerSide) : mc_NumFacetsPerSide(numFacetsPerSide) {  
  const int numFacetVtcs = Facet::NumberOfVertices;
  const int numNodesPerSize = numFacetsPerSide*numFacetVtcs;
  const float X = this->GetX();
  const float Y = this->GetY();
  const float Z = this->GetZ();
  const float dx = this->GetDX();

  this->SetNumberOfNodes(numNodesPerSize*2);
  this->SetNumberOfFacets(numFacetsPerSide*2);
  for (int f = 0; f < numFacetsPerSide; f++) {
    for (int v = 0; v < 3; v++) {
      this->GetAllNodeCoordinates()[3*(3*f+v)] = f*dx;
      this->GetAllNodeCoordinates()[3*(3*f+v)+1] = 0;
      this->GetAllNodeCoordinates()[3*(3*f+v)+2] = 0;
      this->GetFacet(f).NodeIndices[v] = 3*f + v;
    }
    this->GetAllNodeCoordinates()[3*(3*f+2)+0] += dx/2 + 1e-2f;
    this->GetAllNodeCoordinates()[3*(3*f+1)+0] += dx/4;
    this->GetAllNodeCoordinates()[3*(3*f+1)+1] += Y;

    _InsertEdges(f);
  }

  for (int f = numFacetsPerSide; f < 2*numFacetsPerSide; f++) {
    for (int v = 0; v < 3; v++) {
      this->GetAllNodeCoordinates()[3*(3*f+v)] = X - dx/2 - (f - numFacetsPerSide)*dx;
      this->GetAllNodeCoordinates()[3*(3*f+v)+1] = Y + Y/8;
      this->GetAllNodeCoordinates()[3*(3*f+v)+2] = Z;
      this->GetFacet(f).NodeIndices[v] = 3*f + v;
    }
    this->GetAllNodeCoordinates()[3*(3*f+1)+0] -= dx/2 + 1e-2f;
    this->GetAllNodeCoordinates()[3*(3*f+2)+0] -= dx/4;
    this->GetAllNodeCoordinates()[3*(3*f+2)+1] -= Y;

    _InsertEdges(f);
  }

  this->Init();
}

void TestTraverserSurfaceGPU::FakeHistory() {
  std::vector<float3> oldCoords(this->GetNumberOfNodes());
  float3 *dp_cds = static_cast<tledDeformableContactSurfaceT3GPU::GPUSurface&>(this->GetHostGPUSurface()).OldNodeCoordinates;

  tledDeformableContactSurface::Save();
  for (int n = 0; n < 3*this->GetNumberOfFacetsPerSide(); n++) {
    oldCoords[n].x = this->GetNodeCoordinates(n)[0];
    oldCoords[n].y = this->GetNodeCoordinates(n)[1];
    oldCoords[n].z = this->GetNodeCoordinates(n)[2];
  }

  for (int n = 3*this->GetNumberOfFacetsPerSide(); n < 3*this->GetNumberOfFacets(); n++) {
    oldCoords[n].x = this->GetNodeCoordinates(n)[0];
    oldCoords[n].y = this->GetNodeCoordinates(n)[1];
    oldCoords[n].z = -this->GetZ()/8;
  }
  tledCUDAHelpers::CopyToDevice(dp_cds, &oldCoords.front(), this->GetNumberOfNodes());
}

void TestBVHCreator::GenerateMain() {
  assert(BoundingVolume::NumberOfChildBVs == 2);
  
  this->GetOutput().SetMargin(mpc_Surface->GetZ());
  this->GetOutput().SetBVMaxDisplacement(mpc_Surface->GetZ()/4);

  this->AddBV(-1);  
  
  {
    int childInds[2]; // Compiler bug workaround

    childInds[0] = _GenerateRecursive(0, mpc_Surface->GetNumberOfFacetsPerSide(), 0);
    childInds[1] = _GenerateRecursive(mpc_Surface->GetNumberOfFacetsPerSide(), mpc_Surface->GetNumberOfFacets(), 0);

    std::copy(childInds, childInds + 2, this->GetOutput().GetBV(0).ChildIndices);
  }

  this->GetOutput().GetBV(0).SetDisconnectedGeometryFlag();
  this->GetOutput().GetNonAdjacentGeometryNodes().push_back(0);
  this->GetOutput().RefitInteriorBV(0);

  assert(this->GetOutput().GetNumberOfBVs() == 2*(2*mpc_Surface->GetNumberOfFacetsPerSide()) - 1);
}

int TestBVHCreator::_GenerateRecursive(const int startInd, const int endInd, const int parentInd) {
  const int newBVInd = this->AddBV(parentInd);

  if (startInd + 1 >= endInd) {
    this->GetOutput().GetBV(newBVInd).PrimitiveIndex = startInd;
    std::fill(this->GetOutput().GetBV(newBVInd).ChildIndices, this->GetOutput().GetBV(newBVInd).ChildIndices + BVH::BVHOrder, -1);

    this->GetOutput().RefitLeafBV(newBVInd);
    this->GetOutput().GetLeafBVIndices().push_back(newBVInd);
    this->GetOutput().GetAllPrimitiveBVIndices()[this->GetOutput().GetBV(newBVInd).PrimitiveIndex] = newBVInd;
  } else {
    assert((endInd - startInd)%2 == 0);
    this->GetOutput().GetBV(newBVInd).PrimitiveIndex = -1;
    this->GetOutput().GetBV(newBVInd).ChildIndices[0] = _GenerateRecursive(startInd, (startInd + endInd)/2, newBVInd);
    this->GetOutput().GetBV(newBVInd).ChildIndices[1] = _GenerateRecursive((startInd + endInd)/2, endInd, newBVInd);

    this->GetOutput().RefitInteriorBV(newBVInd);
  }

  return newBVInd;
}

void TestBroadPhaseBVHTraverserGPU::PrepareNarrowPhase(const int2 *dpc_inBPList, const int *dpc_actualNumResultItems, const int maxBPResultItems) {
  int numItems;
  tledCUDAHostMemoryBlock &r_hMem = tledCUDAHostMemoryBlock::GetNextFreeBufferWithSize<int2>(maxBPResultItems);

  tledCUDAHelpers::CopyFromDevice(&numItems, dpc_actualNumResultItems);
  tledUnitTestAssert(numItems <= maxBPResultItems);
  tledUnitTestAssert(numItems == mp_Surface->GetNumberOfFacetsPerSide());
  tledCUDAHelpers::CopyFromDevice(r_hMem.GetBuffer<int2>(), dpc_inBPList, numItems);

  for (int2 const *pc_p = r_hMem.GetBuffer<int2>(); pc_p < r_hMem.GetBuffer<int2>() + numItems; pc_p++) {
    tledUnitTestAssert(pc_p->x >= 0 && pc_p->x < mp_Surface->GetNumberOfFacetsPerSide() && pc_p->y >= mp_Surface->GetNumberOfFacetsPerSide() && pc_p->y < mp_Surface->GetNumberOfFacets());
    tledUnitTestAssert(!(this->GetMasterBVH().GetBV(this->GetMasterBVH().GetLeafBVIndices()[pc_p->x]).Bounds[0][1] < this->GetMasterBVH().GetBV(this->GetMasterBVH().GetLeafBVIndices()[pc_p->y]).Bounds[0][0]
    			 || this->GetMasterBVH().GetBV(this->GetMasterBVH().GetLeafBVIndices()[pc_p->y]).Bounds[0][1] < this->GetMasterBVH().GetBV(this->GetMasterBVH().GetLeafBVIndices()[pc_p->x]).Bounds[0][0]));
  }

  m_WasChecked = true;

  r_hMem.ToggleActive();
}

template <class TBVH>
static bool _IsChild(const TBVH &bvh, const int bvInd, const int rootInd) {
  if (bvInd == rootInd) return true;
  else {
    if (!bvh.IsLeaf(rootInd)) {
      for (int const *pc_c = bvh.GetBV(rootInd).ChildIndices; pc_c < bvh.GetBV(rootInd).ChildIndices + TBVH::BVHOrder; pc_c++) if (*pc_c >= 0) {
	  if (_IsChild(bvh, bvInd, *pc_c)) return true;
	}
    }
  }

  return false;
}

void TestBroadPhaseBVHTraverserGPU::RunInitTest() {
  tledCUDADeviceMemoryBlock &r_devMem = tledCUDADeviceMemoryBlock::GetNextFreeBufferWithSize<int2>(mp_Surface->GetNumberOfFacetsPerSide());
  tledCUDAHostMemoryBlock &r_hostMem = tledCUDAHostMemoryBlock::GetNextFreeBufferWithSize<int2>(mp_Surface->GetNumberOfFacetsPerSide());
  int *dp_ctr, counter;

  tledUnitTestAssert(this->GetNumberOfInitialBroadPhaseSearchPairs() == mp_Surface->GetNumberOfFacetsPerSide());
  tledCUDAHelpers::AllocateDeviceMemory(dp_ctr);  
  this->InitialiseBroadPhaseSearch(r_devMem.GetBuffer<int2>(), dp_ctr);
  tledCUDAHelpers::CopyFromDevice(&counter, dp_ctr);
  tledUnitTestAssert(counter == this->GetNumberOfInitialBroadPhaseSearchPairs());
  r_hostMem.CopyFromDevice<int2>(r_devMem, this->GetNumberOfInitialBroadPhaseSearchPairs());
  for (int2 const *pc_pair = r_hostMem.GetBuffer<int2>(); pc_pair < r_hostMem.GetBuffer<int2>() + this->GetNumberOfInitialBroadPhaseSearchPairs(); pc_pair++) {
    tledUnitTestAssert(pc_pair->x == this->GetMasterBVH().GetBV(0).ChildIndices[0]);
    tledUnitTestAssert(_IsChild(this->GetMasterBVH(), pc_pair->y, this->GetMasterBVH().GetBV(0).ChildIndices[1]));
  }
  
  r_devMem.ToggleActive();
  r_hostMem.ToggleActive();
}

void TestBroadPhaseBVHTraverserGPU::RunDetectionTest() {
  m_WasChecked = false;
  this->RunBroadPhase();  
  tledUnitTestAssert(m_WasChecked);
}

static void _RunBroadPhaseTest(const int numFacets) {
  TestTraverserSurfaceGPU surface(numFacets);
  TestBVHGPU bvh(surface);
  tledGreedySelfCollisionBVHUpdaterGPU<TestBVHGPU> *p_updater = new tledGreedySelfCollisionBVHUpdaterGPU<TestBVHGPU>();
  TestBVHCreator builder(surface);
  tledUnstructuredContactManager mgr;

  tledCUDADeviceMemoryBlock::SaveAllocationCounter();

  mgr.SetBVType("AABB2");
  mgr.SetBVMargin(1e-3f);

  bvh.SetUpdater(*p_updater);
  bvh.SetMargin(mgr.GetBVMargin());
  bvh.Init(builder);

  {
    TestBroadPhaseBVHTraverserGPU traverser(bvh, surface);

    traverser.Init(mgr);    
    traverser.RunInitTest();
    traverser.RunDetectionTest();
  }

  tledCUDADeviceMemoryBlock::CheckAllocationCounter();
}

void TestNarrowPhaseBVHTraverserGPU::RunBroadPhase() {
  const int numFacetContacts = mp_Surface->GetNumberOfFacetsPerSide();
   
  this->SetDoMaster(true);

  {
    std::vector<int2> primPairs(numFacetContacts);
    tledCUDADeviceMemoryBlock &r_devPairs = tledCUDADeviceMemoryBlock::GetNextFreeBufferWithSize<int2>(numFacetContacts);
    int *dp_numItems;

    tledCUDAHelpers::AllocateDeviceMemory(dp_numItems);
    for (int n = 0; n < numFacetContacts; n++) {
      primPairs[n].x = numFacetContacts - (n + 1);
      primPairs[n].y = numFacetContacts + n;
      assert(primPairs[n].x >= 0 && primPairs[n].x < mp_Surface->GetNumberOfFacetsPerSide());
      assert(primPairs[n].y >= mp_Surface->GetNumberOfFacetsPerSide() && primPairs[n].y < 2*mp_Surface->GetNumberOfFacetsPerSide());
    }
    tledCUDAHelpers::CopyToDevice(r_devPairs.GetBuffer<int2>(), &primPairs.front(), numFacetContacts);
    tledCUDAHelpers::CopyToDevice(dp_numItems, &numFacetContacts);
    this->PrepareNarrowPhase(r_devPairs.GetBuffer<int2>(), dp_numItems, numFacetContacts);
    tledCheckCUDAErrors(cudaFree(dp_numItems));
    r_devPairs.ToggleActive();
  }

  tledUnitTestAssert(this->GetNumberOfNodeFacetContacts() == 3*numFacetContacts);
  tledUnitTestAssert(this->GetNumberOfEdgeEdgeContacts() == 3*3*numFacetContacts);

  {
    std::vector<int2> items(3*3*numFacetContacts);     

    tledUnitTestAssert(this->GetEdgeEdgeResults().IsActive());
    tledCUDAHelpers::CopyFromDevice(&items.front(), this->GetEdgeEdgeResults().GetBuffer<int2>(), 3*3*numFacetContacts);
    for (std::vector<int2>::const_iterator ic_e = items.begin(); ic_e < items.end(); ic_e++) {
      int mFacetInd, e, sFacetInd;

      for (mFacetInd = 0; mFacetInd < mp_Surface->GetNumberOfFacetsPerSide(); mFacetInd++) {	 
	for (e = 0; e < 3 && mp_Surface->GetFacet(mFacetInd).EdgeIndices[e] != ic_e->y; e++); 
	if (e < 3) break;
      }
      tledUnitTestAssert(mFacetInd < mp_Surface->GetNumberOfFacetsPerSide());

      sFacetInd = 2*numFacetContacts - mFacetInd - 1;
      for (e = 0; e < 3 && mp_Surface->GetFacet(sFacetInd).EdgeIndices[e] != ic_e->x; e++);
      tledUnitTestAssert(e < 3);
    }

    tledUnitTestAssert(this->GetNodeFacetResults().IsActive());
    items.resize(3*numFacetContacts);
    tledCUDAHelpers::CopyFromDevice(&items.front(), this->GetNodeFacetResults().GetBuffer<int2>(), 3*numFacetContacts);
    for (std::vector<int2>::const_iterator ic_nf = items.begin(); ic_nf < items.end(); ic_nf++) {
      int mFacetInd, sFacetInd;

      tledUnitTestAssert(ic_nf->y >= 0 && ic_nf->y < mp_Surface->GetNumberOfFacetsPerSide());
      mFacetInd = ic_nf->y;
      sFacetInd = 2*numFacetContacts - mFacetInd - 1;
      tledUnitTestAssert(ic_nf->x >= 3*sFacetInd && ic_nf->x < 3*(sFacetInd + 1));
    }
  }
}

void TestNarrowPhaseBVHTraverserGPU::CheckNarrowPhaseResults() {  
  tledUnitTestAssert(this->GetNumberOfNodeFacetContacts() == mp_Surface->GetNumberOfFacetsPerSide());
  tledUnitTestAssert(this->GetNumberOfEdgeEdgeContacts() == 2*mp_Surface->GetNumberOfFacetsPerSide());

  {
    typedef tledBVHTraverserGPU::EdgeEdgeNarrowPhaseResult __Result;

    std::vector<__Result> items(2*mp_Surface->GetNumberOfFacetsPerSide());

    tledCUDAHelpers::CopyFromDevice(&items.front(), this->GetEdgeEdgeResults().GetBuffer<__Result>(), 2*mp_Surface->GetNumberOfFacetsPerSide());
    for (std::vector<__Result>::const_iterator ic_i = items.begin(); ic_i < items.end(); ic_i++) {
      int mFacetInd, sFacetInd;

      tledUnitTestAssert(ic_i->MasterEdge.x >= 0 && ic_i->MasterEdge.x < 3*mp_Surface->GetNumberOfFacetsPerSide());
      mFacetInd = ic_i->MasterEdge.x/3;
      tledUnitTestAssert(ic_i->MasterEdge.y >= 3*mFacetInd && ic_i->MasterEdge.y < 3*(mFacetInd + 1));
      sFacetInd = 2*mp_Surface->GetNumberOfFacetsPerSide() - mFacetInd - 1;
      tledUnitTestAssert(ic_i->SlaveEdge.x >= 3*sFacetInd && ic_i->SlaveEdge.x < 3*(sFacetInd + 1));
      tledUnitTestAssert(ic_i->SlaveEdge.y >= 3*sFacetInd && ic_i->SlaveEdge.y < 3*(sFacetInd + 1));
      
      tledUnitTestAssert(ic_i->Xi.y > 0.f && ic_i->Xi.y < 1.f);
      tledUnitTestAssert(ic_i->Xi.z > 0.f && ic_i->Xi.z < 1.f);
      tledUnitTestAssert(ic_i->Xi.x < 0.f);
      
      tledUnitTestAssert(std::fabs(ic_i->Normal.x*ic_i->Normal.x + ic_i->Normal.y*ic_i->Normal.y + ic_i->Normal.z*ic_i->Normal.z - 1.f) < 1e-3f);
      tledUnitTestAssert(std::fabs(ic_i->Normal.z + 1.f) < 1e-2f);
    }
  }

  {
    typedef tledBVHTraverserGPU::NodeFacetNarrowPhaseResult<3> __Result;

    std::vector<__Result> items(mp_Surface->GetNumberOfFacetsPerSide());

    tledCUDAHelpers::CopyFromDevice(&items.front(), this->GetNodeFacetResults().GetBuffer<__Result>(), mp_Surface->GetNumberOfFacetsPerSide());
    for (std::vector<__Result>::const_iterator ic_i = items.begin(); ic_i < items.end(); ic_i++) {
      int mFacetInd, sFacetInd;

      sFacetInd = ic_i->ContactNodeIndices[0]/3;
      tledUnitTestAssert(sFacetInd >= mp_Surface->GetNumberOfFacetsPerSide() && sFacetInd < mp_Surface->GetNumberOfFacets());
      mFacetInd = 2*mp_Surface->GetNumberOfFacetsPerSide() - sFacetInd - 1;
      for (int v = 0; v < 3; v++) {
	tledUnitTestAssert(ic_i->ContactNodeIndices[v+1] >= 3*mFacetInd && ic_i->ContactNodeIndices[v+1] < 3*(mFacetInd + 1));
      }
      
      tledUnitTestAssert(std::fabs(ic_i->ShapeValues[0] + ic_i->ShapeValues[1] + ic_i->ShapeValues[2] - 1.f) < 1e-3f);
      for (int s = 0; s < 3; s++) {
	tledUnitTestAssert(ic_i->ShapeValues[s] > 0.f);
      }

      tledUnitTestAssert(std::fabs(ic_i->Normal.x*ic_i->Normal.x + ic_i->Normal.y*ic_i->Normal.y + ic_i->Normal.z*ic_i->Normal.z - 1.f) < 1e-3f);
      tledUnitTestAssert(std::fabs(ic_i->Normal.z + 1.f) < 1e-2f);
    }
  }

  this->GetEdgeEdgeResults().ToggleActive();
  this->GetNodeFacetResults().ToggleActive();
}

static void _RunNarrowPhaseTest(const int numFacetsPerSide) {
  TestTraverserSurfaceGPU surface(numFacetsPerSide);
  TestBVHGPU bvh(surface);
  tledGreedySelfCollisionBVHUpdaterGPU<TestBVHGPU> *p_updater = new tledGreedySelfCollisionBVHUpdaterGPU<TestBVHGPU>();
  TestBVHCreator builder(surface);
  tledUnstructuredContactManager mgr;
      
  tledCUDADeviceMemoryBlock::SaveAllocationCounter();
  surface.FakeHistory();

  mgr.SetBVType("AABB2");
  mgr.SetBVMargin(1e-3f);

  bvh.SetUpdater(*p_updater);
  bvh.SetMargin(mgr.GetBVMargin());
  bvh.Init(builder);
  
  {
    TestNarrowPhaseBVHTraverserGPU traverser(bvh, surface);

    traverser.Init(mgr);    
    traverser.FindCollisions();
    traverser.CheckNarrowPhaseResults();
  }

  tledCUDADeviceMemoryBlock::CheckAllocationCounter();
}

int main(int argc, char *argv[]) {
  tledUnitTest::InitUnitTest();
  tledCUDAUnitTest::InitCUDATests();

  _RunBroadPhaseTest(8);
  _RunBroadPhaseTest(16);
  _RunBroadPhaseTest(32);

  _RunNarrowPhaseTest(8);
  _RunNarrowPhaseTest(16);
  _RunNarrowPhaseTest(256);

  tledCUDAUnitTest::FinishCUDATests();
  tledUnitTestPrintSuccess;

  return EXIT_SUCCESS;
} /* main */
#else
#include "tledUnitTest.h"

int main(int argc, char *argv[]) {
  tledUnitTestDisabled;
} /* main */
#endif
