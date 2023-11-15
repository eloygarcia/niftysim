// =========================================================================
// File:       tledShellSolver.cu
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    March 2012
// 
// Copyright (c) 2012, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================
#ifndef tledShellSolverGPU_CU
#define tledShellSolverGPU_CU

#include "tledShellSolverGPU.h"
#include "tledShellSolver_kernels.h"
#include "tledElementMembraneSimpleLinear.h"
#include "tledElementMembraneNonLinear.h"
#include "tledElementShellBSTP1.h"
#include "tledSolver.h"
#include "tledHelper.h"

#include <typeinfo>
#include <vector>
#include <cstdlib>
#include <iostream>

#include "tledMembraneMaterialLinear.cu"
#include "tledMembraneMaterialNeoHookean.cu"
#include "tledShellMaterialLinearPlateDecorator.cu"
#include "tledElementMembraneSimpleLinear.cu"
#include "tledElementMembraneNonLinear.cu"
#include "tledElementShellBSTP1.cu"
#include "tledElementMembraneSimpleLinear_kernels.cu"
#include "tledElementMembraneNonLinear_kernels.cu"
#include "tledMembraneMaterialLinear_kernels.cu"
#include "tledMembraneMaterialNeoHookean_kernels.cu"
#include "tledShellMaterialLinearPlateDecorator_kernels.cu"
#include "tledElementShellBSTP1_kernels.cu"
#include "tledShellSolver_kernels.cu"

template <class TShellElement>
void tledShellSolverGPU::_ElementSetImpl<TShellElement>::InitGPU(const tledSolver &mainSolver) {
  std::vector<int> nodeIndexLTB;
  std::vector<typename TShellElement::GPUElement> hostMem;
  std::vector<int> hostNodeElementVtxLTB;
  std::vector<int2> hostNodeElementVtxLTBBase;
  std::vector<std::vector<int> > tmpNodeElementVertexMap;

  hostMem.reserve(this->GetNumberOfElements());
  nodeIndexLTB.reserve(this->GetNumberOfElements()*TShellElement::Facet::NumberOfVertices);
  tmpNodeElementVertexMap.resize(mainSolver.GetMesh()->GetNumNodes());
  for (typename std::vector<TShellElement>::iterator i_el = this->m_Elements.begin(); i_el < this->m_Elements.end(); i_el++) {
    const int elInd = i_el - this->m_Elements.begin();

    typename TShellElement::GPUElement gpuEl;
    int weirdBug; /* work around for NVCC 3.2 bug, do not remove */
    
    i_el->InitGPU(gpuEl);
    hostMem.push_back(gpuEl);
    weirdBug = gpuEl.ElementNodeIndices.x;
    nodeIndexLTB.push_back(weirdBug);
    //nodeIndexLTB.push_back(gpuEl.ElementNodeIndices.x);
    tmpNodeElementVertexMap[gpuEl.ElementNodeIndices.x].push_back(elInd*TShellElement::Facet::NumberOfVertices + 0);
    weirdBug = gpuEl.ElementNodeIndices.y;
    nodeIndexLTB.push_back(weirdBug);
    //nodeIndexLTB.push_back(gpuEl.ElementNodeIndices.y);
    tmpNodeElementVertexMap[gpuEl.ElementNodeIndices.y].push_back(elInd*TShellElement::Facet::NumberOfVertices + 1);
    weirdBug = gpuEl.ElementNodeIndices.z;
    nodeIndexLTB.push_back(weirdBug);
    //nodeIndexLTB.push_back(gpuEl.ElementNodeIndices.z);
    tmpNodeElementVertexMap[gpuEl.ElementNodeIndices.z].push_back(elInd*TShellElement::Facet::NumberOfVertices + 2);
    if (TShellElement::Facet::NumberOfVertices > 3) {
      weirdBug = gpuEl.ElementNodeIndices.z;
      nodeIndexLTB.push_back(weirdBug);
      //nodeIndexLTB.push_back(gpuEl.ElementNodeIndices.w);
      tmpNodeElementVertexMap[gpuEl.ElementNodeIndices.w].push_back(elInd*TShellElement::Facet::NumberOfVertices + 3);
    }
  }
  nodeIndexLTB = tledHelper::MakeSortedUnique(nodeIndexLTB);

  hostNodeElementVtxLTBBase.reserve(nodeIndexLTB.size());
  for (std::vector<int>::const_iterator ic_ind = nodeIndexLTB.begin(); ic_ind < nodeIndexLTB.end(); ic_ind++) {
    int2 base;

    base.x = hostNodeElementVtxLTB.size();
    hostNodeElementVtxLTB.insert(hostNodeElementVtxLTB.end(), tmpNodeElementVertexMap[*ic_ind].begin(), tmpNodeElementVertexMap[*ic_ind].end());
    base.y = hostNodeElementVtxLTB.size();
    hostNodeElementVtxLTBBase.push_back(base);
  }

  if (hostMem.size() > 0) {
    tledCUDAHelpers::AllocateDeviceMemory(mdp_GPUElements, hostMem.size());
    tledCUDAHelpers::CopyToDevice(mdp_GPUElements, &hostMem.front(), hostMem.size());
    mdp_GPUMaterial = this->GetMaterial().InitGPU();
    
    tledCUDAHelpers::AllocateDeviceMemory(mdp_ShellElementForces, this->GetNumberOfElements()*TShellElement::Facet::NumberOfVertices);
    
    m_NumberOfNodes = nodeIndexLTB.size();
    tledCUDAHelpers::AllocateDeviceMemory(mdp_NodeIndexLookupTable, nodeIndexLTB.size());
    tledCUDAHelpers::CopyToDevice(mdp_NodeIndexLookupTable, &nodeIndexLTB.front(), nodeIndexLTB.size());

    m_NodeElementVertexLookupTableSize = hostNodeElementVtxLTB.size();
    tledCUDAHelpers::AllocateDeviceMemory(mdp_NodeElementVertexLookupTable, hostNodeElementVtxLTB.size());
    tledCUDAHelpers::CopyToDevice(mdp_NodeElementVertexLookupTable, &hostNodeElementVtxLTB.front(), hostNodeElementVtxLTB.size());

    tledCUDAHelpers::AllocateDeviceMemory(mdp_NodeElementVertexLookupBaseIndex, hostNodeElementVtxLTBBase.size());
    tledCUDAHelpers::CopyToDevice(mdp_NodeElementVertexLookupBaseIndex, &hostNodeElementVtxLTBBase.front(), hostNodeElementVtxLTBBase.size());

#ifndef NDEBUG
    for (std::vector<int>::const_iterator ic_nInd = nodeIndexLTB.begin(); ic_nInd < nodeIndexLTB.end(); ic_nInd++) {
      for (std::vector<int>::const_iterator ic_vInd = hostNodeElementVtxLTB.begin() + hostNodeElementVtxLTBBase[ic_nInd-nodeIndexLTB.begin()].x; ic_vInd < hostNodeElementVtxLTB.begin() + hostNodeElementVtxLTBBase[ic_nInd-nodeIndexLTB.begin()].y; ic_vInd++) {
	switch ((*ic_vInd)%TShellElement::Facet::NumberOfVertices) {	  
	case 0:
	  assert(hostMem[(*ic_vInd)/TShellElement::Facet::NumberOfVertices].ElementNodeIndices.x == *ic_nInd);
	  continue;

	case 1:
	  assert(hostMem[(*ic_vInd)/TShellElement::Facet::NumberOfVertices].ElementNodeIndices.y == *ic_nInd);
	  continue;

	case 2:
	  assert(hostMem[(*ic_vInd)/TShellElement::Facet::NumberOfVertices].ElementNodeIndices.z == *ic_nInd);
	  continue;

	case 3:
	  assert(hostMem[(*ic_vInd)/TShellElement::Facet::NumberOfVertices].ElementNodeIndices.w == *ic_nInd);
	}
      }
    }
#endif
  } else mdp_GPUElements = NULL; 
}

template <>
tledShellSolver::ElementSetImpl<tledElementShellBSTP1, tledShellSolverGPU::ElementSet>::ElementSetImpl(tledShellMaterial &r_mat, const tledElementShellBSTP1::Surface &surface, const std::vector<int> &elInds) : tledShellSolverGPU::ElementSet(r_mat, elInds) {
  tledSurfaceTopology<tledShellMesh<3> > topology(surface);
  std::vector<int> reverseElSetIndexMap(surface.GetNumberOfFacets(), -1);

  topology.ComputeEdges();
  m_Elements.resize(elInds.size());

  for (std::vector<int>::const_iterator ic_e = elInds.begin(); ic_e < elInds.end(); ic_e++) reverseElSetIndexMap[*ic_e] = ic_e - elInds.begin();
  for (std::vector<tledElementShellBSTP1>::iterator i_el = m_Elements.begin(); i_el < m_Elements.end(); i_el++) {
    i_el->SetMaterial(r_mat);
    i_el->InitialiseElement(topology, elInds[i_el-m_Elements.begin()]);
  }

  for (std::vector<tledElementShellBSTP1>::iterator i_el = m_Elements.begin(); i_el < m_Elements.end(); i_el++) {
    i_el->InitialiseEdgeElements(m_Elements, topology, elInds[i_el-m_Elements.begin()], reverseElSetIndexMap);
  }

  for (std::vector<tledElementShellBSTP1>::iterator i_el = m_Elements.begin(); i_el < m_Elements.end(); i_el++) i_el->FinaliseInitialisation();
}

template <>
void tledShellSolverGPU::_ElementSetImpl<tledElementShellBSTP1>::InitGPU(const tledSolver &mainSolver) {
  std::vector<int> nodeIndexLTB;
  std::vector<tledElementShellBSTP1::GPUElement> hostMem;
  std::vector<int> hostNodeElementVtxLTB;
  std::vector<int2> hostNodeElementVtxLTBBase;
  std::vector<std::vector<int> > tmpNodeElementVertexMap;

  hostMem.reserve(this->GetNumberOfElements());
  nodeIndexLTB.reserve(this->GetNumberOfElements()*3*2);
  tmpNodeElementVertexMap.resize(mainSolver.GetMesh()->GetNumNodes());
  for (std::vector<tledElementShellBSTP1>::iterator i_el = this->m_Elements.begin(); i_el < this->m_Elements.end(); i_el++) {
    tledElementShellBSTP1::GPUElement gpuEl;
    
    i_el->InitGPU(gpuEl);
    hostMem.push_back(gpuEl);
    tmpNodeElementVertexMap[gpuEl.ElementNodeIndices.x].push_back(3*2*(i_el - this->m_Elements.begin()) + 0);
    tmpNodeElementVertexMap[gpuEl.ElementNodeIndices.y].push_back(3*2*(i_el - this->m_Elements.begin()) + 1);
    tmpNodeElementVertexMap[gpuEl.ElementNodeIndices.z].push_back(3*2*(i_el - this->m_Elements.begin()) + 2);

    nodeIndexLTB.push_back(gpuEl.ElementNodeIndices.x);
    nodeIndexLTB.push_back(gpuEl.ElementNodeIndices.y);
    nodeIndexLTB.push_back(gpuEl.ElementNodeIndices.z);

    for (int e = 0; e < 3; e++) {
      if (gpuEl.EdgeElements[e].NeighbourNodeIndex >= 0) {
	tmpNodeElementVertexMap[gpuEl.EdgeElements[e].NeighbourNodeIndex].push_back(3*2*(i_el - this->m_Elements.begin()) + 3 + (e + 2)%3);
	nodeIndexLTB.push_back(gpuEl.EdgeElements[e].NeighbourNodeIndex);
      }
    }
  }
  nodeIndexLTB = tledHelper::MakeSortedUnique(nodeIndexLTB);

  hostNodeElementVtxLTBBase.reserve(nodeIndexLTB.size());
  for (std::vector<int>::const_iterator ic_ind = nodeIndexLTB.begin(); ic_ind < nodeIndexLTB.end(); ic_ind++) {
    int2 base;

    assert(*ic_ind >= 0 && *ic_ind < mainSolver.GetMesh()->GetNumNodes());
    base.x = hostNodeElementVtxLTB.size();
    hostNodeElementVtxLTB.insert(hostNodeElementVtxLTB.end(), tmpNodeElementVertexMap[*ic_ind].begin(), tmpNodeElementVertexMap[*ic_ind].end());
    base.y = hostNodeElementVtxLTB.size();
    hostNodeElementVtxLTBBase.push_back(base);
  }

  if (hostMem.size() > 0) {
    tledCUDAHelpers::AllocateDeviceMemory(mdp_GPUElements, hostMem.size());
    tledCUDAHelpers::CopyToDevice(mdp_GPUElements, &hostMem.front(), hostMem.size());
    mdp_GPUMaterial = this->GetMaterial().InitGPU();
    
    tledCUDAHelpers::AllocateDeviceMemory(mdp_ShellElementForces, this->GetNumberOfElements()*3*2);
    
    m_NumberOfNodes = nodeIndexLTB.size();
    tledCUDAHelpers::AllocateDeviceMemory(mdp_NodeIndexLookupTable, nodeIndexLTB.size());
    tledCUDAHelpers::CopyToDevice(mdp_NodeIndexLookupTable, &nodeIndexLTB.front(), nodeIndexLTB.size());

    m_NodeElementVertexLookupTableSize = hostNodeElementVtxLTB.size();
    tledCUDAHelpers::AllocateDeviceMemory(mdp_NodeElementVertexLookupTable, hostNodeElementVtxLTB.size());
    tledCUDAHelpers::CopyToDevice(mdp_NodeElementVertexLookupTable, &hostNodeElementVtxLTB.front(), hostNodeElementVtxLTB.size());

    tledCUDAHelpers::AllocateDeviceMemory(mdp_NodeElementVertexLookupBaseIndex, hostNodeElementVtxLTBBase.size());
    tledCUDAHelpers::CopyToDevice(mdp_NodeElementVertexLookupBaseIndex, &hostNodeElementVtxLTBBase.front(), hostNodeElementVtxLTBBase.size());

#ifndef NDEBUG
    for (std::vector<int>::const_iterator ic_nInd = nodeIndexLTB.begin(); ic_nInd < nodeIndexLTB.end(); ic_nInd++) {
      for (std::vector<int>::const_iterator ic_vInd = hostNodeElementVtxLTB.begin() + hostNodeElementVtxLTBBase[ic_nInd-nodeIndexLTB.begin()].x; ic_vInd < hostNodeElementVtxLTB.begin() + hostNodeElementVtxLTBBase[ic_nInd-nodeIndexLTB.begin()].y; ic_vInd++) {
	const int elInd = *ic_vInd/6;
	const int lInd = *ic_vInd%6;

	if (lInd < 3) {
	  switch (lInd) {
	  case 0:
	    assert(hostMem[elInd].ElementNodeIndices.x == *ic_nInd);
	    break;

	  case 1:
	    assert(hostMem[elInd].ElementNodeIndices.y == *ic_nInd);
	    break;

	  default:
	    assert(hostMem[elInd].ElementNodeIndices.z == *ic_nInd);
	  }
	} else {
	  assert(hostMem[elInd].EdgeElements[(lInd-3+1)%3].NeighbourNodeIndex == *ic_nInd);
	}
      }
    }
#endif
  } else mdp_GPUElements = NULL; 
}

template <>
void tledShellSolver::ElementSetImpl<tledElementShellBSTP1, tledShellSolverGPU::ElementSet>::ClampNodes(const std::vector<bool> &clampedNodeMask) {
  for (std::vector<tledElementShellBSTP1>::iterator i_el = m_Elements.begin(); i_el < m_Elements.end(); i_el++) {
      i_el->ClampNodes(clampedNodeMask);
    }
}

void tledShellSolverGPU::Init(tledSolver &r_solver, const tledModel &model) {
  Superclass::Init(r_solver, model);

  for (std::vector<Superclass::ElementSet*>::iterator ip_elSet = mvp_ShellElementSets.begin(); ip_elSet < mvp_ShellElementSets.end(); ip_elSet++) {
    dynamic_cast<ElementSet*>(*ip_elSet)->InitGPU(*mp_MainSolver);    
  }
}

template <class TShellElement>
void tledShellSolverGPU::_ElementSetImpl<TShellElement>::FreeGPU() {
  if (mdp_GPUElements != NULL) {
    tledCheckCUDAErrors(cudaFree(mdp_GPUElements));
    tledCheckCUDAErrors(cudaFree(mdp_GPUMaterial));
    tledCheckCUDAErrors(cudaFree(mdp_ShellElementForces));
    tledCheckCUDAErrors(cudaFree(mdp_NodeIndexLookupTable));
    tledCheckCUDAErrors(cudaFree(mdp_NodeElementVertexLookupTable));
    tledCheckCUDAErrors(cudaFree(mdp_NodeElementVertexLookupBaseIndex));

    mdp_GPUElements = NULL;    
  }
}

template <class TShellElement>
void tledShellSolverGPU::_ElementSetImpl<TShellElement>::_InitDisplacementTextures() {
  tledCheckCUDAErrors(cudaBindTexture(0, tledShellSolver_kernels::tx_NodeElementVertexLookupBaseIndex, mdp_NodeElementVertexLookupBaseIndex, sizeof(int2)*m_NumberOfNodes));
  tledCheckCUDAErrors(cudaBindTexture(0, tledShellSolver_kernels::tx_NodeElementVertexLookupTable, mdp_NodeElementVertexLookupTable, sizeof(int)*m_NodeElementVertexLookupTableSize));
  tledCheckCUDAErrors(cudaBindTexture(0, tledShellSolver_kernels::tx_ElementVertexForces, mdp_ShellElementForces, sizeof(float4)*TShellElement::Facet::NumberOfVertices*this->GetNumberOfElements()));
  tledCheckCUDAErrors(cudaBindTexture(0, tledShellSolver_kernels::tx_NodeIndexLookupTable, mdp_NodeIndexLookupTable, sizeof(int)*m_NumberOfNodes));
}

template <>
void tledShellSolverGPU::_ElementSetImpl<tledElementShellBSTP1>::_InitDisplacementTextures() {
  tledCheckCUDAErrors(cudaBindTexture(0, tledShellSolver_kernels::tx_NodeElementVertexLookupBaseIndex, mdp_NodeElementVertexLookupBaseIndex, sizeof(int2)*m_NumberOfNodes));
  tledCheckCUDAErrors(cudaBindTexture(0, tledShellSolver_kernels::tx_NodeElementVertexLookupTable, mdp_NodeElementVertexLookupTable, sizeof(int)*m_NodeElementVertexLookupTableSize));
  tledCheckCUDAErrors(cudaBindTexture(0, tledShellSolver_kernels::tx_ElementVertexForces, mdp_ShellElementForces, sizeof(float4)*3*2*this->GetNumberOfElements()));
  tledCheckCUDAErrors(cudaBindTexture(0, tledShellSolver_kernels::tx_NodeIndexLookupTable, mdp_NodeIndexLookupTable, sizeof(int)*m_NumberOfNodes));
}

template <class TElementSet>
void _ComputeForce(float4 *gdp_f, const tledShellSolverGPU::ElementSet &elSet);

template <>
void tledShellSolverGPU::_ElementSetImpl<tledElementMembraneSimpleLinear<3> >::ComputeForces() {
  static const int szBlock = 128;
  const int szGrid = GetNumberOfElements()/szBlock + (GetNumberOfElements()%szBlock != 0);

  if (typeid(GetMaterial()) == typeid(tledMembraneMaterialLinear)) {    
    tledShellSolver_kernels::ComputeNewForces<tledElementMembraneSimpleLinear<3>, tledMembraneMaterialLinear> <<<szGrid, szBlock>>>(mdp_ShellElementForces, mdp_GPUElements, (const tledMembraneMaterialLinear::GPUMaterial*)mdp_GPUMaterial, this->GetNumberOfElements());
  } else {
    tledFatalError("Element/material combination not supported."); 
  }
}

template <>
void tledShellSolverGPU::_ElementSetImpl<tledElementShellBSTP1>::ComputeForces() {
  static const int szBlock = 128;
  const int szGrid = GetNumberOfElements()/szBlock + (GetNumberOfElements()%szBlock != 0);

  if (typeid(GetMaterial()) == typeid(tledShellMaterialLinearPlateDecorator<tledMembraneMaterialLinear>)) {    
    tledShellSolver_kernels::ComputeNewForces<tledElementShellBSTP1, tledShellMaterialLinearPlateDecorator<tledMembraneMaterialLinear> > <<<szGrid, szBlock>>>(mdp_ShellElementForces, mdp_GPUElements, (const tledShellMaterialLinearPlateDecorator<tledMembraneMaterialLinear>::GPUMaterial*)mdp_GPUMaterial, this->GetNumberOfElements());
  } else {
    tledFatalError("Element/material combination not supported."); 
  }
}

template <>
void tledShellSolverGPU::_ElementSetImpl<tledElementMembraneNonLinear<3> >::ComputeForces(void) {
  static const int szBlock = 128;
  const int szGrid = GetNumberOfElements()/szBlock + (GetNumberOfElements()%szBlock != 0);

  if (typeid(GetMaterial()) == typeid(tledMembraneMaterialNeoHookean)) {    
    tledShellSolver_kernels::ComputeNewForces<tledElementMembraneNonLinear<3>, tledMembraneMaterialNeoHookean> <<<szGrid, szBlock>>>(mdp_ShellElementForces, mdp_GPUElements, (const tledMembraneMaterialNeoHookean::GPUMaterial*)mdp_GPUMaterial, GetNumberOfElements());
  } else {
    tledFatalError("Element/material combination not supported."); 
  }
}

template <class TElement>
void tledShellSolverGPU::_ElementSetImpl<TElement>::ComputeForces(void) {
  tledFatalError("Not implemented.");
}

template <class TElement>
void tledShellSolverGPU::_ElementSetImpl<TElement>::ComputeNodalForces(float4 *dp_f) {
  static const int szBlock = 128;
  const int szGrid = m_NumberOfNodes/szBlock + (m_NumberOfNodes%szBlock != 0);

  _InitDisplacementTextures();
  tledShellSolver_kernels::ComputeNodalForces <<< szGrid, szBlock >>> (dp_f, m_NumberOfNodes);
}

void tledShellSolverGPU::ComputeNewForces() {
  for (std::vector<Superclass::ElementSet*>::iterator ip_elSet = mvp_ShellElementSets.begin(); ip_elSet < mvp_ShellElementSets.end(); ip_elSet++) {
    dynamic_cast<ElementSet*>(*ip_elSet)->ComputeForces();
  }
}

void tledShellSolverGPU::ComputeNodalForces(float4 *dp_f) {
  for (std::vector<Superclass::ElementSet*>::iterator ip_elSet = mvp_ShellElementSets.begin(); ip_elSet < mvp_ShellElementSets.end(); ip_elSet++) {
    dynamic_cast<ElementSet*>(*ip_elSet)->ComputeNodalForces(dp_f);
  }  
}

tledShellSolver::ElementSet* tledShellSolverGPU::CreateElementSetFromIndices3(tledShellMaterial &r_mat, const tledSurface &surface, const std::vector<int> &elInds) {
  const tledShellMesh<3> &mesh = dynamic_cast<const tledShellMesh<3>&>(surface);

  if (r_mat.HasBendingStiffness()) {
    return new _ElementSetImpl<tledElementShellBSTP1>(r_mat, mesh, elInds);
  } else {
    if (!r_mat.IsNonLinear()) return new _ElementSetImpl<tledElementMembraneSimpleLinear<3> >(r_mat, mesh, elInds);
    else return new _ElementSetImpl<tledElementMembraneNonLinear<3> >(r_mat, mesh, elInds);
  }
}

#endif
