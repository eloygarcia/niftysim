// =========================================================================
// File:       tledSubModelManager.cpp
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    July 2011
// 
// Copyright (c) 2011, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================
#include "tledSubModelManager.h"
#include "tledVectorArithmetic.h"
#include "tledSimulator.h"

#include <cassert>
#include <sstream>
#include <iterator>
#include <algorithm>
#include <typeinfo>

using namespace tledVectorArithmetic;
using namespace std;

void tledSubModelManager::_AddSubMesh(const tledMesh &mesh) {
  vector<int> nodeIndMap;

  if (*GetMeshType() == 0) {
    m_MeshType = mesh.GetElType();
    m_NumElementNodes = mesh.GetNodesPerEl();
  }

  if (m_MeshType != mesh.GetElType()) {
    tledLogErrorStream(tledHelper::FatalError() << " sub-mesh type " << mesh.GetElType() << " incompatible with model mesh type " << GetMeshType());
  }

  if (GetMinNodeDistance() >= 0 && m_Nodes.size() > 0) {
    for (float const *pc_meshNode = mesh.GetAllNodeCds(); pc_meshNode < mesh.GetAllNodeCds() + 3*mesh.GetNumNodes(); pc_meshNode += 3) {
      float tmp[3];
      vector<float>::const_iterator ic_nodeCds;

      for (ic_nodeCds = m_Nodes.begin(); ic_nodeCds < m_Nodes.end() && Norm(Sub(tmp, &(*ic_nodeCds), pc_meshNode)) >= GetMinNodeDistance(); ic_nodeCds += 3);
      assert((ic_nodeCds - m_Nodes.begin())%3 == 0);
      if (ic_nodeCds < m_Nodes.end()) {
	assert(Norm(Sub(tmp, &(*ic_nodeCds), pc_meshNode)) < GetMinNodeDistance());
	nodeIndMap.push_back((ic_nodeCds - m_Nodes.begin())/3);
      } else {
	nodeIndMap.push_back(m_Nodes.size()/3);
	m_Nodes.insert(m_Nodes.end(), pc_meshNode, pc_meshNode + 3);
      }
    } /* for mesh points */
    assert((int)nodeIndMap.size() == mesh.GetNumNodes());
  } else {
    const int startInd = m_Nodes.size()/3;

    vector<int>::iterator i_nodeInd;

    m_Nodes.insert(m_Nodes.end(), mesh.GetAllNodeCds(), mesh.GetAllNodeCds() + 3*mesh.GetNumNodes());
    nodeIndMap.resize(mesh.GetNumNodes());
    for (i_nodeInd = nodeIndMap.begin(); i_nodeInd < nodeIndMap.end(); i_nodeInd++) *i_nodeInd = startInd + (i_nodeInd - nodeIndMap.begin());
  } /* if do node merging else ... */

  {
    SubMesh newSubMesh;
    int const *pc_elNodeInd;    

    newSubMesh.ElementStartIndex = m_Elements.size()/m_NumElementNodes;
    newSubMesh.ElementEndIndex = newSubMesh.ElementStartIndex + mesh.GetNumEls();
    newSubMesh.MembraneElementStartIndex = newSubMesh.MembraneElementEndIndex = -1;
    m_Elements.reserve(m_Elements.size() + m_NumElementNodes*mesh.GetNumEls());
    for (pc_elNodeInd = mesh.GetAllElNodeInds(); pc_elNodeInd < mesh.GetAllElNodeInds() + m_NumElementNodes*mesh.GetNumEls(); pc_elNodeInd++) {
      assert((int)nodeIndMap.size() > *pc_elNodeInd && *pc_elNodeInd >= 0);
      m_Elements.push_back(nodeIndMap[*pc_elNodeInd]);
    }
    newSubMesh.SubMeshToMeshNodeIndexMap = nodeIndMap;
    m_SubMeshes.push_back(newSubMesh);
  }
}

void tledSubModelManager::_AddSubMembrane(const tledSurface &membrane, const bool isSurface) {
  SubMesh &r_submesh = m_SubMeshes.back();

  if (this->m_MembraneElements.size() > 0) {
    if (isSurface != this->DoShellUseMeshSurface()) {
      tledFatalError("Either all sub-membranes are of type SURFACE or none.");
    }
  }

  if (membrane.GetNumberOfFacetVertices() == 3) {
    if (m_ShellElementType.length() > 0 && m_ShellElementType != "T3") {
      tledFatalError("Membrane meshes in sub-model cannot be merged due to element type incompatibility.");
    } else {
      m_ShellElementType = "T3";      
    }
  } else {
    tledFatalNotYetImplementedError;
  }

  r_submesh.MembraneElementStartIndex = m_MembraneElements.size();
  m_MembraneElements.reserve(m_MembraneElements.size() + membrane.GetNumberOfFacets()*membrane.GetNumberOfFacetVertices());
  for (int fInd = 0; fInd < membrane.GetNumberOfFacets(); fInd++) {
    for (int const *pc_vInd = membrane.GetFacetNodeIndices(fInd); pc_vInd < membrane.GetFacetNodeIndices(fInd) + membrane.GetNumberOfFacetVertices(); pc_vInd++) m_MembraneElements.push_back(r_submesh.SubMeshToMeshNodeIndexMap[*pc_vInd]);
  }
  r_submesh.MembraneElementEndIndex = m_MembraneElements.size();
  m_IsSurface = isSurface;
}

XMLNode tledSubModelManager::GetSubModelShellElementSet(const int smIndex, const int elSetIndex) const {
  XMLNode smXMLElSet;
  vector<int> smElSet;
  vector<int>::iterator i_elInd;
  ostringstream oss;

  smXMLElSet = m_SubXMLModels[smIndex].getChildNode("ShellElementSet", elSetIndex).deepCopy();

  smElSet = m_SubModels[smIndex].GetShellElementSet(elSetIndex);
  assert(m_SubMeshes[smIndex].MembraneElementStartIndex >= 0);
  for (i_elInd = smElSet.begin(); i_elInd < smElSet.end(); i_elInd++) *i_elInd += m_SubMeshes[smIndex].MembraneElementStartIndex;

  oss << smElSet.size();
  smXMLElSet.updateAttribute(oss.str().c_str(), "Size", "Size");
  
  if (smElSet.size() > 0) {
    oss.str(""), oss.clear();
    copy(smElSet.begin(), smElSet.end() - 1, ostream_iterator<int>(oss, "\n")), oss << smElSet.back();  
    smXMLElSet.updateText(oss.str().c_str());
  }

  assert(equal(smElSet.begin(), smElSet.end(), GetXMLTextAsVector<int>(smXMLElSet).begin()));

  return smXMLElSet;
}

XMLNode tledSubModelManager::GetSubModelElementSet(const int smIndex, const int elSetIndex) const {
  XMLNode smXMLElSet;
  vector<int> smElSet;
  vector<int>::iterator i_elInd;
  ostringstream oss;

  smXMLElSet = m_SubXMLModels[smIndex].getChildNode("ElementSet", elSetIndex).deepCopy();

  smElSet = m_SubModels[smIndex].GetElSet(elSetIndex);
  for (i_elInd = smElSet.begin(); i_elInd < smElSet.end(); i_elInd++) *i_elInd += GetSubMeshElementStartIndex(smIndex);

  oss << smElSet.size();
  smXMLElSet.updateAttribute(oss.str().c_str(), "Size", "Size");
  
  if (smElSet.size() > 0) {
    oss.str(""), oss.clear();
    copy(smElSet.begin(), smElSet.end() - 1, ostream_iterator<int>(oss, "\n")), oss << smElSet.back();  
    smXMLElSet.updateText(oss.str().c_str());
  }

  assert(equal(smElSet.begin(), smElSet.end(), GetXMLTextAsVector<int>(smXMLElSet).begin()));

  return smXMLElSet;
}

void tledSubModelManager::SetNumberOfSubModels(const int numSubModels) {
  m_SubModels.reserve(numSubModels);
  m_SubXMLModels.reserve(numSubModels);
  m_SubMeshes.reserve(numSubModels);
}

void tledSubModelManager::AddSubModel(const XMLNode &subModelRoot) {
  m_SubModels.push_back(tledModel());
  m_SubModels.back().LoadFromXML(subModelRoot);
  m_SubXMLModels.push_back(subModelRoot.deepCopy());
  _AddSubMesh(*m_SubModels.back().GetMesh());

  if (m_SubModels.back().GetNumberOfShellElementSets() > 0) {
    tledSurface *p_shellMesh = m_SubModels.back().GetGenericShellMesh();

    _AddSubMembrane(*p_shellMesh, m_SubModels.back().DoShellUseMeshSurface());

    delete p_shellMesh;
  }  
}

XMLNode tledSubModelManager::CreateShellMeshXML(void) const {
  const int numFctVtcs = m_ShellElementType == "T3"? 3 : 4;

  XMLNode meshXML = XMLNode::createXMLTopNode("ShellElements");
  ostringstream oss;

  meshXML.addAttribute("Type", m_ShellElementType.c_str());

  oss << m_MembraneElements.size()/numFctVtcs;
  meshXML.addAttribute("NumEls", oss.str().c_str());

  oss.str(""), oss.clear();
  copy(m_MembraneElements.begin(), m_MembraneElements.end() - 1, ostream_iterator<int>(oss, "\n")), oss << m_MembraneElements.back();
  meshXML.updateText(oss.str().c_str());

  return meshXML;
}

void tledSubModelManager::WriteMesh(tledMesh &r_mesh) {
  assert(m_SubMeshes.size() > 0);
  r_mesh.SetNumberOfElements(m_Elements.size()/m_NumElementNodes, GetMeshType());
  r_mesh.SetNumberOfNodes(m_Nodes.size()/3, 3);
  
  copy(m_Elements.begin(), m_Elements.end(), r_mesh.GetAllElNodeInds());
  copy(m_Nodes.begin(), m_Nodes.end(), r_mesh.GetAllNodeCds());
}

void tledSubModelManager::ExportSubMeshNodeAttributes(float *p_dst, const float attributes[], const int subMeshInd) const {
  for (std::vector<int>::const_iterator ic_gn = this->GetSubMeshToMeshNodeIndexMap(subMeshInd).begin(); ic_gn < this->GetSubMeshToMeshNodeIndexMap(subMeshInd).end(); ic_gn++, p_dst += 3) {
    std::copy(attributes + 3*(*ic_gn), attributes + 3*(*ic_gn + 1), p_dst);
  }
}

tledMesh* tledSubModelManager::ExportSubMesh(const int subMeshInd, const tledMesh &globalMesh) const {
  tledMesh* p_sub = new tledMesh();
  std::vector<int> reverseNodeMap(globalMesh.GetNumNodes()*3, -1);

  for (int sn = 0; sn < this->GetNumberOfSubMeshNodes(subMeshInd); sn++) reverseNodeMap[this->GetSubMeshToMeshNodeIndexMap(subMeshInd)[sn]] = sn;

  p_sub->SetNumberOfElements(this->GetNumberOfSubMeshElements(subMeshInd), globalMesh.GetElType());
  p_sub->SetNumberOfNodes(this->GetNumberOfSubMeshNodes(subMeshInd), 3);

  this->ExportSubMeshNodeAttributes(p_sub->GetAllNodeCds(), globalMesh.GetAllNodeCds(), subMeshInd);

  for (int e = 0; e < this->GetNumberOfSubMeshElements(subMeshInd); e++) {
    for (int n = 0; n < p_sub->GetNodesPerEl(); n++) {
      p_sub->GetAllElNodeInds()[e*p_sub->GetNodesPerEl()+n] = reverseNodeMap[globalMesh.GetAllElNodeInds()[(this->GetSubMeshElementStartIndex(subMeshInd)+e)*p_sub->GetNodesPerEl()+n]];
    }
  }

  return p_sub;
}

XMLNode tledSubModelManager::GetSubModelConstraint(const int smIndex, const int constraintIndex) const {
  string constType;
  XMLNode xmlConstraint;

  /* 
   * Attn: using fact that xmlParser only does shallow copy.
   * Ideally we should delete "Nodes" child and replace it with a new one 
   */
  xmlConstraint = m_SubXMLModels[smIndex].getChildNode("Constraint", constraintIndex).deepCopy();
  constType = xmlConstraint.getAttribute("Type");
  if (constType == "Pressure" || constType == "Traction") {
    ostringstream oss;
    vector<int> smFaces;
    int numFacetVtcs;

    if (constType == "Pressure") smFaces = m_SubModels[smIndex].GetPressureFaceNodeInds(constraintIndex);  
    else smFaces = m_SubModels[smIndex].GetTractionFaceNodeInds(constraintIndex);  

    for (std::vector<int>::iterator i_smNodeInd = smFaces.begin(); i_smNodeInd < smFaces.end(); i_smNodeInd++) {
      *i_smNodeInd = m_SubMeshes[smIndex].SubMeshToMeshNodeIndexMap[*i_smNodeInd];
      assert(*i_smNodeInd >= 0 && *i_smNodeInd < (int)m_Nodes.size()/3);
    }

    if (constType == "Pressure") {
      numFacetVtcs = std::string(m_SubModels[smIndex].GetPressureFaceType(constraintIndex)) == "Tri"? 3 : 4;
      xmlConstraint.updateAttribute(m_SubModels[smIndex].GetPressureFaceType(constraintIndex), "FaceType", "FaceType");
    } else {
      numFacetVtcs = std::string(m_SubModels[smIndex].GetTractionFaceType(constraintIndex)) == "Tri"? 3 : 4;
      xmlConstraint.updateAttribute(m_SubModels[smIndex].GetTractionFaceType(constraintIndex), "FaceType", "FaceType");
    }

    xmlConstraint.updateAttribute("FACES", "SpecType", "SpecType");    
    oss << smFaces.size()/numFacetVtcs;
    xmlConstraint.updateAttribute(oss.str().c_str(), "NumFaces", "NumFaces");

    oss.str(""), oss.clear();
    copy(smFaces.begin(), smFaces.end() - 1, ostream_iterator<int>(oss, "\n")), oss << smFaces.back();
    if (xmlConstraint.nChildNode("Faces") == 0) xmlConstraint.addChild("Faces");
    xmlConstraint.getChildNode("Faces").updateText(oss.str().c_str());
    assert(xmlConstraint.nChildNode("Faces") == 0 || equal(smFaces.begin(), smFaces.end(), GetXMLTextAsVector<int>(xmlConstraint.getChildNode("Faces")).begin()));
  } else {
    ostringstream oss;
    vector<int> smNodeInds;

    smNodeInds = m_SubModels[smIndex].GetConstraintInd(constraintIndex);  
    for (vector<int>::iterator i_smNodeInd = smNodeInds.begin(); i_smNodeInd < smNodeInds.end(); i_smNodeInd++) {
      *i_smNodeInd = m_SubMeshes[smIndex].SubMeshToMeshNodeIndexMap[*i_smNodeInd];
      assert(*i_smNodeInd >= 0 && *i_smNodeInd < (int)m_Nodes.size()/3);
    }

    if (smNodeInds.size() == 0) {
      tledLogErrorStream(tledHelper::FatalError() << "Boundary " << constraintIndex << " of sub-model " << smIndex << " has no nodes. Wrong normal?");
    }

    oss << smNodeInds.size();
    xmlConstraint.updateAttribute("NODES", "SpecType", "SpecType");
    xmlConstraint.updateAttribute(oss.str().c_str(), "NumNodes", "NumNodes");

    oss.str(""), oss.clear();
    copy(smNodeInds.begin(), smNodeInds.end() - 1, ostream_iterator<int>(oss, "\n")), oss << smNodeInds.back();
    if (xmlConstraint.nChildNode("Nodes") == 0) xmlConstraint.addChild("Nodes");
    xmlConstraint.getChildNode("Nodes").updateText(oss.str().c_str());
    assert(xmlConstraint.nChildNode("Nodes") == 0 || equal(smNodeInds.begin(), smNodeInds.end(), GetXMLTextAsVector<int>(xmlConstraint.getChildNode("Nodes")).begin()));
  }

  return xmlConstraint;
}

tledSubModelManager::~tledSubModelManager() {  
  /* sub-models destroy their respective XML trees */
}
