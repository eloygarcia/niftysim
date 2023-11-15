// =========================================================================
// File:       tledShellSolver.h
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    February 2012
// 
// Copyright (c) 2012, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================
#ifndef tledShellSolver_H
#define tledShellSolver_H

#include "tledModel.h"
#include "tledHelper.h"
#include "tledElementMembrane.h"
#include "tledSolver.h"

#include <vector>
#include <string>

/**
 * \brief tledSolver backend for shell/membrane problems
 * \ingroup shell
 * \ingroup solver
 */
class tledShellSolver {
  /**
   * \name Solution Computation
   * @{
   */
private:
  std::vector<float> m_MShell;

public:
  /** Mass vector for shell elements, has same dimensions as the tledSolver one (1 component/solid mesh node) */
  const std::vector<float>& GetMass(void) const { return m_MShell; }
  /** @} */

  /**
   * \name Setup 
   * @{
   */
public:
  class ElementSet {
  private:
    tledShellMaterial *mp_Material;
    std::vector<int> m_ElementIndices;

  public:
    virtual tledElementMembrane& GetElement(const int elInd) = 0;
    virtual const tledElementMembrane& GetElement(const int elInd) const = 0;
    virtual int GetNumberOfElements(void) const = 0;
    const tledShellMaterial& GetMaterial(void) const { return *mp_Material; }
    const std::vector<int>& GetElementIndices(void) const { return m_ElementIndices; }

  public:
    virtual void ComputeElementThicknesses(float *p_ts, const float U[]) const = 0;
    //virtual void ComputeElementNormal(float *p_ns, const float U[]) const = 0;

  public:
    virtual void ComputeMass(float *p_mass) = 0;
    virtual void ClampNodes(const std::vector<bool> &clampedNodeMask) = 0;

  public:
    ElementSet(tledShellMaterial &r_mat, const std::vector<int> &elInds) : mp_Material(&r_mat), m_ElementIndices(elInds) {}
    virtual ~ElementSet(void) { delete mp_Material; }
  };

protected:
  template <class TShellElement, class TElementSetInterface>
  class ElementSetImpl : public TElementSetInterface {
  protected:
    std::vector<TShellElement> m_Elements;

  public:
    virtual tledElementMembrane& GetElement(const int elInd) { return m_Elements[elInd]; }
    virtual const tledElementMembrane& GetElement(const int elInd) const { return m_Elements[elInd]; }
    virtual int GetNumberOfElements(void) const { return m_Elements.size(); }

  public:
    virtual void ComputeElementThicknesses(float *p_ts, const float U[]) const;

  public:
    virtual void ComputeMass(float *p_mass);
    virtual void ClampNodes(const std::vector<bool> &clampedNodeMask) { /* Only relevant for plate elements; default: do nothing */ }

  protected:
    ElementSetImpl(tledShellMaterial &r_mat) : TElementSetInterface(r_mat) {}

  public:
    /* Work-around for bizarre MSVC++ 2010 bug where constructor gets defined multiple times when contained in tpp file. */
    ElementSetImpl(tledShellMaterial &r_mat, const typename TShellElement::Surface &surface, const std::vector<int> &elInds) : TElementSetInterface(r_mat, elInds) {
      std::vector<int> actualElInds = elInds;

      if (actualElInds.size() == 1 && actualElInds.front() == -1) actualElInds = tledSequenceGenerator::MakeSequence(0, surface.GetNumberOfFacets());
      m_Elements.reserve(actualElInds.size());
      for (std::vector<int>::const_iterator ic_elInd = actualElInds.begin(); ic_elInd < actualElInds.end(); ic_elInd++) {
	TShellElement el;
    
	el.InitialiseElement(surface, *ic_elInd);
	el.SetMaterial(r_mat);
	m_Elements.push_back(el);
      }
    }

    virtual ~ElementSetImpl() {}
  };

protected:
  tledSurface *mp_Surface;
  std::vector<ElementSet*> mvp_ShellElementSets;
  tledSolver *mp_MainSolver;
  float m_Dt;
  std::vector<std::pair<int, int> > m_ElementElementSetIndices;

protected:
  /** Element set (triangular elements) factory, CPU/GPU implementations have to provide functions instantiating element sets suitable for their respective platforms. */
  virtual tledShellSolver::ElementSet* CreateElementSetFromIndices3(tledShellMaterial &r_material, const tledSurface &surface, const std::vector<int> &elSetInds) = 0;

  ElementSet* CreateElementSet(const tledSurface &surface, const tledModel &model, const int elSetInd);

  template <class TSurface>
  TSurface& GetSurface(void) { return *dynamic_cast<TSurface*>(mp_Surface); }

  template <class TSurface>
  const TSurface& GetSurface(void) const { return *dynamic_cast<const TSurface*>(mp_Surface); }

  float GetTimeStep(void) const { return m_Dt; }

public:
  /** Does the precomputation, has to be called before first call to GetMass or any solution computation member functions. */
  virtual void Init(tledSolver &r_solver, const tledModel &model);  

  int GetNumberOfElementSets(void) const { return mvp_ShellElementSets.size(); }
  int GetNumberOfElementSetElements(const int elSetIndex) const { return mvp_ShellElementSets[elSetIndex]->GetNumberOfElements(); }

  /** Direct access to element set, only intended for debugging and testing. */
  ElementSet& GetElementSet(const int elSetIndex) { return *mvp_ShellElementSets[elSetIndex]; }

  tledSurface& GetSurface(void) { return GetSurface<tledSurface>(); }
  const tledSurface& GetSurface(void) const { return GetSurface<tledSurface>(); }

  /** 
   * \brief Returns the element-set index (first component) and index of the element in the element-set (second component) of a membrane mesh element.
   *
   * Elements that are not contained in any element set can be identified by a (-1, -1) entry.
   */
  const std::pair<int, int>& GetElementSetIndexForElement(const int elIndex) const { return m_ElementElementSetIndices[elIndex]; }

  /** 
   * \brief Returns a computation element for a mesh element index. 
   *
   * A NULL reference is returned if the element is not part of any element set. 
   */
  const tledElementMembrane& GetElement(const int elIndex) const;
  /** @} */

  /**
   * \name Post-Processing
   * @{
   */
public:
  /** Computes the deformed-state element thickness for an element set. */
  void ComputeElementThicknesses(float *p_dst, const int elSetInd, const float U[]) const { mvp_ShellElementSets[elSetInd]->ComputeElementThicknesses(p_dst, U); }

  /** Computes the thickness of all elements in all element sets, indexing based on mesh element indices */
  virtual void ComputeAllElementThicknesses(float *p_dst, const float U[]) const;

  /** Copies the appropriate initial thickness values for all elements to a buffer */
  void GetAllInitialElementThicknesses(float *p_dst) const;
  /** @} */

  /**
   * \name Construction, Destruction
   * @{
   */
public:
  tledShellSolver(void);
  virtual ~tledShellSolver(void);
  /** @} */
};

#include "tledShellSolver.tpp"
#endif
