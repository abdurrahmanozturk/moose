//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef ADKERNELSUPG_H
#define ADKERNELSUPG_H

#include "ADKernel.h"

#define usingTemplKernelSUPGMembers(type)                                                          \
  usingTemplKernelMembers(type);                                                                   \
  using ADKernelSUPGTempl<type, compute_stage>::_velocity;                                         \
  using ADKernelSUPGTempl<type, compute_stage>::_tau
#define usingKernelSUPGMembers usingTemplKernelSUPGMembers(Real)
#define usingVectorKernelSUPGMembers usingTemplKernelSUPGMembers(RealVectorValue)

template <typename, ComputeStage>
class ADKernelSUPGTempl;

template <ComputeStage compute_stage>
using ADKernelSUPG = ADKernelSUPGTempl<Real, compute_stage>;
template <ComputeStage compute_stage>
using ADVectorKernelSUPG = ADKernelSUPGTempl<RealVectorValue, compute_stage>;

declareADValidParams(ADKernelSUPG);
declareADValidParams(ADVectorKernelSUPG);

template <typename T, ComputeStage compute_stage>
class ADKernelSUPGTempl : public ADKernelTempl<T, compute_stage>
{
public:
  ADKernelSUPGTempl(const InputParameters & parameters);

  virtual void computeResidual() override;
  virtual void computeJacobian() override;
  virtual void computeADOffDiagJacobian() override;

protected:
  /**
   * Called before forming the residual for an element
   */
  virtual typename OutputTools<typename Moose::ValueType<T, compute_stage>::type>::OutputValue
  precomputeQpStrongResidual() = 0;

  virtual ADResidual computeQpResidual() override final;

  const ADMaterialProperty(Real) & _tau;
  const ADVectorVariableValue & _velocity;

  usingTemplKernelMembers(T);
};

#endif /* ADKERNELSUPG_H */
