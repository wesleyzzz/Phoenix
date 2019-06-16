#pragma once

#include "InterfaceKernel.h"
#include "PostprocessorInterface.h"
#include "VectorPostprocessorInterface.h"

class InterfaceForce_sd_2_ld_Noflux;

template<>
InputParameters validParams<InterfaceForce_sd_2_ld_Noflux>();

class InterfaceForce_sd_2_ld_Noflux :
  public InterfaceKernel,
  public PostprocessorInterface,
  protected VectorPostprocessorInterface
{
public:
  InterfaceForce_sd_2_ld_Noflux(const InputParameters & parameters);
  virtual const PostprocessorValue & getPostprocessorValue(const std::string & name);
  virtual const PostprocessorValue & getPostprocessorValueByName(const PostprocessorName & name);


protected:
  virtual Real computeQpResidual(Moose::DGResidualType type);
  virtual Real computeQpJacobian(Moose::DGJacobianType type);

  /// Depend UserObjects
  std::set<std::string> _depend_uo;

  const MaterialProperty<Real> & _driving_rate;

  const MaterialProperty<Real> & _solubility_in_liquid;
  const MaterialProperty<Real> & _solubility_in_solid;

  const MaterialProperty<Real> & _diffusivity_in_liquid;
  const MaterialProperty<Real> & _diffusivity_in_solid;

};
