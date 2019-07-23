#pragma once

#include "InterfaceKernel.h"

class InterfaceForce_sd_2_ld_flux;

template<>
InputParameters validParams<InterfaceForce_sd_2_ld_flux>();

class InterfaceForce_sd_2_ld_flux :
  public InterfaceKernel
{
public:
  InterfaceForce_sd_2_ld_flux(const InputParameters & parameters);

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
