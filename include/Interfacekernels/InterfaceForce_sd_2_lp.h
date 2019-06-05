
#ifndef INTERFACEFORCE_SD_2_LP_H
#define INTERFACEFORCE_SD_2_LP_H

#include "InterfaceKernel.h"
#include "PostprocessorInterface.h"
#include "VectorPostprocessorInterface.h"

class InterfaceForce_sd_2_lp;

template<>
InputParameters validParams<InterfaceForce_sd_2_lp>();

class InterfaceForce_sd_2_lp :
  public InterfaceKernel,
  public PostprocessorInterface,
  protected VectorPostprocessorInterface
{
public:
  InterfaceForce_sd_2_lp(const InputParameters & parameters);
  virtual const PostprocessorValue & getPostprocessorValue(const std::string & name);
  virtual const PostprocessorValue & getPostprocessorValueByName(const PostprocessorName & name);


protected:
  virtual Real computeQpResidual(Moose::DGResidualType type);
  virtual Real computeQpJacobian(Moose::DGJacobianType type);

  /// Depend UserObjects
  std::set<std::string> _depend_uo;
  const Real _unit_scalor;
  const Real _stop_threshold;
  const MaterialProperty<Real> & _driving_rate;
  const VariableValue & _coupled_s_p;
  // const Real _driving_rate;
  // const MaterialProperty<Real> & _a_0;
  // const Real _a_0;
  const MaterialProperty<Real> & _solubility_in_solid;
  const MaterialProperty<Real> & _diffusivity_in_solid_sd;
  const MaterialProperty<Real> & _diffusivity_in_solid_lp;
  // const Real & _Ave_s_d;
  // const Real & _Ave_l_p;
};

#endif
