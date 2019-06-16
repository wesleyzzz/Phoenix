// Describe s_d --> l_d, coupled neighbor _coupled_s_p, consider flux term, neglect metal reaction, no flux term (J=D*grad(u))

#include "InterfaceForce_sd_2_ld_Noflux.h"
#include <cmath>

registerMooseObject("PhoenixApp", InterfaceForce_sd_2_ld_Noflux);

template<>
InputParameters validParams<InterfaceForce_sd_2_ld_Noflux>()
{
  InputParameters params = validParams<InterfaceKernel>();
  params.addParam<MaterialPropertyName>("driving_rate", 1.0, "driving force rate between the interface");
  params.addCoupledVar("coupled_l_p", "Ln precipitation in liquid");
  params.addCoupledVar("coupled_s_p", "Ln precipitation in solid from neighbor domain");
  params.addParam<MaterialPropertyName>("solubility_in_liquid", 0.0, "Name of solubility material in the liquid");
  params.addParam<MaterialPropertyName>("solubility_in_solid", 0.0, "Name of solubility material in the solid");
  params.addParam<MaterialPropertyName>("diffusivity_in_liquid", 0.0, "Name of diffusivity material in the liquid");
  params.addParam<MaterialPropertyName>("diffusivity_in_solid", 0.0, "Name of diffusivity material in the solid");
  return params;
}

InterfaceForce_sd_2_ld_Noflux::InterfaceForce_sd_2_ld_Noflux(const InputParameters & parameters) :
    InterfaceKernel(parameters),
    PostprocessorInterface(this),
    VectorPostprocessorInterface(this),
    _driving_rate(getMaterialProperty<Real>(getParam<MaterialPropertyName>("driving_rate"))),
    _solubility_in_liquid(getMaterialProperty<Real>(getParam<MaterialPropertyName>("solubility_in_liquid"))),
    _solubility_in_solid(getMaterialProperty<Real>(getParam<MaterialPropertyName>("solubility_in_solid"))),
    _diffusivity_in_liquid(getMaterialProperty<Real>(getParam<MaterialPropertyName>("diffusivity_in_liquid"))),
    _diffusivity_in_solid(getMaterialProperty<Real>(getParam<MaterialPropertyName>("diffusivity_in_solid")))
//why not getNeighborMaterialProperty

{
  if (!parameters.isParamValid("boundary"))
  {
    mooseError("In order to use the InterfaceForce_ld dgkernel, you must specify a boundary where it will live.");
  }
}

const PostprocessorValue &
InterfaceForce_sd_2_ld_Noflux::getPostprocessorValue(const std::string & name)
{
  _depend_uo.insert(_pars.get<PostprocessorName>(name));
  return PostprocessorInterface::getPostprocessorValue(name);
}

const PostprocessorValue &
InterfaceForce_sd_2_ld_Noflux::getPostprocessorValueByName(const PostprocessorName & name)
{
  _depend_uo.insert(name);
  return PostprocessorInterface::getPostprocessorValueByName(name);
}

Real
InterfaceForce_sd_2_ld_Noflux::computeQpResidual(Moose::DGResidualType type)
{
  Real r = 0.0;

  switch (type)
  {
  case Moose::Element:
    r = -_test[_i][_qp]* (_driving_rate[_qp] * _solubility_in_solid[_qp] * ( _neighbor_value[_qp]/_solubility_in_solid[_qp] - _u[_qp]/_solubility_in_liquid[_qp]));
    break;

  case Moose::Neighbor:
    r = _test_neighbor[_i][_qp]* (_driving_rate[_qp] * _solubility_in_solid[_qp] * ( _neighbor_value[_qp]/_solubility_in_solid[_qp] - _u[_qp]/_solubility_in_liquid[_qp]));
    break;
  }

  return r;
}

Real
InterfaceForce_sd_2_ld_Noflux::computeQpJacobian(Moose::DGJacobianType type)
{
  Real jac = 0.0;

  switch (type)
  {

    case Moose::ElementElement:
      jac = -_test[_i][_qp] * (_driving_rate[_qp] * _solubility_in_solid[_qp] * ( - _phi[_j][_qp]/_solubility_in_liquid[_qp]));
      break;

    case Moose::NeighborNeighbor:
      jac = _test_neighbor[_i][_qp]* (_driving_rate[_qp] * _solubility_in_solid[_qp] * ( _phi_neighbor[_j][_qp]/_solubility_in_solid[_qp]));;
      break;

    case Moose::NeighborElement:
      jac = _test_neighbor[_i][_qp]* (_driving_rate[_qp] * _solubility_in_solid[_qp] * ( - _phi[_j][_qp]/_solubility_in_liquid[_qp]));
      break;

    case Moose::ElementNeighbor:
      jac = -_test[_i][_qp]* (_driving_rate[_qp] * _solubility_in_solid[_qp] * ( _phi_neighbor[_j][_qp]/_solubility_in_solid[_qp]));
      break;
  }

  return jac;
}
