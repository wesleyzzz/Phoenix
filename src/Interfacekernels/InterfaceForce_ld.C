// Describe s_d --> l_d

#include "InterfaceForce_ld.h"
#include <cmath>

registerMooseObject("PhoenixApp", InterfaceForce_ld);

template<>
InputParameters validParams<InterfaceForce_ld>()
{
  InputParameters params = validParams<InterfaceKernel>();
  params.addRequiredParam<Real>("unit_scalor", "unit_scalor for the precipitate");
  params.addParam<MaterialPropertyName>("driving_rate", 1.0, "driving force rate between the interface");
  params.addCoupledVar("coupled_l_p", "Ln precipitation in liquid");
  params.addCoupledVar("coupled_s_p", "Ln precipitation in solid");
  params.addParam<MaterialPropertyName>("solubility_in_liquid", 0.0, "Name of solubility material in the liquid");
  params.addParam<MaterialPropertyName>("solubility_in_solid", 0.0, "Name of solubility material in the solid");
  params.addParam<MaterialPropertyName>("diffusivity_in_liquid", 0.0, "Name of diffusivity material in the liquid");
  params.addParam<MaterialPropertyName>("diffusivity_in_solid", 0.0, "Name of diffusivity material in the solid");
  return params;
}

InterfaceForce_ld::InterfaceForce_ld(const InputParameters & parameters) :
    InterfaceKernel(parameters),
    PostprocessorInterface(this),
    VectorPostprocessorInterface(this),
    _unit_scalor(getParam<Real>("unit_scalor")),
    _driving_rate(getMaterialProperty<Real>(getParam<MaterialPropertyName>("driving_rate"))),
    _coupled_l_p(coupledValue("coupled_l_p")),
    _coupled_s_p(coupledValue("coupled_s_p")),
    _solubility_in_liquid(getMaterialProperty<Real>(getParam<MaterialPropertyName>("solubility_in_liquid"))),
    _solubility_in_solid(getMaterialProperty<Real>(getParam<MaterialPropertyName>("solubility_in_solid"))),
    _diffusivity_in_liquid(getMaterialProperty<Real>(getParam<MaterialPropertyName>("diffusivity_in_liquid"))),
    _diffusivity_in_solid(getMaterialProperty<Real>(getParam<MaterialPropertyName>("diffusivity_in_solid")))

{
  if (!parameters.isParamValid("boundary"))
  {
    mooseError("In order to use the InterfaceForce_ld dgkernel, you must specify a boundary where it will live.");
  }
}

const PostprocessorValue &
InterfaceForce_ld::getPostprocessorValue(const std::string & name)
{
  _depend_uo.insert(_pars.get<PostprocessorName>(name));
  return PostprocessorInterface::getPostprocessorValue(name);
}

const PostprocessorValue &
InterfaceForce_ld::getPostprocessorValueByName(const PostprocessorName & name)
{
  _depend_uo.insert(name);
  return PostprocessorInterface::getPostprocessorValueByName(name);
}

Real
InterfaceForce_ld::computeQpResidual(Moose::DGResidualType type)
{
  Real r = 0.0;
  Real Diff_eff = _diffusivity_in_solid[_qp];
  Real Force_eff = _driving_rate[_qp] * _solubility_in_solid[_qp] * (_neighbor_value[_qp] / _solubility_in_solid[_qp] - _u[_qp] / _solubility_in_liquid[_qp]);

  r = Diff_eff * Force_eff * (1.0 - _coupled_l_p[_qp] / _unit_scalor) * (1.0 - _coupled_s_p[_qp] / _unit_scalor);

  switch (type)
  {
  case Moose::Element:
    r *= -_test[_i][_qp];
    break;

  case Moose::Neighbor:
    r *= _test_neighbor[_i][_qp];
    break;
  }

  return r;
}

Real
InterfaceForce_ld::computeQpJacobian(Moose::DGJacobianType type)
{
  Real jac = 0.0;
  Real Diff_eff = _diffusivity_in_solid[_qp];
  Real prefactor = Diff_eff * _driving_rate[_qp] * _solubility_in_solid[_qp] * (1.0 - _coupled_l_p[_qp] / _unit_scalor) * (1.0 - _coupled_s_p[_qp] / _unit_scalor);

  switch (type)
  {

    case Moose::ElementElement:
      jac += prefactor * _phi[_i][_qp] / _solubility_in_liquid[_qp] * _test[_i][_qp];
      break;

    case Moose::NeighborNeighbor:
      jac += prefactor * _phi_neighbor[_j][_qp] / _solubility_in_solid[_qp] * _test_neighbor[_i][_qp];
      break;

    case Moose::NeighborElement:
      jac -= prefactor * _phi[_i][_qp] / _solubility_in_liquid[_qp] * _test_neighbor[_i][_qp];
      break;

    case Moose::ElementNeighbor:
      jac -= prefactor * _phi_neighbor[_j][_qp] / _solubility_in_solid[_qp] * _test[_i][_qp];
      break;
  }

  return jac;
}
