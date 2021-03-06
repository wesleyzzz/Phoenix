// Describe s_d --> l_p

#include "InterfaceForce_lp.h"
#include <cmath>

registerMooseObject("PhoenixApp", InterfaceForce_lp);

template<>
InputParameters validParams<InterfaceForce_lp>()
{
  InputParameters params = validParams<InterfaceKernel>();
  params.addParam<Real>("unit_scalor", 1e6, "atomic density for the pure Ln");
  params.addParam<Real>("stop_threshold", 9.9e5, "Ln concentration stoping threshold for such flux");
  params.addParam<MaterialPropertyName>("driving_rate", 1.0, "driving force rate between the interface");
  params.addCoupledVar("coupled_s_p", "Ln precipitation in solid");

  params.addParam<MaterialPropertyName>("solubility_in_solid", 0.0, "Name of solubility material in the solid");
  params.addParam<MaterialPropertyName>("diffusivity_in_solid_sd", 1.0, "Name of diffusivity material in the solid fuel");
  params.addParam<MaterialPropertyName>("diffusivity_in_solid_lp", 1.0, "Name of diffusivity material in the Ln precipitate in liquid");
  return params;
}



InterfaceForce_lp::InterfaceForce_lp(const InputParameters & parameters) :
    InterfaceKernel(parameters),
    PostprocessorInterface(this),
    VectorPostprocessorInterface(this),
    _unit_scalor(getParam<Real>("unit_scalor")),
    _stop_threshold(getParam<Real>("stop_threshold")),
    _driving_rate(getMaterialProperty<Real>(getParam<MaterialPropertyName>("driving_rate"))),
    _coupled_s_p(coupledValue("coupled_s_p")),
    _solubility_in_solid(getMaterialProperty<Real>(getParam<MaterialPropertyName>("solubility_in_solid"))),
    _diffusivity_in_solid_sd(getMaterialProperty<Real>(getParam<MaterialPropertyName>("diffusivity_in_solid_sd"))),
    _diffusivity_in_solid_lp(getMaterialProperty<Real>(getParam<MaterialPropertyName>("diffusivity_in_solid_lp")))

{
  if (!parameters.isParamValid("boundary"))
  {
    mooseError("In order to use the InterfaceForce_lp dgkernel, you must specify a boundary where it will live.");
  }

}

const PostprocessorValue &
InterfaceForce_lp::getPostprocessorValue(const std::string & name)
{
  _depend_uo.insert(_pars.get<PostprocessorName>(name));
  return PostprocessorInterface::getPostprocessorValue(name);
}

const PostprocessorValue &
InterfaceForce_lp::getPostprocessorValueByName(const PostprocessorName & name)
{
  _depend_uo.insert(name);
  return PostprocessorInterface::getPostprocessorValueByName(name);
}

Real
InterfaceForce_lp::computeQpResidual(Moose::DGResidualType type)
{
  Real r = 0.0;
  Real Diff_eff = _diffusivity_in_solid_sd[_qp];
  Real Force_eff = _driving_rate[_qp] * (_neighbor_value[_qp] - _solubility_in_solid[_qp]);

  if (_u[_qp] < _stop_threshold)
  {
      r = Diff_eff * Force_eff * _u[_qp] / _unit_scalor * (1.0 - _coupled_s_p[_qp] / _unit_scalor);    // The precipitation in the pore serves as a Ln pump
  }

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
InterfaceForce_lp::computeQpJacobian(Moose::DGJacobianType type)
{
  Real jac = 0.0;
  // effective diffusivity for the s_d to l_p
  Real Diff_eff = _diffusivity_in_solid_sd[_qp];
  Real prefactor = Diff_eff * _driving_rate[_qp] * (1.0 - _coupled_s_p[_qp] / _unit_scalor);

  if (_u[_qp] < _stop_threshold)
  {
    switch (type)
    {
      case Moose::ElementElement:
        jac -= prefactor * (_neighbor_value[_qp] - _solubility_in_solid[_qp]) * _phi[_i][_qp] / _unit_scalor * _test[_i][_qp];
        break;

      case Moose::NeighborNeighbor:
        jac += prefactor * _phi_neighbor[_j][_qp] * _u[_qp] / _unit_scalor * _test_neighbor[_i][_qp];
        break;

      case Moose::NeighborElement:
        jac += prefactor * (_neighbor_value[_qp] - _solubility_in_solid[_qp]) * _phi[_i][_qp] / _unit_scalor * _test_neighbor[_i][_qp];
        break;

      case Moose::ElementNeighbor:
        jac -= prefactor * _phi_neighbor[_j][_qp] * _u[_qp] / _unit_scalor * _test[_i][_qp];
        break;
    }
  }

  return jac;
}
