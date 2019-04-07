#include "Precipitation_liquid.h"
#include <cmath>

registerMooseObject("PhoenixApp", Precipitation_liquid);

template<>
InputParameters validParams<Precipitation_liquid>()
{
  InputParameters params = validParams<Kernel>();
  params.addCoupledVar("precipitate_variable", "Ln solute in liquid precipitate.");
  params.addCoupledVar("dissolve_variable", "Ln precipitation in liquid redissolve.");
  params.addParam<MaterialPropertyName>("scale_factor", 1.0, "scale factor of the precipitation");
  params.addParam<MaterialPropertyName>("solubility", 0.0, "Solubility of the Ln solute in corresponding media");

  return params;
}

Precipitation_liquid::Precipitation_liquid(const InputParameters & parameters) :
  Kernel(parameters),
  _scale_factor(getMaterialProperty<Real>(getParam<MaterialPropertyName>("scale_factor"))),
  _has_coupled_precipitate(isCoupled("precipitate_variable")),
  _has_coupled_dissolve(isCoupled("dissolve_variable")),
  _sign((_has_coupled_precipitate) ? -1 : 1),
  _coupled_variable((_has_coupled_precipitate) ? coupledValue("precipitate_variable") : coupledValue("dissolve_variable")),
  _solubility(getMaterialProperty<Real>(getParam<MaterialPropertyName>("solubility"))),
  _coupled_variable_offJac((_has_coupled_precipitate) ? coupled("precipitate_variable") : coupled("dissolve_variable"))

{}

Real
Precipitation_liquid::computeQpResidual()
{
  Real r = 0.0;
  Real prefactor = _scale_factor[_qp] * _sign;

  if (_has_coupled_precipitate)
  {
    if (_u[_qp] > _solubility[_qp])      
      r = prefactor * (_u[_qp] - _solubility[_qp]);
    else if (_u[_qp] < _solubility[_qp] && _coupled_variable[_qp] > 0.0)
      r = prefactor * (_u[_qp] - _solubility[_qp]) * _coupled_variable[_qp];
  }
  else if (_has_coupled_dissolve)
  {
    if (_coupled_variable[_qp] > _solubility[_qp])
      r = prefactor * (_coupled_variable[_qp] - _solubility[_qp]);
    else if (_coupled_variable[_qp] < _solubility[_qp] && _u[_qp] > 0.0)
      r = prefactor * (_coupled_variable[_qp] - _solubility[_qp]) * _u[_qp];
  }
  else
    mooseError("Either precipitate_variable or dissolve_variable must be provided");

  return -_test[_i][_qp] * r;
}

Real
Precipitation_liquid::computeQpJacobian()
{
  Real jac = 0.0;
  Real prefactor = _scale_factor[_qp] * _sign;

  if (_has_coupled_precipitate)
  {
    if (_u[_qp] > _solubility[_qp])
      jac = prefactor * _phi[_j][_qp];
    else if (_u[_qp] < _solubility[_qp] && _coupled_variable[_qp] > 0.0)
      jac = prefactor * _phi[_j][_qp] * _coupled_variable[_qp];
  }
  else if (_has_coupled_dissolve)
  {
    if (_coupled_variable[_qp] > _solubility[_qp])
      jac = 0.0;
    else if (_coupled_variable[_qp] < _solubility[_qp] && _u[_qp] > 0.0)
      jac = prefactor * (_coupled_variable[_qp] - _solubility[_qp]) * _phi[_j][_qp];
  }
  else
    mooseError("Either precipitate_variable or dissolve_variable must be provided");

  return -_test[_i][_qp] * jac;
}


Real
Precipitation_liquid::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (jvar == _coupled_variable_offJac)
  {
  Real offjac = 0.0;
  Real prefactor = _scale_factor[_qp] * _sign;

  if (_has_coupled_precipitate)
  {  
    if (_u[_qp] > _solubility[_qp])
      offjac = 0.0;
    else if (_u[_qp] < _solubility[_qp] && _coupled_variable[_qp] > 0.0)
      offjac = prefactor * _phi[_j][_qp] * ( _u[_qp] - _solubility[_qp]);
  }

// this is for precipitation eqution
  else if (_has_coupled_dissolve)
  {
    if (_coupled_variable[_qp] > _solubility[_qp])
      offjac = prefactor * _phi[_j][_qp];  
    else if (_coupled_variable[_qp] < _solubility[_qp] && _u[_qp] > 0.0)
      offjac = prefactor * _phi[_j][_qp]* _u[_qp];
  }

  return _test[_i][_qp] * offjac;
  
  }

  else return 0.0;
}

