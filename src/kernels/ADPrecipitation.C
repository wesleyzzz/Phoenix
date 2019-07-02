// This kernel include precipitation and redissolution USING AD

#include "ADPrecipitation.h"
#include <cmath>

registerADMooseObject("PhoenixApp", ADPrecipitation);

defineADValidParams(ADPrecipitation,
                    ADKernel,
                    params.addCoupledVar("precipitate_variable", "Ln solute in liquid precipitate.");
                    params.addCoupledVar("dissolve_variable", "Ln precipitation in liquid redissolve.");
                    params.addParam<MaterialPropertyName>("scale_factor", 1.0, "scale factor of the precipitation");
                    params.addParam<MaterialPropertyName>("diffusivity", 1.0, "Name of diffusivity material properties");
                    params.addParam<MaterialPropertyName>("solubility", 0.0, "Solubility of the Ln solute in corresponding media");
                    params.addRequiredParam<Real>("unit_scalor", "unit_scalor for the precipitate, C_Ln_pure");
                    );

template <ComputeStage compute_stage>
ADPrecipitation<compute_stage>::ADPrecipitation(const InputParameters & parameters) :
  ADKernel<compute_stage>(parameters),
  PI(3.1415926),
  _scale_factor(getADMaterialProperty<Real>("scale_factor")),
  _has_coupled_precipitate(isCoupled("precipitate_variable")),
  _has_coupled_dissolve(isCoupled("dissolve_variable")),
  _sign((_has_coupled_precipitate) ? -1 : 1),
  _coupled_variable((_has_coupled_precipitate) ? adCoupledValue("precipitate_variable") : adCoupledValue("dissolve_variable")),
  _diffusivity(getADMaterialProperty<Real>("diffusivity")),
  _solubility(getADMaterialProperty<Real>("solubility")),
  _unit_scalor(getParam<Real>("unit_scalor"))
{
}

template <ComputeStage compute_stage>
ADResidual
ADPrecipitation<compute_stage>::computeQpResidual()
{
  ADReal r = 0.0;
  ADReal prefactor = 4.0 * PI * _scale_factor[_qp] * _sign * _diffusivity[_qp];

  if (_has_coupled_precipitate)
  {
    if (_u[_qp] > _solubility[_qp])
      r = prefactor * (_u[_qp] - _solubility[_qp]);
    else if (_u[_qp] < _solubility[_qp] && _coupled_variable[_qp] > 0.0)
      r = prefactor * (_u[_qp] - _solubility[_qp]) * _coupled_variable[_qp]/_unit_scalor;
  }
  else if (_has_coupled_dissolve)
  {
    if (_coupled_variable[_qp] > _solubility[_qp])
      r = prefactor * (_coupled_variable[_qp] - _solubility[_qp]);
    else if (_coupled_variable[_qp] < _solubility[_qp] && _u[_qp] > 0.0)
      r = prefactor * (_coupled_variable[_qp] - _solubility[_qp]) * _u[_qp]/_unit_scalor;
  }
  else
    mooseError("Either precipitate_variable or dissolve_variable must be provided");

  return -_test[_i][_qp] * r;
}
