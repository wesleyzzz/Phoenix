#ifndef PRECIPITATION_H
#define PRECIPITATION_H

#include "Kernel.h"

class Precipitation;

template<>
InputParameters validParams<Precipitation>();

class Precipitation : public Kernel
{
public:
  Precipitation(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;

  const Real PI;
  const MaterialProperty<Real> & _scale_factor;
  const bool _has_coupled_precipitate;
  const bool _has_coupled_dissolve;
  const int _sign;
  const VariableValue & _coupled_variable;
  const MaterialProperty<Real> & _diffusivity;
  const MaterialProperty<Real> & _solubility;
  unsigned int _coupled_variable_offJac;
};
#endif //PRECIPITATION_H
