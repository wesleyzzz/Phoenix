
#ifndef PRECIPITATION_LIQUID_H
#define PRECIPITATION_LIQUID_H

#include "Kernel.h"

// Forward Declaration
class Precipitation_liquid;

template<>
InputParameters validParams<Precipitation_liquid>();

class Precipitation_liquid : public Kernel
{
public:
  Precipitation_liquid(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;

  const MaterialProperty<Real> & _scale_factor;
  const bool _has_coupled_precipitate;
  const bool _has_coupled_dissolve;
  const int _sign;
  const VariableValue & _coupled_variable;
  const MaterialProperty<Real> & _solubility;
  unsigned int _coupled_variable_offJac;

};
#endif //PRECIPITATION_LIQUID_H
