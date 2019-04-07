#ifndef ADPRECIPITATION_H
#define ADPRECIPITATION_H

#include "ADKernel.h"

template <ComputeStage compute_stage>
class ADPrecipitation;

declareADValidParams(ADPrecipitation);

template <ComputeStage compute_stage>
class ADPrecipitation : public ADKernel<compute_stage>
{
public:
  ADPrecipitation(const InputParameters & parameters);

protected:
  virtual ADResidual computeQpResidual() override;

  const ADReal PI;
  const ADMaterialProperty(Real) & _scale_factor;
  const bool _has_coupled_precipitate;
  const bool _has_coupled_dissolve;
  const int _sign;
  const ADVariableValue & _coupled_variable;
  const ADMaterialProperty(Real) & _diffusivity;
  const ADMaterialProperty(Real) & _solubility;
  const Real _unit_scalor;

  usingKernelMembers;


};
#endif //ADPRECIPITATION_H
