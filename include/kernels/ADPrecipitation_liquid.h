#ifndef ADPRECIPITATION_LIQUID_H
#define ADPRECIPITATION_LIQUID_H

#include "ADKernel.h"

template <ComputeStage compute_stage>
class ADPrecipitation_liquid;

declareADValidParams(ADPrecipitation_liquid);

template <ComputeStage compute_stage>
class ADPrecipitation_liquid : public ADKernel<compute_stage>
{
public:
  ADPrecipitation_liquid(const InputParameters & parameters);

protected:
  virtual ADResidual computeQpResidual() override;

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
