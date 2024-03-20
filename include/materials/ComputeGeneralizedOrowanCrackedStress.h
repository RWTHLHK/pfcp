#pragma once

#include "Material.h"
#include "RankTwoTensor.h"
#include "RankFourTensor.h"
#include "RotationTensor.h"
#include "DerivativeMaterialInterface.h"

/**
 * Computes energy and modifies the stress for phase field fracture. Can be used with any
 * constitutive model or elastic symmetry.
 */
class ComputeGeneralizedOrowanCrackedStress : public Material
{
public:
  static InputParameters validParams();

  ComputeGeneralizedOrowanCrackedStress(const InputParameters & parameters);

protected:
  virtual void computeQpProperties();
  virtual void initQpStatefulProperties();
  /// Base name of the stress after being modified to include cracks
  const std::string _base_name;

  /// Base name of the uncracked stress and strain
  const std::string _uncracked_base_name;
  //neo Hookean constants mu/2 and lambda/2
  const Real _nH1;
  const Real _nH2;
  /// deformation gradient
  const MaterialProperty<RankTwoTensor> &_deformation_gradient;
  /// plastic deformation gradient
  const MaterialProperty<RankTwoTensor> &_plastic_deformation_gradient;
  /// micro crack formation degradation
  const MaterialProperty<Real> &_crack_formation_degradation;
  /// crack opening degratdaion
  const MaterialProperty<Real> &_crack_open_degradation;
  /// uncracked cauchy stress
  MaterialProperty<RankTwoTensor> &_uncracked_cauchy_stress;
  /// cracked cauchy stress
  MaterialProperty<RankTwoTensor> &_cracked_cauchy_stress;
};
