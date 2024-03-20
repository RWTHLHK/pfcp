#pragma once
#include "Material.h"
#include "RankTwoTensor.h"
class ComputeMicroCrackFormation : public Material
{
public:
  static InputParameters validParams();
  ComputeMicroCrackFormation(const InputParameters & parameters);

protected:
  virtual void initQpStatefulProperties() override;
  virtual void computeQpProperties() override;
  // member property holds base name
  const std::string _base_name;
  // number of slip systems
  const unsigned int _ns;
  //member property to hold total micro crack formation
  MaterialProperty<Real> & _micro_crack_formation;
  const MaterialProperty<Real> & _micro_crack_formation_old;
  // member property holds initial micro crack formation rate
  const Real & _dot_micro_crack_formation_zero;
  // micro formation degradation coefficient
  const Real & _alpha;
  // micro crack free energy coeffiecient
  const Real & _cm;
  // critical resolved shear stress to trigger micro crack formation
  const Real & _tau_d;
  // degree of tau/taud
  const Real & _pm;
  // tensile part of strain energy
  const MaterialProperty<Real> & _cp_phi_pos;
  // resolved shear stress
  const MaterialProperty<std::vector<Real>> & _tau;
  // damage value
  const VariableValue & _dm;
  // crtitical damage value
  const Real _d_duc;
};
