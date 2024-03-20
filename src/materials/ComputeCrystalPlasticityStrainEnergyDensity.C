#include "ComputeCrystalPlasticityStrainEnergyDensity.h"

registerMooseObject("PFCPApp", ComputeCrystalPlasticityStrainEnergyDensity);

InputParameters
ComputeCrystalPlasticityStrainEnergyDensity ::validParams()
{
  InputParameters params = Material::validParams();
  params.addClassDescription("compute tensile part of NeoHookean strain energy");
  params.addParam<std::string>("base_name", "base name of deformation gradient");
  return params;
}

ComputeCrystalPlasticityStrainEnergyDensity ::ComputeCrystalPlasticityStrainEnergyDensity(
    const InputParameters & parameters)
  : Material(parameters),
    _base_name(isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : ""),
    _cp_phi(declareProperty<Real>(_base_name + "cp_phi")),
    _cp_phi_old(getMaterialPropertyOld<Real>(_base_name + "cp_phi")),
    _pk2(getMaterialProperty<RankTwoTensor>("second_piola_kirchhoff_stress")),
    _pk2_old(getMaterialPropertyOld<RankTwoTensor>("second_piola_kirchhoff_stress")),
    _total_lagrangian_strain(getMaterialProperty<RankTwoTensor>("total_lagrangian_strain")),
    _total_lagrangian_strain_old(getMaterialPropertyOld<RankTwoTensor>("total_lagrangian_strain"))
{
}

void
ComputeCrystalPlasticityStrainEnergyDensity ::initQpStatefulProperties()
{
  _cp_phi[_qp] = 0.0;
}

void
ComputeCrystalPlasticityStrainEnergyDensity ::computeQpProperties()
{
  // compute strain increment
  RankTwoTensor strain_increment = _total_lagrangian_strain[_qp] - _total_lagrangian_strain_old[_qp];
  _cp_phi[_qp] =
        _cp_phi_old[_qp] +
        MetaPhysicL::raw_value(_pk2[_qp])
                .doubleContraction(MetaPhysicL::raw_value(strain_increment)) / 2.0 +
        _pk2_old[_qp].doubleContraction(MetaPhysicL::raw_value(strain_increment)) / 2.0;
}
