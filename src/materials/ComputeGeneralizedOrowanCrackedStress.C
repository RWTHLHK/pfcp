#include "ComputeGeneralizedOrowanCrackedStress.h"

registerMooseObject("PFCPApp", ComputeGeneralizedOrowanCrackedStress);

InputParameters
ComputeGeneralizedOrowanCrackedStress::validParams()
{
  InputParameters params = Material::validParams();
  params.addClassDescription("compute cracked cauchy stress based on neo Hookean Law");
  params.addRequiredParam<Real>("nH1", "the 1st neo Hookean constant");
  params.addRequiredParam<Real>("nH2", "the 2nd neo Hookean constant");
  params.addParam<std::string>("base_name", "The base name used to save the cracked stress");
  params.addRequiredParam<std::string>("uncracked_base_name",
                                       "The base name used to calculate the original stress");
  return params;
}

ComputeGeneralizedOrowanCrackedStress::ComputeGeneralizedOrowanCrackedStress(
    const InputParameters & parameters)
  : Material(parameters),
    _base_name(isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : ""),
    _uncracked_base_name(isParamValid("uncracked_base_name")
                             ? getParam<std::string>("uncracked_base_name") + "_"
                             : ""),
    _nH1(getParam<Real>("nH1")),
    _nH2(getParam<Real>("nH2")),
    _deformation_gradient(
        getMaterialProperty<RankTwoTensor>(_uncracked_base_name + "deformation_gradient")),
    _plastic_deformation_gradient(
        getMaterialProperty<RankTwoTensor>("plastic_deformation_gradient")),
    _crack_formation_degradation(getMaterialProperty<Real>("crack_formation_degradation")),
    _crack_open_degradation(getMaterialProperty<Real>("crack_open_degradation")),
    _uncracked_cauchy_stress(declareProperty<RankTwoTensor>("uncracked_cauchy_stress")),
    _cracked_cauchy_stress(declareProperty<RankTwoTensor>("cracked_cauchy_stress"))

{
}

void
ComputeGeneralizedOrowanCrackedStress::initQpStatefulProperties()
{
  _uncracked_cauchy_stress[_qp].zero();
  _cracked_cauchy_stress[_qp].zero();
}

void
ComputeGeneralizedOrowanCrackedStress::computeQpProperties()
{
  // compute elastic deformtation gradient
  RankTwoTensor Fe = _deformation_gradient[_qp] * _plastic_deformation_gradient[_qp].inverse();
  Real J = Fe.det();
  _uncracked_cauchy_stress[_qp] =
      _nH1 * 1.0 / J * (Fe * Fe.transpose() - RankTwoTensor::Identity()) +
      _nH2 * (J - 1.0) * RankTwoTensor::Identity();
  if (J < 1)
  {
    _cracked_cauchy_stress[_qp] =
        (_crack_formation_degradation[_qp] + _crack_open_degradation[_qp]) * _nH1 * 1.0 / J *
        (Fe * Fe.transpose() - RankTwoTensor::Identity());
  }
  else
  {
    _cracked_cauchy_stress[_qp] =
        (_crack_formation_degradation[_qp] + _crack_open_degradation[_qp]) *
        _uncracked_cauchy_stress[_qp];
  }
}
