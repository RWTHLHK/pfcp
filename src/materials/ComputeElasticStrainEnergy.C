#include "ComputeElasticStrainEnergy.h"

registerMooseObject("PFCPApp", ComputeElasticStrainEnergy);

InputParameters
ComputeElasticStrainEnergy ::validParams()
{
  InputParameters params = Material::validParams();
  params.addClassDescription("compute tensile part of NeoHookean strain energy");
  params.addRequiredParam<unsigned int>("dimension", "dimension of the problem");
  params.addRequiredParam<Real>("nH1", "neo Hookean constant mu/2");
  params.addRequiredParam<Real>("nH2", "neo Hookean constant lambda/2");
  params.addParam<std::string>("base_name", "base name of deformation gradient");
  return params;
}

ComputeElasticStrainEnergy ::ComputeElasticStrainEnergy(
    const InputParameters & parameters)
  : Material(parameters),
    _base_name(isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : ""),
    _d(getParam<unsigned int>("dimension")),
    _nH1(getParam<Real>("nH1")),
    _nH2(getParam<Real>("nH2")),
    _phi_pos(declareProperty<Real>(_base_name + "neo_Hookean_pos")),
    _phi_pos_old(getMaterialPropertyOld<Real>(_base_name + "neo_Hookean_pos")),
    _deformation_gradient(getMaterialProperty<RankTwoTensor>(_base_name + "deformation_gradient")),
    _plastic_deformation_gradient(
        getMaterialProperty<RankTwoTensor>("plastic_deformation_gradient"))
{
}

void
ComputeElasticStrainEnergy ::initQpStatefulProperties()
{
  _phi_pos[_qp] = 0.0;
}

void
ComputeElasticStrainEnergy ::computeQpProperties()
{
  // compute elastic deformtation gradient
  RankTwoTensor Fe = _deformation_gradient[_qp] * _plastic_deformation_gradient[_qp].inverse();
  // compute elastic part of right cauchy strain tensor
  Real phi_pos = 0.0;
  RankTwoTensor c = Fe.transpose() * Fe;
  // compute principal stretches of right cauchy strain tensor
  std::vector<Real> lamda;
  c.symmetricEigenvalues(lamda);

  // determinant of deformation gradient
  Real J = 1.0;
  // compute tensile part of neo Hookean strain energy
  for (unsigned int i = 0; i < _d; ++i)
  {
    J *= std::sqrt(lamda[i]);
    if (lamda[i] >= 1)
    {
      phi_pos += lamda[i];
    }
  }
  if (J > 1)
  {
    phi_pos -= 2.0 * std::log(J);
    phi_pos *= _nH1;
    phi_pos += _nH2 * (J - 1.0) * (J - 1.0);
  }
  else
  {
    phi_pos *= _nH1; // ?
  }
  _phi_pos[_qp] = phi_pos;
}
