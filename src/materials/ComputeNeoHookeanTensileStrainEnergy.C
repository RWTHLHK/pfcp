#include "ComputeNeoHookeanTensileStrainEnergy.h"

registerMooseObject("PFCPApp", ComputeNeoHookeanTensileStrainEnergy);

InputParameters
ComputeNeoHookeanTensileStrainEnergy :: validParams(){
    InputParameters params = Material::validParams();
    params.addClassDescription("compute tensile part of NeoHookean strain energy");
    params.addRequiredParam<unsigned int>("dimension","dimension of the problem");
    params.addRequiredParam<Real>("nH1","neo Hookean constant mu/2");
    params.addRequiredParam<Real>("nH2","neo Hookean constant lambda/2");
    params.addParam<std::string>("base_name","base name of deformation gradient");
    return params;
}

ComputeNeoHookeanTensileStrainEnergy :: ComputeNeoHookeanTensileStrainEnergy(const InputParameters & parameters):
Material(parameters),
_base_name(getParam<std::string>("base_name")),
_d(getParam<unsigned int>("dimension")),
_nH1(getParam<Real>("nH1")),
_nH2(getParam<Real>("nH2")),
_phi_pos(declareProperty<Real>(_base_name + "_neo_Hookean_pos")),
_phi_pos_old(getMaterialPropertyOld<Real>(_base_name + "_neo_Hookean_pos")),
_deformation_gradient(getMaterialProperty<RankTwoTensor>(_base_name + "_deformation_gradient"))
{}

void
ComputeNeoHookeanTensileStrainEnergy :: initQpStatefulProperties()
{
    _phi_pos[_qp] = 0.0;
}

void
ComputeNeoHookeanTensileStrainEnergy :: computeQpProperties()
{
    // compute right cauchy strain tensor 
    RankTwoTensor c = _deformation_gradient[_qp].transpose() * _deformation_gradient[_qp];
    // compute principal stretches of right cauchy strain tensor
    std::vector<Real> lamda;
    c.symmetricEigenvalues(lamda);
    
    // determinant of deformation gradient
    Real J = 1.0;
    //compute tensile part of neo Hookean strain energy
    for(unsigned int i=0;i<_d;++i){
        J *= std::sqrt(lamda[i]);
        if (lamda[i] >= 1){
            _phi_pos[_qp] += lamda[i];
        }
    }
    if(J>1){
        _phi_pos[_qp] -= 2.0 * std::log(J);
        _phi_pos[_qp] *= _nH1;
        _phi_pos[_qp] += _nH2*(J - 1.0) * (J - 1.0);
    }
    else{
        _phi_pos[_qp] *= _nH1; // ?
    }
}
