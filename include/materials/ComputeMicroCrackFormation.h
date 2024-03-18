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
    //member property holds base name
    const std::string _base_name;
    // number of slip systems
    const unsigned int _ns;
    // member property holds micro crack formation
    MaterialProperty<std::vector<Real>> &_micro_crack_formation;
    const MaterialProperty<std::vector<Real>> &_micro_crack_formation_old;
    // member property holds initial micro crack formation rate
    const Real &_dot_micro_crack_formation_zero;
    // micro formation degradation coefficient
    const Real &_alpha;
    // micro crack free energy coeffiecient
    const Real &_cm;
    // critical resolved shear stress to trigger micro crack formation
    const Real &_tau_d;
    // degree of tau/taud
    const Real &_pm;
    // tensile part of neo Hookean strain energy
    const  MaterialProperty<Real> &_phi_pos; 
    // resolved shear stress
    const MaterialProperty<std::vector<Real>> &_tau;
};