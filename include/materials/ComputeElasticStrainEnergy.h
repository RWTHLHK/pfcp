#pragma once
#include "Material.h"
#include "RankTwoTensor.h"

class ComputeElasticStrainEnergy : public Material {
    public:
        static InputParameters validParams();
        ComputeElasticStrainEnergy(const InputParameters & parameters);
    protected:
        virtual void computeQpProperties() override;
        virtual void initQpStatefulProperties() override;
        //member property holds base name
        const std::string _base_name;
        //dimension of the problem
        const unsigned int _d;
        //neo Hookean constants mu/2 and lambda/2
        const Real _nH1;
        const Real _nH2;
        //tensile part of neo Hookean strain energy
        MaterialProperty<Real> &_phi_pos;
        const MaterialProperty<Real> &_phi_pos_old;
        // deformation gradient to compute neo Hookean strain energy
        const MaterialProperty<RankTwoTensor> &_deformation_gradient;
        // deformation gradient to compute neo Hookean strain energy
        const MaterialProperty<RankTwoTensor> &_plastic_deformation_gradient;
};
