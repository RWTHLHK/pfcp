#pragma once
#include "Material.h"
#include "RankTwoTensor.h"
#include "DerivativeMaterialInterface.h"

class ComputeCrystalPlasticityStrainEnergyDensity : public Material {
    public:
        static InputParameters validParams();
        ComputeCrystalPlasticityStrainEnergyDensity(const InputParameters & parameters);
    protected:
        virtual void computeQpProperties() override;
        virtual void initQpStatefulProperties() override;
        //member property holds base name
        const std::string _base_name;
        //strain energy
        MaterialProperty<Real> &_cp_phi;
        const MaterialProperty<Real> &_cp_phi_old;
        // 2nd piola kirchhoff stress
        const MaterialProperty<RankTwoTensor> &_pk2;
        const MaterialProperty<RankTwoTensor> &_pk2_old;
        // strain tensor to compute strain energy
        const MaterialProperty<RankTwoTensor> &_total_lagrangian_strain;
        // old strain tensor to compute strain energy
        const MaterialProperty<RankTwoTensor> &_total_lagrangian_strain_old;
};
