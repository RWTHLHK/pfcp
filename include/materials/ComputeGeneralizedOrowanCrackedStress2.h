#pragma once

#include "Material.h"
#include "RankTwoTensor.h"
#include "RankFourTensor.h"
#include "RotationTensor.h"
#include "DerivativeMaterialInterface.h"
#include "DelimitedFileReader.h"

/**
 * Computes energy and modifies the stress for phase field fracture. Can be used with any
 * constitutive model or elastic symmetry.
 */
class ComputeGeneralizedOrowanCrackedStress2 : public Material
{
public:
  static InputParameters validParams();

  ComputeGeneralizedOrowanCrackedStress2(const InputParameters & parameters);

protected:
  virtual void computeQpProperties();
  virtual void initQpStatefulProperties();
  void computeFractureToughness();
  void computeFractureEnergy();
  void computeCrackFormation();
  void computeCrackOpen();
  void getSlipSystems();
  void transformHexagonalMillerBravaisSlipSystems(const MooseUtils::DelimitedFileReader & reader);
  /// Base name of the stress after being modified to include cracks
  const std::string _uncracked_base_name;
  // number of slip systems
  const unsigned int _ns;
  // initial fracture toughness
  const Real _g0;
  // damaged fracture toughness
  const Real _gf;
  // varied fracture toughness
  MaterialProperty<Real> &_gc;
  // length parameter
  const MaterialProperty<Real> &_l;
  // fracture energy
  MaterialProperty<Real> &_phi_f;
  //member property to hold total micro crack formation
  MaterialProperty<Real> & _micro_crack_formation;
  const MaterialProperty<Real> & _micro_crack_formation_old;
  // member property holds initial micro crack formation rate
  const Real & _dot_micro_crack_formation_zero;
  // micro formation degradation coefficient
  const Real _alpha;
  // micro crack free energy coeffiecient
  const Real _cm;
  // degree of tau/tau_d
  const Real _pm;
  // critical resolved shear stress to trigger micro crack formation
  const Real _tau_d;
  // resolved shear stress
  const MaterialProperty<std::vector<Real>> & _tau;
  // damage value
  const VariableValue & _d;
  // damage gradient 
  const MaterialProperty<Real> & _grad_d;
  /// micro crack formation degradation
  MaterialProperty<Real> &_crack_formation_degradation;
  
  // member property holds micro crack formation
  MaterialProperty<Real> & _crack_open;
  const MaterialProperty<Real> & _crack_open_old;
  // member property holds initial micro crack formation rate
  const Real & _dot_crack_open_zero;
  // micro formation degradation coefficient
  const Real _beta;
  // micro crack free energy coeffiecient
  const Real _co;
  // degree of sigma / sigma_d
  const Real _po;
  // critical resolved normal stress to trigger micro crack formation
  // resolved normal stress
  MaterialProperty<std::vector<Real>> &_sigma;
  const Real _sigma_d;
  // cleavage plane normal tensor
  MaterialProperty<std::vector<RankTwoTensor>> & _normal_tensor;
  // crysrot
  const MaterialProperty<RankTwoTensor> & _crysrot;
  // slip plane normal
  std::vector<RealVectorValue> _slip_plane_normal;
  // file name
  std::string _slip_sys_file_name;
  const enum class CrystalLatticeType { BCC, FCC, HCP } _crystal_lattice_type;
  const std::vector<Real> _unit_cell_dimension;
  Real _zero_tol;
  /// crack opening degratdaion
  MaterialProperty<Real> &_crack_open_degradation;
  
  /// uncracked pk2 stress
  const MaterialProperty<RankTwoTensor> &_uncracked_pk2;
  MaterialProperty<RankTwoTensor> &_uncracked_pk2_pos;
  const MaterialProperty<RankTwoTensor> &_uncracked_pk2_pos_old;
  /// cracked pk2 stress
  MaterialProperty<RankTwoTensor> &_cracked_pk2;
  // tensile part strain energy
  MaterialProperty<Real> &_cp_phi_pos;
  const MaterialProperty<Real> &_cp_phi_pos_old;
  // Green-Lagrangian strain 
  const MaterialProperty<RankTwoTensor> &_total_lagrangian_strain;
  const MaterialProperty<RankTwoTensor> &_total_lagrangian_strain_old;
};
