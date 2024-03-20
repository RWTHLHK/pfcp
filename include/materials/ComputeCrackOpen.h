#pragma once
#include "Material.h"
#include "RankTwoTensor.h"
#include "DelimitedFileReader.h"
class ComputeCrackOpen : public Material
{
public:
  static InputParameters validParams();
  ComputeCrackOpen(const InputParameters & parameters);

protected:
  virtual void initQpStatefulProperties() override;
  virtual void computeQpProperties() override;
  void getSlipSystems();
  void transformHexagonalMillerBravaisSlipSystems(const MooseUtils::DelimitedFileReader & reader);
  // member property holds base name
  const std::string _base_name;
  // number of slip systems
  const unsigned int _ns;
  // member property holds micro crack formation
  MaterialProperty<Real> & _crack_open;
  const MaterialProperty<Real> & _crack_open_old;
  // member property holds initial micro crack formation rate
  const Real & _dot_crack_open_zero;
  // micro formation degradation coefficient
  const Real _beta;
  // micro crack free energy coeffiecient
  const Real _co;
  // critical resolved normal stress to trigger micro crack formation
  // resolved normal stress
  MaterialProperty<std::vector<Real>> &_sigma;
  const Real _sigma_d;
  // degree of tau/taud
  const Real _po;
  // tensile part of neo Hookean strain energy
  const MaterialProperty<Real> & _cp_phi_pos;
  // resolved shear stress
  const MaterialProperty<RankTwoTensor> & _pk2;
  // damage value
  const VariableValue & _do;
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
};
