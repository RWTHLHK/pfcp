#include "ComputeGeneralizedOrowanCrackedStress2.h"
#include "RankFourTensor.h"

registerMooseObject("PFCPApp", ComputeGeneralizedOrowanCrackedStress2);

InputParameters
ComputeGeneralizedOrowanCrackedStress2::validParams()
{
  InputParameters params = Material::validParams();
  params.addClassDescription("compute cracked pk2 stress and tensile strain energy");
  params.addParam<std::string>("uncracked_base_name",
                               "The base name used to save the cracked stress");
  params.addRequiredParam<Real>("dot_m0", "initial formation rate");
  params.addRequiredParam<Real>("alpha", "micro formation degradation coefficient");
  params.addRequiredParam<Real>("cm", "micro crack free energy coeffiecient");
  params.addRequiredParam<Real>("pm", "degree of tau/tau_d");
  params.addRequiredParam<Real>("tau_d",
                                "critical resolved shear stress to trigger micro crack formation");
  params.addRequiredParam<unsigned int>("number_slip_systems",
                                        "number of slip systems of the crystal");
  params.addRequiredCoupledVar("d", "damage variable");
  ///crack opening part
  params.addRequiredParam<Real>("dot_o0", "initial crack rate");
  params.addRequiredParam<Real>("beta", "crack opening coefficient");
  params.addRequiredParam<Real>("co", "crack opening free energy coeffiecient");
  params.addRequiredParam<Real>("po", "degree of sigma/sigma_d");
  params.addRequiredParam<Real>("sigma_d",
                                "critical resolved normal stress to trigger crack opening");
  params.addRequiredParam<Real>("g0", "initial fracture toughness");
  params.addRequiredParam<Real>("gf", "final fracture toughness");
  params.addRequiredParam<unsigned int>("number_slip_systems",
                                        "number of slip systems of the crystal");
  params.addRequiredParam<FileName>(
      "slip_sys_file_name",
      "Name of the file containing the slip systems, one slip system per row, with the slip plane "
      "normal given before the slip plane direction.");
  params.addParam<MooseEnum>(
      "crystal_lattice_type",
      MooseEnum("BCC FCC HCP", "FCC"),
      "Crystal lattice type or representative unit cell, i.e., BCC, FCC, HCP, etc.");
  params.addRangeCheckedParam<std::vector<Real>>(
      "unit_cell_dimension",
      std::vector<Real>{1.0, 1.0, 1.0},
      "unit_cell_dimension_size = 3",
      "The dimension of the unit cell along three directions, where a cubic unit cell is assumed "
      "for cubic crystals and a hexagonal unit cell (a, a, c) is assumed for HCP crystals. These "
      "dimensions will be taken into account while computing the slip systems."
      " Default size is 1.0 along all three directions.");
  params.addParam<Real>("zero_tol",
                        1e-12,
                        "Tolerance for residual check when variable value is zero for each "
                        "individual constitutive model");
  return params;
}

ComputeGeneralizedOrowanCrackedStress2::ComputeGeneralizedOrowanCrackedStress2(
    const InputParameters & parameters)
  : Material(parameters),
    _uncracked_base_name(isParamValid("uncracked_base_name")
                             ? getParam<std::string>("uncracked_base_name") + "_"
                             : ""),
    _ns(getParam<unsigned int>("number_slip_systems")),
    _g0(getParam<Real>("g0")),
    _gf(getParam<Real>("gf")),
    _gc(declareProperty<Real>("gc_prop")),
    _l(getMaterialProperty<Real>("l")),
    _phi_f(declareProperty<Real>("fracture_energy")),
    _micro_crack_formation(declareProperty<Real>("mf")),
    _micro_crack_formation_old(getMaterialPropertyOld<Real>("mf")),
    _dot_micro_crack_formation_zero(getParam<Real>("dot_m0")),
    _alpha(getParam<Real>("alpha")),
    _cm(getParam<Real>("cm")),
    _pm(getParam<Real>("pm")),
    _tau_d(getParam<Real>("tau_d")),
    _tau(getMaterialProperty<std::vector<Real>>(_uncracked_base_name + "applied_shear_stress")),
    _d(coupledValue("d")),
    _grad_d(getMaterialProperty<Real>("grad_d")),
    _crack_formation_degradation(
        declareProperty<Real>("crack_formation_degradation")), /// crack formation part

    _crack_open(declareProperty<Real>("ro")),
    _crack_open_old(getMaterialPropertyOld<Real>("ro")),
    _dot_crack_open_zero(getParam<Real>("dot_o0")),
    _beta(getParam<Real>("beta")),
    _co(getParam<Real>("co")),
    _po(getParam<Real>("po")),
    _sigma(declareProperty<std::vector<Real>>("resolved_normal_stress")),
    _sigma_d(getParam<Real>("sigma_d")),
    _normal_tensor(declareProperty<std::vector<RankTwoTensor>>("normal_tensor")),
    _crysrot(getMaterialProperty<RankTwoTensor>(_uncracked_base_name + "crysrot")),
    _slip_plane_normal(_ns),
    _slip_sys_file_name(getParam<FileName>("slip_sys_file_name")),
    _crystal_lattice_type(
        getParam<MooseEnum>("crystal_lattice_type").getEnum<CrystalLatticeType>()),
    _unit_cell_dimension(getParam<std::vector<Real>>("unit_cell_dimension")),
    _zero_tol(getParam<Real>("zero_tol")),
    _crack_open_degradation(declareProperty<Real>("crack_open_degradation")), ///crack opening part

    _uncracked_pk2(getMaterialProperty<RankTwoTensor>("second_piola_kirchhoff_stress")),
    _uncracked_pk2_pos(declareProperty<RankTwoTensor>("pk2_pos")),
    _uncracked_pk2_pos_old(getMaterialPropertyOld<RankTwoTensor>("pk2_pos")),
    _cracked_pk2(declareProperty<RankTwoTensor>("cracked_pk2")),
    _cp_phi_pos(declareProperty<Real>("cp_phi_pos")),
    _cp_phi_pos_old(getMaterialPropertyOld<Real>("cp_phi_pos")),
    _total_lagrangian_strain(getMaterialProperty<RankTwoTensor>("total_lagrangian_strain")),
    _total_lagrangian_strain_old(getMaterialPropertyOld<RankTwoTensor>("total_lagrangian_strain"))

{
}

void 
ComputeGeneralizedOrowanCrackedStress2::computeFractureToughness()
{
  _gc[_qp] = _g0 + (_gf - _g0) * 
  (1.0 - std::exp(-_alpha*_micro_crack_formation_old[_qp] 
                  -_beta*_crack_open_old[_qp]));
}
void 
ComputeGeneralizedOrowanCrackedStress2::computeFractureEnergy()
{
  _phi_f[_qp] = _gc[_qp] / (2.0 * _l[_qp]) * _d[_qp] * _d[_qp] + 
                _gc[_qp] * _l[_qp] / 2.0 * _grad_d[_qp] * _grad_d[_qp];
}

void
ComputeGeneralizedOrowanCrackedStress2::initQpStatefulProperties()
{
  /// set fracture toughness
  _gc[_qp] = _g0;
  /// crack formation part
  _cp_phi_pos[_qp] = 0.0;
  _uncracked_pk2_pos[_qp].zero();
  _micro_crack_formation[_qp] = 0.0;
  /// crack open part
  _crack_open[_qp] = 0.0;
  _normal_tensor[_qp].resize(_ns);
  _sigma[_qp].resize(_ns);
  for (const auto i : make_range(_ns))
  {
    _normal_tensor[_qp][i].zero();
    _sigma[_qp][i] = 0.0;
  }
}

void
ComputeGeneralizedOrowanCrackedStress2::computeCrackFormation()
{
  // compute degration function of damage
  const Real damage = _d[_qp];
  Real hm = std::pow(1.0 - damage, 2);
  // compute fracture energy
  //  compute generalized energetic force(phim = cm * mf * mf)
  Real fm = - (hm * _cm * 2.0 * _micro_crack_formation_old[_qp] + _alpha * (_gf - _g0) * 
                std::exp(-_alpha * _micro_crack_formation_old[_qp] - _beta * _crack_open_old[_qp]) *
                _phi_f[_qp]);
  if (fm < 0)
  {
    fm = 0.0;
  }
  // compute micro cracks formation rate
  std::vector<Real> dot_m(_ns, 0);
  for (const auto i : make_range(_ns))
  {
    dot_m[i] = _dot_micro_crack_formation_zero * fm * std::pow(abs(_tau[_qp][i]/ _tau_d), _pm);
  }
  _micro_crack_formation[_qp] = _micro_crack_formation_old[_qp];
  // update micro cracks formation
  for (const auto i : make_range(_ns))
  {
    _micro_crack_formation[_qp] += dot_m[i] * _dt;
  }
}

void
ComputeGeneralizedOrowanCrackedStress2::transformHexagonalMillerBravaisSlipSystems(
    const MooseUtils::DelimitedFileReader & reader)
{
  const unsigned int miller_bravais_indices = 4;
  RealVectorValue temporary_slip_direction, temporary_slip_plane;
  // temporary_slip_plane.resize(LIBMESH_DIM);
  // temporary_slip_direction.resize(LIBMESH_DIM);

  if (_unit_cell_dimension[0] != _unit_cell_dimension[1] ||
      _unit_cell_dimension[0] == _unit_cell_dimension[2])
    mooseError("ComputeCrackOpen Error: The specified unit cell dimensions are "
               "not consistent with expectations for "
               "HCP crystal hexagonal lattices.");
  else if (reader.getData(0).size() != miller_bravais_indices * 2)
    mooseError("ComputeCrackOpen Error: The number of entries in the first row of "
               "the slip system file is not consistent with the expectations for the 4-index "
               "Miller-Bravais assumption for HCP crystals. This file should represent both the "
               "slip plane normal and the slip direction with 4-indices each.");

  // set up the tranformation matrices
  RankTwoTensor transform_matrix;
  transform_matrix.zero();
  transform_matrix(0, 0) = 1.0 / _unit_cell_dimension[0];
  transform_matrix(1, 0) = 1.0 / (_unit_cell_dimension[0] * std::sqrt(3.0));
  transform_matrix(1, 1) = 2.0 / (_unit_cell_dimension[0] * std::sqrt(3.0));
  transform_matrix(2, 2) = 1.0 / (_unit_cell_dimension[2]);

  for (const auto i : make_range(_ns))
  {
    // read in raw data from file and store in the temporary vectors
    for (const auto j : index_range(reader.getData(i)))
    {
      // Check that the slip plane normal indices of the basal plane sum to zero for consistency
      Real basal_pl_sum = 0.0;
      for (const auto k : make_range(LIBMESH_DIM))
        basal_pl_sum += reader.getData(i)[k];

      if (basal_pl_sum > _zero_tol)
        mooseError(
            "ComputeCrackOpen Error: The specified HCP basal plane Miller-Bravais "
            "indices do not sum to zero. Check the values supplied in the associated text file.");

      // Check that the slip direction indices of the basal plane sum to zero for consistency
      Real basal_dir_sum = 0.0;
      for (const auto k : make_range(miller_bravais_indices, miller_bravais_indices + LIBMESH_DIM))
        basal_dir_sum += reader.getData(i)[k];

      if (basal_dir_sum > _zero_tol)
        mooseError("ComputeCrackOpen Error: The specified HCP slip direction "
                   "Miller-Bravais indices in the basal plane (U, V, and T) do not sum to zero "
                   "within the user specified tolerance (try loosing zero_tol if using the default "
                   "value). Check the values supplied in the associated text file.");

      if (j < miller_bravais_indices)
      {
        // Planes are directly copied over, per a_1 = x convention used here:
        // Store the first two indices for the basal plane, (h and k), and drop
        // the redundant third basal plane index (i)
        if (j < 2)
          temporary_slip_plane(j) = reader.getData(i)[j];
        // Store the c-axis index as the third entry in the orthorombic index convention
        else if (j == 3)
          temporary_slip_plane(j - 1) = reader.getData(i)[j];
      }
      else
      {
        const auto direction_j = j - miller_bravais_indices;
        // Store the first two indices for the slip direction in the basal plane,
        //(U, V), and drop the redundant third basal plane index (T)
        if (direction_j < 2)
          temporary_slip_direction(direction_j) = reader.getData(i)[j];
        // Store the c-axis index as the third entry in the orthorombic index convention
        else if (direction_j == 3)
          temporary_slip_direction(direction_j - 1) = reader.getData(i)[j];
      }
    }

    // perform transformation calculation
    _slip_plane_normal[i] = transform_matrix * temporary_slip_plane;
  }
}

void
ComputeGeneralizedOrowanCrackedStress2::getSlipSystems()
{
  // read in the slip system data from auxiliary text file
  MooseUtils::DelimitedFileReader _reader(_slip_sys_file_name);
  _reader.setFormatFlag(MooseUtils::DelimitedFileReader::FormatFlag::ROWS);
  _reader.read();

  // check the size of the input
  if (_reader.getData().size() != _ns)
    paramError(
        "ns", "The number of rows in the slip system file should match the number of slip system.");

  for (const auto i : make_range(_ns))
  {
    // initialize to zero
    _slip_plane_normal[i].zero();
  }

  if (_crystal_lattice_type == CrystalLatticeType::HCP)
    transformHexagonalMillerBravaisSlipSystems(_reader);
  else if (_crystal_lattice_type == CrystalLatticeType::BCC ||
           _crystal_lattice_type == CrystalLatticeType::FCC)
  {
    for (const auto i : make_range(_ns))
    {
      // directly grab the raw data and scale it by the unit cell dimension
      for (const auto j : index_range(_reader.getData(i)))
      {
        if (j < LIBMESH_DIM)
          _slip_plane_normal[i](j) = _reader.getData(i)[j] / _unit_cell_dimension[j];
      }
    }
  }

  for (const auto i : make_range(_ns))
  {
    // normalize
    _slip_plane_normal[i] /= _slip_plane_normal[i].norm();
  }
}

void
ComputeGeneralizedOrowanCrackedStress2::computeCrackOpen()
{
  getSlipSystems();
  // compute resolved normal stress
  std::vector<RealVectorValue> local_plane_normal;
  local_plane_normal.resize(_ns);
  for (const auto i : make_range(_ns))
  {
    local_plane_normal[i].zero();

    for (const auto j : make_range(LIBMESH_DIM))
      for (const auto k : make_range(LIBMESH_DIM))
      {
        local_plane_normal[i](j) =
            local_plane_normal[i](j) + _crysrot[_qp](j, k) * _slip_plane_normal[i](k);
      }

    // Calculate normal tensor
    for (const auto j : make_range(LIBMESH_DIM))
      for (const auto k : make_range(LIBMESH_DIM))
      {
        _normal_tensor[_qp][i](j, k) = local_plane_normal[i](j) * local_plane_normal[i](k);
      }
  }

  std::vector<Real> sigma(_ns, 0.0);
  for (const auto i : make_range(_ns))
  {
    _sigma[_qp][i] = _uncracked_pk2[_qp].doubleContraction(_normal_tensor[_qp][i]);
    sigma[i] = _sigma[_qp][i];
    if(sigma[i] < 0){
      sigma[i] = 0.0;
    }
  }
  // compute crack open
  // compute degration function of damage
  const Real damage = _d[_qp];
  Real ho = std::pow(1.0 - damage, 2);
  //  compute generalized energetic force
  Real fo = - (ho * _co * 2.0 * _crack_open_old[_qp] + _beta * (_gf - _g0) * 
                std::exp(-_alpha * _micro_crack_formation_old[_qp] - _beta * _crack_open_old[_qp]) *
                _phi_f[_qp]);
  if (fo < 0)
  {
    fo = 0.0;
  }
  // compute cracks open rate
  std::vector<Real> dot_o(_ns, 0);
  for (const auto i : make_range(_ns))
  {
    dot_o[i] = _dot_crack_open_zero * fo * std::pow(sigma[i] / _sigma_d, _po);
  }
  _crack_open[_qp] = _crack_open_old[_qp];
  // update micro cracks formation
  for (const auto i : make_range(_ns))
  {
    _crack_open[_qp] += dot_o[i] * _dt;
  }
}

void
ComputeGeneralizedOrowanCrackedStress2::computeQpProperties()
{ 
  // new gc -> new fracture energy -> new mf and ro -> new damage
  // compute fracture toughness
  computeFractureToughness();
  // compute fracture energy 
  computeFractureEnergy();
  // compute micro crack formation using old strain energy
  computeCrackFormation();
  // compute crack opening using old strain energy
  computeCrackOpen();
  // compute strain increment
  RankTwoTensor strain_increment =
      _total_lagrangian_strain[_qp] - _total_lagrangian_strain_old[_qp];
  // Create the positive and negative projection tensors
  RankFourTensor I4sym(RankFourTensor::initIdentitySymmetricFour);
  std::vector<Real> eigval;
  RankTwoTensor eigvec;
  RankFourTensor Ppos = _uncracked_pk2[_qp].positiveProjectionEigenDecomposition(eigval, eigvec);
  RankFourTensor Pneg = I4sym - Ppos;

  // Project the positive and negative stresses
  _uncracked_pk2_pos[_qp] = Ppos * _uncracked_pk2[_qp];
  RankTwoTensor uncracked_pk2_neg = Pneg * _uncracked_pk2[_qp];

  // compute cracked stress
  _cracked_pk2[_qp] =
      _uncracked_pk2_pos[_qp] * (1.0 - _d[_qp]) * ( 1.0 - _d[_qp]) +
      uncracked_pk2_neg;
  // compute tensile part strain energy
  _cp_phi_pos[_qp] =
      _cp_phi_pos_old[_qp] +
      MetaPhysicL::raw_value(_uncracked_pk2_pos[_qp])
              .doubleContraction(MetaPhysicL::raw_value(strain_increment)) /
          2.0 +
      _uncracked_pk2_pos_old[_qp].doubleContraction(MetaPhysicL::raw_value(strain_increment)) / 2.0;
}
