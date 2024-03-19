#include "ComputeCrackOpen.h"
#include "libmesh/utility.h"
#include "libmesh/int_range.h"
#include "MooseException.h"
#include "Conversion.h"

registerMooseObject("PFCPApp", ComputeCrackOpen);

InputParameters
ComputeCrackOpen ::validParams()
{
  InputParameters params = Material::validParams();
  params.addClassDescription("update the internal variable associated with "
                             "crack opening");
  params.addParam<std::string>(
      "base_name", "Optional parameter that distincs between cracked and uncracked stress");
  params.addRequiredParam<Real>("dot_o0", "initial crack rate");
  params.addRequiredParam<Real>("beta", "crack opening coefficient");
  params.addRequiredParam<Real>("co", "crack opening free energy coeffiecient");
  params.addRequiredParam<Real>("sigma_d",
                                "critical resolved normal stress to trigger crack opening");
  params.addRequiredParam<Real>("po", "degree of sigma/sigma_d");
  params.addRequiredParam<unsigned int>("number_slip_systems",
                                        "number of slip systems of the crystal");
  params.addRequiredCoupledVar("do", "damage variable");
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

ComputeCrackOpen ::ComputeCrackOpen(const InputParameters & parameters)
  : Material(parameters),
    _base_name(isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : ""),
    _ns(getParam<unsigned int>("number_slip_systems")),
    _crack_open(declareProperty<std::vector<Real>>("ro")),
    _crack_open_old(getMaterialPropertyOld<std::vector<Real>>("ro")),
    _dot_crack_open_zero(getParam<Real>("dot_o0")),
    _beta(getParam<Real>("beta")),
    _co(getParam<Real>("co")),
    _sigma_d(getParam<Real>("sigma_d")),
    _po(getParam<Real>("po")),
    _phi_pos(getMaterialProperty<Real>(_base_name + "neo_Hookean_pos")),
    _pk2(getMaterialProperty<RankTwoTensor>("second_piola_kirchhoff_stress")),
    _do(coupledValue("do")),
    _normal_tensor(declareProperty<std::vector<RankTwoTensor>>(_base_name + "normal_tensor")),
    _crysrot(getMaterialProperty<RankTwoTensor>(_base_name + "crysrot")),
    _slip_plane_normal(_ns),
    _slip_sys_file_name(getParam<FileName>("slip_sys_file_name")),
    _crystal_lattice_type(
        getParam<MooseEnum>("crystal_lattice_type").getEnum<CrystalLatticeType>()),
    _unit_cell_dimension(getParam<std::vector<Real>>("unit_cell_dimension")),
    _zero_tol(getParam<Real>("zero_tol"))
{
}

void
ComputeCrackOpen::initQpStatefulProperties()
{
  _crack_open[_qp].resize(_ns, 0.0);
  _normal_tensor[_qp].resize(_ns);
  for (const auto i : make_range(_ns))
  {
    _normal_tensor[_qp][i].zero();
  }
}

void
ComputeCrackOpen::transformHexagonalMillerBravaisSlipSystems(
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
ComputeCrackOpen::getSlipSystems()
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
ComputeCrackOpen::computeQpProperties()
{
  const Real d = _do[_qp];
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
    sigma[i] = _pk2[_qp].doubleContraction(_normal_tensor[_qp][i]);
  }
  // compute degration function of damage
  Real ho = std::pow(1.0 - d, 2);
  //  compute generalized energetic force
  std::vector<Real> fo(_ns, 0.0);
  for (const auto i : make_range(_ns))
  {
    fo[i] = 0.5 * _beta * std::exp(-_beta * _crack_open[_qp][i]) * _phi_pos[_qp] - ho * _co;
    if (fo[i] < 0)
    {
      fo[i] = 0.0;
    }
  }

  // compute micro cracks formation rate
  std::vector<Real> dot_o(_ns, 0);
  for (const auto i : make_range(_ns))
  {
    dot_o[i] = _dot_crack_open_zero * fo[i] * std::pow(abs(sigma[i] / _sigma_d), _po);
  }

  // update micro cracks formation
  for (const auto i : make_range(_ns))
  {
    _crack_open[_qp][i] += dot_o[i] * _dt;
  }
}
