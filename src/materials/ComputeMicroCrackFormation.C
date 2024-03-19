#include "ComputeMicroCrackFormation.h"
#include "libmesh/utility.h"

registerMooseObject("PFCPApp", ComputeMicroCrackFormation);

InputParameters
ComputeMicroCrackFormation ::validParams()
{
  InputParameters params = Material::validParams();
  params.addClassDescription("update the internal variable associated with "
                             "micro cracks formation");
  params.addParam<std::string>(
      "base_name", "Optional parameter that distincs between cracked and uncracked stress");
  params.addRequiredParam<Real>("dot_m0", "initial formation rate");
  params.addRequiredParam<Real>("alpha", "micro formation degradation coefficient");
  params.addRequiredParam<Real>("cm", "micro crack free energy coeffiecient");
  params.addRequiredParam<Real>("tau_d",
                                "critical resolved shear stress to trigger micro crack formation");
  params.addRequiredParam<Real>("pm", "degree of tau/taud");
  params.addRequiredParam<unsigned int>("number_slip_systems",
                                        "number of slip systems of the crystal");
  params.addRequiredCoupledVar("dm", "damage variable");
  params.addRequiredParam<Real>("d_duc", "critical damage value");
  return params;
}

ComputeMicroCrackFormation ::ComputeMicroCrackFormation(const InputParameters & parameters)
  : Material(parameters),
    _base_name(isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : ""),
    _ns(getParam<unsigned int>("number_slip_systems")),
    _micro_crack_formation(declareProperty<Real>("mf")),
    _micro_crack_formation_old(getMaterialPropertyOld<Real>("mf")),
    _dot_micro_crack_formation_zero(getParam<Real>("dot_m0")),
    _alpha(getParam<Real>("alpha")),
    _cm(getParam<Real>("cm")),
    _tau_d(getParam<Real>("tau_d")),
    _pm(getParam<Real>("pm")),
    _phi_pos(getMaterialProperty<Real>(_base_name + "neo_Hookean_pos")),
    _tau(getMaterialProperty<std::vector<Real>>(_base_name + "applied_shear_stress")),
    _dm(coupledValue("dm")),
    _d_duc(getParam<Real>("d_duc"))

{
}

void
ComputeMicroCrackFormation::initQpStatefulProperties()
{
  _micro_crack_formation[_qp] = 0.0;
}

void
ComputeMicroCrackFormation::computeQpProperties()
{
  const Real d = _dm[_qp];
  if (d < _d_duc)
  {
    _micro_crack_formation[_qp] = _micro_crack_formation_old[_qp];
    // compute degration function of damage
    Real hm = 1.0 / std::pow(_d_duc, 2) * std::pow(_d_duc - d, 2);
    //  compute generalized energetic force
    Real fm = 0.5 * _alpha * std::exp(-_alpha * _micro_crack_formation[_qp]) * _phi_pos[_qp] - hm * _cm;
    if (fm < 0){
      fm = 0.0;
    }

    // compute micro cracks formation rate
    std::vector<Real> dot_m(_ns, 0);
    for (const auto i : make_range(_ns))
    {
      dot_m[i] =
          _dot_micro_crack_formation_zero * fm * std::pow(abs(_tau[_qp][i] / _tau_d), _pm);
    }

    // update micro cracks formation
    for (const auto i : make_range(_ns))
    {
      _micro_crack_formation[_qp] += dot_m[i] * _dt;
    }
  }
}
