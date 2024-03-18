#include "VectorComponentAux.h"

// MOOSE includes
#include "MooseVariableInterface.h"
registerMooseObject("PFCPApp", VectorComponentAux);
InputParameters VectorComponentAux::validParams()
{
  InputParameters params = AuxKernel::validParams();
  params.addClassDescription("An AuxKernel to access a specific component of a std::vector<Real> material property.");
  params.addRequiredParam<MaterialPropertyName>("property", "The name of the vector material property.");
  params.addRequiredParam<unsigned int>("component", "The index of the component to access.");
  return params;
}

VectorComponentAux::VectorComponentAux(const InputParameters & parameters) :
  AuxKernel(parameters),
  _vectorProperty(getMaterialProperty<std::vector<Real>>("property")),
  _componentIndex(getParam<unsigned int>("component"))
{
}

Real VectorComponentAux::computeValue()
{
  if (_componentIndex >= _vectorProperty.size())
  {
    mooseError("VectorComponentAux: component_index is out of range for the vector property.");
  }
  return _vectorProperty[_qp][_componentIndex];
}
