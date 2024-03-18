 #pragma once 
// MOOSE includes
#include "AuxKernel.h"
// Forward declarations
class VectorComponentAux : public AuxKernel
{
public:
  static InputParameters validParams();
  VectorComponentAux(const InputParameters & parameters);

protected:
  virtual Real computeValue() override;

private:
  const MaterialProperty<std::vector<Real>> & _vectorProperty;
  unsigned int _componentIndex;
};
