#include "PFCPApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "ModulesApp.h"
#include "MooseSyntax.h"

InputParameters
PFCPApp::validParams()
{
  InputParameters params = MooseApp::validParams();
  params.set<bool>("use_legacy_material_output") = false;
  return params;
}

PFCPApp::PFCPApp(InputParameters parameters) : MooseApp(parameters)
{
  PFCPApp::registerAll(_factory, _action_factory, _syntax);
}

PFCPApp::~PFCPApp() {}

void 
PFCPApp::registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  ModulesApp::registerAllObjects<PFCPApp>(f, af, s);
  Registry::registerObjectsTo(f, {"PFCPApp"});
  Registry::registerActionsTo(af, {"PFCPApp"});

  /* register custom execute flags, action syntax, etc. here */
}

void
PFCPApp::registerApps()
{
  registerApp(PFCPApp);
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
extern "C" void
PFCPApp__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  PFCPApp::registerAll(f, af, s);
}
extern "C" void
PFCPApp__registerApps()
{
  PFCPApp::registerApps();
}
