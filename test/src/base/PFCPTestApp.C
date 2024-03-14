//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
#include "PFCPTestApp.h"
#include "PFCPApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "MooseSyntax.h"

InputParameters
PFCPTestApp::validParams()
{
  InputParameters params = PFCPApp::validParams();
  params.set<bool>("use_legacy_material_output") = false;
  return params;
}

PFCPTestApp::PFCPTestApp(InputParameters parameters) : MooseApp(parameters)
{
  PFCPTestApp::registerAll(
      _factory, _action_factory, _syntax, getParam<bool>("allow_test_objects"));
}

PFCPTestApp::~PFCPTestApp() {}

void
PFCPTestApp::registerAll(Factory & f, ActionFactory & af, Syntax & s, bool use_test_objs)
{
  PFCPApp::registerAll(f, af, s);
  if (use_test_objs)
  {
    Registry::registerObjectsTo(f, {"PFCPTestApp"});
    Registry::registerActionsTo(af, {"PFCPTestApp"});
  }
}

void
PFCPTestApp::registerApps()
{
  registerApp(PFCPApp);
  registerApp(PFCPTestApp);
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
// External entry point for dynamic application loading
extern "C" void
PFCPTestApp__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  PFCPTestApp::registerAll(f, af, s);
}
extern "C" void
PFCPTestApp__registerApps()
{
  PFCPTestApp::registerApps();
}
