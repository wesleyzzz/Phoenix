//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
#ifndef PHOENIXAPP_H
#define PHOENIXAPP_H

#include "MooseApp.h"

class PhoenixApp;

template <>
InputParameters validParams<PhoenixApp>();

class PhoenixApp : public MooseApp
{
public:
  PhoenixApp(InputParameters parameters);
  virtual ~PhoenixApp();

  static void registerApps();
  static void registerAll(Factory & f, ActionFactory & af, Syntax & s);
};

#endif /* PHOENIXAPP_H */
