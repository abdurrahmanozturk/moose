//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "VariableWarehouse.h"
#include "MooseVariableFE.h"
#include "MooseVariableScalar.h"
#include "MooseTypes.h"

VariableWarehouse::VariableWarehouse() {}

VariableWarehouse::~VariableWarehouse()
{
  for (auto & var : _all_objects)
    delete var;
}

void
VariableWarehouse::add(const std::string & var_name, MooseVariableBase * var)
{
  _names.push_back(var_name);
  _var_name[var_name] = var;
  _all_objects.push_back(var);

  if (auto * tmp_var = dynamic_cast<MooseVariableFEBase *>(var))
  {
    _vars.push_back(tmp_var);
    if (auto * tmp_var = dynamic_cast<MooseVariable *>(var))
    {
      _regular_vars_by_number[tmp_var->number()] = tmp_var;
      _regular_vars_by_name[var_name] = tmp_var;
    }
    else if (auto * tmp_var = dynamic_cast<VectorMooseVariable *>(var))
    {
      _vector_vars_by_number[tmp_var->number()] = tmp_var;
      _vector_vars_by_name[var_name] = tmp_var;
    }
    else
      mooseError("Unknown variable class passed into VariableWarehouse. Attempt to hack us?");
  }
  else if (auto * tmp_var = dynamic_cast<MooseVariableScalar *>(var))
    _scalar_vars.push_back(tmp_var);
  else
    mooseError("Unknown variable class passed into VariableWarehouse. Attempt to hack us?");
}

void
VariableWarehouse::addBoundaryVar(BoundaryID bnd, MooseVariableFEBase * var)
{
  _boundary_vars[bnd].insert(var);
}

void
VariableWarehouse::addBoundaryVar(const std::set<BoundaryID> & boundary_ids,
                                  MooseVariableFEBase * var)
{
  for (const auto & bid : boundary_ids)
    addBoundaryVar(bid, var);
}

void
VariableWarehouse::addBoundaryVars(
    const std::set<BoundaryID> & boundary_ids,
    const std::map<std::string, std::vector<MooseVariableFEBase *>> & vars)
{
  for (const auto & bid : boundary_ids)
    for (const auto & it : vars)
      for (const auto & var : it.second)
        addBoundaryVar(bid, var);
}

MooseVariableBase *
VariableWarehouse::getVariable(const std::string & var_name)
{
  return _var_name[var_name];
}

MooseVariableBase *
VariableWarehouse::getVariable(unsigned int var_number)
{
  if (var_number < _all_objects.size())
    return _all_objects[var_number];
  else
    return NULL;
}

const std::vector<VariableName> &
VariableWarehouse::names() const
{
  return _names;
}

const std::vector<MooseVariableFEBase *> &
VariableWarehouse::fieldVariables()
{
  return _vars;
}

const std::vector<MooseVariableScalar *> &
VariableWarehouse::scalars()
{
  return _scalar_vars;
}

const std::set<MooseVariableFEBase *> &
VariableWarehouse::boundaryVars(BoundaryID bnd)
{
  return _boundary_vars[bnd];
}

template <typename T>
MooseVariableFE<T> *
VariableWarehouse::getFieldVariable(const std::string & var_name)
{
  return _regular_vars_by_name.at(var_name);
}

template <typename T>
MooseVariableFE<T> *
VariableWarehouse::getFieldVariable(unsigned int var_number)
{
  return _regular_vars_by_number.at(var_number);
}

template <>
MooseVariableFE<RealVectorValue> *
VariableWarehouse::getFieldVariable<RealVectorValue>(const std::string & var_name)
{
  return _vector_vars_by_name.at(var_name);
}

template <>
MooseVariableFE<RealVectorValue> *
VariableWarehouse::getFieldVariable<RealVectorValue>(unsigned int var_number)
{
  return _vector_vars_by_number.at(var_number);
}

template MooseVariableFE<Real> *
VariableWarehouse::getFieldVariable<Real>(const std::string & var_name);
template MooseVariableFE<Real> * VariableWarehouse::getFieldVariable<Real>(unsigned int var_number);
