#include "RadiativeHeatFluxBC.h"
#include "PenetrationLocator.h"

registerMooseObject("HeatConductionApp", RadiativeHeatFluxBC);

template <>
InputParameters
validParams<RadiativeHeatFluxBC>()
{
  InputParameters params = validParams<IntegratedBC>();
  params.addClassDescription("Radiation Heat Transfer BC Model with view factors.");
  params.addRequiredParam<UserObjectName>("viewfactor_userobject",
                                          "The name of the UserObject that is "
                                          "used for view factor calculations.");
  params.addParam<Real>("stefan_boltzmann",
                        5.670367e-8,
                        "The Stefan-Boltzmann constant [kg.s^3.K^4].");
  params.addParam<std::vector<Real>>("emissivity","The emissivities of boundaries (sort by ids).");
  params.addParam<Real>("ambient_temperature",
                        300,
                        "Ambient Temperature(in K) to calculate "
                        "radiative heat loss from a surface.");
  return params;
}

RadiativeHeatFluxBC::RadiativeHeatFluxBC(const InputParameters & parameters)
  : IntegratedBC(parameters),
    _var_number(_subproblem.getVariable(_tid,
                                 parameters.get<NonlinearVariableName>("variable"),
                                 Moose::VarKindType::VAR_ANY,
                                 Moose::VarFieldType::VAR_FIELD_STANDARD).number()),
    _system(_subproblem.getSystem(getParam<NonlinearVariableName>("variable"))),
    _viewfactor(getUserObject<ViewFactor>("viewfactor_userobject")),
    _boundary_ids(boundaryIDs()),
    _master_boundary_ids(_viewfactor.getMasterBoundaries()),
    _slave_boundary_ids(_viewfactor.getSlaveBoundaries()),
    _stefan_boltzmann(getParam<Real>("stefan_boltzmann")),
    _ambient_temp(getParam<Real>("ambient_temperature")),
    _heat_loss(isParamValid("ambient_temperature") && (_boundary_ids.size() == 1))
{
  std::vector<Real> emissivity = getParam<std::vector<Real>>("emissivity");
  if (emissivity.size()!=_boundary_ids.size())
    mooseError("The number of emissivities does not match the number of boundaries.");
  unsigned int i{0};
  for (const auto bid : _boundary_ids)
    _emissivity[bid]=emissivity[i++];

  for (const auto & t : _mesh.buildSideList())    //buildSideList(el,side,bnd)
  {
    auto elem_id = std::get<0>(t);
    auto side_id = std::get<1>(t);
    auto bnd_id = std::get<2>(t);
    if (_boundary_ids.find(bnd_id)!=_boundary_ids.end())
      _elem_side_map[elem_id]=side_id;
  }
}

Real
RadiativeHeatFluxBC::computeQpResidual()
{
  Real _u_slave, q_ms, q_sm, q_net{0};
  Real temp_func_master = _u[_qp] * _u[_qp] * _u[_qp] * _u[_qp];
  BoundaryID current_boundary_id = _mesh.getBoundaryIDs(_current_elem, _current_side)[0];
  if (_boundary_ids.find(current_boundary_id)!=_boundary_ids.end())
  {
    _master_elem_id = _current_elem->id();
    if (_heat_loss)
    {
      Real temp_func_ambient = _ambient_temp * _ambient_temp * _ambient_temp * _ambient_temp;
      q_net = _emissivity[current_boundary_id] * _stefan_boltzmann *
              (temp_func_master - temp_func_ambient); // master-ambient
    }
    else
    {
      for (const auto & elem : _elem_side_map)
      {
        Elem * el = _mesh.elemPtr(elem.first);
        unsigned int side = elem.second;
        BoundaryID bnd_id = _mesh.getBoundaryIDs(el,side)[0];
        _slave_elem_id = el->id();   //elem.first
        if (_master_elem_id == _slave_elem_id)
          continue;
        std::map<unsigned int, std::vector<Point>> side_map{_viewfactor.getSideMap(el,side)};
        const Point center = _viewfactor.getCenterPoint(side_map);
        _u_slave = _system.point_value(_var_number, center, false);
        Real temp_func_slave = _u_slave * _u_slave * _u_slave * _u_slave;
        Real f_ms = _viewfactor.getViewFactor(current_boundary_id,_master_elem_id, bnd_id, _slave_elem_id);
        q_ms = _emissivity[current_boundary_id] * _stefan_boltzmann * f_ms * temp_func_master; // master-slave : from Modest
        q_sm = _emissivity[bnd_id] * _stefan_boltzmann * f_ms * temp_func_slave;  // slave-master :  from Modest
        q_net += q_ms - q_sm;
      }
    }
  }
  return _test[_i][_qp] * q_net;
}

Real
RadiativeHeatFluxBC::computeQpJacobian()
{
  // Real dq_net{0.0};
  // BoundaryID current_boundary_id = _mesh.getBoundaryIDs(_current_elem, _current_side)[0];
  // if (_boundary_ids.find(current_boundary_id)!=_boundary_ids.end())
  // {
  //   Real dtemp_func_master = 4 * _phi[_j][_qp] * _u[_qp] * _u[_qp] * _u[_qp];
  //   dq_net = _stefan_boltzmann * dtemp_func_master; // black body
  // }
  // return _test[_i][_qp] * dq_net;
  return 0;
}
