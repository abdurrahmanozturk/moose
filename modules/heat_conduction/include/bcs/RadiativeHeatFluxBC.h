#ifndef RADIATIVEHEATFLUXBC_H
#define RADIATIVEHEATFLUXBC_H

#include "IntegratedBC.h"
#include "ViewFactor.h"

// Forward declarations
class RadiativeHeatFluxBC;

template <>
InputParameters validParams<RadiativeHeatFluxBC>();

/**
 * Base class for deriving any boundary condition of a integrated type
 */
class RadiativeHeatFluxBC : public IntegratedBC
{
public:
  RadiativeHeatFluxBC(const InputParameters & parameters);
  const std::map<unsigned int, std::vector<Real>> getSideMap(const Elem * elem,const unsigned int side);
  const Real getArea(const Elem * elem,const unsigned int side);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  const unsigned int _var_number;
  const System & _system;
  // Point _point;
  // Real _value;

private:
  const ViewFactor & _viewfactor;
  const std::set<BoundaryID> & _boundary_ids;
  const std::set<BoundaryID> & _master_boundary_ids;
  const std::set<BoundaryID> & _slave_boundary_ids;
  const Real _stefan_boltzmann;
  std::map<BoundaryID, Real> _emissivity;
  std::map<unsigned int, unsigned int> _elem_side_map;
  unsigned int _master_elem_id,_slave_elem_id;
};

#endif /* RADIATIVEHEATFLUXBC_H */
