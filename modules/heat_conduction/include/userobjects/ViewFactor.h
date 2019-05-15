#ifndef VIEWFACTOR_H
#define VIEWFACTOR_H

#include "ViewFactorBase.h"

// Forward Declarations
class ViewFactor;

template <>
InputParameters validParams<ViewFactor>();

class ViewFactor : public ViewFactorBase
{
public:
  ViewFactor(const InputParameters & parameters);
  virtual void initialize() override;
  virtual void execute() override;
  virtual void finalize() override;
  virtual void threadJoin(const UserObject & y) override {}

  virtual Real getViewFactor(BoundaryID master_bnd, unsigned int master_elem,
                             BoundaryID slave_bnd, unsigned int slave_elem) const;

protected:
  const unsigned int _samplingNumber,_sourceNumber;
  const std::string _method;
  //
  std::map<unsigned int, unsigned int> _elem_side_map;
  std::map<unsigned int, std::vector<Point>> _master_side_map;
  std::map<unsigned int, std::vector<Point>> _slave_side_map;
};

#endif // VIEWFACTOR_H
