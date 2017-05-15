#pragma once

#include "MechanicsBase.h"
#include "ProcessLib/Parameter/SpatialPosition.h"

namespace MaterialLib
{
namespace Solids
{
template <int DisplacementDim>
struct DamageBase : MechanicsBase<DisplacementDim>
{
    virtual double updateDamage(
        double const t, ProcessLib::SpatialPosition const& x,
        double const kappa_damage,
        typename MechanicsBase<DisplacementDim>::MaterialStateVariables&
            material_state_variables) = 0;
};
}
}
