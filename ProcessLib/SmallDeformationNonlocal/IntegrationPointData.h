/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "LocalAssemblerInterface.h"

namespace ProcessLib
{
namespace SmallDeformationNonlocal
{

template <typename BMatricesType, int DisplacementDim>
struct IntegrationPointData final
{
    explicit IntegrationPointData(
        MaterialLib::Solids::MechanicsBase<DisplacementDim>& solid_material)
        : _solid_material(solid_material),
          _material_state_variables(
              _solid_material.createMaterialStateVariables())
    {
        if (auto const msv =
                dynamic_cast<typename MaterialLib::Solids::Ehlers::SolidEhlers<
                    DisplacementDim>::MaterialStateVariables*>(
                    _material_state_variables.get()))
        {
            _eps_p_V = &msv->eps_p_V;
            _eps_p_D_xx = &(msv->eps_p_D[0]);
        }
    }

#if defined(_MSC_VER) && _MSC_VER < 1900
    // The default generated move-ctor is correctly generated for other
    // compilers.
    explicit IntegrationPointData(IntegrationPointData&& other)
        : _b_matrices(std::move(other._b_matrices)),
          _sigma(std::move(other._sigma)),
          _sigma_prev(std::move(other._sigma_prev)),
          _eps(std::move(other._eps)),
          _eps_prev(std::move(other._eps_prev)),
          _solid_material(other._solid_material),
          _material_state_variables(std::move(other._material_state_variables)),
          _C(std::move(other._C)),
          _detJ(std::move(other._detJ)),
          _integralMeasure(other._integralMeasure)
    {
    }
#endif  // _MSC_VER

    typename BMatricesType::BMatrixType _b_matrices;
    typename BMatricesType::KelvinVectorType _sigma, _sigma_prev;
    typename BMatricesType::KelvinVectorType _eps, _eps_prev;
    double _damage = 0;

    MaterialLib::Solids::MechanicsBase<DisplacementDim>& _solid_material;
    std::unique_ptr<typename MaterialLib::Solids::MechanicsBase<
        DisplacementDim>::MaterialStateVariables>
        _material_state_variables;

    typename BMatricesType::KelvinMatrixType _C;
    double _detJ;
    double _integralMeasure;

    double const* _eps_p_V;
    double const* _eps_p_D_xx;

    void pushBackState()
    {
        _eps_prev = _eps;
        _sigma_prev = _sigma;
        _material_state_variables->pushBackState();
    }

    double getLocalVariable() const
    {
        return _material_state_variables->getLocalVariable();
    }

    double updateDamage(double const t, SpatialPosition const& x_position,
                      double const kappa_d)
    {
        return static_cast<
                   MaterialLib::Solids::Ehlers::SolidEhlers<DisplacementDim>&>(
                   _solid_material)
            .updateDamage(t, x_position, kappa_d, *_material_state_variables);
    }

    std::vector<std::tuple<
        // element's local assembler
        SmallDeformationNonlocalLocalAssemblerInterface const* const,
        int,     // integration point id,
        double,  // squared distance to current integration point
        double   // alpha_kl
        >>
        non_local_assemblers;
};

}  // namespace SmallDeformationNonlocal
}  // namespace ProcessLib