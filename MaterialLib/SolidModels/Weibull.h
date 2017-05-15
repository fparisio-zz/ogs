/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

/*
 * Implementation of Weibull's isotropic damage model.
 *
 * Refer to "Single-surface benchmark of OpenGeoSys documentation
 * (https://docs.opengeosys.org/docs/benchmarks/small-deformations/mechanics-plasticity-single-surface)"
 * for more details for the tests.
 */

#pragma once

#include <cfloat>
#include <memory>
#ifndef NDEBUG
#include <ostream>
#endif

#include <Eigen/Dense>
#include <logog/include/logog.hpp>
#include <utility>

#include "BaseLib/Error.h"
#include "NumLib/NewtonRaphson.h"

#include "KelvinVector.h"
#include "DamageBase.h"

namespace MaterialLib
{
namespace Solids
{
namespace Weibull
{
struct WeibullDamageProperties
{
    using P = ProcessLib::Parameter<double>;
    P const& alpha_d;
    P const& beta_d;
    P const& h_d;
};

template <int DisplacementDim>
class SolidWeibull final : public DamageBase<DisplacementDim>
{
public:
    static int const KelvinVectorSize =
        ProcessLib::KelvinVectorDimensions<DisplacementDim>::value;
    static int const JacobianResidualSize =
        2 * KelvinVectorSize + 3;  // 2 is the number of components in the
                                   // jacobian/residual, not the space
                                   // dimension. And 3 is for additional
                                   // variables.
    using ResidualVectorType = Eigen::Matrix<double, JacobianResidualSize, 1>;
    using JacobianMatrix = Eigen::Matrix<double, JacobianResidualSize,
                                         JacobianResidualSize, Eigen::RowMajor>;

public:
    //
    // Variables specific to the material model.
    //
    struct MaterialProperties
    {
        using P = ProcessLib::Parameter<double>;

        using KelvinVector = ProcessLib::KelvinVectorType<DisplacementDim>;
        using KelvinMatrix = ProcessLib::KelvinMatrixType<DisplacementDim>;

        MaterialProperties(P const& G_, P const& K_) : G(G_), K(K_) {}
        // basic material parameters
        P const& G;  ///< shear modulus
        P const& K;  ///< bulk modulus
    };

    struct MaterialStateVariables
        : public MechanicsBase<DisplacementDim>::MaterialStateVariables
    {
        MaterialStateVariables()
            : eps_p_D(KelvinVector::Zero()), eps_p_D_prev(KelvinVector::Zero())
        {
        }

        void setInitialConditions()
        {
            eps_p_D = eps_p_D_prev;
            eps_p_V = eps_p_V_prev;
            eps_p_eff = eps_p_eff_prev;
            kappa_d = kappa_d_prev;
            damage = damage_prev;
            lambda = 0;
        }

        void pushBackState() override
        {
            eps_p_D_prev = eps_p_D;
            eps_p_V_prev = eps_p_V;
            eps_p_eff_prev = eps_p_eff;  // effective part of trace(eps_p)
            kappa_d_prev = kappa_d;
            damage_prev = damage;
            lambda = 0;
        }

        double getLocalVariable() const override { return kappa_d; }
        using KelvinVector = ProcessLib::KelvinVectorType<DisplacementDim>;

        KelvinVector eps_p_D;  ///< deviatoric plastic strain
        double eps_p_V = 0;    ///< volumetric strain
        double eps_p_eff = 0;  ///< effective plastic strain
        double kappa_d = 0;    ///< damage driving variable
        double damage = 0;     ///< isotropic damage variable

        // Initial values from previous timestep
        KelvinVector eps_p_D_prev;  ///< \copydoc eps_p_D
        double eps_p_V_prev = 0;    ///< \copydoc eps_p_V
        double eps_p_eff_prev = 0;  ///< \copydoc eps_p_eff
        double kappa_d_prev = 0;    ///< \copydoc kappa_d
        double damage_prev = 0;     ///< \copydoc damage
        double lambda = 0;          ///< plastic multiplier

#ifndef NDEBUG
        friend std::ostream& operator<<(std::ostream& os,
                                        MaterialStateVariables const& m)
        {
            os << "State:\n"
               << "eps_p_D: " << m.eps_p_D << "\n"
               << "eps_p_eff: " << m.eps_p_eff << "\n"
               << "kappa_d: " << m.kappa_d << "\n"
               << "damage: " << m.damage << "\n"
               << "eps_p_D_prev: " << m.eps_p_D_prev << "\n"
               << "eps_p_eff_prev: " << m.eps_p_eff_prev << "\n"
               << "kappa_d_prev: " << m.kappa_d_prev << "\n"
               << "damage_prev: " << m.damage_prev << "\n"
               << "lambda: " << m.lambda << "\n";
            return os;
        }
#endif  // NDEBUG
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
    };

    std::unique_ptr<
        typename MechanicsBase<DisplacementDim>::MaterialStateVariables>
    createMaterialStateVariables() override
    {
        return std::unique_ptr<
            typename MechanicsBase<DisplacementDim>::MaterialStateVariables>{
            new MaterialStateVariables};
    }

public:
    using KelvinVector = ProcessLib::KelvinVectorType<DisplacementDim>;
    using KelvinMatrix = ProcessLib::KelvinMatrixType<DisplacementDim>;

public:
    explicit SolidWeibull(
        NumLib::NewtonRaphsonSolverParameters nonlinear_solver_parameters,
        MaterialProperties material_properties,
        std::unique_ptr<WeibullDamageProperties>&& damage_properties,
        bool const compute_local_damage)
        : _nonlinear_solver_parameters(std::move(nonlinear_solver_parameters)),
          _mp(std::move(material_properties)),
          _damage_properties(std::move(damage_properties)),
          _compute_local_damage(compute_local_damage)

    {
    }

    bool computeConstitutiveRelation(
        double const t,
        ProcessLib::SpatialPosition const& x,
        double const dt,
        KelvinVector const& eps_prev,
        KelvinVector const& eps,
        KelvinVector const& sigma_prev,
        KelvinVector& sigma,
        KelvinMatrix& C,
        typename MechanicsBase<DisplacementDim>::MaterialStateVariables&
            material_state_variables) override;

    /// Updates the internal damage based on the given kappa_damage value.
    /// \returns the new updated damage value.
    double updateDamage(
        double const t, ProcessLib::SpatialPosition const& x,
        double const kappa_damage,
        typename MechanicsBase<DisplacementDim>::MaterialStateVariables&
            material_state_variables) override;

#ifdef PROTOBUF_FOUND
    OGS::MaterialState writeMaterialState(
        typename MechanicsBase<DisplacementDim>::MaterialStateVariables const&
            material_state_variables) const override
    {
        assert(dynamic_cast<MaterialStateVariables const*>(
                   &material_state_variables) != nullptr);
        MaterialStateVariables const& state =
            static_cast<MaterialStateVariables const&>(
                material_state_variables);

        OGS::MaterialState material_state;
        auto weibull_material_state = material_state.mutable_weibull();

        auto eps_p_D = weibull_material_state->mutable_eps_p_d();
        eps_p_D->set_dimension(DisplacementDim);
        for (int i = 0; i < state.eps_p_D.size(); ++i)
            eps_p_D->add_value(state.eps_p_D[i]);

        weibull_material_state->set_eps_p_v(state.eps_p_V);
        weibull_material_state->set_eps_p_eff(state.eps_p_eff);
        weibull_material_state->set_kappa_d(state.kappa_d);
        weibull_material_state->set_damage(state.damage);

        return material_state;
    };
#endif  // PROTOBUF_FOUND

private:
    void calculateLocalKappaD(
        double const t, ProcessLib::SpatialPosition const& x, double const dt,
        KelvinVector const& sigma, KelvinVector const& eps_prev,
        KelvinVector const& eps,
        typename MechanicsBase<DisplacementDim>::MaterialStateVariables&
            material_state_variables);

private:
    NumLib::NewtonRaphsonSolverParameters const _nonlinear_solver_parameters;

    MaterialProperties _mp;
    std::unique_ptr<WeibullDamageProperties> _damage_properties;
    bool const _compute_local_damage;
};

}  // namespace Weibull
}  // namespace Solids
}  // namespace MaterialLib
#include "Weibull-impl.h"
