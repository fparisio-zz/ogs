/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

/*
 * Implementation of Ehler's single-surface model.
 * see Ehler's paper "A single-surface yield function for geomaterials" for more
 * details. \cite{Ehlers1995}
 *
 * Refer to "Single-surface benchmark of OpenGeoSys documentation
 * (https://docs.opengeosys.org/docs/benchmarks/small-deformations/mechanics-plasticity-single-surface)"
 * for more details for the tests.
 */

#pragma once

#ifndef NDEBUG
#include <ostream>
#endif

#include "BaseLib/Error.h"
#include "NumLib/NewtonRaphson.h"
#include "ProcessLib/Parameter/Parameter.h"

#include "KelvinVector.h"
#include "MechanicsBase.h"

namespace MaterialLib
{
namespace Solids
{
namespace Ehlers
{
//
// Variables specific to the material model.
//
struct MaterialProperties
{
    using P = ProcessLib::Parameter<double>;

    /// material parameters in relation to Ehler's single-surface model
    /// see Ehler's paper "A single-surface yield function for geomaterials"
    /// for more details.
    MaterialProperties(P const& G_, P const& K_, P const& alpha_,
                       P const& beta_, P const& gamma_, P const& delta_,
                       P const& epsilon_, P const& m_, P const& alpha_p_,
                       P const& beta_p_, P const& gamma_p_, P const& delta_p_,
                       P const& epsilon_p_, P const& m_p_, P const& kappa_,
                       P const& hardening_coefficient_)
        : G(G_),
          K(K_),
          alpha(alpha_),
          beta(beta_),
          gamma(gamma_),
          delta(delta_),
          epsilon(epsilon_),
          m(m_),
          alpha_p(alpha_p_),
          beta_p(beta_p_),
          gamma_p(gamma_p_),
          delta_p(delta_p_),
          epsilon_p(epsilon_p_),
          m_p(m_p_),
          kappa(kappa_),
          hardening_coefficient(hardening_coefficient_)
    {
    }
    // basic material parameters
    P const& G;  ///< shear modulus
    P const& K;  ///< bulk modulus

    P const& alpha;
    P const& beta;
    P const& gamma;
    P const& delta;
    P const& epsilon;
    P const& m;

    P const& alpha_p;
    P const& beta_p;
    P const& gamma_p;
    P const& delta_p;
    P const& epsilon_p;
    P const& m_p;

    P const& kappa;
    P const& hardening_coefficient;
};

struct DamageProperties
{
    using P = ProcessLib::Parameter<double>;
    P const& alpha_d;
    P const& beta_d;
    P const& h_d;
};

/// Evaluated MaterialProperties container.
struct MaterialPropertiesV final
{
    MaterialPropertiesV(double const t, ProcessLib::SpatialPosition const& x,
                        MaterialProperties const& mp)
        : G(mp.G(t, x)[0]),
          K(mp.K(t, x)[0]),
          alpha(mp.alpha(t, x)[0]),
          beta(mp.beta(t, x)[0]),
          gamma(mp.gamma(t, x)[0]),
          delta(mp.delta(t, x)[0]),
          epsilon(mp.epsilon(t, x)[0]),
          m(mp.m(t, x)[0]),
          alpha_p(mp.alpha_p(t, x)[0]),
          beta_p(mp.beta_p(t, x)[0]),
          gamma_p(mp.gamma_p(t, x)[0]),
          delta_p(mp.delta_p(t, x)[0]),
          epsilon_p(mp.epsilon_p(t, x)[0]),
          m_p(mp.m_p(t, x)[0]),
          kappa(mp.kappa(t, x)[0]),
          hardening_coefficient(mp.hardening_coefficient(t, x)[0])
    {
    }
    // basic material parameters
    double const G;  ///< shear modulus
    double const K;  ///< bulk modulus

    double const alpha;
    double const beta;
    double const gamma;
    double const delta;
    double const epsilon;
    double const m;

    double const alpha_p;
    double const beta_p;
    double const gamma_p;
    double const delta_p;
    double const epsilon_p;
    double const m_p;

    double const kappa;
    double const hardening_coefficient;
};

/// Evaluated DamageProperties container.
struct DamagePropertiesV
{
    DamagePropertiesV(double const t,
                      ProcessLib::SpatialPosition const& x,
                      DamageProperties const& dp)
        : alpha_d(dp.alpha_d(t, x)[0]),
          beta_d(dp.beta_d(t, x)[0]),
          h_d(dp.h_d(t, x)[0])
    {
    }
    double const alpha_d;
    double const beta_d;
    double const h_d;
};

template <typename KelvinVector>
struct PlasticStrain final
{
    PlasticStrain() : D(KelvinVector::Zero()) {}
    PlasticStrain(KelvinVector const& eps_p_D_, double const eps_p_V_,
                  double const eps_p_eff_)
        : D(eps_p_D_), V(eps_p_V_), eff(eps_p_eff_){};

    KelvinVector D;  ///< deviatoric plastic strain
    double V = 0;    ///< volumetric strain
    double eff = 0;  ///< effective plastic strain
};

class Damage final
{
public:
    Damage(){};
    Damage(double const kappa_d, double const value)
        : _kappa_d(kappa_d), _value(value)
    {
    }

    double kappa_d() const { return _kappa_d; }
    double value() const { return _value; }

private:
    double _kappa_d = 0;  ///< damage driving variable
    double _value = 0;    ///< isotropic damage variable
};

template <int DisplacementDim>
struct StateVariables
    : public MechanicsBase<DisplacementDim>::MaterialStateVariables
{
    StateVariables& operator=(StateVariables const&) = default;
    typename MechanicsBase<DisplacementDim>::MaterialStateVariables& operator=(
        typename MechanicsBase<DisplacementDim>::MaterialStateVariables const&
            state) noexcept override
    {
        return operator=(static_cast<StateVariables const&>(state));
    }

    void setInitialConditions()
    {
        eps_p = eps_p_prev;
        damage = damage_prev;
    }

    void pushBackState() override
    {
        eps_p_prev = eps_p;
        damage_prev = damage;
    }

    double getLocalVariable() const override { return damage.kappa_d(); }

    using KelvinVector = ProcessLib::KelvinVectorType<DisplacementDim>;

    PlasticStrain<KelvinVector> eps_p;  ///< plastic part of the state.
    Damage damage;                      ///< damage part of the state.

    // Initial values from previous timestep
    PlasticStrain<KelvinVector> eps_p_prev;  ///< \copydoc eps_p
    Damage damage_prev;                      ///< \copydoc damage

#ifndef NDEBUG
    friend std::ostream& operator<<(
        std::ostream& os, StateVariables<DisplacementDim> const& m)
    {
        os << "State:\n"
           << "eps_p_D: " << m.eps_p.D << "\n"
           << "eps_p_eff: " << m.eps_p.eff << "\n"
           << "kappa_d: " << m.damage.kappa_d << "\n"
           << "damage: " << m.damage.damage << "\n"
           << "eps_p_D_prev: " << m.eps_p_prev.D << "\n"
           << "eps_p_eff_prev: " << m.eps_p_prev.eff << "\n"
           << "kappa_d_prev: " << m.damage_prev.kappa_d << "\n"
           << "damage_prev: " << m.damage_prev.damage << "\n"
           << "lambda: " << m.lambda << "\n";
        return os;
    }
#endif  // NDEBUG
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
};

template <int DisplacementDim>
class SolidEhlers final : public MechanicsBase<DisplacementDim>
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
    std::unique_ptr<
        typename MechanicsBase<DisplacementDim>::MaterialStateVariables>
    createMaterialStateVariables() override
    {
        return std::unique_ptr<StateVariables<DisplacementDim>>{
            new StateVariables<DisplacementDim>};
    }

public:
    using KelvinVector = ProcessLib::KelvinVectorType<DisplacementDim>;
    using KelvinMatrix = ProcessLib::KelvinMatrixType<DisplacementDim>;

public:
    explicit SolidEhlers(
        NumLib::NewtonRaphsonSolverParameters nonlinear_solver_parameters,
        MaterialProperties material_properties,
        std::unique_ptr<DamageProperties>&& damage_properties,
        bool const compute_local_damage)
        : _nonlinear_solver_parameters(std::move(nonlinear_solver_parameters)),
          _mp(std::move(material_properties)),
          _damage_properties(std::move(damage_properties)),
          _compute_local_damage(compute_local_damage)
    {
    }

    std::tuple<KelvinVector, std::unique_ptr<typename MechanicsBase<
                                 DisplacementDim>::MaterialStateVariables>,
               std::unique_ptr<KelvinMatrix>>
    integrateStress(
        double const t,
        ProcessLib::SpatialPosition const& x,
        double const dt,
        KelvinVector const& eps_prev,
        KelvinVector const& eps,
        KelvinVector const& sigma_prev,
        typename MechanicsBase<DisplacementDim>::MaterialStateVariables const&
            material_state_variables) override;

#ifdef PROTOBUF_FOUND
    OGS::MaterialState writeMaterialState(
        typename MechanicsBase<DisplacementDim>::MaterialStateVariables const&
            material_state_variables) const override
    {
        assert(dynamic_cast<StateVariables<DisplacementDim> const*>(
                   &material_state_variables) != nullptr);
        auto const& state =
            static_cast<StateVariables<DisplacementDim> const&>(
                material_state_variables);

        OGS::MaterialState material_state;
        auto ehlers_material_state = material_state.mutable_ehlers();

        auto eps_p_D = ehlers_material_state->mutable_eps_p_d();
        eps_p_D->set_dimension(DisplacementDim);
        for (int i = 0; i < state.eps_p.D.size(); ++i)
            eps_p_D->add_value(state.eps_p.D[i]);

        ehlers_material_state->set_eps_p_v(state.eps_p.V);
        ehlers_material_state->set_eps_p_eff(state.eps_p.eff);
        ehlers_material_state->set_kappa_d(state.damage.kappa_d());
        ehlers_material_state->set_damage(state.damage.value());

        return material_state;
    };
#endif  // PROTOBUF_FOUND

private:
    NumLib::NewtonRaphsonSolverParameters const _nonlinear_solver_parameters;

    MaterialProperties _mp;
    std::unique_ptr<DamageProperties> _damage_properties;
    bool const _compute_local_damage;
};

}  // namespace Ehlers
}  // namespace Solids
}  // namespace MaterialLib
#include "Ehlers-impl.h"
