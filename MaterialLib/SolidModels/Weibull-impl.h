/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

/**
 * Common convenitions for naming:
 * x_D              - deviatoric part of tensor x
 * x_V              - volumetric part of tensor x
 * x_p              - a variable related to plastic potential
 * x_prev           - value of x in previous time step
 *e
 * Variables used in the code:
 * eps_D            - deviatoric strain
 * eps_p_D_dot      - deviatoric increment of plastic strain
 * eps_p_eff_dot    - increment of effective plastic strain
 * eps_p_V_dot      - volumetric increment of plastic strain
 * sigma_D_inverse_D - deviatoric part of sigma_D_inverse
 *
 */
#pragma once

#include <boost/math/special_functions/pow.hpp>
#include <logog/include/logog.hpp>
#include "MaterialLib/SolidModels/KelvinVector.h"

namespace MaterialLib
{
namespace Solids
{
namespace Weibull
{
/// Special product of \c v with itself: \f$v \odot v\f$.
/// The tensor \c v is given in Kelvin mapping.
/// \note Implementation only for 2 and 3 dimensions.
/// \attention Pay attention to the sign of the result, which normally would be
/// negative, but the returned value is not negated. This has to do with \f$
/// d(A^{-1})/dA = -A^{-1} \odot A^{-1} \f$.
template <int DisplacementDim>
ProcessLib::KelvinMatrixType<DisplacementDim> sOdotS(
    ProcessLib::KelvinVectorType<DisplacementDim> const& v);

template <int DisplacementDim>
struct PhysicalStressWithInvariants final
{
    static int const KelvinVectorSize =
        ProcessLib::KelvinVectorDimensions<DisplacementDim>::value;
    using Invariants = MaterialLib::SolidModels::Invariants<KelvinVectorSize>;
    using KelvinVector = ProcessLib::KelvinVectorType<DisplacementDim>;

    explicit PhysicalStressWithInvariants(KelvinVector const& stress)
        : value{stress},
          D{Invariants::deviatoric_projection * stress},
          I_1{Invariants::trace(stress)},
          J_2{Invariants::J2(D)},
          J_3{Invariants::J3(D)}
    {
    }

    PhysicalStressWithInvariants(PhysicalStressWithInvariants const&) = default;
    PhysicalStressWithInvariants& operator=(
        PhysicalStressWithInvariants const&) = default;
#if defined(_MSC_VER) && (_MSC_VER >= 1900)
    PhysicalStressWithInvariants(PhysicalStressWithInvariants&&) = default;
    PhysicalStressWithInvariants& operator=(PhysicalStressWithInvariants&&) =
        default;
#endif  // _MSC_VER

    KelvinVector value;
    KelvinVector D;
    double I_1;
    double J_2;
    double J_3;
};

template <int DisplacementDim>
double SolidWeibull<DisplacementDim>::vectorNorm(KelvinVector const& eps,
                                                 normType const flag)
{
    Eigen::Matrix<double, DisplacementDim, DisplacementDim> strain_mat =
        Eigen::Matrix<double, DisplacementDim, DisplacementDim>::Zero();

    for (int i = 0; i < DisplacementDim; ++i)
        for (int j = 0; j < DisplacementDim; ++j)
        {
            if (i == j)
            {
                strain_mat(i, j) = eps(i);
            }
            else
            {
                strain_mat(i, j) = eps(i + j + 2);
            };
        };

    Eigen::EigenSolver<decltype(strain_mat)> eigen_solver(strain_mat);
    auto const& principal_strain = eigen_solver.eigenvalues();
    // building kappa_d (damage driving variable)
    double prod_strain_tens = 0.;
    double prod_strain_comp = 0.;
    switch (flag)
    {
        case normType::positive:
            for (int i = 0; i < DisplacementDim; ++i)
            {
                Eigen::Matrix<double, DisplacementDim, 1> eig_val_strain;
                eig_val_strain(i, 0) = real(principal_strain(i, 0));
                if (eig_val_strain(i, 0) > 0)
                    prod_strain_tens =
                        prod_strain_tens +
                        eig_val_strain(i, 0) * eig_val_strain(i, 0);
            }
            return std::sqrt(prod_strain_tens);
            break;

        case normType::negative:
            for (int i = 0; i < DisplacementDim; ++i)
            {
                Eigen::Matrix<double, DisplacementDim, 1> eig_val_strain;
                eig_val_strain(i, 0) = real(principal_strain(i, 0));
                if (eig_val_strain(i, 0) < 0)
                    prod_strain_tens =
                        prod_strain_tens +
                        eig_val_strain(i, 0) * eig_val_strain(i, 0);
            }
            return std::sqrt(prod_strain_tens);
            break;

        case normType::total:
            for (int i = 0; i < DisplacementDim; ++i)
            {
                Eigen::Matrix<double, DisplacementDim, 1> eig_val_strain;
                eig_val_strain(i, 0) = real(principal_strain(i, 0));
                prod_strain_tens = prod_strain_tens +
                                   eig_val_strain(i, 0) * eig_val_strain(i, 0);
            }
            return std::sqrt(prod_strain_tens);
            break;
    }
};

template <int DisplacementDim>
void SolidWeibull<DisplacementDim>::calculateLocalKappaD(
    double const t, ProcessLib::SpatialPosition const& x, double const dt,
    KelvinVector const& sigma, KelvinVector const& eps_prev,
    KelvinVector const& eps,
    typename MechanicsBase<DisplacementDim>::MaterialStateVariables&
        material_state_variables)
{
    assert(dynamic_cast<MaterialStateVariables*>(&material_state_variables) !=
           nullptr);
    auto& _state =
        static_cast<MaterialStateVariables&>(material_state_variables);

    // Default case of the rate problem. Updated below if volumetric plastic
    // strain rate is positive (dilatancy).
    _state.kappa_d = _state.kappa_d_prev;

    // Compute damage current step
    double const eps_p_V_dot = _state.eps_p_V - _state.eps_p_V_prev;

    // Compute elastic-damage simple model
    // compute norm in principal strain space

    double const r_d = _damage_properties->r_d(t, x)[0];
    double const alpha_d = _damage_properties->alpha_d(t, x)[0];
    double const m_w = _damage_properties->m_w(t, x)[0];
    double const G = _mp.G(t, x)[0];

    double prod_strain_tens = vectorNorm(eps, normType::positive);
    double d_prod_strain_tens = vectorNorm(eps - eps_prev, normType::positive);
    double prod_strain_comp = vectorNorm(eps, normType::negative);
    double d_prod_strain_comp = vectorNorm(eps - eps_prev, normType::negative);

    std::cout << "tensile\n" << prod_strain_tens << std::endl;
    std::cout << "compressive\n" << prod_strain_comp << std::endl;
    std::cout << "check\n"
              << prod_strain_comp + prod_strain_tens -
                     vectorNorm(eps, normType::total)
              << std::endl;

    std::cout << "tensile rate\n" << d_prod_strain_tens << std::endl;
    std::cout << "compressive rate\n" << d_prod_strain_comp << std::endl;

    double strain_norm_0 =
        alpha_d *
        std::pow(std::log(1. / (1. - _state.damage / ((1. - 1e-8)))), 1. / m_w);
    if ((std::sqrt(prod_strain_tens + r_d * prod_strain_comp) -
         strain_norm_0) >= 0)
        if (std::sqrt(d_prod_strain_tens + r_d * d_prod_strain_comp) > 0)
            _state.kappa_d +=
                std::sqrt(d_prod_strain_tens + r_d * d_prod_strain_comp);

    assert(_state.kappa_d >= 0.);
}

template <int DisplacementDim>
double SolidWeibull<DisplacementDim>::updateDamage(
    double const t, ProcessLib::SpatialPosition const& x,
    double const kappa_damage,
    typename MechanicsBase<DisplacementDim>::MaterialStateVariables&
        material_state_variables)
{
    assert(dynamic_cast<MaterialStateVariables*>(&material_state_variables) !=
           nullptr);
    MaterialStateVariables& _state =
        static_cast<MaterialStateVariables&>(material_state_variables);

    double const alpha_d = _damage_properties->alpha_d(t, x)[0];
    double const m_w = _damage_properties->m_w(t, x)[0];

    // Update internal damage variable.
    _state.damage =
        (1. - 1e-8) * (1 - std::exp(-std::pow(kappa_damage / alpha_d, m_w)));

    return _state.damage;
}

template <int DisplacementDim>
typename SolidWeibull<DisplacementDim>::KelvinVector predict_sigma(
    double const G, double const K,
    typename SolidWeibull<DisplacementDim>::KelvinVector const& sigma_prev,
    typename SolidWeibull<DisplacementDim>::KelvinVector const& eps,
    typename SolidWeibull<DisplacementDim>::KelvinVector const& eps_prev,
    double const eps_V)
{
    static int const KelvinVectorSize =
        ProcessLib::KelvinVectorDimensions<DisplacementDim>::value;
    using Invariants = MaterialLib::SolidModels::Invariants<KelvinVectorSize>;
    auto const& P_dev = Invariants::deviatoric_projection;

    // dimensionless initial hydrostatic pressure
    double const pressure_prev = Invariants::trace(sigma_prev) / (-3. * G);
    // initial strain invariant
    double const e_prev = Invariants::trace(eps_prev);
    // dimensioness hydrostatic stress increment
    double const pressure = pressure_prev - K / G * (eps_V - e_prev);
    // dimensionless deviatoric initial stress
    typename SolidWeibull<DisplacementDim>::KelvinVector const sigma_D_prev =
        P_dev * sigma_prev / G;
    // dimensionless deviatoric stress
    typename SolidWeibull<DisplacementDim>::KelvinVector const sigma_D =
        sigma_D_prev + 2 * P_dev * (eps - eps_prev);
    return sigma_D - pressure * Invariants::identity2;
}

template <int DisplacementDim>
bool SolidWeibull<DisplacementDim>::computeConstitutiveRelation(
    double const t,
    ProcessLib::SpatialPosition const& x,
    double const dt,
    KelvinVector const& eps_prev,
    KelvinVector const& eps,
    KelvinVector const& sigma_prev,
    KelvinVector& sigma_final,
    KelvinMatrix& C,
    typename MechanicsBase<DisplacementDim>::MaterialStateVariables&
        material_state_variables)
{
    assert(dynamic_cast<MaterialStateVariables*>(&material_state_variables) !=
           nullptr);
    auto& _state =
        static_cast<MaterialStateVariables&>(material_state_variables);
    _state.setInitialConditions();

    using Invariants = MaterialLib::SolidModels::Invariants<KelvinVectorSize>;
    C.setZero();

    // volumetric strain
    double const eps_V = Invariants::trace(eps);

    auto const& P_dev = Invariants::deviatoric_projection;
    // deviatoric strain
    KelvinVector const eps_D = P_dev * eps;

    // dimensionless stress/hydrostatic pressure
    double const G = _mp.G(t, x)[0];
    double const K = _mp.K(t, x)[0];

    KelvinVector sigma_eff_prev = sigma_prev;  // In case without damage the
                                               // effective value is same as the
                                               // previous one.
    if (_damage_properties)  // Valid for local and non-local damage
    {
        // Compute sigma_eff from damage total stress sigma, which is given by
        // sigma_eff=sigma_prev / (1-damage)
        sigma_eff_prev = sigma_prev / (1 - _state.damage_prev);
    }
    KelvinVector sigma = predict_sigma<DisplacementDim>(G, K, sigma_eff_prev,
                                                        eps, eps_prev, eps_V);

    PhysicalStressWithInvariants<DisplacementDim> s{G * sigma};

    if (_damage_properties)
    {
        calculateLocalKappaD(t, x, dt, sigma, eps, eps_prev, _state);

        if (_compute_local_damage)  // The non-local damage update is called
                                    // from the FEM.
            updateDamage(t, x, _state.kappa_d, material_state_variables);

        // Compute damage update of consistent tangent matrix
        // compute principal strains etc.
        double prod_strain_tens = vectorNorm(eps, normType::positive);
        double prod_strain_tens_prev = vectorNorm(eps_prev, normType::positive);
        double prod_strain_comp = vectorNorm(eps, normType::negative);
        double prod_strain_comp_prev = vectorNorm(eps_prev, normType::negative);

        double d_prod_strain_tens = prod_strain_tens - prod_strain_tens_prev;
        double d_prod_strain_comp = prod_strain_comp - prod_strain_comp_prev;

        double const r_d = _damage_properties->r_d(t, x)[0];
        double const alpha_d = _damage_properties->alpha_d(t, x)[0];
        double const m_w = _damage_properties->m_w(t, x)[0];
        double g_prime = 0;
        KelvinVector eta_dam = eps;

        if (_state.kappa_d > 0)
        {
            g_prime = (1. - 1e-8) * m_w *
                      std::exp(-std::pow(_state.kappa_d / alpha_d, m_w)) *
                      std::pow(_state.kappa_d / alpha_d, m_w) / _state.kappa_d;
            eta_dam =
                eta_dam / std::sqrt(prod_strain_tens + r_d * prod_strain_comp);
        }
        else
        {
            g_prime = 0;
        }

        Eigen::Matrix<double, KelvinVectorSize, KelvinVectorSize,
                      Eigen::RowMajor>
            dam_matrix =
                g_prime * (1 - _state.damage) * sigma * eta_dam.transpose();

        C.setZero();
        C.template topLeftCorner<3, 3>().setConstant(K - 2. / 3 * G);
        C.noalias() +=
            (1 - _state.damage) * 2 * G * KelvinMatrix::Identity() - dam_matrix;
    }

    // Update sigma.
    if (_damage_properties && _compute_local_damage)
        sigma_final = G * sigma * (1. - _state.damage);
    else
        sigma_final.noalias() = G * sigma;  // Plastic part only

    return true;
}

}  // namespace Weibull
}  // namespace Solids
}  // namespace MaterialLib
