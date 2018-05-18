/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <memory>
#include <utility>

namespace MeshLib
{
class Element;
}

namespace ProcessLib
{
namespace SmallDeformationNonlocalHydroMechanics
{
template <int DisplacementDim>
struct SmallDeformationNonlocalHydroMechanicsProcessData
{
    SmallDeformationNonlocalHydroMechanicsProcessData(
        std::unique_ptr<MaterialLib::Solids::MechanicsBase<DisplacementDim>>&&
            material_,
        Parameter<double> const& intrinsic_permeability_,
        Parameter<double> const& specific_storage_,
        Parameter<double> const& fluid_viscosity_,
        Parameter<double> const& fluid_density_,
        Parameter<double> const& biot_coefficient_,
        Parameter<double> const& porosity_,
        Parameter<double> const& solid_density_,
        Eigen::Matrix<double, DisplacementDim, 1>
            specific_body_force_,
            double const internal_length_)
        : material{std::move(material_)},
          intrinsic_permeability(intrinsic_permeability_),
          specific_storage(specific_storage_),
          fluid_viscosity(fluid_viscosity_),
          fluid_density(fluid_density_),
          biot_coefficient(biot_coefficient_),
          porosity(porosity_),
          solid_density(solid_density_),
          specific_body_force(std::move(specific_body_force_)),
          internal_length_squared(internal_length_ * internal_length_)
    {
    }

    SmallDeformationNonlocalHydroMechanicsProcessData(SmallDeformationNonlocalHydroMechanicsProcessData&& other)
        : material{std::move(other.material)},
          intrinsic_permeability(other.intrinsic_permeability),
          specific_storage(other.specific_storage),
          fluid_viscosity(other.fluid_viscosity),
          fluid_density(other.fluid_density),
          biot_coefficient(other.biot_coefficient),
          porosity(other.porosity),
          solid_density(other.solid_density),
          specific_body_force(other.specific_body_force),
          dt(other.dt),
          t(other.t),
          internal_length_squared{other.internal_length_squared}
    {
    }

    //! Copies are forbidden.
    SmallDeformationNonlocalHydroMechanicsProcessData(
        SmallDeformationNonlocalHydroMechanicsProcessData const&) = delete;

    //! Assignments are not needed.
    void operator=(SmallDeformationNonlocalHydroMechanicsProcessData const&) = delete;

    //! Assignments are not needed.
    void operator=(SmallDeformationNonlocalHydroMechanicsProcessData&&) = delete;

    /// The constitutive relation for the mechanical part.
    /// \note Linear elasticity is the only supported one in the moment.
    std::unique_ptr<MaterialLib::Solids::MechanicsBase<DisplacementDim>>
        material;
    /// Permeability of the solid. A scalar quantity, Parameter<double>.
    Parameter<double> const& intrinsic_permeability;
    /// Volumetric average specific storage of the solid and fluid phases.
    /// A scalar quantity, Parameter<double>.
    Parameter<double> const& specific_storage;
    /// Fluid's viscosity. A scalar quantity, Parameter<double>.
    Parameter<double> const& fluid_viscosity;
    /// Fluid's density. A scalar quantity, Parameter<double>.
    Parameter<double> const& fluid_density;
    /// Biot coefficient. A scalar quantity, Parameter<double>.
    Parameter<double> const& biot_coefficient;
    /// Porosity of the solid. A scalar quantity, Parameter<double>.
    Parameter<double> const& porosity;
    /// Solid's density. A scalar quantity, Parameter<double>.
    Parameter<double> const& solid_density;
    /// Specific body forces applied to solid and fluid.
    /// It is usually used to apply gravitational forces.
    /// A vector of displacement dimension's length.
    Eigen::Matrix<double, DisplacementDim, 1> const specific_body_force;
    double dt = 0.0;
    double t = 0.0;
    double const internal_length_squared;

    double injected_volume = 0.0;
    double crack_volume_old = 0.0;
    double crack_volume = 0.0;
    bool propagating_crack = true;

    double stiffness = 2.15e11;

    double pressure = 0.0;
    double pressure_old = 0.0;
    double pressure_error = 0.0;
};

}  // namespace SmallDeformationNonlocalHydroMechanics
}  // namespace ProcessLib
