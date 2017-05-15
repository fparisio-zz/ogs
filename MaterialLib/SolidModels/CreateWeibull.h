/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <logog/include/logog.hpp>

#include "ProcessLib/Utils/ProcessUtils.h"  // required for findParameter

#include "Weibull.h"
#include "MechanicsBase.h"

namespace MaterialLib
{
namespace Solids
{
namespace Weibull
{
inline NumLib::NewtonRaphsonSolverParameters
createNewtonRaphsonSolverParameters(BaseLib::ConfigTree const& config)
{
    DBUG("Create local nonlinear solver parameters.");
    auto const maximum_iterations =
        //! \ogs_file_param{material__solid__constitutive_relation__Weibull__nonlinear_solver__maximum_iterations}
        config.getConfigParameter<int>("maximum_iterations");

    DBUG("\tmaximum_iterations: %d.", maximum_iterations);

    auto const error_tolerance =
        //! \ogs_file_param{material__solid__constitutive_relation__Weibull__nonlinear_solver__error_tolerance}
        config.getConfigParameter<double>("error_tolerance");

    DBUG("\terror_tolerance: %g.", error_tolerance);

    return {maximum_iterations, error_tolerance};
}

inline std::unique_ptr<WeibullDamageProperties> createDamageProperties(
    std::vector<std::unique_ptr<ProcessLib::ParameterBase>> const& parameters,
    BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param_special{material__solid__constitutive_relation__Weibull__damage_properties__alpha_d}
    auto& alpha_d =
        ProcessLib::findParameter<double>(config, "alpha_d", parameters, 1);

    DBUG("Use \'%s\' as alpha_d.", alpha_d.name.c_str());

    //! \ogs_file_param_special{material__solid__constitutive_relation__Weibull__damage_properties__beta_d}
    auto& beta_d =
        ProcessLib::findParameter<double>(config, "beta_d", parameters, 1);

    DBUG("Use \'%s\' as beta_d.", beta_d.name.c_str());

    //! \ogs_file_param_special{material__solid__constitutive_relation__Weibull__damage_properties__h_d}
    auto& h_d = ProcessLib::findParameter<double>(config, "h_d", parameters, 1);

    DBUG("Use \'%s\' as h_d.", h_d.name.c_str());

    return std::unique_ptr<WeibullDamageProperties>{
        new WeibullDamageProperties{alpha_d, beta_d, h_d}};
}

template <int DisplacementDim>
std::unique_ptr<MechanicsBase<DisplacementDim>> createWeibull(
    std::vector<std::unique_ptr<ProcessLib::ParameterBase>> const& parameters,
    BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{material__solid__constitutive_relation__type}
    config.checkConfigParameter("type", "Weibull");
    DBUG("Create Weibull material");

    //! \ogs_file_param_special{material__solid__constitutive_relation__Weibull__shear_modulus}
    auto& shear_modulus = ProcessLib::findParameter<double>(
        config, "shear_modulus", parameters, 1);

    DBUG("Use \'%s\' as shear modulus parameter.", shear_modulus.name.c_str());

    //! \ogs_file_param_special{material__solid__constitutive_relation__Weibull__bulk_modulus}
    auto& bulk_modulus = ProcessLib::findParameter<double>(
        config, "bulk_modulus", parameters, 1);

    DBUG("Use \'%s\' as bulk modulus parameter.", bulk_modulus.name.c_str());

    typename SolidWeibull<DisplacementDim>::MaterialProperties mp{
        shear_modulus, bulk_modulus};

    // Damage properties.
    std::unique_ptr<WeibullDamageProperties> weibull_damage_properties;

    auto const& weibull_damage_config =
        //! \ogs_file_param{material__solid__constitutive_relation__Weibull__damage_properties}
        config.getConfigSubtreeOptional("damage_properties");
    if (weibull_damage_config)
    {
        weibull_damage_properties =
            createDamageProperties(parameters, *weibull_damage_config);
    }

    auto const& nonlinear_solver_config =
        //! \ogs_file_param{material__solid__constitutive_relation__Weibull__nonlinear_solver}
        config.getConfigSubtree("nonlinear_solver");
    auto const nonlinear_solver_parameters =
        createNewtonRaphsonSolverParameters(nonlinear_solver_config);

    auto const compute_local_damage =
        config.getConfigParameter<bool>("compute_local_damage");

    return std::unique_ptr<MechanicsBase<DisplacementDim>>{
        new SolidWeibull<DisplacementDim>{nonlinear_solver_parameters, mp,
                                         std::move(weibull_damage_properties),
                                         compute_local_damage}};
}

}  // namespace Weibull
}  // namespace Solids
}  // namespace MaterialLib
