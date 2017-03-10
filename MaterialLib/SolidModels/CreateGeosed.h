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

#include "Geosed.h"
#include "MechanicsBase.h"

namespace MaterialLib
{
namespace Solids
{
namespace Geosed
{
inline std::unique_ptr<GeosedDamageProperties> createDamageProperties(
    std::vector<std::unique_ptr<ProcessLib::ParameterBase>> const& parameters,
    BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param_special{material__solid__constitutive_relation__Geosed__damage_properties__alpha_d}
    auto& alpha_d =
        ProcessLib::findParameter<double>(config, "alpha_d", parameters, 1);

    DBUG("Use \'%s\' as alpha_d.", alpha_d.name.c_str());

    //! \ogs_file_param_special{material__solid__constitutive_relation__Geosed__damage_properties__beta_d}
    auto& beta_d =
        ProcessLib::findParameter<double>(config, "beta_d", parameters, 1);

    DBUG("Use \'%s\' as beta_d.", beta_d.name.c_str());

    //! \ogs_file_param_special{material__solid__constitutive_relation__Geosed__damage_properties__h_d}
    auto& h_d = ProcessLib::findParameter<double>(config, "h_d", parameters, 1);

    DBUG("Use \'%s\' as h_d.", h_d.name.c_str());

    return std::unique_ptr<GeosedDamageProperties>{
        new GeosedDamageProperties{alpha_d, beta_d, h_d}};
}

template <int DisplacementDim>
std::unique_ptr<MechanicsBase<DisplacementDim>> createGeosed(
    std::vector<std::unique_ptr<ProcessLib::ParameterBase>> const& parameters,
    BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{material__solid__constitutive_relation__type}
    config.checkConfigParameter("type", "Geosed");
    DBUG("Create Geosed material");

    //! \ogs_file_param_special{material__solid__constitutive_relation__Geosed__shear_modulus}
    auto& shear_modulus = ProcessLib::findParameter<double>(
        config, "shear_modulus", parameters, 1);

    DBUG("Use \'%s\' as shear modulus parameter.", shear_modulus.name.c_str());

    //! \ogs_file_param_special{material__solid__constitutive_relation__Geosed__bulk_modulus}
    auto& bulk_modulus = ProcessLib::findParameter<double>(
        config, "bulk_modulus", parameters, 1);

    DBUG("Use \'%s\' as bulk modulus parameter.", bulk_modulus.name.c_str());

    //! \ogs_file_param_special{material__solid__constitutive_relation__Geosed__kappa}
    auto& kappa =
        ProcessLib::findParameter<double>(config, "kappa", parameters, 1);

    DBUG("Use \'%s\' as kappa.", kappa.name.c_str());

    //! \ogs_file_param_special{material__solid__constitutive_relation__Geosed__beta}
    auto& beta =
        ProcessLib::findParameter<double>(config, "beta", parameters, 1);

    DBUG("Use \'%s\' as beta.", beta.name.c_str());

    //! \ogs_file_param_special{material__solid__constitutive_relation__Geosed__gamma}
    auto& gamma =
        ProcessLib::findParameter<double>(config, "gamma", parameters, 1);

    DBUG("Use \'%s\' as gamma.", gamma.name.c_str());

    //! \ogs_file_param_special{material__solid__constitutive_relation__Geosed__hardening_modulus}
    auto& hardening_modulus = ProcessLib::findParameter<double>(
        config, "hardening_modulus", parameters, 1);

    DBUG("Use \'%s\' as hardening modulus parameter.",
         hardening_modulus.name.c_str());

    //! \ogs_file_param_special{material__solid__constitutive_relation__Geosed__alpha}
    auto& alpha =
        ProcessLib::findParameter<double>(config, "alpha", parameters, 1);

    DBUG("Use \'%s\' as alpha.", alpha.name.c_str());

    //! \ogs_file_param_special{material__solid__constitutive_relation__Geosed__delta}
    auto& delta =
        ProcessLib::findParameter<double>(config, "delta", parameters, 1);

    DBUG("Use \'%s\' as delta.", delta.name.c_str());

    //! \ogs_file_param_special{material__solid__constitutive_relation__Geosed__eps}
    auto& eps = ProcessLib::findParameter<double>(config, "eps", parameters, 1);

    DBUG("Use \'%s\' as eps.", eps.name.c_str());

    //! \ogs_file_param_special{material__solid__constitutive_relation__Geosed__m}
    auto& m = ProcessLib::findParameter<double>(config, "m", parameters, 1);

    DBUG("Use \'%s\' as m.", m.name.c_str());

    //! \ogs_file_param_special{material__solid__constitutive_relation__Geosed__alphap}
    auto& alphap =
        ProcessLib::findParameter<double>(config, "alphap", parameters, 1);

    DBUG("Use \'%s\' as alphap.", alphap.name.c_str());

    //! \ogs_file_param_special{material__solid__constitutive_relation__Geosed__deltap}
    auto& deltap =
        ProcessLib::findParameter<double>(config, "deltap", parameters, 1);

    DBUG("Use \'%s\' as deltap.", deltap.name.c_str());

    //! \ogs_file_param_special{material__solid__constitutive_relation__Geosed__epsp}
    auto& epsp =
        ProcessLib::findParameter<double>(config, "epsp", parameters, 1);

    DBUG("Use \'%s\' as epsp.", epsp.name.c_str());

    //! \ogs_file_param_special{material__solid__constitutive_relation__Geosed__mp}
    auto& paremeter_mp =
        ProcessLib::findParameter<double>(config, "mp", parameters, 1);

    DBUG("Use \'%s\' as mp.", paremeter_mp.name.c_str());

    //! \ogs_file_param_special{material__solid__constitutive_relation__Geosed__betap}
    auto& betap =
        ProcessLib::findParameter<double>(config, "betap", parameters, 1);

    DBUG("Use \'%s\' as betap.", betap.name.c_str());

    //! \ogs_file_param_special{material__solid__constitutive_relation__Geosed__gammap}
    auto& gammap =
        ProcessLib::findParameter<double>(config, "gammap", parameters, 1);

    DBUG("Use \'%s\' as gammap.", gammap.name.c_str());

    typename SolidGeosed<DisplacementDim>::MaterialProperties mp{
        shear_modulus, bulk_modulus, alpha,  beta,
        gamma,         delta,        eps,    m,
        alphap,        betap,        gammap, deltap,
        epsp,          paremeter_mp, kappa,  hardening_modulus};

    // Damage properties.
    std::unique_ptr<GeosedDamageProperties> geosed_damage_properties;

    auto const& geosed_damage_config =
        //! \ogs_file_param{material__solid__constitutive_relation__Geosed__damage_properties}
        config.getConfigSubtreeOptional("damage_properties");
    if (geosed_damage_config)
    {
        geosed_damage_properties =
            createDamageProperties(parameters, *geosed_damage_config);
    }

    return std::unique_ptr<MechanicsBase<DisplacementDim>>{
        new SolidGeosed<DisplacementDim>{mp,
                                         std::move(geosed_damage_properties)}};
}

}  // namespace Geosed
}  // namespace Solids
}  // namespace MaterialLib
