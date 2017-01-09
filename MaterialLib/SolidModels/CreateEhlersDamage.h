/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef MATERIALLIB_SOLIDMODELS_CREATEEHLERSDAMAGE_H_
#define MATERIALLIB_SOLIDMODELS_CREATEEHLERSDAMAGE_H_

#include <logog/include/logog.hpp>

#include "MechanicsBase.h"
#include "ProcessLib/Utils/ProcessUtils.h"  // required for findParameter
#include "EhlersDamage.h"

namespace MaterialLib
{
namespace Solids
{
namespace EhlersDamage
{
template <int DisplacementDim>
std::unique_ptr<MechanicsBase<DisplacementDim>> createEhlersDamage(
    std::vector<std::unique_ptr<ProcessLib::ParameterBase>> const& parameters,
    BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{material__solid__constitutive_relation__type}
    config.checkConfigParameter("type", "EhlersDamage");
    DBUG("Create EhlersDamage material");

    //! \ogs_file_param_special{material__solid__constitutive_relation__EhlersDamage__shear_modulus}
    auto& shear_modulus = ProcessLib::findParameter<double>(
        config, "shear_modulus", parameters, 1);

    DBUG("Use \'%s\' as shear modulus parameter.", shear_modulus.name.c_str());

    //! \ogs_file_param_special{material__solid__constitutive_relation__EhlersDamage__bulk_modulus}
    auto& bulk_modulus = ProcessLib::findParameter<double>(
        config, "bulk_modulus", parameters, 1);

    DBUG("Use \'%s\' as bulk modulus parameter.", bulk_modulus.name.c_str());

    //! \ogs_file_param_special{material__solid__constitutive_relation__EhlersDamage__kappa}
    auto& kappa = ProcessLib::findParameter<double>(
        config, "kappa", parameters, 1);

    DBUG("Use \'%s\' as kappa.",
         kappa.name.c_str());

    //! \ogs_file_param_special{material__solid__constitutive_relation__EhlersDamage__beta}
    auto& beta = ProcessLib::findParameter<double>(
        config, "beta", parameters, 1);

    DBUG("Use \'%s\' as beta.",
         beta.name.c_str());

    //! \ogs_file_param_special{material__solid__constitutive_relation__EhlersDamage__gamma}
    auto& gamma = ProcessLib::findParameter<double>(
        config, "gamma", parameters, 1);

    DBUG("Use \'%s\' as gamma.",
         gamma.name.c_str());

    //! \ogs_file_param_special{material__solid__constitutive_relation__EhlersDamage__hardening_modulus}
    auto& hardening_modulus = ProcessLib::findParameter<double>(
        config, "hardening_modulus", parameters, 1);

    DBUG("Use \'%s\' as hardening modulus parameter.",
         hardening_modulus.name.c_str());

    //! \ogs_file_param_special{material__solid__constitutive_relation__EhlersDamage__alpha}
    auto& alpha = ProcessLib::findParameter<double>(
        config, "alpha", parameters, 1);

    DBUG("Use \'%s\' as alpha.",
         alpha.name.c_str());

    //! \ogs_file_param_special{material__solid__constitutive_relation__EhlersDamage__delta}
    auto& delta = ProcessLib::findParameter<double>(
        config, "delta", parameters, 1);

    DBUG("Use \'%s\' as delta.",
         delta.name.c_str());

    //! \ogs_file_param_special{material__solid__constitutive_relation__Ehlers__eps}
    auto& eps = ProcessLib::findParameter<double>(
        config, "eps", parameters, 1);

    DBUG("Use \'%s\' as eps.",
         eps.name.c_str());

    //! \ogs_file_param_special{material__solid__constitutive_relation__Ehlers__m}
    auto& m = ProcessLib::findParameter<double>(
        config, "m", parameters, 1);

    DBUG("Use \'%s\' as m.",
         m.name.c_str());

    //! \ogs_file_param_special{material__solid__constitutive_relation__Ehlers__alphap}
    auto& alphap = ProcessLib::findParameter<double>(
        config, "alphap", parameters, 1);

    DBUG("Use \'%s\' as alphap.",
         alphap.name.c_str());

    //! \ogs_file_param_special{material__solid__constitutive_relation__Ehlers__deltap}
    auto& deltap = ProcessLib::findParameter<double>(
        config, "deltap", parameters, 1);

    DBUG("Use \'%s\' as deltap.",
         deltap.name.c_str());

    //! \ogs_file_param_special{material__solid__constitutive_relation__Ehlers__epsp}
    auto& epsp = ProcessLib::findParameter<double>(
        config, "epsp", parameters, 1);

    DBUG("Use \'%s\' as epsp.",
         epsp.name.c_str());

    //! \ogs_file_param_special{material__solid__constitutive_relation__Ehlers__mp}
    auto& paremeter_mp = ProcessLib::findParameter<double>(
        config, "mp", parameters, 1);

    DBUG("Use \'%s\' as mp.",
         paremeter_mp.name.c_str());

    //! \ogs_file_param_special{material__solid__constitutive_relation__EhlersDamage__betap}
    auto& betap = ProcessLib::findParameter<double>(
        config, "betap", parameters, 1);

    DBUG("Use \'%s\' as betap.",
         betap.name.c_str());

    //! \ogs_file_param_special{material__solid__constitutive_relation__EhlersDamage__gammap}
    auto& gammap = ProcessLib::findParameter<double>(
        config, "gammap", parameters, 1);

    DBUG("Use \'%s\' as gammap.",
         gammap.name.c_str());

    //! \ogs_file_param_special{material__solid__constitutive_relation__EhlersDamage__alpha_d}
    auto& alpha_d = ProcessLib::findParameter<double>(
        config, "alpha_d", parameters, 1);

    DBUG("Use \'%s\' as alpha_d parameter.", alpha_d.name.c_str());

    //! \ogs_file_param_special{material__solid__constitutive_relation__EhlersDamage__beta_d}
    auto& beta_d = ProcessLib::findParameter<double>(
    config, "beta_d", parameters, 1);

    DBUG("Use \'%s\' as beta_d parameter.", beta_d.name.c_str());

    typename SolidEhlersDamage<DisplacementDim>::MaterialProperties mp{
        shear_modulus, bulk_modulus, alpha,  beta,
        gamma,         delta,        eps,    m,
        alphap,        betap,        gammap, deltap,
        epsp,          paremeter_mp, kappa,  hardening_modulus, alpha_d, beta_d};

    return std::unique_ptr<MechanicsBase<DisplacementDim>>{
        new SolidEhlersDamage<DisplacementDim>{mp}};
}

}  // namespace EhlersDamage
}  // namespace Solids
}  // namespace MaterialLib

#endif  // MATERIALLIB_SOLIDMODELS_CREATEEHLERS_H_
