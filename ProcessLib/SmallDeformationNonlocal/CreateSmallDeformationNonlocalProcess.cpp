/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "CreateSmallDeformationNonlocalProcess.h"

#include <cassert>

#include "MaterialLib/SolidModels/CreateEhlers.h"
#include "MaterialLib/SolidModels/CreateLinearElasticIsotropic.h"
#include "MaterialLib/SolidModels/CreateLubby2.h"
#include "ProcessLib/Utils/ParseSecondaryVariables.h"

#include "SmallDeformationNonlocalProcess.h"
#include "SmallDeformationNonlocalProcessData.h"

namespace ProcessLib
{
namespace SmallDeformationNonlocal
{
template <int DisplacementDim>
class SmallDeformationNonlocalProcess;

extern template class SmallDeformationNonlocalProcess<2>;
extern template class SmallDeformationNonlocalProcess<3>;

template <int DisplacementDim>
std::unique_ptr<Process> createSmallDeformationNonlocalProcess(
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<ProcessVariable> const& variables,
    std::vector<std::unique_ptr<ParameterBase>> const& parameters,
    unsigned const integration_order,
    BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{prj__processes__process__type}
    config.checkConfigParameter("type", "SMALL_DEFORMATION_NONLOCAL");
    DBUG("Create SmallDeformationNonlocalProcess.");

    // Process variable.

    //! \ogs_file_param{prj__processes__process__SMALL_DEFORMATION__process_variables}
    auto const pv_config = config.getConfigSubtree("process_variables");

    auto process_variables = findProcessVariables(
        variables, pv_config,
        {//! \ogs_file_param_special{prj__processes__process__SMALL_DEFORMATION__process_variables__process_variable}
         "process_variable"});

    DBUG("Associate displacement with process variable \'%s\'.",
         process_variables.back().get().getName().c_str());

    if (process_variables.back().get().getNumberOfComponents() !=
        DisplacementDim)
    {
        OGS_FATAL(
            "Number of components of the process variable '%s' is different "
            "from the displacement dimension: got %d, expected %d",
            process_variables.back().get().getName().c_str(),
            process_variables.back().get().getNumberOfComponents(),
            DisplacementDim);
    }

    // Constitutive relation.
    // read type;
    auto const constitutive_relation_config =
        //! \ogs_file_param{prj__processes__process__SMALL_DEFORMATION__constitutive_relation}
        config.getConfigSubtree("constitutive_relation");

    auto const type =
        //! \ogs_file_param{prj__processes__process__SMALL_DEFORMATION__constitutive_relation__type}
        constitutive_relation_config.peekConfigParameter<std::string>("type");

    std::unique_ptr<MaterialLib::Solids::MechanicsBase<DisplacementDim>>
        material = nullptr;
    if (type == "Ehlers")
    {
        material = MaterialLib::Solids::Ehlers::createEhlers<DisplacementDim>(
            parameters, constitutive_relation_config);
    }
    else if (type == "LinearElasticIsotropic")
    {
        material =
            MaterialLib::Solids::createLinearElasticIsotropic<DisplacementDim>(
                parameters, constitutive_relation_config);
    }
    else if (type == "Lubby2")
    {
        material = MaterialLib::Solids::Lubby2::createLubby2<DisplacementDim>(
            parameters, constitutive_relation_config);
    }
    else
    {
        OGS_FATAL(
            "Cannot construct constitutive relation of given type \'%s\'.",
            type.c_str());
    }

    auto const internal_length =
        //! \ogs_file_param{prj__processes__process__SMALL_DEFORMATION__internal_length}
        config.getConfigParameter<double>("internal_length");

    SmallDeformationNonlocalProcessData<DisplacementDim> process_data{
        std::move(material), internal_length};

    SecondaryVariableCollection secondary_variables;

    NumLib::NamedFunctionCaller named_function_caller(
        {"SmallDeformationNonlocal_displacement"});

    ProcessLib::parseSecondaryVariables(config, secondary_variables,
                                        named_function_caller);

    return std::make_unique<SmallDeformationNonlocalProcess<DisplacementDim>>(
        mesh, std::move(jacobian_assembler), parameters, integration_order,
        std::move(process_variables), std::move(process_data),
        std::move(secondary_variables), std::move(named_function_caller));
}

template std::unique_ptr<Process> createSmallDeformationNonlocalProcess<2>(
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<ProcessVariable> const& variables,
    std::vector<std::unique_ptr<ParameterBase>> const& parameters,
    unsigned const integration_order,
    BaseLib::ConfigTree const& config);

template std::unique_ptr<Process> createSmallDeformationNonlocalProcess<3>(
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<ProcessVariable> const& variables,
    std::vector<std::unique_ptr<ParameterBase>> const& parameters,
    unsigned const integration_order,
    BaseLib::ConfigTree const& config);

}  // namespace SmallDeformationNonlocal
}  // namespace ProcessLib
