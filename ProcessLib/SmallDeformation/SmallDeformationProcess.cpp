/**
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "SmallDeformationProcess.h"

#include <cassert>
#include <nlohmann/json.hpp>

#include "ProcessLib/Process.h"
#include "ProcessLib/SmallDeformation/CreateLocalAssemblers.h"
#include "ProcessLib/Utils/RegisterInternalVariables.h"

#include "SmallDeformationFEM.h"

namespace ProcessLib
{
namespace SmallDeformation
{
template <int DisplacementDim>
SmallDeformationProcess<DisplacementDim>::SmallDeformationProcess(
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    unsigned const integration_order,
    std::vector<std::vector<std::reference_wrapper<ProcessVariable>>>&&
        process_variables,
    SmallDeformationProcessData<DisplacementDim>&& process_data,
    SecondaryVariableCollection&& secondary_variables,
    NumLib::NamedFunctionCaller&& named_function_caller)
    : Process(mesh, std::move(jacobian_assembler), parameters,
              integration_order, std::move(process_variables),
              std::move(secondary_variables), std::move(named_function_caller)),
      _process_data(std::move(process_data))
{
    _nodal_forces = MeshLib::getOrCreateMeshProperty<double>(
        mesh, "NodalForces", MeshLib::MeshItemType::Node, DisplacementDim);

    _material_forces = MeshLib::getOrCreateMeshProperty<double>(
        mesh, "MaterialForces", MeshLib::MeshItemType::Node, DisplacementDim);

    _integration_point_writer.emplace_back(
        std::make_unique<SigmaIntegrationPointWriter>(
            static_cast<int>(mesh.getDimension() == 2 ? 4 : 6) /*n components*/,
            2 /*integration order*/, [this]() {
                // Result containing integration point data for each local
                // assembler.
                std::vector<std::vector<double>> result;
                result.resize(_local_assemblers.size());

                for (std::size_t i = 0; i < _local_assemblers.size(); ++i)
                {
                    auto const& local_asm = *_local_assemblers[i];

                    result[i] = local_asm.getSigma();
                }

                return result;
            }));
}

template <int DisplacementDim>
bool SmallDeformationProcess<DisplacementDim>::isLinear() const
{
    return false;
}

template <int DisplacementDim>
void SmallDeformationProcess<DisplacementDim>::initializeConcreteProcess(
    NumLib::LocalToGlobalIndexMap const& dof_table,
    MeshLib::Mesh const& mesh,
    unsigned const integration_order)
{
    using nlohmann::json;

    ProcessLib::SmallDeformation::createLocalAssemblers<
        DisplacementDim, SmallDeformationLocalAssembler>(
        mesh.getElements(), dof_table, _local_assemblers,
        mesh.isAxiallySymmetric(), integration_order, _process_data);

    // TODO move the two data members somewhere else.
    // for extrapolation of secondary variables
    std::vector<MeshLib::MeshSubset> all_mesh_subsets_single_component{
        *_mesh_subset_all_nodes};
    _local_to_global_index_map_single_component =
        std::make_unique<NumLib::LocalToGlobalIndexMap>(
            std::move(all_mesh_subsets_single_component),
            // by location order is needed for output
            NumLib::ComponentOrder::BY_LOCATION);

    auto add_secondary_variable = [this](std::string const& name,
                                         int const num_components,
                                         auto function) {
        _secondary_variables.addSecondaryVariable(
            name,
            makeExtrapolator(num_components, getExtrapolator(),
                             _local_assemblers, function));
    };

    add_secondary_variable("free_energy_density", 1,
                           &LocalAssemblerInterface::getIntPtFreeEnergyDensity);

    add_secondary_variable("sigma",
                           MathLib::KelvinVector::KelvinVectorType<
                               DisplacementDim>::RowsAtCompileTime,
                           &LocalAssemblerInterface::getIntPtSigma);

    add_secondary_variable("epsilon",
                           MathLib::KelvinVector::KelvinVectorType<
                               DisplacementDim>::RowsAtCompileTime,
                           &LocalAssemblerInterface::getIntPtEpsilon);

    //
    // enable output of internal variables defined by material models
    //

    // Collect the internal variables for all solid materials.
    std::vector<typename MaterialLib::Solids::MechanicsBase<
        DisplacementDim>::InternalVariable>
        internal_variables;
    for (auto const& material_id__solid_material :
         _process_data.solid_materials)
    {
        auto const variables =
            material_id__solid_material.second->getInternalVariables();
        copy(begin(variables), end(variables),
             back_inserter(internal_variables));
    }
    Utils::registerInternalVariables<LocalAssemblerInterface>(
        internal_variables, add_secondary_variable);

    // Set initial conditions for integration point data.
    for (auto const& ip_writer : _integration_point_writer)
    {
        // Find the mesh property with integration point writer's name.
        auto const& name = ip_writer->name();
        if (!mesh.getProperties().existsPropertyVector<double>(name))
        {
            continue;
        }
        auto const& mesh_property =
            *mesh.getProperties().template getPropertyVector<double>(name);

        // The mesh property must be defined on integration points.
        if (mesh_property.getMeshItemType() !=
            MeshLib::MeshItemType::IntegrationPoint)
        {
            continue;
        }

        auto const ip_meta_data = getIntegrationPointMetaData(mesh, name);

        // Check the number of components.
        if (ip_meta_data.n_components != mesh_property.getNumberOfComponents())
        {
            OGS_FATAL(
                "Different number of components in meta data (%d) than in "
                "the integration point field data for '%s': %d.",
                ip_meta_data.n_components, name.c_str(),
                mesh_property.getNumberOfComponents());
        }

        // Now we have a properly named vtk's field data array and the
        // corresponding meta data.
        std::size_t position = 0;
        for (auto& local_asm : _local_assemblers)
        {
            std::size_t const integration_points_read =
                local_asm->setIPDataInitialConditions(
                    name, &mesh_property[position],
                    ip_meta_data.integration_order);
            if (integration_points_read == 0)
            {
                OGS_FATAL(
                    "No integration points read in the integration point "
                    "initial conditions set function.");
            }
            position += integration_points_read * ip_meta_data.n_components;
        }
    }
}

template <int DisplacementDim>
void SmallDeformationProcess<DisplacementDim>::assembleConcreteProcess(
    const double t, GlobalVector const& x, GlobalMatrix& M, GlobalMatrix& K,
    GlobalVector& b)
{
    DBUG("Assemble SmallDeformationProcess.");

    const int process_id = 0;
    ProcessLib::ProcessVariable const& pv = getProcessVariables(process_id)[0];

    std::vector<std::reference_wrapper<NumLib::LocalToGlobalIndexMap>>
        dof_table = {std::ref(*_local_to_global_index_map)};
    // Call global assembler for each local assembly item.
    GlobalExecutor::executeSelectedMemberDereferenced(
        _global_assembler, &VectorMatrixAssembler::assemble, _local_assemblers,
        pv.getActiveElementIDs(), dof_table, t, x, M, K, b,
        _coupled_solutions);
}

template <int DisplacementDim>
void SmallDeformationProcess<DisplacementDim>::
    assembleWithJacobianConcreteProcess(const double t, GlobalVector const& x,
                                        GlobalVector const& xdot,
                                        const double dxdot_dx,
                                        const double dx_dx, GlobalMatrix& M,
                                        GlobalMatrix& K, GlobalVector& b,
                                        GlobalMatrix& Jac)
{
    DBUG("AssembleWithJacobian SmallDeformationProcess.");

    const int process_id = 0;
    ProcessLib::ProcessVariable const& pv = getProcessVariables(process_id)[0];

    std::vector<std::reference_wrapper<NumLib::LocalToGlobalIndexMap>>
        dof_table = {std::ref(*_local_to_global_index_map)};
    // Call global assembler for each local assembly item.
    GlobalExecutor::executeSelectedMemberDereferenced(
        _global_assembler, &VectorMatrixAssembler::assembleWithJacobian,
        _local_assemblers, pv.getActiveElementIDs(), dof_table, t, x,
        xdot, dxdot_dx, dx_dx, M, K, b, Jac, _coupled_solutions);

    transformVariableFromGlobalVector(b, 0, *_local_to_global_index_map,
                                      *_nodal_forces, std::negate<double>());
}

template <int DisplacementDim>
void SmallDeformationProcess<DisplacementDim>::preTimestepConcreteProcess(
    GlobalVector const& x, double const t, double const dt,
    const int process_id)
{
    DBUG("PreTimestep SmallDeformationProcess.");

    _process_data.dt = dt;
    _process_data.t = t;

    ProcessLib::ProcessVariable const& pv = getProcessVariables(process_id)[0];

    GlobalExecutor::executeSelectedMemberOnDereferenced(
        &LocalAssemblerInterface::preTimestep, _local_assemblers,
        pv.getActiveElementIDs(), *_local_to_global_index_map,
        x, t, dt);
}

template <int DisplacementDim>
void SmallDeformationProcess<DisplacementDim>::postTimestepConcreteProcess(
    GlobalVector const& x, const double /*t*/, const double /*delta_t*/,
    int const process_id)
{
    DBUG("PostTimestep SmallDeformationProcess.");

    ProcessLib::ProcessVariable const& pv = getProcessVariables(process_id)[0];

    GlobalExecutor::executeSelectedMemberOnDereferenced(
        &LocalAssemblerInterface::postTimestep, _local_assemblers,
        pv.getActiveElementIDs(), *_local_to_global_index_map, x);

    std::unique_ptr<GlobalVector> material_forces;
    ProcessLib::SmallDeformation::writeMaterialForces(
        material_forces, _local_assemblers, *_local_to_global_index_map, x);

    material_forces->copyValues(*_material_forces);
}

template class SmallDeformationProcess<2>;
template class SmallDeformationProcess<3>;

}  // namespace SmallDeformation
}  // namespace ProcessLib
