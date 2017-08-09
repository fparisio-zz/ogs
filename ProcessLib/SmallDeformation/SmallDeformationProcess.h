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
#include "ProcessLib/Process.h"
#include "SmallDeformationProcessData.h"
#include "CreateLocalAssemblers.h"
#include "ProcessLib/SmallDeformationCommon/Common.h"

namespace ProcessLib
{
namespace SmallDeformation
{
template <int DisplacementDim>
class SmallDeformationProcess final : public Process
{
    using Base = Process;

public:
    SmallDeformationProcess(
        MeshLib::Mesh& mesh,
        std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&&
            jacobian_assembler,
        std::vector<std::unique_ptr<ParameterBase>> const& parameters,
        unsigned const integration_order,
        std::vector<std::reference_wrapper<ProcessVariable>>&&
            process_variables,
        SmallDeformationProcessData<DisplacementDim>&& process_data,
        SecondaryVariableCollection&& secondary_variables,
        NumLib::NamedFunctionCaller&& named_function_caller);

    //! \name ODESystem interface
    //! @{
    bool isLinear() const override;
    //! @}

private:
    using LocalAssemblerInterface =
        SmallDeformationLocalAssemblerInterface<DisplacementDim>;

    void initializeConcreteProcess(
        NumLib::LocalToGlobalIndexMap const& dof_table,
        MeshLib::Mesh const& mesh,
        unsigned const integration_order) override;

    std::size_t writeIntegrationPointData(MeshLib::PropertyVector<char>& output,
            MeshLib::PropertyVector<std::size_t>& offsets)
    {
        output.clear();
        offsets.clear();
        std::vector<char> local_data;
        std::size_t offset = 0;
        for (auto& la : _local_assemblers)
        {
            offsets.push_back(offset);
            std::size_t const local_offset =
                la->writeIntegrationPointData(local_data);
            std::copy_n(std::begin(local_data), local_offset,
                        std::back_inserter(output));
            offset += local_offset;
        }
        return offset;
    }

    void assembleConcreteProcess(const double t, GlobalVector const& x,
                                         GlobalMatrix& M, GlobalMatrix& K,
                                         GlobalVector& b,
                                         StaggeredCouplingTerm const&
                                         coupling_term) override;

    void assembleWithJacobianConcreteProcess(
        const double t, GlobalVector const& x, GlobalVector const& xdot,
        const double dxdot_dx, const double dx_dx, GlobalMatrix& M,
        GlobalMatrix& K, GlobalVector& b, GlobalMatrix& Jac,
        StaggeredCouplingTerm const& coupling_term) override;

    void setInitialConditionsConcreteProcess(double const t,
                                             GlobalVector const& x) override
    {
        DBUG("SetInitialConditions SmallDeformationProcess.");

        if (!_mesh.getProperties().hasPropertyVector("integration_point_data"))
            return;
        if (!_mesh.getProperties().hasPropertyVector("integration_point_offsets"))
            OGS_FATAL(
                "integration_point_data field exists in the input but there is "
                "no integration_point_offsets cell data.");

        auto const& data =
            *_mesh.getProperties().template getPropertyVector<char>(
                "integration_point_data");
        assert(data.getMeshItemType() ==
               MeshLib::MeshItemType::IntegrationPoint);

        auto const& offsets =
            *_mesh.getProperties().template getPropertyVector<std::size_t>(
                "integration_point_offsets");
        assert(offsets.getMeshItemType() == MeshLib::MeshItemType::Cell);

        std::vector<char> local_data;
        assert(_local_assemblers.size() == offsets.size());
        // Starting counting from one; the last cell is handled after the loop.
        std::size_t i = 0;
        for (; i < _local_assemblers.size() - 1; ++i)
        {
            std::size_t const size = offsets[i + 1] - offsets[i];
            local_data.resize(size);
            std::memcpy(local_data.data(), &data[offsets[i]], size);
            _local_assemblers[i]->readIntegrationPointData(local_data);
        }
        {   // last cell
            std::size_t const size = data.size() - offsets[i];
            local_data.resize(size);
            std::memcpy(local_data.data(), &data[offsets[i]], size);
            _local_assemblers[i]->readIntegrationPointData(local_data);
        }
    }

    void preTimestepConcreteProcess(GlobalVector const& x, double const t,
                                    double const dt) override;

    void postTimestepConcreteProcess(GlobalVector const& x) override
    {
        DBUG("PostTimestep SmallDeformationProcess.");

        ProcessLib::SmallDeformation::writeNodalForces(
            *_nodal_forces, _local_assemblers, *_local_to_global_index_map);

        ProcessLib::SmallDeformation::writeMaterialForces(
            *_material_forces, _local_assemblers, *_local_to_global_index_map,
            x);
    }

private:
    SmallDeformationProcessData<DisplacementDim> _process_data;

    std::vector<std::unique_ptr<LocalAssemblerInterface>> _local_assemblers;

    std::unique_ptr<NumLib::LocalToGlobalIndexMap>
        _local_to_global_index_map_single_component;
    MeshLib::PropertyVector<double>* _nodal_forces = nullptr;
    MeshLib::PropertyVector<double>* _material_forces = nullptr;
};

extern template class SmallDeformationProcess<2>;
extern template class SmallDeformationProcess<3>;

}  // namespace SmallDeformation
}  // namespace ProcessLib
