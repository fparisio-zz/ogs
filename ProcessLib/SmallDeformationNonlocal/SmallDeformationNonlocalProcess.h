/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <cassert>

#include "BaseLib/Functional.h"
#include "NumLib/DOF/DOFTableUtil.h"
#include "ProcessLib/Process.h"
#include "ProcessLib/SmallDeformationNonlocal/CreateLocalAssemblers.h"

#include "SmallDeformationNonlocalFEM.h"
#include "SmallDeformationNonlocalProcessData.h"

namespace ProcessLib
{
namespace SmallDeformationNonlocal
{
template <int DisplacementDim>
class SmallDeformationNonlocalProcess final : public Process
{
public:
    SmallDeformationNonlocalProcess(
        MeshLib::Mesh& mesh,
        std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&&
            jacobian_assembler,
        std::vector<std::unique_ptr<ParameterBase>> const& parameters,
        unsigned const integration_order,
        std::vector<std::vector<std::reference_wrapper<ProcessVariable>>>&&
            process_variables,
        SmallDeformationNonlocalProcessData<DisplacementDim>&& process_data,
        SecondaryVariableCollection&& secondary_variables,
        NumLib::NamedFunctionCaller&& named_function_caller);

    //! \name ODESystem interface
    //! @{

    bool isLinear() const override;
    //! @}

private:
    void initializeConcreteProcess(
        NumLib::LocalToGlobalIndexMap const& dof_table,
        MeshLib::Mesh const& mesh,
        unsigned const integration_order) override;

    void assembleConcreteProcess(const double t, GlobalVector const& x,
                                 GlobalMatrix& M, GlobalMatrix& K,
                                 GlobalVector& b) override;

    void preAssembleConcreteProcess(const double t,
                                    GlobalVector const& x) override;

    void assembleWithJacobianConcreteProcess(
        const double t, GlobalVector const& x, GlobalVector const& xdot,
        const double dxdot_dx, const double dx_dx, GlobalMatrix& M,
        GlobalMatrix& K, GlobalVector& b, GlobalMatrix& Jac) override;

    void preTimestepConcreteProcess(GlobalVector const& x, double const t,
                                    double const dt,
                                    int const /*process_id*/) override;

    void postTimestepConcreteProcess(GlobalVector const& x,
                                     int const /*process_id*/) override;

private:
    SmallDeformationNonlocalProcessData<DisplacementDim> _process_data;

    using LocalAssemblerInterface =
        SmallDeformationNonlocalLocalAssemblerInterface<DisplacementDim>;
    std::vector<std::unique_ptr<LocalAssemblerInterface>> _local_assemblers;

    std::unique_ptr<NumLib::LocalToGlobalIndexMap>
        _local_to_global_index_map_single_component;

    MeshLib::PropertyVector<double>* _nodal_forces = nullptr;
    MeshLib::PropertyVector<double>* _material_forces_property = nullptr;
    std::unique_ptr<GlobalVector> _material_forces;
};

extern template class ProcessLib::SmallDeformationNonlocal::
    SmallDeformationNonlocalProcess<2>;
extern template class ProcessLib::SmallDeformationNonlocal::
    SmallDeformationNonlocalProcess<3>;

}  // namespace SmallDeformationNonlocal
}  // namespace ProcessLib
