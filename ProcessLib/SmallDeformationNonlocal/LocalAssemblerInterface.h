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
#include <vector>

#include "NumLib/Extrapolation/ExtrapolatableElement.h"
#include "ProcessLib/IntegrationPointSerialization.h"
#include "ProcessLib/LocalAssemblerInterface.h"

namespace ProcessLib
{
namespace SmallDeformationNonlocal
{
struct SmallDeformationNonlocalLocalAssemblerInterface
    : public ProcessLib::LocalAssemblerInterface,
      public ProcessLib::SmallDeformation::NodalForceCalculationInterface,
      public NumLib::ExtrapolatableElement,
      public ProcessLib::IntegrationPointSerialization
{
    virtual std::vector<double> const& getIntPtEpsPV(
            const double /*t*/,
            GlobalVector const& /*current_solution*/,
            NumLib::LocalToGlobalIndexMap const& /*dof_table*/,
        std::vector<double>& cache) const = 0;
    virtual std::vector<double> const& getIntPtEpsPDXX(
            const double /*t*/,
            GlobalVector const& /*current_solution*/,
            NumLib::LocalToGlobalIndexMap const& /*dof_table*/,
        std::vector<double>& cache) const = 0;

    virtual std::vector<double> const& getIntPtFreeEnergyDensity(
            const double /*t*/,
            GlobalVector const& /*current_solution*/,
            NumLib::LocalToGlobalIndexMap const& /*dof_table*/,
        std::vector<double>& cache) const = 0;

    virtual std::vector<double> const& getIntPtDamage(
            const double /*t*/,
            GlobalVector const& /*current_solution*/,
            NumLib::LocalToGlobalIndexMap const& /*dof_table*/,
        std::vector<double>& cache) const = 0;

    virtual std::vector<double> const& getIntPtSigmaXX(
            const double /*t*/,
            GlobalVector const& /*current_solution*/,
            NumLib::LocalToGlobalIndexMap const& /*dof_table*/,
        std::vector<double>& cache) const = 0;

    virtual std::vector<double> const& getIntPtSigmaYY(
            const double /*t*/,
            GlobalVector const& /*current_solution*/,
            NumLib::LocalToGlobalIndexMap const& /*dof_table*/,
        std::vector<double>& cache) const = 0;

    virtual std::vector<double> const& getIntPtSigmaZZ(
            const double /*t*/,
            GlobalVector const& /*current_solution*/,
            NumLib::LocalToGlobalIndexMap const& /*dof_table*/,
        std::vector<double>& cache) const = 0;

    virtual std::vector<double> const& getIntPtSigmaXY(
            const double /*t*/,
            GlobalVector const& /*current_solution*/,
            NumLib::LocalToGlobalIndexMap const& /*dof_table*/,
        std::vector<double>& cache) const = 0;

    virtual std::vector<double> const& getIntPtSigmaXZ(
            const double /*t*/,
            GlobalVector const& /*current_solution*/,
            NumLib::LocalToGlobalIndexMap const& /*dof_table*/,
        std::vector<double>& cache) const = 0;

    virtual std::vector<double> const& getIntPtSigmaYZ(
            const double /*t*/,
            GlobalVector const& /*current_solution*/,
            NumLib::LocalToGlobalIndexMap const& /*dof_table*/,
        std::vector<double>& cache) const = 0;

    virtual std::vector<double> const& getIntPtEpsilonXX(
            const double /*t*/,
            GlobalVector const& /*current_solution*/,
            NumLib::LocalToGlobalIndexMap const& /*dof_table*/,
        std::vector<double>& cache) const = 0;

    virtual std::vector<double> const& getIntPtEpsilonYY(
            const double /*t*/,
            GlobalVector const& /*current_solution*/,
            NumLib::LocalToGlobalIndexMap const& /*dof_table*/,
        std::vector<double>& cache) const = 0;

    virtual std::vector<double> const& getIntPtEpsilonZZ(
            const double /*t*/,
            GlobalVector const& /*current_solution*/,
            NumLib::LocalToGlobalIndexMap const& /*dof_table*/,
        std::vector<double>& cache) const = 0;

    virtual std::vector<double> const& getIntPtEpsilonXY(
            const double /*t*/,
            GlobalVector const& /*current_solution*/,
            NumLib::LocalToGlobalIndexMap const& /*dof_table*/,
        std::vector<double>& cache) const = 0;

    virtual std::vector<double> const& getIntPtEpsilonXZ(
            const double /*t*/,
            GlobalVector const& /*current_solution*/,
            NumLib::LocalToGlobalIndexMap const& /*dof_table*/,
        std::vector<double>& cache) const = 0;

    virtual std::vector<double> const& getIntPtEpsilonYZ(
            const double /*t*/,
            GlobalVector const& /*current_solution*/,
            NumLib::LocalToGlobalIndexMap const& /*dof_table*/,
        std::vector<double>& cache) const = 0;

    virtual std::vector<double> const& getNodalValues(
        std::vector<double>& nodal_values) const = 0;

    virtual void nonlocal(
        std::size_t const mesh_item_id,
        std::vector<std::unique_ptr<
            SmallDeformationNonlocalLocalAssemblerInterface>> const&
            local_assemblers) = 0;

    virtual std::vector<std::tuple<int, int, Eigen::Vector3d, double>>
    getIntegrationPointCoordinates(Eigen::Vector3d const& coords) const = 0;
};

}  // namespace SmallDeformationNonlocal
}  // namespace ProcessLib
