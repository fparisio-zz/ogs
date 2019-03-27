/**
 * \file
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "BaseLib/Functional.h"
#include "MathLib/LinAlg/Eigen/EigenMapTools.h"

namespace ProcessLib::Utils
{
template <typename LocalAssemblerInterface, typename InternalVariable>
std::function<const std::vector<double>&(const LocalAssemblerInterface&,
                                         double,
                                         const MathLib::EigenVector&,
                                         const NumLib::LocalToGlobalIndexMap&,
                                         std::vector<double>&)>
registerInternalVariable(InternalVariable const& internal_variable)
{
    auto const& fct = internal_variable.getter;
    auto const num_components = internal_variable.num_components;

    return BaseLib::easyBind(
        [fct, num_components](
            LocalAssemblerInterface const& loc_asm,
            const double /*t*/,
            GlobalVector const& /*current_solution*/,
            NumLib::LocalToGlobalIndexMap const& /*dof_table*/,
            std::vector<double>& cache) -> std::vector<double> const& {
            const unsigned num_int_pts = loc_asm.getNumberOfIntegrationPoints();

            cache.clear();
            auto cache_mat =
                MathLib::createZeroedMatrix<Eigen::Matrix<double,
                                                          Eigen::Dynamic,
                                                          Eigen::Dynamic,
                                                          Eigen::RowMajor>>(
                    cache, num_components, num_int_pts);

            // TODO avoid the heap allocation (one per finite element)
            std::vector<double> cache_column(num_int_pts);

            for (unsigned i = 0; i < num_int_pts; ++i)
            {
                auto const& state = loc_asm.getMaterialStateVariablesAt(i);

                auto const& int_pt_values = fct(state, cache_column);
                assert(int_pt_values.size() == num_components);
                auto const int_pt_values_vec = MathLib::toVector(int_pt_values);

                cache_mat.col(i).noalias() = int_pt_values_vec;
            }

            return cache;
        });
}

template <typename LocalAssemblerInterface, typename InternalVariable,
          typename RegisterFunction>
void registerInternalVariables(
    std::vector<InternalVariable> const& internal_variables,
    RegisterFunction register_function)
{
    // Register the internal variables.
    for (auto const& internal_variable : internal_variables)
    {
        auto get_integration_point_values =
            registerInternalVariable<LocalAssemblerInterface>(
                internal_variable);

        DBUG("Registering internal variable %s.",
             internal_variable.name.c_str());
        register_function(internal_variable.name,
                          internal_variable.num_components,
                          std::move(get_integration_point_values));
    }
}
}  // namespace ProcessLib::Utils
