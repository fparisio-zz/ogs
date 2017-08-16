/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <iostream>
#include <limits>
#include <memory>
#include <vector>

#include "MaterialLib/SolidModels/Ehlers.h"
#include "MaterialLib/SolidModels/LinearElasticIsotropic.h"
#include "MaterialLib/SolidModels/Lubby2.h"
#include "MathLib/LinAlg/Eigen/EigenMapTools.h"
#include "NumLib/Fem/FiniteElement/TemplateIsoparametric.h"
#include "NumLib/Fem/ShapeMatrixPolicy.h"
#include "NumLib/Function/Interpolation.h"
#include "ProcessLib/Deformation/BMatrixPolicy.h"
#include "ProcessLib/Deformation/LinearBMatrix.h"
#include "ProcessLib/LocalAssemblerTraits.h"
#include "ProcessLib/Parameter/Parameter.h"
#include "ProcessLib/Utils/InitShapeMatrices.h"

#include "ProcessLib/SmallDeformationCommon/integration_point_data.h"

#include "IntegrationPointData.h"
#include "LocalAssemblerInterface.h"
#include "SmallDeformationNonlocalProcessData.h"

namespace ProcessLib
{
namespace SmallDeformationNonlocal
{
/// Used by for extrapolation of the integration point values. It is ordered
/// (and stored) by integration points.
template <typename ShapeMatrixType>
struct SecondaryData
{
    std::vector<ShapeMatrixType, Eigen::aligned_allocator<ShapeMatrixType>> N;
};

template <typename ShapeFunction, typename IntegrationMethod,
          int DisplacementDim_>
class SmallDeformationNonlocalLocalAssembler
    : public SmallDeformationNonlocalLocalAssemblerInterface<DisplacementDim_>
{
public:
    static int const DisplacementDim = DisplacementDim_;
    using ShapeMatricesType =
        ShapeMatrixPolicyType<ShapeFunction, DisplacementDim>;
    using NodalMatrixType = typename ShapeMatricesType::NodalMatrixType;
    using NodalVectorType = typename ShapeMatricesType::NodalVectorType;
    using GlobalDimVectorType = typename ShapeMatricesType::GlobalDimVectorType;
    using ShapeMatrices = typename ShapeMatricesType::ShapeMatrices;
    using BMatricesType = BMatrixPolicyType<ShapeFunction, DisplacementDim>;

    using BMatrixType = typename BMatricesType::BMatrixType;
    using GradientVectorType = typename BMatricesType::GradientVectorType;
    using GradientMatrixType = typename BMatricesType::GradientMatrixType;
    using StiffnessMatrixType = typename BMatricesType::StiffnessMatrixType;
    using NodalForceVectorType = typename BMatricesType::NodalForceVectorType;
    using NodalDisplacementVectorType =
        typename BMatricesType::NodalForceVectorType;

    SmallDeformationNonlocalLocalAssembler(
        SmallDeformationNonlocalLocalAssembler const&) = delete;
    SmallDeformationNonlocalLocalAssembler(
        SmallDeformationNonlocalLocalAssembler&&) = delete;

    SmallDeformationNonlocalLocalAssembler(
        MeshLib::Element const& e,
        std::size_t const /*local_matrix_size*/,
        bool const is_axially_symmetric,
        unsigned const integration_order,
        SmallDeformationNonlocalProcessData<DisplacementDim>& process_data)
        : _process_data(process_data),
          _integration_method(integration_order),
          _element(e),
          _is_axially_symmetric(is_axially_symmetric)
    {
        unsigned const n_integration_points =
            _integration_method.getNumberOfPoints();

        _ip_data.reserve(n_integration_points);
        _secondary_data.N.resize(n_integration_points);
        // TODO (naumov) This is compile-time size. Use eigen matrix.
        _material_forces.resize(DisplacementDim * ShapeFunction::NPOINTS);

        auto const shape_matrices =
            initShapeMatrices<ShapeFunction, ShapeMatricesType,
                              IntegrationMethod, DisplacementDim>(
                e, is_axially_symmetric, _integration_method);

        for (unsigned ip = 0; ip < n_integration_points; ip++)
        {
            _ip_data.emplace_back(*_process_data.material);
            auto& ip_data = _ip_data[ip];
            auto const& sm = shape_matrices[ip];
            _ip_data[ip].integration_weight =
                _integration_method.getWeightedPoint(ip).getWeight() *
                sm.integralMeasure * sm.detJ;

            ip_data.N = sm.N;
            ip_data.dNdx = sm.dNdx;

            // Initialize current time step values
            ip_data.sigma.setZero(
                KelvinVectorDimensions<DisplacementDim>::value);
            ip_data.eps.setZero(KelvinVectorDimensions<DisplacementDim>::value);

            // Previous time step values are not initialized and are set later.
            ip_data.sigma_prev.resize(
                KelvinVectorDimensions<DisplacementDim>::value);
            ip_data.eps_prev.resize(
                KelvinVectorDimensions<DisplacementDim>::value);

            _secondary_data.N[ip] = shape_matrices[ip].N;
        }
    }

    double alpha_0(double const distance2) const
    {
        double const internal_length2 =
            _process_data.internal_length * _process_data.internal_length;
        return (distance2 > internal_length2)
                   ? 0
                   : (1 - distance2 / (internal_length2)) *
                         (1 - distance2 / (internal_length2));
    }

    void nonlocal(
        std::size_t const /*mesh_item_id*/,
        std::vector<
            std::unique_ptr<SmallDeformationNonlocalLocalAssemblerInterface<
                DisplacementDim>>> const& local_assemblers) override
    {
        // std::cout << "\nXXX nonlocal in element " << _element.getID() <<
        // "\n";

        unsigned const n_integration_points =
            _integration_method.getNumberOfPoints();

        //
        // For every integration point in this element collect the neighbouring
        // integration points falling in given radius (internal length) and
        // compute the alpha_kl weight.
        //
        for (unsigned k = 0; k < n_integration_points; k++)
        {
            //
            // Collect the integration points.
            //
            // std::cout << "\n\tip = " << k << "\n";

            auto const xyz = getSingleIntegrationPointCoordinates(k);
            // std::cout << "\tCurrent ip_k coords : " << xyz << "\n";

            // For all neighbors of element
            for (auto const& la : local_assemblers)
            {
                auto const neighbor_ip_coords =
                    la->getIntegrationPointCoordinates(xyz);
                for (auto const& n : neighbor_ip_coords)
                {
                    // output
                    // std::cout << "\t[" << std::get<0>(n) << ", "
                    //          << std::get<1>(n) << ", (";
                    // for (int i = 0; i < std::get<2>(n).size(); ++i)
                    //    std::cout << std::get<2>(n)[i] << ", ";
                    // std::cout << "), " << std::get<3>(n) << "]\n";

                    // save into current ip_k
                    _ip_data[k].non_local_assemblers.push_back(std::make_tuple(
                        la.get(), std::get<1>(n), std::get<3>(n),
                        std::numeric_limits<double>::quiet_NaN()));
                    // std::cout << "\tadd "
                    //          << _ip_data[k].non_local_assemblers.size()
                    //          << " points.\n";
                }
            }
            if (_ip_data[k].non_local_assemblers.size() == 0) {
                OGS_FATAL("no neighbours found!");
            }

            double a_k_sum_m = 0;
            for (auto const& tuple_m : _ip_data[k].non_local_assemblers)
            {
                auto const& la_m =
                    *static_cast<SmallDeformationNonlocalLocalAssembler<
                        ShapeFunction, IntegrationMethod,
                        DisplacementDim> const* const>(std::get<0>(tuple_m));

                int const m = std::get<1>(tuple_m);
                double const distance2_m = std::get<2>(tuple_m);

                auto const& w_m = la_m._ip_data[m].integration_weight;

                a_k_sum_m += w_m * alpha_0(distance2_m);

                // int const m_ele = la_m._element.getID();
                // std::cout
                //    << "\tCompute sum_a_km for k = " << k << " and m = ("
                //    << m_ele << ", " << m
                //    << "); distance^2_m = " << distance2_m
                //    << "alpha_0(d^2_m) = " << alpha_0(distance2_m)
                //    << "; sum_alpha_km = " << a_k_sum_m << "\n";
            }


            //
            // Calculate alpha_kl =
            //       alpha_0(|x_k - x_l|) / int_{m \in ip} alpha_0(|x_k - x_m|)
            //
            for (auto& tuple : _ip_data[k].non_local_assemblers)
            {
                // auto const& la_l =
                //    *static_cast<SmallDeformationNonlocalLocalAssembler<
                //        ShapeFunction, IntegrationMethod,
                //        DisplacementDim> const* const>(std::get<0>(tuple));

                double const distance2_l = std::get<2>(tuple);

                // int const l_ele = la_l._element.getID();
                // int const l = std::get<1>(tuple);
                // std::cout << "Compute a_kl for k = " << k << " and l = ("
                //          << l_ele << ", " << l
                //          << "); distance^2_l = " << distance2_l << "\n";
                double const a_kl = alpha_0(distance2_l) / a_k_sum_m;

                // std::cout << "alpha_0(d^2_l) = " << alpha_0(distance2_l)
                //          << "\n";
                // std::cout << "alpha_kl = " << a_kl << "done\n";
                std::get<3>(tuple) = a_kl;
            }
        }
    }

    // TODO (naumov) is this the same as the
    // interpolateXCoordinate<ShapeFunction, ShapeMatricesType>(_element, N) ???
    // NO, it is indeed only X-coordinate
    Eigen::Vector3d getSingleIntegrationPointCoordinates(
        int integration_point) const
    {
        auto const& N = _secondary_data.N[integration_point];

        Eigen::Vector3d xyz = Eigen::Vector3d::Zero();  // Resulting coordinates
        auto* nodes = _element.getNodes();
        for (int i = 0; i < N.size(); ++i)
        {
            Eigen::Map<Eigen::Vector3d const> node_coordinates{
                nodes[i]->getCoords(), 3};
            xyz += node_coordinates * N[i];
        }

        // std::cout << "\t\t singleIPcoords: xyz = " << xyz[0] << " " << xyz[1]
        //          << " " << xyz[2] << "\n";
        return xyz;
    }

    /// \returns for each of the current element's integration points the
    /// element's id, the integration point number, its coordinates, and the
    /// squared distance from the current integration point.
    std::vector<std::tuple<int, int, Eigen::Vector3d, double>>
    getIntegrationPointCoordinates(Eigen::Vector3d const& coords) const override
    {
        unsigned const n_integration_points =
            _integration_method.getNumberOfPoints();

        std::vector<std::tuple<int, int, Eigen::Vector3d, double>> result;
        result.reserve(n_integration_points);

        for (unsigned ip = 0; ip < n_integration_points; ip++)
        {
            // std::cout << _element.getID() << ", " << ip << "\n";

            auto const xyz = getSingleIntegrationPointCoordinates(ip);
            double const distance2 = (xyz - coords).squaredNorm();
            if (distance2 <
                _process_data.internal_length * _process_data.internal_length)
                result.emplace_back(_element.getID(), ip, xyz, distance2);
        }
        // std::cout << "\tfor element " << _element.getID() << " got "
        //          << result.size() << " point in internal_length\n";
        return result;
    }

    void assemble(double const /*t*/, std::vector<double> const& /*local_x*/,
                  std::vector<double>& /*local_M_data*/,
                  std::vector<double>& /*local_K_data*/,
                  std::vector<double>& /*local_b_data*/) override
    {
        OGS_FATAL(
            "SmallDeformationNonlocalLocalAssembler: assembly without jacobian "
            "is not "
            "implemented.");
    }

    void preAssemble(double const t,
                     std::vector<double> const& local_x) override
    {
        // auto const local_matrix_size = local_x.size();

        // auto local_Jac = MathLib::createZeroedMatrix<StiffnessMatrixType>(
        //    local_Jac_data, local_matrix_size, local_matrix_size);

        unsigned const n_integration_points =
            _integration_method.getNumberOfPoints();

        SpatialPosition x_position;
        x_position.setElementID(_element.getID());

        //std::cout << "\nXXX nonlocal in element " << _element.getID() <<
        //"\n";

        for (unsigned ip = 0; ip < n_integration_points; ip++)
        {
            //std::cout << "\n\tip = " << ip << "\n";

            auto const xyz = getSingleIntegrationPointCoordinates(ip);
            //std::cout << "\tCurrent ip_k coords : " << xyz << "\n";

            x_position.setIntegrationPoint(ip);

            auto const& N = _ip_data[ip].N;
            auto const& dNdx = _ip_data[ip].dNdx;

            auto& x = _ip_data[ip].coordinates;
            // Take only the first DisplacementDim values, because the function
            // alway returns 3D vector.
            x = getSingleIntegrationPointCoordinates(ip).head(DisplacementDim);

            auto const x_coord =
                interpolateXCoordinate<ShapeFunction, ShapeMatricesType>(
                    _element, N);
            auto const B = LinearBMatrix::computeBMatrix<
                DisplacementDim, ShapeFunction::NPOINTS,
                typename BMatricesType::BMatrixType>(dNdx, N, x_coord,
                                                     _is_axially_symmetric);
            auto const& eps_prev = _ip_data[ip].eps_prev;
            auto const& sigma_prev = _ip_data[ip].sigma_prev;

            auto& eps = _ip_data[ip].eps;
            auto& sigma = _ip_data[ip].sigma;
            auto& C = _ip_data[ip].C;
            auto& state = _ip_data[ip].material_state_variables;

            eps.noalias() =
                B *
                Eigen::Map<typename BMatricesType::NodalForceVectorType const>(
                    local_x.data(), ShapeFunction::NPOINTS * DisplacementDim);

            // sigma is for plastic part only.
            std::unique_ptr<KelvinMatrixType<DisplacementDim>> new_C;
            std::unique_ptr<typename MaterialLib::Solids::MechanicsBase<
                DisplacementDim>::MaterialStateVariables>
                new_state;
            auto&& solution = _ip_data[ip].solid_material.integrateStress(
                t, x_position, _process_data.dt, eps_prev, eps, sigma_prev,
                *state);

            if (!solution)
                OGS_FATAL("Computation of local constitutive relation failed.");

            std::tie(sigma, state, C) = std::move(*solution);
        }

        // Compute material forces, needed in the non-local assembly, storing
        // them locally and interpolating them to integration points.
        getMaterialForces(local_x, _material_forces);
        for (unsigned ip = 0; ip < n_integration_points; ip++)
        {
            x_position.setIntegrationPoint(ip);

            auto const& N = _ip_data[ip].N;
            auto& g = _ip_data[ip].material_force;
            if (DisplacementDim == 2)
                NumLib::shapeFunctionInterpolate(_material_forces, N, g[0],
                                                 g[1]);
            if (DisplacementDim == 3)
                NumLib::shapeFunctionInterpolate(_material_forces, N, g[0],
                                                 g[1], g[2]);
        }
    }

    void assembleWithJacobian(double const t,
                              std::vector<double> const& local_x,
                              std::vector<double> const& /*local_xdot*/,
                              const double /*dxdot_dx*/, const double /*dx_dx*/,
                              std::vector<double>& /*local_M_data*/,
                              std::vector<double>& /*local_K_data*/,
                              std::vector<double>& local_b_data,
                              std::vector<double>& local_Jac_data) override
    {
        auto const local_matrix_size = local_x.size();

        auto local_Jac = MathLib::createZeroedMatrix<StiffnessMatrixType>(
            local_Jac_data, local_matrix_size, local_matrix_size);

        auto local_b = MathLib::createZeroedVector<NodalDisplacementVectorType>(
            local_b_data, local_matrix_size);

        unsigned const n_integration_points =
            _integration_method.getNumberOfPoints();

        SpatialPosition x_position;
        x_position.setElementID(_element.getID());

        // Compute the non-local internal length depending on the material
        // forces and directions dir := x_ip - \xi:
        // length(x_ip) := 1/Vol \int_{Vol} <g(\xi), dir> / scaling * beta  d\xi,
        // where:
        // beta(x_ip, \xi) := exp(-|dir|) * (1 - exp(-|g(\xi)|)),
        // scaling(x_ip, \xi) := |<g(\xi), dir>|
        for (unsigned ip = 0; ip < n_integration_points; ip++)
        {
            auto const& N = _ip_data[ip].N;
            auto const& x = _ip_data[ip].coordinates;
            auto& l = _ip_data[ip].nonlocal_internal_length;
            l = 0;

            for (auto const& tuple : _ip_data[ip].non_local_assemblers)
            {
                // TODO (naumov) FIXME FIXME FIXME
                // Static cast is not possible because the neighbouring element
                // might have a different ShapeFunction!!!!
                auto const& la_l =
                    *static_cast<SmallDeformationNonlocalLocalAssembler<
                        ShapeFunction, IntegrationMethod,
                        DisplacementDim> const* const>(std::get<0>(tuple));
                int const& ip_l = std::get<1>(tuple);
                double const elements_volume = la_l._element.getContent();

                auto const& w_l = la_l._ip_data[ip_l].integration_weight;
                auto const& g_l = la_l._ip_data[ip_l].material_force;
                auto const& x_l = la_l._ip_data[ip].coordinates;
                GlobalDimVectorType d_l = x - x_l;
                l +=  g_l.dot(d_l) * w_l / elements_volume;
            }
            //INFO("local_length(%d) %g", ip, l)
        }

        // Non-local integration.
        for (unsigned ip = 0; ip < n_integration_points; ip++)
        {
            x_position.setIntegrationPoint(ip);
            auto const& w = _ip_data[ip].integration_weight;

            auto const& N = _ip_data[ip].N;
            auto const& dNdx = _ip_data[ip].dNdx;

            auto const x_coord =
                interpolateXCoordinate<ShapeFunction, ShapeMatricesType>(
                    _element, N);
            auto const B = LinearBMatrix::computeBMatrix<
                DisplacementDim, ShapeFunction::NPOINTS,
                typename BMatricesType::BMatrixType>(dNdx, N, x_coord,
                                                     _is_axially_symmetric);
            // auto const& eps_prev = _ip_data[ip].eps_prev;
            // auto const& sigma_prev = _ip_data[ip].sigma_prev;

            auto& eps = _ip_data[ip].eps;
            auto& sigma = _ip_data[ip].sigma;
            auto& C = _ip_data[ip].C;
            // auto& material_state_variables =
            //    *_ip_data[ip].material_state_variables;

            eps.noalias() =
                B *
                Eigen::Map<typename BMatricesType::NodalForceVectorType const>(
                    local_x.data(), ShapeFunction::NPOINTS * DisplacementDim);

            /*
            if (!_ip_data[ip].solid_material.updateNonlocalDamage(
                    t, x_position, _process_data.dt, eps_prev, eps, sigma_prev,
                    sigma, C, material_state_variables))
                OGS_FATAL("Computation of non-local damage update failed.");
            */

            {
                double test_alpha = 0;  // Integration of one-function.
                double nonlocal_kappa_d = 0;

                for (auto const& tuple : _ip_data[ip].non_local_assemblers)
                {
                    auto const& la_l =
                        *static_cast<SmallDeformationNonlocalLocalAssembler<
                            ShapeFunction, IntegrationMethod,
                            DisplacementDim> const* const>(std::get<0>(tuple));
                    // Get damage from the local assembler and its corresponding
                    // integration point l.
                    int const& l = std::get<1>(tuple);
                    double const kappa_d =
                        static_cast<SmallDeformationNonlocalLocalAssembler<
                            ShapeFunction, IntegrationMethod,
                            DisplacementDim> const* const>(std::get<0>(tuple))
                            ->_ip_data[l]
                            .getLocalVariable();
                    // std::cerr << kappa_d << "\n";
                    double const a_kl = std::get<3>(tuple);

                    auto const& w_l = la_l._ip_data[l].integration_weight;

                    test_alpha += a_kl * w_l;
                    nonlocal_kappa_d += a_kl * kappa_d * w_l;
                }
                if (std::abs(test_alpha - 1) >= 1e-6)
                    OGS_FATAL(
                        "One-function integration failed. v: %f, diff: %f",
                        test_alpha, test_alpha - 1);

                double const _gamma_nonlocal = 1.0;
                // === Overnonlocal formulation ===
                // Update nonlocal damage with local damage (scaled with 1 -
                // \gamma_nonlocal) for the current integration point and the
                // nonlocal integral part.
                nonlocal_kappa_d =
                    (1. - _gamma_nonlocal) * _ip_data[ip].getLocalVariable() +
                    _gamma_nonlocal * nonlocal_kappa_d;

                if (nonlocal_kappa_d < 0.)
                {
                    ERR("set kappa_d zero %g", nonlocal_kappa_d);
                    nonlocal_kappa_d = 0;
                }

                _ip_data[ip].nonlocal_kappa_d = nonlocal_kappa_d;
                _ip_data[ip].damage =
                    _ip_data[ip].updateDamage(t, x_position, nonlocal_kappa_d);
                if (_ip_data[ip].damage < 0. || _ip_data[ip].damage > 1.)
                    std::cerr << "DD " << _ip_data[ip].damage << "\n\n";

                sigma = sigma * (1 - _ip_data[ip].damage);
            }

            local_b.noalias() -= B.transpose() * sigma * w;
            local_Jac.noalias() += B.transpose() * C * B * w;
        }
    }

    void preTimestepConcrete(std::vector<double> const& /*local_x*/,
                             double const /*t*/,
                             double const /*delta_t*/) override
    {
        unsigned const n_integration_points =
            _integration_method.getNumberOfPoints();

        for (unsigned ip = 0; ip < n_integration_points; ip++)
        {
            _ip_data[ip].pushBackState();
        }
    }

    std::vector<double> const& getMaterialForces(
        std::vector<double> const& local_x,
        std::vector<double>& nodal_values) override
    {
        return ProcessLib::SmallDeformation::getMaterialForces<
            DisplacementDim, ShapeFunction, ShapeMatricesType,
            typename BMatricesType::NodalForceVectorType,
            NodalDisplacementVectorType, GradientVectorType,
            GradientMatrixType>(local_x, nodal_values, _integration_method,
                                _ip_data, _element, _is_axially_symmetric);
    }

    friend void readSmallDeformationIntegrationPointData<
        SmallDeformationNonlocalLocalAssembler>(
        std::vector<char> const& data,
        SmallDeformationNonlocalLocalAssembler& local_assembler);
    void readIntegrationPointData(std::vector<char> const& data) override
    {
        readSmallDeformationIntegrationPointData(data, *this);
    }
#ifdef PROTOBUF_FOUND
    friend OGS::SmallDeformationCommon
    getSmallDeformationCommonIntegrationPointData<
        SmallDeformationNonlocalLocalAssembler>(
        SmallDeformationNonlocalLocalAssembler const& local_assembler);
    std::size_t writeIntegrationPointData(std::vector<char>& data) override
    {
        unsigned const n_integration_points =
            _integration_method.getNumberOfPoints();

        OGS::ElementData element_data;
        element_data.set_element_id(_element.getID());
        element_data.set_n_integration_points(n_integration_points);

        auto small_deformation_nonlocal =
            element_data.mutable_small_deformation_nonlocal();
        auto common = small_deformation_nonlocal->mutable_common();
        common->CopyFrom(getSmallDeformationCommonIntegrationPointData(*this));

        {  // SmallDeformationNonlocal specific output.
            unsigned const n_integration_points =
                _integration_method.getNumberOfPoints();

            for (unsigned ip = 0; ip < n_integration_points; ip++)
            {
                small_deformation_nonlocal->add_nonlocal_damage(
                    _ip_data[ip].nonlocal_kappa_d);
                small_deformation_nonlocal->add_nonlocal_length(
                    _ip_data[ip].nonlocal_internal_length);
            }
        }

        data.resize(element_data.ByteSize());
        element_data.SerializeToArray(data.data(), element_data.ByteSize());

        return element_data.ByteSize();
    };
#else
    std::size_t writeIntegrationPointData(std::vector<char>& data) override
    {
        return 0;
    }
#endif

    Eigen::Map<const Eigen::RowVectorXd> getShapeMatrix(
        const unsigned integration_point) const override
    {
        auto const& N = _secondary_data.N[integration_point];

        // assumes N is stored contiguously in memory
        return Eigen::Map<const Eigen::RowVectorXd>(N.data(), N.size());
    }

    std::vector<double> const& getNodalValues(
        std::vector<double>& nodal_values) const override
    {
        nodal_values.clear();
        auto local_b = MathLib::createZeroedVector<NodalDisplacementVectorType>(
            nodal_values, ShapeFunction::NPOINTS * DisplacementDim);

        unsigned const n_integration_points =
            _integration_method.getNumberOfPoints();

        SpatialPosition x_position;
        x_position.setElementID(_element.getID());

        for (unsigned ip = 0; ip < n_integration_points; ip++)
        {
            x_position.setIntegrationPoint(ip);
            auto const& w = _ip_data[ip].integration_weight;

            auto const& N = _ip_data[ip].N;
            auto const& dNdx = _ip_data[ip].dNdx;

            auto const x_coord =
                interpolateXCoordinate<ShapeFunction, ShapeMatricesType>(
                    _element, N);
            auto const B = LinearBMatrix::computeBMatrix<
                DisplacementDim, ShapeFunction::NPOINTS,
                typename BMatricesType::BMatrixType>(dNdx, N, x_coord,
                                                     _is_axially_symmetric);
            auto& sigma = _ip_data[ip].sigma;

            local_b.noalias() += B.transpose() * sigma * w;
        }

        return nodal_values;
    }

    std::vector<double> const& getIntPtFreeEnergyDensity(
        const double /*t*/,
        GlobalVector const& /*current_solution*/,
        NumLib::LocalToGlobalIndexMap const& /*dof_table*/,
        std::vector<double>& cache) const override
    {
        cache.clear();
        cache.reserve(_ip_data.size());

        for (auto const& ip_data : _ip_data)
        {
            cache.push_back(ip_data.free_energy_density);
        }

        return cache;
    }

    std::vector<double> const& getIntPtEpsPV(
        const double /*t*/,
        GlobalVector const& /*current_solution*/,
        NumLib::LocalToGlobalIndexMap const& /*dof_table*/,
        std::vector<double>& cache) const override
    {
        cache.clear();
        cache.reserve(_ip_data.size());

        for (auto const& ip_data : _ip_data)
        {
            cache.push_back(*ip_data.eps_p_V);
        }

        return cache;
    }

    std::vector<double> const& getIntPtEpsPDXX(
        const double /*t*/,
        GlobalVector const& /*current_solution*/,
        NumLib::LocalToGlobalIndexMap const& /*dof_table*/,
        std::vector<double>& cache) const override
    {
        cache.clear();
        cache.reserve(_ip_data.size());

        for (auto const& ip_data : _ip_data)
        {
            cache.push_back(*ip_data.eps_p_D_xx);
        }

        return cache;
    }

    std::vector<double> const& getIntPtSigma(
        const double /*t*/,
        GlobalVector const& /*current_solution*/,
        NumLib::LocalToGlobalIndexMap const& /*dof_table*/,
        std::vector<double>& cache) const override
    {
        using KelvinVectorType = typename BMatricesType::KelvinVectorType;
        auto const kelvin_vector_size =
            KelvinVectorDimensions<DisplacementDim>::value;
        auto const num_intpts = _ip_data.size();

        cache.clear();
        auto cache_mat = MathLib::createZeroedMatrix<Eigen::Matrix<
            double, kelvin_vector_size, Eigen::Dynamic, Eigen::RowMajor>>(
            cache, kelvin_vector_size, num_intpts);

        // TODO make a general implementation for converting KelvinVectors
        // back to symmetric rank-2 tensors.
        for (unsigned ip = 0; ip < num_intpts; ++ip)
        {
            auto const& sigma = _ip_data[ip].sigma;

            for (typename KelvinVectorType::Index component = 0;
                 component < kelvin_vector_size && component < 3;
                 ++component)
            {  // xx, yy, zz components
                cache_mat(component, ip) = sigma[component];
            }
            for (typename KelvinVectorType::Index component = 3;
                 component < kelvin_vector_size;
                 ++component)
            {  // mixed xy, yz, xz components
                cache_mat(component, ip) = sigma[component] / std::sqrt(2);
            }
        }

        return cache;
    }

    std::vector<double> const& getIntPtDamage(
        const double /*t*/,
        GlobalVector const& /*current_solution*/,
        NumLib::LocalToGlobalIndexMap const& /*dof_table*/,
        std::vector<double>& cache) const override
    {
        cache.clear();
        cache.reserve(_ip_data.size());

        for (auto const& ip_data : _ip_data)
        {
            auto const& state =
                static_cast<MaterialLib::Solids::Ehlers::StateVariables<
                    DisplacementDim> const&>(*ip_data.material_state_variables);
            cache.push_back(state.damage.value());
        }

        return cache;
    }

    // TODO remove the component-wise methods.
    std::vector<double> const& getIntPtSigmaXX(
        const double /*t*/,
        GlobalVector const& /*current_solution*/,
        NumLib::LocalToGlobalIndexMap const& /*dof_table*/,
        std::vector<double>& cache) const override
    {
        return getIntPtSigma(cache, 0);
    }

    std::vector<double> const& getIntPtSigmaYY(
        const double /*t*/,
        GlobalVector const& /*current_solution*/,
        NumLib::LocalToGlobalIndexMap const& /*dof_table*/,
        std::vector<double>& cache) const override
    {
        return getIntPtSigma(cache, 1);
    }

    std::vector<double> const& getIntPtSigmaZZ(
        const double /*t*/,
        GlobalVector const& /*current_solution*/,
        NumLib::LocalToGlobalIndexMap const& /*dof_table*/,
        std::vector<double>& cache) const override
    {
        return getIntPtSigma(cache, 2);
    }

    std::vector<double> const& getIntPtSigmaXY(
        const double /*t*/,
        GlobalVector const& /*current_solution*/,
        NumLib::LocalToGlobalIndexMap const& /*dof_table*/,
        std::vector<double>& cache) const override
    {
        return getIntPtSigma(cache, 3);
    }

    std::vector<double> const& getIntPtSigmaXZ(
        const double /*t*/,
        GlobalVector const& /*current_solution*/,
        NumLib::LocalToGlobalIndexMap const& /*dof_table*/,
        std::vector<double>& cache) const override
    {
        assert(DisplacementDim == 3);
        return getIntPtSigma(cache, 4);
    }

    std::vector<double> const& getIntPtSigmaYZ(
        const double /*t*/,
        GlobalVector const& /*current_solution*/,
        NumLib::LocalToGlobalIndexMap const& /*dof_table*/,
        std::vector<double>& cache) const override
    {
        assert(DisplacementDim == 3);
        return getIntPtSigma(cache, 5);
    }

    virtual std::vector<double> const& getIntPtEpsilon(
        const double /*t*/,
        GlobalVector const& /*current_solution*/,
        NumLib::LocalToGlobalIndexMap const& /*dof_table*/,
        std::vector<double>& cache) const override
    {
        using KelvinVectorType = typename BMatricesType::KelvinVectorType;
        auto const kelvin_vector_size =
            KelvinVectorDimensions<DisplacementDim>::value;
        auto const num_intpts = _ip_data.size();

        cache.clear();
        auto cache_mat = MathLib::createZeroedMatrix<Eigen::Matrix<
            double, kelvin_vector_size, Eigen::Dynamic, Eigen::RowMajor>>(
            cache, kelvin_vector_size, num_intpts);

        // TODO make a general implementation for converting KelvinVectors
        // back to symmetric rank-2 tensors.
        for (unsigned ip = 0; ip < num_intpts; ++ip)
        {
            auto const& eps = _ip_data[ip].eps;

            for (typename KelvinVectorType::Index component = 0;
                 component < kelvin_vector_size && component < 3;
                 ++component)
            {  // xx, yy, zz components
                cache_mat(component, ip) = eps[component];
            }
            for (typename KelvinVectorType::Index component = 3;
                 component < kelvin_vector_size;
                 ++component)
            {  // mixed xy, yz, xz components
                cache_mat(component, ip) = eps[component] / std::sqrt(2);
            }
        }

        return cache;
    }

    // TODO remove the component-wise methods
    std::vector<double> const& getIntPtEpsilonXX(
        const double /*t*/,
        GlobalVector const& /*current_solution*/,
        NumLib::LocalToGlobalIndexMap const& /*dof_table*/,
        std::vector<double>& cache) const override
    {
        return getIntPtEpsilon(cache, 0);
    }

    std::vector<double> const& getIntPtEpsilonYY(
        const double /*t*/,
        GlobalVector const& /*current_solution*/,
        NumLib::LocalToGlobalIndexMap const& /*dof_table*/,
        std::vector<double>& cache) const override
    {
        return getIntPtEpsilon(cache, 1);
    }

    std::vector<double> const& getIntPtEpsilonZZ(
        const double /*t*/,
        GlobalVector const& /*current_solution*/,
        NumLib::LocalToGlobalIndexMap const& /*dof_table*/,
        std::vector<double>& cache) const override
    {
        return getIntPtEpsilon(cache, 2);
    }

    std::vector<double> const& getIntPtEpsilonXY(
        const double /*t*/,
        GlobalVector const& /*current_solution*/,
        NumLib::LocalToGlobalIndexMap const& /*dof_table*/,
        std::vector<double>& cache) const override
    {
        return getIntPtEpsilon(cache, 3);
    }

    std::vector<double> const& getIntPtEpsilonYZ(
        const double /*t*/,
        GlobalVector const& /*current_solution*/,
        NumLib::LocalToGlobalIndexMap const& /*dof_table*/,
        std::vector<double>& cache) const override
    {
        assert(DisplacementDim == 3);
        return getIntPtEpsilon(cache, 4);
    }

    std::vector<double> const& getIntPtEpsilonXZ(
        const double /*t*/,
        GlobalVector const& /*current_solution*/,
        NumLib::LocalToGlobalIndexMap const& /*dof_table*/,
        std::vector<double>& cache) const override
    {
        assert(DisplacementDim == 3);
        return getIntPtEpsilon(cache, 5);
    }

    unsigned getNumberOfIntegrationPoints() const override
    {
        return _integration_method.getNumberOfPoints();
    }

    typename MaterialLib::Solids::MechanicsBase<
        DisplacementDim>::MaterialStateVariables const&
    getMaterialStateVariablesAt(unsigned integration_point) const override
    {
        return *_ip_data[integration_point].material_state_variables;
    }

private:
    std::vector<double> const& getIntPtSigma(std::vector<double>& cache,
                                             std::size_t const component) const
    {
        cache.clear();
        cache.reserve(_ip_data.size());

        for (auto const& ip_data : _ip_data)
        {
            if (component < 3)  // xx, yy, zz components
                cache.push_back(ip_data.sigma[component]);
            else  // mixed xy, yz, xz components
                cache.push_back(ip_data.sigma[component] / std::sqrt(2));
        }

        return cache;
    }

    std::vector<double> const& getIntPtEpsilon(
        std::vector<double>& cache, std::size_t const component) const
    {
        cache.clear();
        cache.reserve(_ip_data.size());

        for (auto const& ip_data : _ip_data)
        {
            if (component < 3)  // xx, yy, zz components
                cache.push_back(ip_data.eps[component]);
            else  // mixed xy, yz, xz components
                cache.push_back(ip_data.eps[component] / std::sqrt(2));
        }

        return cache;
    }

    SmallDeformationNonlocalProcessData<DisplacementDim>& _process_data;

    std::vector<
        IntegrationPointData<BMatricesType, ShapeMatricesType, DisplacementDim>,
        Eigen::aligned_allocator<IntegrationPointData<
            BMatricesType, ShapeMatricesType, DisplacementDim>>>
        _ip_data;

    std::vector<double> _material_forces;

    IntegrationMethod _integration_method;
    MeshLib::Element const& _element;
    SecondaryData<typename ShapeMatrices::ShapeType> _secondary_data;
    bool const _is_axially_symmetric;
};

}  // namespace SmallDeformationNonlocal
}  // namespace ProcessLib
