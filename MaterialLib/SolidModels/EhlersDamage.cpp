/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "EhlersDamage.h"

namespace MaterialLib
{
namespace Solids
{
namespace EhlersDamage
{
template class SolidEhlersDamage<2>;
template class SolidEhlersDamage<3>;

}  // namespace EhlersDamage
}  // namespace Solids
}  // namespace MaterialLib

