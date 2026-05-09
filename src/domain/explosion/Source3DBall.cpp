/**
 * @file Source3DBall.cpp
 * @brief Pass-11 (axis 1b) scaffolded factory for the 3D source-ball
 *        solver. Pass-11 throws on construction; pass-12 replaces
 *        the throw with a concrete implementation.
 */

#include "domain/explosion/Source3DBall.hpp"

#include <stdexcept>

namespace FSRM {

std::unique_ptr<Source3DBall> makeSource3DBall(const Source3DBallConfig& cfg)
{
    (void)cfg;
    throw std::runtime_error(
        "Source3DBall: axis-1b 3D source ball is pass-12 work, not yet "
        "implemented. See docs/AXIS_1B_DESIGN.md for the design stub. "
        "Use cavity_geometry=SPHERICAL for the pass-10 1D radial "
        "Lagrangian path.");
}

}  // namespace FSRM
