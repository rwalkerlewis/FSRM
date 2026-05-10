/**
 * @file Source3DBall.cpp
 * @brief Pass-13a (axis-1b foundation) factory wiring for the 3D
 *        source-ball solver.
 *
 *  Pass-11 throw replaced with construction of Source3DBallImpl. The
 *  factory returns an uninitialized Source3DBallImpl; callers must
 *  call initialize(cfg) to actually load the mesh. The factory itself
 *  performs no I/O so the foundation regression test
 *  (BackwardCompat.Source3DBallScaffoldThrowGoneOnInstantiation)
 *  succeeds even on hosts without a mesh fixture.
 */

#include "domain/explosion/Source3DBall.hpp"

#include "domain/explosion/Source3DBallImpl.hpp"

namespace FSRM {

std::unique_ptr<Source3DBall> makeSource3DBall(const Source3DBallConfig& cfg)
{
    (void)cfg;
    return std::unique_ptr<Source3DBall>(new Source3DBallImpl());
}

}  // namespace FSRM
