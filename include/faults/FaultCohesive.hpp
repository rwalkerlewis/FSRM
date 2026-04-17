#ifndef FSRM_FAULTS_FAULT_COHESIVE_HPP
#define FSRM_FAULTS_FAULT_COHESIVE_HPP

#include <petsc.h>
#include <petscdm.h>

namespace FSRM {

// Forward declaration: Phase 2 subclasses will take a populated FaultConfig
// for initialize(). Declared here so the pure-virtual signature is
// stable from Phase 1 onward.
struct FaultConfig;

namespace faults {

/**
 * @brief Abstract base for fault interface objects (PyLith-style).
 *
 * Phase 1 scope: provides the cohesive interface DMLabel helper and an
 * accessor. The initialize() hook is pure virtual; concrete subclasses
 * (FaultCohesiveKin in Phase 3, FaultCohesiveImpulses in a future phase)
 * will be added when the kernel and aux field infrastructure lands.
 *
 * Design mirrors PyLith `pylith::faults::FaultCohesive`
 * (`libsrc/pylith/faults/FaultCohesive.hh` @ geodynamics/pylith
 * 6accf7e7395aa792444ad1b489d455e41ece34bc).
 */
class FaultCohesive {
public:
    FaultCohesive() = default;
    virtual ~FaultCohesive() = default;

    // Non-copyable (owns PETSc state in concrete subclasses).
    FaultCohesive(const FaultCohesive &) = delete;
    FaultCohesive &operator=(const FaultCohesive &) = delete;

    /**
     * @brief Subclass-specific setup. Implemented in Phase 2+.
     */
    virtual PetscErrorCode initialize(DM dm, const FaultConfig &config) = 0;

    /**
     * @brief Create (or fetch) the "cohesive interface" DMLabel on @p dm.
     *
     * Thin wrapper over
     * @ref FSRM::faults::topology::createInterfacesLabel . Caches the
     * returned label pointer in @ref cohesive_interface_label_ so that
     * @ref getCohesiveInterfaceLabel returns a non-NULL value without
     * re-querying the DM.
     *
     * @param dm           DMPlex after DMPlexConstructCohesiveCells.
     * @param fault_label  Reserved for Phase 2. Pass NULL in Phase 1.
     */
    PetscErrorCode createCohesiveInterfaceLabel(DM dm, DMLabel fault_label);

    /**
     * @brief Accessor for the cached DMLabel pointer.
     *
     * @return DMLabel or nullptr if @ref createCohesiveInterfaceLabel has
     *         not been called successfully yet.
     */
    DMLabel getCohesiveInterfaceLabel() const { return cohesive_interface_label_; }

protected:
    DMLabel cohesive_interface_label_ = nullptr;
};

}  // namespace faults
}  // namespace FSRM

#endif  // FSRM_FAULTS_FAULT_COHESIVE_HPP
