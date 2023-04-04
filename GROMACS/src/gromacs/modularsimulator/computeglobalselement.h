/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2019,2020, by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed in the AUTHORS file in the
 * top-level source directory and at http://www.gromacs.org.
 *
 * GROMACS is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * GROMACS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GROMACS; if not, see
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
 */
/*! \internal \file
 * \brief Declares the global reduction element for the modular simulator
 *
 * \author Pascal Merz <pascal.merz@me.com>
 * \ingroup module_modularsimulator
 *
 * This header is only used within the modular simulator module
 */

#ifndef GMX_MODULARSIMULATOR_COMPUTEGLOBALSELEMENT_H
#define GMX_MODULARSIMULATOR_COMPUTEGLOBALSELEMENT_H

#include "gromacs/mdlib/simulationsignal.h"
#include "gromacs/mdlib/vcm.h"

#include "energydata.h"
#include "modularsimulatorinterfaces.h"
#include "statepropagatordata.h"
#include "topologyholder.h"

struct gmx_global_stat;
struct gmx_wallcycle;
struct t_nrnb;

namespace gmx
{
class FreeEnergyPerturbationData;
class LegacySimulatorData;
class MDAtoms;
class MDLogger;

//! \addtogroup module_modularsimulator
//! \{

//! The different global reduction schemes we know about
enum class ComputeGlobalsAlgorithm
{
    LeapFrog,
    VelocityVerlet
};

//! The function type allowing to request a check of the number of bonded interactions
typedef std::function<void()> CheckBondedInteractionsCallback;

/*! \internal
 * \brief Encapsulate the calls to `compute_globals`
 *
 * This element aims at offering an interface to the legacy
 * implementation which is compatible with the new simulator approach.
 *
 * The element comes in 3 (templated) flavors: the leap-frog case, the first
 * call during a velocity-verlet integrator, and the second call during a
 * velocity-verlet integrator. In velocity verlet, the state at the beginning
 * of the step corresponds to
 *     positions at time t
 *     velocities at time t - dt/2
 * The first velocity propagation (+dt/2) therefore actually corresponds to the
 * previous step, bringing the state to the full timestep at time t. Most global
 * reductions are made at this point. The second call is needed to correct the
 * constraint virial after the second propagation of velocities (+dt/2) and of
 * the positions (+dt).
 *
 * \tparam algorithm  The global reduction scheme
 */
template<ComputeGlobalsAlgorithm algorithm>
class ComputeGlobalsElement final :
    public ISimulatorElement,
    public IEnergySignallerClient,
    public ITrajectorySignallerClient,
    public ITopologyHolderClient
{
public:
    //! Constructor
    ComputeGlobalsElement(StatePropagatorData*        statePropagatorData,
                          EnergyData*                 energyData,
                          FreeEnergyPerturbationData* freeEnergyPerturbationData,
                          SimulationSignals*          signals,
                          int                         nstglobalcomm,
                          FILE*                       fplog,
                          const MDLogger&             mdlog,
                          t_commrec*                  cr,
                          const t_inputrec*           inputrec,
                          const MDAtoms*              mdAtoms,
                          t_nrnb*                     nrnb,
                          gmx_wallcycle*              wcycle,
                          t_forcerec*                 fr,
                          const gmx_mtop_t*           global_top,
                          Constraints*                constr);

    //! Destructor
    ~ComputeGlobalsElement() override;

    /*! \brief Element setup - first call to compute_globals
     *
     */
    void elementSetup() override;

    /*! \brief Register run function for step / time
     *
     * This registers the call to compute_globals when needed.
     *
     * \param step                 The step number
     * \param time                 The time
     * \param registerRunFunction  Function allowing to register a run function
     */
    void scheduleTask(Step step, Time time, const RegisterRunFunction& registerRunFunction) override;

    //! Get callback to request checking of bonded interactions
    CheckBondedInteractionsCallback getCheckNumberOfBondedInteractionsCallback();

    //! No element teardown needed
    void elementTeardown() override {}

    /*! \brief Factory method implementation
     *
     * \param legacySimulatorData  Pointer allowing access to simulator level data
     * \param builderHelper  ModularSimulatorAlgorithmBuilder helper object
     * \param statePropagatorData  Pointer to the \c StatePropagatorData object
     * \param energyData  Pointer to the \c EnergyData object
     * \param freeEnergyPerturbationData  Pointer to the \c FreeEnergyPerturbationData object
     * \param globalCommunicationHelper  Pointer to the \c GlobalCommunicationHelper object
     *
     * \throws std::bad_any_cast  on internal error in VelocityVerlet algorithm builder.
     * \throws std::bad_alloc  when out of memory.
     *
     * \return  Pointer to the element to be added. Element needs to have been stored using \c storeElement
     */
    static ISimulatorElement* getElementPointerImpl(LegacySimulatorData* legacySimulatorData,
                                                    ModularSimulatorAlgorithmBuilderHelper* builderHelper,
                                                    StatePropagatorData*        statePropagatorData,
                                                    EnergyData*                 energyData,
                                                    FreeEnergyPerturbationData* freeEnergyPerturbationData,
                                                    GlobalCommunicationHelper* globalCommunicationHelper);

private:
    //! ITopologyClient implementation
    void setTopology(const gmx_localtop_t* top) override;
    //! IEnergySignallerClient implementation
    std::optional<SignallerCallback> registerEnergyCallback(EnergySignallerEvent event) override;
    //! ITrajectorySignallerClient implementation
    std::optional<SignallerCallback> registerTrajectorySignallerCallback(TrajectoryEvent event) override;
    //! The compute_globals call
    void compute(Step step, unsigned int flags, SimulationSignaller* signaller, bool useLastBox, bool isInit = false);

    //! Next step at which energy needs to be reduced
    Step energyReductionStep_;
    //! Next step at which virial needs to be reduced
    Step virialReductionStep_;

    //! For VV only, we need to schedule twice per step. This keeps track of the scheduling stage.
    Step vvSchedulingStep_;

    //! Whether center of mass motion stopping is enabled
    const bool doStopCM_;
    //! Number of steps after which center of mass motion is removed
    int nstcomm_;
    //! Compute globals communication period
    int nstglobalcomm_;
    //! The last (planned) step (only used for LF)
    const Step lastStep_;
    //! The initial step (only used for VV)
    const Step initStep_;
    //! A dummy signaller (used for setup and VV)
    std::unique_ptr<SimulationSignaller> nullSignaller_;

    /*! \brief Check that DD doesn't miss bonded interactions
     *
     * Domain decomposition could incorrectly miss a bonded
     * interaction, but checking for that requires a global
     * communication stage, which does not otherwise happen in DD
     * code. So we do that alongside the first global energy reduction
     * after a new DD is made. These variables handle whether the
     * check happens, and the result it returns.
     */
    //! \{
    int  totalNumberOfBondedInteractions_;
    bool shouldCheckNumberOfBondedInteractions_;
    //! \}

    /*! \brief Signal to ComputeGlobalsElement that it should check for DD errors
     *
     * Note that this should really be the responsibility of the DD element.
     * MDLogger, global and local topology are only needed due to the call to
     * checkNumberOfBondedInteractions(...).
     *
     * The DD element should have a single variable which gets reduced, and then
     * be responsible for the checking after a global reduction has happened.
     * This would, however, require a new approach for the compute_globals calls,
     * which is not yet implemented. So for now, we're leaving this here.
     */
    void needToCheckNumberOfBondedInteractions();

    //! Global reduction struct
    gmx_global_stat* gstat_;

    // TODO: Clarify relationship to data objects and find a more robust alternative to raw pointers (#3583)
    //! Pointer to the microstate
    StatePropagatorData* statePropagatorData_;
    //! Pointer to the energy data (needed for the tensors and mu_tot)
    EnergyData* energyData_;
    //! Pointer to the local topology (only needed for checkNumberOfBondedInteractions)
    const gmx_localtop_t* localTopology_;
    //! Pointer to the free energy perturbation data
    FreeEnergyPerturbationData* freeEnergyPerturbationData_;

    //! Center of mass motion removal
    t_vcm vcm_;
    //! Signals
    SimulationSignals* signals_;

    // Access to ISimulator data
    //! Handles logging.
    FILE* fplog_;
    //! Handles logging.
    const MDLogger& mdlog_;
    //! Handles communication.
    t_commrec* cr_;
    //! Contains user input mdp options.
    const t_inputrec* inputrec_;
    //! Full system topology - only needed for checkNumberOfBondedInteractions.
    const gmx_mtop_t* top_global_;
    //! Atom parameters for this domain.
    const MDAtoms* mdAtoms_;
    //! Handles constraints.
    Constraints* constr_;
    //! Manages flop accounting.
    t_nrnb* nrnb_;
    //! Manages wall cycle accounting.
    gmx_wallcycle* wcycle_;
    //! Parameters for force calculations.
    t_forcerec* fr_;
};

//! \}
} // namespace gmx

#endif // GMX_MODULARSIMULATOR_COMPUTEGLOBALSELEMENT_H
