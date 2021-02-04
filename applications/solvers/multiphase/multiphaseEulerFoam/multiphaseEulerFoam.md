# multiphaseEulerFoam

- [multiphaseEulerFoam](#multiphaseeulerfoam)
  - [phaseSystems](#phasesystems)
    - [phaseSystems](#phasesystems-1)
      - [MomentumTransferPhaseSystem](#momentumtransferphasesystem)
        - [MomentumTransferPhaseSystem.H](#momentumtransferphasesystemh)
        - [MomentumTransferPhaseSystem.C](#momentumtransferphasesystemc)
          - [DByAfs](#dbyafs)
    - [phaseModel](#phasemodel)
      - [Summary](#summary)
      - [phaseModel](#phasemodel-1)
        - [phaseModel.H](#phasemodelh)
        - [phaseModel.C](#phasemodelc)
          - [constructor](#constructor)
        - [phaseModelNew.C](#phasemodelnewc)
        - [phaseModels.C](#phasemodelsc)
          - [purePhaseModel](#purephasemodel)
          - [pureStationaryPhaseModel](#purestationaryphasemodel)
          - [pureIsothermalPhaseModel](#pureisothermalphasemodel)
          - [pureStationaryIsothermalPhaseModel](#purestationaryisothermalphasemodel)
          - [multiComponentPhaseModel](#multicomponentphasemodel)
          - [multiComponentIsothermalPhaseModel](#multicomponentisothermalphasemodel)
          - [reactingPhaseModel](#reactingphasemodel)
      - [AnisothermalPhaseModel](#anisothermalphasemodel)
        - [AnisothermalPhaseModel.C](#anisothermalphasemodelc)
        - [AnisothermalPhaseModel.H](#anisothermalphasemodelh)
      - [InertPhaseModel](#inertphasemodel)
        - [InertPhaseModel.C](#inertphasemodelc)
        - [InertPhaseModel.H](#inertphasemodelh)
      - [IsothermalPhaseModel](#isothermalphasemodel)
        - [IsothermalPhaseModel.C](#isothermalphasemodelc)
        - [IsothermalPhaseModel.H](#isothermalphasemodelh)
      - [MovingPhaseModel](#movingphasemodel)
        - [MovingPhaseModel.C](#movingphasemodelc)
        - [MovingPhaseModel.H](#movingphasemodelh)
      - [MultiComponentPhaseModel](#multicomponentphasemodel-1)
        - [MultiComponentPhaseModel.C](#multicomponentphasemodelc)
        - [MultiComponentPhaseModel.H](#multicomponentphasemodelh)
      - [PurePhaseModel](#purephasemodel-1)
        - [PurePhaseModel.C](#purephasemodelc)
        - [PurePhaseModel.H](#purephasemodelh)
      - [ReactingPhaseModel](#reactingphasemodel-1)
        - [ReactingPhaseModel.C](#reactingphasemodelc)
        - [ReactingPhaseModel.H](#reactingphasemodelh)
      - [StationaryPhaseModel](#stationaryphasemodel)
        - [StationaryPhaseModel.C](#stationaryphasemodelc)
        - [StationaryPhaseModel.H](#stationaryphasemodelh)
      - [ThermoPhaseModel](#thermophasemodel)
        - [ThermoPhaseModel.C](#thermophasemodelc)
        - [ThermoPhaseModel.H](#thermophasemodelh)
  - [multiPhaseEulerFoam](#multiphaseeulerfoam-1)
    - [createFieldRefs.H](#createfieldrefsh)
    - [createFields.H](#createfieldsh)
    - [CourantNo.H](#courantnoh)
    - [setDeltaT.H](#setdeltath)
    - [YEqns.H](#yeqnsh)
    - [EEqns.H](#eeqnsh)
    - [PU](#pu)
      - [UEqns.H](#ueqnsh)
      - [pEqn.H](#peqnh)
        - [Face volume fractions](#face-volume-fractions)
        - [Diagnocal coefficients](#diagnocal-coefficients)
        - [Phase diagonal coefficients](#phase-diagonal-coefficients)
        - [Explicit force fluxes](#explicit-force-fluxes)
        - [Pressure corrector loop](#pressure-corrector-loop)
          - [Correct fixed-flux BCs to be consistent with the velocity BCs](#correct-fixed-flux-bcs-to-be-consistent-with-the-velocity-bcs)
          - [Combined buoyancy and force fluxes](#combined-buoyancy-and-force-fluxes)
          - [Predicted velocities and fluxes for each phase](#predicted-velocities-and-fluxes-for-each-phase)
          - [Add explicit drag forces and fluxes](#add-explicit-drag-forces-and-fluxes)
          - [Total predicted flux](#total-predicted-flux)
          - [Construct pressure "diffusivity"](#construct-pressure-diffusivity)
          - [Update the fixedFluxPressure BCs to ensure flux consistency](#update-the-fixedfluxpressure-bcs-to-ensure-flux-consistency)
          - [Compressible pressure equations](#compressible-pressure-equations)
          - [Cache p prior to solve for density update](#cache-p-prior-to-solve-for-density-update)
          - [Iterate over the pressure equation to correct for non-orthogonality](#iterate-over-the-pressure-equation-to-correct-for-non-orthogonality)
          - [Update and limit the static pressure](#update-and-limit-the-static-pressure)
          - [Account for static pressure reference](#account-for-static-pressure-reference)
          - [Limit p_rgh](#limit-p_rgh)
          - [Update densities from change in p_rgh](#update-densities-from-change-in-p_rgh)
          - [Correct p_rgh for consistency with p and the updated densities](#correct-p_rgh-for-consistency-with-p-and-the-updated-densities)
        - [clear rAUs](#clear-raus)
    - [PUf](#puf)
    - [multiphaseSystems](#multiphasesystems)
      - [make](#make)
        - [files](#files)
        - [options](#options)
      - [multiphaseSystems.C](#multiphasesystemsc)
        - [basicMultiphaseSystem](#basicmultiphasesystem)
        - [interfaceCompositionPhaseChangeMultiphaseSystem](#interfacecompositionphasechangemultiphasesystem)
        - [thermalPhaseChangeMultiphaseSystem](#thermalphasechangemultiphasesystem)
        - [populationBalanceMultiphaseSystem](#populationbalancemultiphasesystem)
        - [thermalPhaseChangePopulationBalanceMultiphaseSystem](#thermalphasechangepopulationbalancemultiphasesystem)

## phaseSystems

### phaseSystems

#### MomentumTransferPhaseSystem

##### MomentumTransferPhaseSystem.H

Class which models interfacial momentum transfer between a number of phases. Drag, virtual mass, lift, wall lubrication and turbulent dispersion are all modelled. The explicit contribution from the drag is omitted from the transfer matrices, as this forms part of the solution of the pressure equation.

##### MomentumTransferPhaseSystem.C

###### DByAfs

```cpp
template<class BasePhaseSystem>
Foam::PtrList<Foam::surfaceScalarField>
Foam::MomentumTransferPhaseSystem<BasePhaseSystem>::DByAfs
(
    const PtrList<volScalarField>& rAUs,
    const PtrList<surfaceScalarField>& rAUfs
) const
{
    ...
}
```

definition of `DByAfs`

```cpp
    PtrList<surfaceScalarField> DByAfs(this->phaseModels_.size());

    if (rAUfs.size())
    {
        ...
    }
```

* define the list of `DByAfs`
* if `rAUfs` exist, start the calculation

**Add the phase pressure**

```cpp
        // Add the phase pressure
        forAll(this->phaseModels_, phasei)
        {
            const phaseModel& phase = this->phaseModels_[phasei];

            const surfaceScalarField pPrimeByAf
            (
                rAUfs[phasei]*fvc::interpolate(phase.pPrime())
            );

            addField(phase, "DByAf", pPrimeByAf, DByAfs);

            forAll(this->phaseModels_, phasej)
            {
                if (phasej != phasei)
                {
                    const phaseModel& phase2 = this->phaseModels_[phasej];

                    addField
                    (
                        phase2,
                        "DByAf",
                        fvc::interpolate
                        (
                            phase2
                           /max(1 - phase, phase2.residualAlpha())
                        )*pPrimeByAf,
                        DByAfs
                    );
                }
            }
        }
```

Add the phase pressure:
* for every phase
  * get current `phase`, or $\alpha^k$
  * define `pPrimeByAf` as
    * $$pPrimeByAf^k = rAUf^k \left(\frac{\partial p^k}{\partial \alpha^k}\right)_f = \|\sum_{l=1}^N \alpha_f^l (\frac{\alpha}{A_p})^l_f\| \left(\frac{\partial p^k}{\partial \alpha^k}\right)_f$$
  * add `DByAf` to field, and `DByAf` is initialized as `pPrimeByAf`
    * $$DByAf^k = pPrimeByAf^k = \|\sum_{l=1}^N \alpha_f^l (\frac{\alpha}{A_p})^l_f\| \left(\frac{\partial p^k}{\partial \alpha^k}\right)_f$$
  * for every phasej
    * if `phasej` is not current `phase`
      * get `phase2` or $\alpha^l$
      * addField to or create `DByAf`, namely
        * $$DbyAf^k = pPrimeByAf^k + \sum_{l = 1, l \neq k}^N \left[\frac{\alpha^l}{\max(1-\alpha^k, residualAlpha()^l)}\right]_f \cdot pPrimeByAf^k$$

`pPrime()` can be found in `applications\solvers\multiphase\multiphaseEulerFoam\phaseSystems\phaseModel\MovingPhaseModel\MovingPhaseModel.C` is to return the phase-pressure' (derivative of phase-pressure w.r.t. phase-fraction), namely,

$$
pPrime() = \frac{\partial p^k}{\partial \alpha^k}
$$

**Add the turbulent dispersion**

```cpp
        // Add the turbulent dispersion
        forAllConstIter
        (
            turbulentDispersionModelTable,
            turbulentDispersionModels_,
            turbulentDispersionModelIter
        )
        {
            const phasePair&
                pair(this->phasePairs_[turbulentDispersionModelIter.key()]);

            const surfaceScalarField Df
            (
                fvc::interpolate(turbulentDispersionModelIter()->D())
            );

            const surfaceScalarField alpha12f
            (
                fvc::interpolate(pair.phase1() + pair.phase2())
            );

            addField
            (
                pair.phase1(),
                "DByAf",
                rAUfs[pair.phase1().index()]*Df
               /max(alpha12f, pair.phase1().residualAlpha()),
                DByAfs
            );
            addField
            (
                pair.phase2(),
                "DByAf",
                rAUfs[pair.phase2().index()]*Df
               /max(alpha12f, pair.phase2().residualAlpha()),
                DByAfs
            );
        }
```

Add the turbulent dispersion

for every turbulent dispersion needs to be added
  * set `pair`
  * get turbulent diffusivity `D` at surface
    * $$Df = (D)_f$$
  * get the sum of phase fraction of phase 1 and phase 2
    * $$alpah12f^m = (\alpha_1^m + \alpha_2^m)_f$$
  * add `phase1` field to `DByAf`
    * $$DByAf^k = DByAf^k + \sum_{m}^M \frac{rAUf^m_1 (D)_f}{\max(alpha12f^m, residualAlpha()^m_1)}$$
  * add `phase2` field to `DByAf`
    * $$DByAf^k = DByAf^k + \sum_{m}^M \frac{rAUf^m_2 (D)_f}{\max(alpha12f^m, residualAlpha()^m_2)}$$
  * namely, 
  * $$DByAf^k = DByAf^k + \sum_{m}^M \frac{rAUf^m_1 (D)_f}{\max(alpha12f^m, residualAlpha()^m_1)} + \sum_{m}^M \frac{rAUf^m_2 (D)_f}{\max(alpha12f^m, residualAlpha()^m_2)}$$
    * where $M$ is the number of pairs of turbulent dispersion

```cpp
    else
    {
        ...
    }
```

else, `rAUfs` does not exist, then

**Add the phase pressure**

```cpp
        // Add the phase pressure
        forAll(this->phaseModels_, phasei)
        {
            const phaseModel& phase = this->phaseModels_[phasei];

            const surfaceScalarField pPrimeByAf
            (
                fvc::interpolate(rAUs[phasei]*phase.pPrime())
            );

            addField(phase, "DByAf", pPrimeByAf, DByAfs);

            forAll(this->phaseModels_, phasej)
            {
                if (phasej != phasei)
                {
                    const phaseModel& phase2 = this->phaseModels_[phasej];

                    addField
                    (
                        phase2,
                        "DByAf",
                        fvc::interpolate
                        (
                            phase2
                           /max(1 - phase, phase2.residualAlpha())
                        )*pPrimeByAf,
                        DByAfs
                    );
                }
            }
        }
```

the difference between this and above is that, since `rAUfs` is not deifne here, in deifinition of `pPrimeByAf`, `rAUs` is also interpolated to obtain the on surface values

* $$pPrimeByAf = \left(rAU^k \frac{\partial p^k}{\partial \alpha^k} \right)_f$$

others are all same

**Add the turbulent dispersion**

```cpp
        // Add the turbulent dispersion
        forAllConstIter
        (
            turbulentDispersionModelTable,
            turbulentDispersionModels_,
            turbulentDispersionModelIter
        )
        {
            const phasePair&
                pair(this->phasePairs_[turbulentDispersionModelIter.key()]);

            const volScalarField D(turbulentDispersionModelIter()->D());

            const volScalarField alpha12(pair.phase1() + pair.phase2());

            addField
            (
                pair.phase1(),
                "DByAf",
                fvc::interpolate
                (
                    rAUs[pair.phase1().index()]*D
                   /max(alpha12, pair.phase1().residualAlpha())
                ),
                DByAfs
            );
            addField
            (
                pair.phase2(),
                "DByAf",
                fvc::interpolate
                (
                    rAUs[pair.phase2().index()]*D
                   /max(alpha12, pair.phase2().residualAlpha())
                ),
                DByAfs
            );
        }
```

similarly, interpolation is performed in addition to `DByAf`, like

* $$ \left(\frac{rAU^m_1 D}{\max(alpha12^m, residualAlpha()^m_1)}\right)_f $$

**return**

```cpp
    return DByAfs;
```

finally return `DByAfs`



            
### phaseModel

#### Summary

* phaseModel: typedef phase models, details below
* AnisothermalPhaseModel: Class which represents a phase for which the **temperature (strictly energy) varies**. Returns the energy equation and corrects the thermodynamic model.
* InertPhaseModel: Class which represents an **inert phase, with no reactions**. Returns zero reaction rate and heat.
* IsothermalPhaseModel: Class which represents a phase for which the **temperature (strictly energy) remains constant**. Returns an empty energy equation and does nothing when correctThermo is called.
* MovingPhaseModel: Class which represents **a moving fluid phase**. Holds the velocity, fluxes and turbulence model and can generate the momentum equation. The interface is quite restrictive as it also has to support an equivalent stationary model, which does not store motion fields or a turbulence model. Possible future extensions include separating the turbulent functionality into another layer.
* MultiComponentPhaseModel: Class which represents **a phase with multiple species**. Returns the species' mass fractions, and their governing equations.
* PurePhaseModel: Class which represents **pure phases, i.e. without any species**. Returns an empty list of mass fractions.
* ReactingPhaseModel: Class which represents **phases with volumetric reactions**. Returns the reaction rate and heat.
* StationaryPhaseModel: Class which represents **a stationary (and therefore probably solid) phase**. Generates, but does not store, zero velocity and flux field and turbulent quantities. Throws an error when non-const access is requested to the motion fields or when the momentum equation is requested. Usage must must protect against such calls.
* ThermoPhaseModel: Class which represents **a phase with a thermodynamic model**. Provides access to the thermodynamic variables. Note that the thermo model itself is not returned as this class could be substituted in the hierarchy for one which mirrors the functionality, but does not include a thermo model; an incompressible phase model, for example.

#### phaseModel

##### phaseModel.H

##### phaseModel.C

###### constructor

```cpp
Foam::phaseModel::phaseModel
(
    const phaseSystem& fluid,
    const word& phaseName,
    const bool referencePhase,
    const label index
)
:
    // initialize
    // initialize parent class volScalarField, with alpha
    volScalarField
    (
        referencePhase
      ? volScalarField
        (
            IOobject
            (
                IOobject::groupName("alpha", phaseName),
                fluid.mesh().time().timeName(),
                fluid.mesh(),
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            fluid.mesh(),
            dimensionedScalar(dimless, 0)
        )
      : volScalarField
        (
            IOobject
            (
                IOobject::groupName("alpha", phaseName),
                fluid.mesh().time().timeName(),
                fluid.mesh(),
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            fluid.mesh()
        )
    ),

    // initialize private data
    fluid_(fluid),
    name_(phaseName),
    index_(index),
    residualAlpha_
    (
        "residualAlpha",
        dimless,
        fluid.subDict(phaseName).lookup("residualAlpha")
    ),
    alphaMax_(fluid.subDict(phaseName).lookupOrDefault("alphaMax", 1.0))
{
    diameterModel_ = diameterModel::New(fluid.subDict(phaseName), *this);
}
```

initialize parent class and private data.

it should be noted that there is only one definition in class `phaseSystem` owns the type of `volScalarField`, which is alos the parent class.

##### phaseModelNew.C

##### phaseModels.C

###### purePhaseModel

```cpp
    typedef
        AnisothermalPhaseModel
        <
            PurePhaseModel
            <
                InertPhaseModel
                <
                    MovingPhaseModel
                    <
                        ThermoPhaseModel<phaseModel, rhoThermo>
                    >
                >
            >
        >
        purePhaseModel;
```

purePhaseModel:

* temperature changes
* pure specie
* no reaction
* **moving**
* thermodynamic model

###### pureStationaryPhaseModel

```cpp:
    typedef
        AnisothermalPhaseModel
        <
            PurePhaseModel
            <
                InertPhaseModel
                <
                    StationaryPhaseModel
                    <
                        ThermoPhaseModel<phaseModel, rhoThermo>
                    >
                >
            >
        >
        pureStationaryPhaseModel;
```

pureStationaryPhaseModel:

* temperature changes
* pure specie
* no reaction
* **stationary**
* thermodynamic model

###### pureIsothermalPhaseModel

```cpp
    typedef
        IsothermalPhaseModel
        <
            PurePhaseModel
            <
                InertPhaseModel
                <
                    MovingPhaseModel
                    <
                        ThermoPhaseModel<phaseModel, rhoThermo>
                    >
                >
            >
        >
        pureIsothermalPhaseModel;
```

pureIsothermalPhaseModel:

* constant temperature
* pure specie
* no reaction
* **moving**
* thermodynamic model

###### pureStationaryIsothermalPhaseModel

```cpp
    typedef
        IsothermalPhaseModel
        <
            PurePhaseModel
            <
                InertPhaseModel
                <
                    StationaryPhaseModel
                    <
                        ThermoPhaseModel<phaseModel, rhoThermo>
                    >
                >
            >
        >
        pureStationaryIsothermalPhaseModel;
```

pureStationaryIsothermalPhaseModel:

* constant temperature
* pure specie
* no reaction
* **stationary**
* thermodynamic model

###### multiComponentPhaseModel

```cpp
    typedef
        AnisothermalPhaseModel
        <
            MultiComponentPhaseModel
            <
                InertPhaseModel
                <
                    MovingPhaseModel
                    <
                        ThermoPhaseModel<phaseModel, rhoReactionThermo>
                    >
                >
            >
        >
        multiComponentPhaseModel;
```

multiComponentPhaseModel:

* temperature changes
* **multi-component**
* no reaction
* **moving**
* thermodynamic model

###### multiComponentIsothermalPhaseModel

```cpp
    typedef
        IsothermalPhaseModel
        <
            MultiComponentPhaseModel
            <
                InertPhaseModel
                <
                    MovingPhaseModel
                    <
                        ThermoPhaseModel<phaseModel, rhoReactionThermo>
                    >
                >
            >
        >
        multiComponentIsothermalPhaseModel;
```

multiComponentIsothermalPhaseModel:

* constant temperature
* **multi-component**
* no reaction
* **moving**
* thermodynamic model

###### reactingPhaseModel

```cpp
    typedef
        AnisothermalPhaseModel
        <
            MultiComponentPhaseModel
            <
                ReactingPhaseModel
                <
                    MovingPhaseModel
                    <
                        ThermoPhaseModel<phaseModel, rhoReactionThermo>
                    >,
                    CombustionModel<rhoReactionThermo>
                >
            >
        >
        reactingPhaseModel;
```

reactingPhaseModel:

* temperature changes
* **multi-component**
* **reacting**
* **moving**
* thermodynamic model
* **combustion**

#### AnisothermalPhaseModel

##### AnisothermalPhaseModel.C

##### AnisothermalPhaseModel.H

#### InertPhaseModel

##### InertPhaseModel.C

##### InertPhaseModel.H

#### IsothermalPhaseModel

##### IsothermalPhaseModel.C

##### IsothermalPhaseModel.H

#### MovingPhaseModel

##### MovingPhaseModel.C

##### MovingPhaseModel.H

#### MultiComponentPhaseModel

##### MultiComponentPhaseModel.C

##### MultiComponentPhaseModel.H

#### PurePhaseModel

##### PurePhaseModel.C

##### PurePhaseModel.H

#### ReactingPhaseModel

##### ReactingPhaseModel.C

##### ReactingPhaseModel.H

#### StationaryPhaseModel

##### StationaryPhaseModel.C

##### StationaryPhaseModel.H

#### ThermoPhaseModel

##### ThermoPhaseModel.C

##### ThermoPhaseModel.H

## multiPhaseEulerFoam

### createFieldRefs.H

```cpp
surfaceScalarField& phi = fluid.phi();

const IOMRFZoneList& MRF = fluid.MRF();
fv::options& fvOptions = fluid.fvOptions();
```

define `phi`, `MRF` and `fvOptions`

### createFields.H

```cpp
#include "createRDeltaT.H"
#include "readGravitationalAcceleration.H"
#include "readhRef.H"
```

include some required materials

```cpp
Info<< "Creating phaseSystem\n" << endl;

autoPtr<phaseSystem> fluidPtr
(
    phaseSystem::New(mesh)
);
phaseSystem& fluid = fluidPtr(); // define fluid
phaseSystem::phaseModelList& phases = fluid.phases(); // define phase
```

define a pointer for phase as `fluidPtr`, then define `fluid` as `fluidPtr` and define a phase list s`phases` of the phase sytem of `fluid`

```cpp
dimensionedScalar pMin
(
    "pMin",
    dimPressure,
    fluid
);

#include "gh.H"

volScalarField& p = phases[0].thermoRef().p();

Info<< "Reading field p_rgh\n" << endl;
volScalarField p_rgh // using p_rgh
(
    IOobject
    (
        "p_rgh",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

label pRefCell = 0;
scalar pRefValue = 0.0;
if (fluid.incompressible())
{
    p = max(p_rgh + fluid.rho()*gh, pMin);

    if (p_rgh.needReference())
    {
        setRefCell
        (
            p,
            p_rgh,
            pimple.dict(),
            pRefCell,
            pRefValue
        );

        p += dimensionedScalar
        (
            "p",
            p.dimensions(),
            pRefValue - getRefCellValue(p, pRefCell)
        );
        p_rgh = p - fluid.rho()*gh;
    }
}
mesh.setFluxRequired(p_rgh.name());
```

define pMin, g, p, p_rgh, pRefCell, pRefValue etc.

if the fluid is incompressible, then

$$
p = \min(p_{rgh} + \rho gh, pMin)
$$

set refference for p, and 

$$
p_{rgh} = p - \rho gh
$$

```cpp
PtrList<volScalarField> rAUs;
PtrList<surfaceScalarField> rAUfs;
```

Define lists for `rAUs` and `rAUfs`

### CourantNo.H

```cpp
scalar CoNum = 0.0;
scalar meanCoNum = 0.0;
```

initialize `CoNum` and `meanCoNum`

```cpp
if (mesh.nInternalFaces())
{
    scalarField sumPhi
    (
        fvc::surfaceSum(mag(phi))().primitiveField()
    );

    forAll(phases, phasei)
    {
        sumPhi = max
        (
            sumPhi,
            fvc::surfaceSum(mag(phases[phasei].phi()))().primitiveField()
        );
    }

    CoNum = 0.5*gMax(sumPhi/mesh.V().field())*runTime.deltaTValue();

    meanCoNum =
        0.5*(gSum(sumPhi)/gSum(mesh.V().field()))*runTime.deltaTValue();
}
```

$$
CoNum = 0.5 \left(\frac{sunPhi}{V}\right)_{\max} \Delta t
$$

$$
meanCoNum = 0.5 \frac{\sum sumPhi}{\sum V} \Delta t
$$

```cpp
Info<< "Courant Number mean: " << meanCoNum
    << " max: " << CoNum << endl;
```

output mean and max Courant number

### setDeltaT.H

```cpp
volScalarField& rDeltaT = trDeltaT.ref();

const dictionary& pimpleDict = pimple.dict();
```

initialize `rDeltaT` and `pimpleDict`

```cpp
scalar maxCo
(
    pimpleDict.lookupOrDefault<scalar>("maxCo", 0.2)
);

scalar maxDeltaT
(
    pimpleDict.lookupOrDefault<scalar>("maxDeltaT", great)
);

scalar rDeltaTSmoothingCoeff
(
    pimpleDict.lookupOrDefault<scalar>("rDeltaTSmoothingCoeff", 0.02)
);
```

define lookup and initialize `maxCo`, `maxDeltaT` and `rDeltaTSmoothingCoeff`

```cpp
surfaceScalarField maxPhi("maxPhi", phi);

forAll(phases, phasei)
{
    maxPhi = max(maxPhi, mag(phases[phasei].phi()));
}
```

define `maxPhi`

$$
maxPhi = \max(phi, mag(phases[phasei].phi()))
$$

```cpp
// Set the reciprocal time-step from the local Courant number
rDeltaT.ref() = max
(
    1/dimensionedScalar(dimTime, maxDeltaT),
    fvc::surfaceSum(maxPhi)()()
    /((2*maxCo)*mesh.V())
);

// Update the boundary values of the reciprocal time-step
rDeltaT.correctBoundaryConditions();

fvc::smooth(rDeltaT, rDeltaTSmoothingCoeff);

Info<< "Flow time scale min/max = "
    << gMin(1/rDeltaT.primitiveField())
    << ", " << gMax(1/rDeltaT.primitiveField()) << endl;
```

$$
rDeltaT.ref() = \max(1, \frac{\sum maxPhi}{2 \cdot maxCo \cdot V})
$$

$$
rDeltaT = \frac{1}{\Delta t}
$$

correct boundaries of `rDeltaT` and smooth `rDeltaT`

### YEqns.H

This should be the specie transport equations for different phases of the phase system.

```cpp
autoPtr<phaseSystem::specieTransferTable>
    specieTransferPtr(fluid.specieTransfer());

phaseSystem::specieTransferTable&
    specieTransfer(specieTransferPtr());

fluid.correctReactions();
```

define the list of species and rename it as `specieTransfer`, correct reactions

```cpp
forAll(fluid.multiComponentPhases(), multiComponentPhasei)
{
    phaseModel& phase = fluid.multiComponentPhases()[multiComponentPhasei];

    UPtrList<volScalarField>& Y = phase.YActiveRef();
    const volScalarField& alpha = phase;
    const volScalarField& rho = phase.rho();

    forAll(Y, i)
    {
        fvScalarMatrix YiEqn
        (
            phase.YiEqn(Y[i])
            ==
            *specieTransfer[Y[i].name()]
            + fvOptions(alpha, rho, Y[i])
        );

        YiEqn.relax();
        YiEqn.solve("Yi");
    }
}

fluid.correctSpecies();
```

first, get one of the phase of phases with multiple species;

second, obetain the list of species within the phase i

third, get $\alpha$ ($\alpha$ is phase) and $\rho$

fourth, define $Y$ equation for every specie, as

`YiEqn` can be found in `applications\solvers\multiphase\multiphaseEulerFoam\phaseSystems\phaseModel\MultiComponentPhaseModel\MultiComponentPhaseModel.C`

```cpp
template<class BasePhaseModel>
Foam::tmp<Foam::fvScalarMatrix>
Foam::MultiComponentPhaseModel<BasePhaseModel>::YiEqn(volScalarField& Yi)
{
    const volScalarField& alpha = *this;
    const surfaceScalarField alphaRhoPhi(this->alphaRhoPhi());
    const volScalarField& rho = this->thermo().rho();

    return
    (
        fvm::ddt(alpha, rho, Yi)
      + fvm::div(alphaRhoPhi, Yi, "div(" + alphaRhoPhi.name() + ",Yi)")
      + this->divj(Yi)
     ==
        alpha*this->R(Yi)

      + fvc::ddt(residualAlpha_*rho, Yi)
      - fvm::ddt(residualAlpha_*rho, Yi)
    );
}
```

$$
\frac{\partial \alpha \rho Y_i}{\partial t} + \nabla \cdot (\alpha \rho \phi Y_i) - \nabla \cdot (\alpha \alpha_{thermo} \nabla Y_i) = \alpha R(Y_i) + \left(\frac{\partial residualAlpha\_\rho Y_i}{\partial t}\right)_{explicit} - \left(\frac{\partial residualAlpha\_\rho Y_i}{\partial t}\right)_{implicit}
$$

`alphaRhoPhi` is the Mass flux, in other words, $\rho \mathbf{U}$, and $\phi$ is the flux, $\alpha \phi$ is the Volumetric flux

`divj(Yi)` can be found in `src\ThermophysicalTransportModels\laminar\Fourier\Fourier.C`, it's about Fourier's gradient heat flux model for laminar flow.


```cpp
template<class BasicThermophysicalTransportModel>
tmp<fvScalarMatrix>
Fourier<BasicThermophysicalTransportModel>::divj(volScalarField& Yi) const
{
    return -fvm::laplacian(this->alpha()*this->thermo().alpha(), Yi);
}
```

$\alpha$ is the phase fraction, and the $\alpha$ in thermo models is the Laminar thermal diffusivity [kg/m/s], for clarity, it represented by $\alpha_{thermo}$

So,

$$
divj(Y_i) = - \nabla \cdot (\alpha \alpha_{thermo} \nabla Y_i)
$$

`R(Yi)` represents the fuel consumption rate matrix, so it's about reaction.

`residualAlpha_` is the residual phase-fraction for given phase, which is used to stabilize the phase momentum as the phase-fraction -> 0

`fvm::` represents implicit while `fvc::` represents explicit.

The last two terms are added to improve numerical stability of phase momentum as the phase-fraction -> 0. 

fifth, relax and solve $Y$ equations

finally, correct species

### EEqns.H

```cpp
for (int Ecorr=0; Ecorr<nEnergyCorrectors; Ecorr++)
{
    fluid.correctEnergyTransport();

    autoPtr<phaseSystem::heatTransferTable>
        heatTransferPtr(fluid.heatTransfer());

    phaseSystem::heatTransferTable& heatTransfer = heatTransferPtr();

    forAll(fluid.anisothermalPhases(), anisothermalPhasei)
    {
        phaseModel& phase = fluid.anisothermalPhases()[anisothermalPhasei];

        const volScalarField& alpha = phase;
        tmp<volScalarField> tRho = phase.rho();
        const volScalarField& rho = tRho();
        tmp<volVectorField> tU = phase.U();
        const volVectorField& U = tU();

        fvScalarMatrix EEqn
        (
            phase.heEqn()
         ==
           *heatTransfer[phase.name()]
          + alpha*rho*(U&g)
          + fvOptions(alpha, rho, phase.thermoRef().he())
        );

        EEqn.relax();
        fvOptions.constrain(EEqn);
        EEqn.solve();
        fvOptions.correct(phase.thermoRef().he());
    }

    fluid.correctThermo();
    fluid.correctContinuityError();
}
```

first, it's a loop to corrent energy enquations.

In this loop:

* correct energy transport equations
* define `heatTransfer` as the list of heat transfer matrices
* for every anisothermal phases:
  * define `phase` as current anisothermal phase
  * get $\alpha$, $\rho$ and $\mathbf{U}$
  * define energy equation `EEqn` as:
  * relax, constrain and solve `EEqn`
  * correct `he`
* correct `thermo` and `continuityError`

Output the minimum and maximum temperature for every phases


$$
phase.heEqn() = *heatTransfer[phase.name()] + \alpha \rho \mathbf{U} \cdot \mathbf{g} + fvOptions
$$

`phase.heEqn()` can be found in `applications\solvers\multiphase\multiphaseEulerFoam\phaseSystems\phaseModel\AnisothermalPhaseModel\AnisothermalPhaseModel.C`

```cpp
template<class BasePhaseModel>
Foam::tmp<Foam::fvScalarMatrix>
Foam::AnisothermalPhaseModel<BasePhaseModel>::heEqn()
{
    const volScalarField& alpha = *this;

    const volVectorField U(this->U());
    const surfaceScalarField alphaPhi(this->alphaPhi());
    const surfaceScalarField alphaRhoPhi(this->alphaRhoPhi());

    const volScalarField contErr(this->continuityError());
    const volScalarField K(this->K());

    volScalarField& he = this->thermo_->he();

    tmp<fvScalarMatrix> tEEqn
    (
        fvm::ddt(alpha, this->rho(), he)
      + fvm::div(alphaRhoPhi, he)
      - fvm::Sp(contErr, he)

      + fvc::ddt(alpha, this->rho(), K) + fvc::div(alphaRhoPhi, K)
      - contErr*K
      + this->divq(he)
     ==
        alpha*this->Qdot()
    );

    // Add the appropriate pressure-work term
    if (he.name() == this->thermo_->phasePropertyName("e"))
    {
        tEEqn.ref() += filterPressureWork
        (
            fvc::div(fvc::absolute(alphaPhi, alpha, U), this->thermo().p())
          + (fvc::ddt(alpha) - contErr/this->rho())*this->thermo().p()
        );
    }
    else if (this->thermo_->dpdt())
    {
        tEEqn.ref() -= filterPressureWork(alpha*this->fluid().dpdt());
    }

    return tEEqn;
}
```

* obtain $\alpha$, $\mathbf{U}$, volumetric flux $\alpha \phi$, mass flux $\alpha \rho \phi$, the phase kinetic energy $K$, the continuity error $contErr$, and enthalpy or internal energy $he$
* define `tEEqn` as

$$
\frac{\partial \alpha \rho he}{\partial t} + \nabla \cdot (\alpha \rho \phi he) - Sp(contErr, he) + \frac{\partial \alpha \rho K}{\partial t} + \nabla \cdot (\alpha \rho \phi K) - contErr \cdot K - \nabla \cdot (\alpha \alpha_{thermo} \nabla he) = \dot{Q}
$$

`divq` is generally defined as:

$$
divq = \nabla \cdot (\alpha \alpha_{thermo} \nabla he) 
$$

`Qdot()` or $\dot{Q}$ is the heat release rate

* according to energy equation or ethalpy equation, adding required terms
* for energy equation:

$$
tEEqn.ref()  + \nabla \cdot ((\alpha \phi + \alpha \phi) p) + (\frac{\partial \alpha}{\partial t} - \frac{contErr}{\rho})p 
$$

`filterPressureWork` is to optionally filter the pressure work term as the phase-fraction -> 0

```cpp
template<class BasePhaseModel>
Foam::tmp<Foam::volScalarField>
Foam::AnisothermalPhaseModel<BasePhaseModel>::filterPressureWork
(
    const tmp<volScalarField>& pressureWork
) const
{
    const volScalarField& alpha = *this;

    scalar pressureWorkAlphaLimit =
        this->thermo_->lookupOrDefault("pressureWorkAlphaLimit", 0.0);

    if (pressureWorkAlphaLimit > 0)
    {
        return
        (
            max(alpha - pressureWorkAlphaLimit, scalar(0))
           /max(alpha - pressureWorkAlphaLimit, pressureWorkAlphaLimit)
        )*pressureWork;
    }
    else
    {
        return pressureWork;
    }
}
```

* obtain $\alpha$
* define and look up for `pressureWorkAlphaLimit`, whose default value is 0
* if `pressureWorkAlphaLimit` > 0 return

$$
\frac{\max(\alpha-pressureWorkAlphaLimit, 0)}{\max(\alpha - pressureWorkAlphaLimit, pressureWorkAlphaLimit)} \cdot pressureWork
$$

When $\alpha \rArr 0$, it returns 0.

* else, return `pressureWork`

`absolute()` is to return the given relative flux in absolute form as in `src\finiteVolume\finiteVolume\fvc\fvcMeshPhi.C`

```cpp
Foam::tmp<Foam::surfaceScalarField> Foam::fvc::absolute
(
    const tmp<surfaceScalarField>& tphi,
    const volScalarField& rho,
    const volVectorField& U
)
{
    if (tphi().mesh().moving())
    {
        return tphi + fvc::interpolate(rho)*fvc::meshPhi(rho, U);
    }
    else
    {
        return tmp<surfaceScalarField>(tphi, true);
    }
}
```
 if moving, then return

 $$
tphi + \rho phi
 $$

 so

 $$
fvc::div(fvc::absolute(alphaPhi, alpha, U), this->thermo().p()) = \nabla \cdot ((\alpha \phi + \alpha \phi) p)
 $$

 * else for ethalpy equation

$$
tEEqn.ref() - \alpha \frac{\partial p}{\partial t}
$$
 
So, for energy equation:

$$
\frac{\partial \rho e}{\partial t}+\nabla \cdot (\rho \mathbf{U} e) + \frac{\partial \rho K}{\partial t}+\nabla \cdot (\rho \mathbf{U} K)= -\nabla\cdot(p\mathbf{U})+\rho r -\nabla\cdot\mathbf{q} + \rho \mathbf{g} \cdot \mathbf{U}+\nabla \cdot(\tau \cdot \mathbf{U})
$$

$$
\frac{\partial \rho e}{\partial t}+\nabla \cdot (\rho \mathbf{U} e) + \frac{\partial \rho K}{\partial t}+\nabla \cdot (\rho \mathbf{U} K)- \nabla \cdot (\alpha_\mathrm{eff}\nabla e)= -\nabla\cdot(p\mathbf{U}) 
$$

$$
\frac{\partial \alpha \rho he}{\partial t} + \nabla \cdot (\alpha \rho \phi he) - Sp(contErr, he) + \frac{\partial \alpha \rho K}{\partial t} + \nabla \cdot (\alpha \rho \phi K) - contErr \cdot K - \nabla \cdot (\alpha \alpha_{thermo} \nabla he) + \nabla \cdot ((\alpha \phi + \alpha \phi) p) + (\frac{\partial \alpha}{\partial t} - \frac{contErr}{\rho})p = \dot{Q} 
$$

for ethalpy equation:

$$
\frac{\partial \rho h}{\partial t}+\nabla \cdot (\rho \mathbf{U} h) + \frac{\partial \rho K}{\partial t}+\nabla \cdot (\rho \mathbf{U} K) =\frac{\partial p}{\partial t}+ \rho r -\nabla\cdot\mathbf{q} + \rho \mathbf{g} \cdot \mathbf{U}+\nabla \cdot(\tau \cdot \mathbf{U})
$$

$$
\frac{\partial \rho h}{\partial t}+\nabla \cdot (\rho \mathbf{U} h) + \frac{\partial \rho K}{\partial t}+\nabla \cdot (\rho \mathbf{U} K) - \nabla \cdot (\alpha_\mathrm{eff}\nabla h) =\frac{\partial p}{\partial t}
$$

$$
\frac{\partial \alpha \rho he}{\partial t} + \nabla \cdot (\alpha \rho \phi he) - Sp(contErr, he) + \frac{\partial \alpha \rho K}{\partial t} + \nabla \cdot (\alpha \rho \phi K) - contErr \cdot K - \nabla \cdot (\alpha \alpha_{thermo} \nabla he) - \alpha \frac{\partial p}{\partial t} = \dot{Q}
$$

also the same

finally, the `EEqn` is:

$$
\frac{\partial \alpha \rho he}{\partial t} + \nabla \cdot (\alpha \rho \phi he) - Sp(contErr, he) + \frac{\partial \alpha \rho K}{\partial t} + \nabla \cdot (\alpha \rho \phi K) - contErr \cdot K - \nabla \cdot (\alpha \alpha_{thermo} \nabla he) + \nabla \cdot ((\alpha \phi + \alpha \phi) p) + (\frac{\partial \alpha}{\partial t} - \frac{contErr}{\rho})p = *heatTransfer[phase.name()] + \alpha \rho \mathbf{U} \cdot \mathbf{g} + fvOptions
$$

$$
\frac{\partial \alpha \rho he}{\partial t} + \nabla \cdot (\alpha \rho \phi he) - Sp(contErr, he) + \frac{\partial \alpha \rho K}{\partial t} + \nabla \cdot (\alpha \rho \phi K) - contErr \cdot K - \nabla \cdot (\alpha \alpha_{thermo} \nabla he) - \alpha \frac{\partial p}{\partial t} = *heatTransfer[phase.name()] + \alpha \rho \mathbf{U} \cdot \mathbf{g} + fvOptions
$$

### PU

#### UEqns.H

```cpp
Info<< "Constructing momentum equations" << endl;

PtrList<fvVectorMatrix> UEqns(phases.size());

{
    autoPtr<phaseSystem::momentumTransferTable>
        momentumTransferPtr(fluid.momentumTransfer());

    phaseSystem::momentumTransferTable&
        momentumTransfer(momentumTransferPtr());

    forAll(fluid.movingPhases(), movingPhasei)
    {
        phaseModel& phase = fluid.movingPhases()[movingPhasei];

        const volScalarField& alpha = phase;
        const volScalarField& rho = phase.rho();
        volVectorField& U = phase.URef();

        UEqns.set
        (
            phase.index(),
            new fvVectorMatrix
            (
                phase.UEqn()
             ==
               *momentumTransfer[phase.name()]
              + fvOptions(alpha, rho, U)
            )
        );

        UEqns[phase.index()].relax();
        fvOptions.constrain(UEqns[phase.index()]);
        fvOptions.correct(U);
    }
}
```

* output information about constructing momentum equation
* define a list for velocity equations with the size of phase number 
* define `momentumTransfer` as a list for terms in momentrum equations
* then start a loop for every moving phase, in which velocity equation for every moving phase is defined 
  * define $\alpha$, $\rho$, $\mathbf{U}$
  * define `UEqn` in the list `UEqns` as
  * relax `UEqn`
  * constrain and correct fvOpentions

`phase.UEqn()` is defined in `applications\solvers\multiphase\multiphaseEulerFoam\phaseSystems\phaseModel\MovingPhaseModel\MovingPhaseModel.C`

```cpp
template<class BasePhaseModel>
Foam::tmp<Foam::fvVectorMatrix>
Foam::MovingPhaseModel<BasePhaseModel>::UEqn()
{
    const volScalarField& alpha = *this;
    const volScalarField& rho = this->thermo().rho();

    return
    (
        fvm::ddt(alpha, rho, U_)
      + fvm::div(alphaRhoPhi_, U_)
      + fvm::SuSp(-this->continuityError(), U_)
      + this->fluid().MRF().DDt(alpha*rho, U_)
      + turbulence_->divDevTau(U_)
    );
}
```

`turbulence_->divDevTau(U_)` can be found in `D:\Documents\Git\OpenFOAM-8\src\MomentumTransportModels\momentumTransportModels\linearViscousStress\linearViscousStress.C`

```cpp
template<class BasicMomentumTransportModel>
Foam::tmp<Foam::fvVectorMatrix>
Foam::linearViscousStress<BasicMomentumTransportModel>::divDevTau
(
    volVectorField& U
) const
{
    return
    (
      - fvc::div((this->alpha_*this->rho_*this->nuEff())*dev2(T(fvc::grad(U))))
      - fvm::laplacian(this->alpha_*this->rho_*this->nuEff(), U)
    );
}
```

`dev2` can be found in `src\OpenFOAM\primitives\Tensor\TensorI.H`

```cpp
//- Return the deviatoric part of a tensor
template<class Cmpt>
inline Tensor<Cmpt> dev2(const Tensor<Cmpt>& t)
{
    return t - SphericalTensor<Cmpt>::twoThirdsI*tr(t);
}
```

$$
dev2(\mathbf{T}) = \mathbf{T} - \frac{2}{3} tr(\mathbf{T}) \mathbf{I}
$$

`T()` can be found in `D:\Documents\Git\OpenFOAM-8\src\OpenFOAM\primitives\Tensor\Tensor.H` or in `src\OpenFOAM\primitives\SphericalTensor\SphericalTensor.H`, is to return the transpose

$$
T(\mathbf{T}) = \mathbf{T}^T
$$

So

$$
T(\nabla \mathbf{U}) = (\nabla \mathbf{U})^T
$$

$$
dev2((\nabla \mathbf{U})^T) = (\nabla \mathbf{U})^T - \frac{2}{3} (\nabla \cdot \mathbf{U}) \mathbf{I}
$$

So `divDevTau` is

$$
divDevTau = -\nabla \cdot (\alpha \rho \nu_{Eff} ((\nabla \mathbf{U})^T - \frac{2}{3} (\nabla \cdot \mathbf{U}) \mathbf{I})) - \nabla \cdot (\alpha \rho \nu_{Eff} \nabla \mathbf{U}) \\= -\nabla \cdot \left[\alpha \rho \nu_{Eff} \left((\nabla \mathbf{U}+(\nabla \mathbf{U})^T) - \frac{2}{3} (\nabla \cdot \mathbf{U}) \mathbf{I}\right)\right]
$$

It should be noted that 

$$
\nabla \cdot \tau = \nabla \cdot (-\frac{2}{3} \mu (\nabla \cdot \mathbf{U})\mathbf{I} + 2 \mu (\frac{1}{2}(\nabla \mathbf{U} + (\nabla \mathbf{U})^T)))
$$

which is the same with `divDevTau`

The `UEqn` is defined as

$$
\frac{\partial \alpha^k \rho^k \mathbf{U}^k}{\partial t} + \nabla \cdot (\alpha^k \rho^k \phi \mathbf{U}^k) + SuSp(continuityError_, \mathbf{U}^k) + MRF(\alpha^k \rho^k \mathbf{U}^k) -\nabla \cdot \left[\alpha^k \rho^k \nu_{Eff}^k \left((\nabla \mathbf{U}^k+(\nabla \mathbf{U}^k)^T) - \frac{2}{3} (\nabla \cdot \mathbf{U}^k) \mathbf{I}\right)\right]
$$

or

$$
\frac{\partial \alpha^k \rho^k \mathbf{U}^k}{\partial t} + \nabla \cdot (\alpha^k \rho^k \phi \mathbf{U}^k) + SuSp(continuityError_, \mathbf{U}^k) + MRF(\alpha^k \rho^k \mathbf{U}^k) -\nabla \cdot \tau^k
$$

So

$$
\frac{\partial \alpha^k \rho^k \mathbf{U}^k}{\partial t} + \nabla \cdot (\alpha^k \rho^k \phi \mathbf{U}^k) + SuSp(continuityError_, \mathbf{U}^k) + MRF(\alpha^k \rho^k \mathbf{U}^k) -\nabla \cdot \tau^k = momentumTransfer + fvOptions(\alpha, \rho, \mathbf{U})
$$

#### pEqn.H

##### Face volume fractions

```cpp
// Face volume fractions
PtrList<surfaceScalarField> alphafs(phases.size());
forAll(phases, phasei)
{
    phaseModel& phase = phases[phasei];
    const volScalarField& alpha = phase;

    alphafs.set(phasei, fvc::interpolate(alpha).ptr());
    alphafs[phasei].rename("pEqn" + alphafs[phasei].name());
}
```

define pointer list of $\alpha_f$, $\alpha$ at surface, for all phases as `alphafs`

##### Diagnocal coefficients

```cpp
// Diagonal coefficients
rAUs.clear();
rAUs.setSize(phases.size());

forAll(fluid.movingPhases(), movingPhasei)
{
    phaseModel& phase = fluid.movingPhases()[movingPhasei];
    const volScalarField& alpha = phase;

    rAUs.set
    (
        phase.index(),
        new volScalarField
        (
            IOobject::groupName("rAU", phase.name()),
            1.0
           /(
               UEqns[phase.index()].A()
             + byDt(max(phase.residualAlpha() - alpha, scalar(0))*phase.rho())
            )
        )
    );
}
fluid.fillFields("rAU", dimTime/dimDensity, rAUs);
```

* clear `rAUs`, it's a pointer list of `volScalarField` variables defined in `createFields.H`
* set the size of `rAUs` as the number of phases
* for every moving phases:
  * define and get $\alpha$
  * then, set `rAUs`
  * fill up gaps in a phase-indexed table of fields with zeros

`byDt` is defined in `D:\Documents\Git\OpenFOAM-8\applications\solvers\multiphase\multiphaseEulerFoam\phaseSystems\phaseSystem\phaseSystem.C` for `volScalarField` and `surfaceScalarField`

```cpp
Foam::tmp<Foam::volScalarField> Foam::byDt(const volScalarField& vf)
{
    if (fv::localEulerDdt::enabled(vf.mesh()))
    {
        return fv::localEulerDdt::localRDeltaT(vf.mesh())*vf;
    }
    else
    {
        return vf/vf.mesh().time().deltaT();
    }
}


Foam::tmp<Foam::surfaceScalarField> Foam::byDt(const surfaceScalarField& sf)
{
    if (fv::localEulerDdt::enabled(sf.mesh()))
    {
        return fv::localEulerDdt::localRDeltaTf(sf.mesh())*sf;
    }
    else
    {
        return sf/sf.mesh().time().deltaT();
    }
}
```

`UEqns[phase.index()].A()` is to return the central coefficient, namely $A_P$

$$
byDt = \frac{v_f}{\Delta t}
$$

$$
rAU = \frac{1}{A_P + \frac{\max(residualAlpha() - \alpha, 0)\rho}{\Delta t}}
$$

in general $\alpha$ is not too small, then $rAU = \frac{1}{A_P}$

##### Phase diagonal coefficients

```cpp
// Phase diagonal coefficients
PtrList<surfaceScalarField> alpharAUfs(phases.size());
forAll(phases, phasei)
{
    phaseModel& phase = phases[phasei];
    const volScalarField& alpha = phase;

    alpharAUfs.set
    (
        phasei,
        (
            fvc::interpolate(max(alpha, phase.residualAlpha())*rAUs[phasei])
        ).ptr()
    );
}
```

* define a pointer list of `surfaceScalarField` variables `alpharAUfs` with size of phase number
* for every phases:
  * define and get $\alpha$
  * set `alpharAUfs`

$$
alpharAUfs[i] = \left(\max(\alpha, residualAlpha()) \cdot rAUs[i] \right)_f
$$

In general, $alpharAUf = \left(\frac{\alpha}{A_P}\right)_f$

##### Explicit force fluxes

```cpp
// Explicit force fluxes
PtrList<surfaceScalarField> phiFs(fluid.phiFs(rAUs));
```

`phiFs()` can be found in `applications\solvers\multiphase\multiphaseEulerFoam\phaseSystems\PhaseSystems\MomentumTransferPhaseSystem\MomentumTransferPhaseSystem.C`, is to return the explicit force fluxes for the cell-based algorithm, that do not depend on phase mass/volume fluxes, and can therefore be evaluated outside the corrector loop. This includes things like lift, turbulent dispersion, and wall lubrication.

##### Pressure corrector loop

```cpp
// --- Pressure corrector loop
while (pimple.correct())
{
    volScalarField rho("rho", fluid.rho());

    // Correct p_rgh for consistency with p and the updated densities
    p_rgh = p - rho*gh;

    ...

}
```

this is the pressure corrector loop, obtain $\rho$ and correct $p_{rgh}$ with updated $\rho$

###### Correct fixed-flux BCs to be consistent with the velocity BCs

```cpp
    // Correct fixed-flux BCs to be consistent with the velocity BCs
    forAll(fluid.movingPhases(), movingPhasei)
    {
        phaseModel& phase = fluid.movingPhases()[movingPhasei];
        MRF.correctBoundaryFlux(phase.U(), phase.phiRef());
    }
```

* for every moving phase:
  * get phase
  * correct boundary flux

###### Combined buoyancy and force fluxes

```cpp
    // Combined buoyancy and force fluxes
    PtrList<surfaceScalarField> phigFs(phases.size());
    {
        const surfaceScalarField ghSnGradRho
        (
            "ghSnGradRho",
            ghf*fvc::snGrad(rho)*mesh.magSf()
        );

        forAll(phases, phasei)
        {
            phaseModel& phase = phases[phasei];

            phigFs.set
            (
                phasei,
                (
                    alpharAUfs[phasei]
                   *(
                       ghSnGradRho
                     - (fvc::interpolate(phase.rho() - rho))*(g & mesh.Sf())
                     - fluid.surfaceTension(phase)*mesh.magSf()
                    )
                ).ptr()
            );

            if (phiFs.set(phasei))
            {
                phigFs[phasei] += phiFs[phasei];
            }
        }
    }
```

* define pointer list of `surfaceScalarField` variables `phigFs` with size of phase number
* define `ghSnGradRho` as:

$$
ghSnGradRho = (\mathbf{g}h)_f \cdot [\mathbf{n} \cdot (\nabla \rho)_f] \|\mathbf{S}_f\|
$$

* for every phase:
  * get phase or $\alpha$
  * set `phigFs`:

$$
phigFs[i] = phigF = \left(\frac{\alpha}{A_P}\right)_f \left[ (\mathbf{g}h)_f \cdot [\mathbf{n} \cdot (\nabla \rho)_f] \|\mathbf{S}_f\| - ((\rho_k - \rho)_f \mathbf{g} \cdot \mathbf{S}_f) - (F_{st}\|\mathbf{S}_f\|) \right]
$$

where $\rho_k$ is the density of phase $k$, $F_{st}$ is the surface tension

  * if `phiFs` is set, then

$$
phigFs[i] = phigFs[i] + phiFs[i]
$$

###### Predicted velocities and fluxes for each phase

```cpp
    // Predicted velocities and fluxes for each phase
    PtrList<volVectorField> HbyAs(phases.size());
    PtrList<surfaceScalarField> phiHbyAs(phases.size());
```

Define pointer list of `volVectorField` $HbyAs$ and of `surfaceScalarField` $phiHbyAs$

```cpp
    {
        // Correction force fluxes
        PtrList<surfaceScalarField> ddtCorrByAs(fluid.ddtCorrByAs(rAUs));

        forAll(fluid.movingPhases(), movingPhasei)
        {
            phaseModel& phase = fluid.movingPhases()[movingPhasei];
            const volScalarField& alpha = phase;

            HbyAs.set
            (
                phase.index(),
                new volVectorField
                (
                    IOobject::groupName("HbyA", phase.name()),
                    phase.U()
                )
            );

            HbyAs[phase.index()] =
                rAUs[phase.index()]
               *(
                    UEqns[phase.index()].H()
                  + byDt
                    (
                        max(phase.residualAlpha() - alpha, scalar(0))
                       *phase.rho()
                    )
                   *phase.U()().oldTime()
                );

            phiHbyAs.set
            (
                phase.index(),
                new surfaceScalarField
                (
                    IOobject::groupName("phiHbyA", phase.name()),
                    fvc::flux(HbyAs[phase.index()])
                  - phigFs[phase.index()]
                  - ddtCorrByAs[phase.index()]
                )
            );
        }
    }
    fluid.fillFields("HbyA", dimVelocity, HbyAs);
    fluid.fillFields("phiHbyA", dimForce/dimDensity/dimVelocity, phiHbyAs);
```

* correct force fluxes
* for every moving phase:
  * define and get $\alpha$
  * set $HbyAs$:

$$
HbyA = \mathbf{U}
$$

  * further initialize $HbyAs$:

$$
HbyA = \frac{1}{rAU} \left( (-\sum A_N \mathbf{U}_N^r + S_P^n) + \frac{\max(residualAlpha - \alpha, 0) \rho}{\Delta t} \mathbf{U}^{n-1}\right)
$$

in general, $HbyA^n = \frac{1}{rAU} \left( -\sum A_N \mathbf{U}_N^r + S_P^n \right) = \frac{1}{A_P} \left( -\sum A_N \mathbf{U}_N^r + S_P^n \right)$

  * set $phiHbyAs$:

$$
phiHbyA = \int \mathbf{S}_f \cdot HbyA - phigF - ddtCorrByA
$$

* fill empty lists so that null pointer will not be returned

###### Add explicit drag forces and fluxes

```cpp
    // Add explicit drag forces and fluxes
    PtrList<volVectorField> KdUByAs(fluid.KdUByAs(rAUs));
    PtrList<surfaceScalarField> phiKdPhis(fluid.phiKdPhis(rAUs));

    forAll(phases, phasei)
    {
        if (KdUByAs.set(phasei))
        {
            HbyAs[phasei] -= KdUByAs[phasei];
        }

        if (phiKdPhis.set(phasei))
        {
            phiHbyAs[phasei] -= phiKdPhis[phasei];
        }
    }
```

* define and get explicit drag forces and fluxes `KdUByAs` and `phiKdPhis` 
* add them to `HbyAs` and `phiHbyAs`

$$
HbyA = HbyA - KdUByA
$$

$$
phiHbyA = phiHbyA - phiKdUByA
$$

###### Total predicted flux

```cpp
    // Total predicted flux
    surfaceScalarField phiHbyA
    (
        IOobject
        (
            "phiHbyA",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar(dimArea*dimVelocity, 0)
    );

    forAll(phases, phasei)
    {
        phiHbyA += alphafs[phasei]*phiHbyAs[phasei];
    }

    MRF.makeRelative(phiHbyA);
```

* define total phiHbyA, namely phiHbyA for all phases, obtained by average:

$$
phiHbyA = \sum_{k=1}^N \alpha^k phiHbyA^k
$$

where $^k$ means $k$th phase

###### Construct pressure "diffusivity"

```cpp
    // Construct pressure "diffusivity"
    surfaceScalarField rAUf
    (
        IOobject
        (
            "rAUf",
            runTime.timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar(dimensionSet(-1, 3, 1, 0, 0), 0)
    );

    forAll(phases, phasei)
    {
        rAUf += alphafs[phasei]*alpharAUfs[phasei];
    }

    rAUf = mag(rAUf);
```

* define total $rAU$ at surface, $rAUf$

$$
rAUf = \sum_{k=1}^N \alpha_f^k (\frac{\alpha}{A_p})^k_f
$$

$$
rAUf = \|\sum_{k=1}^N \alpha_f^k (\frac{\alpha}{A_p})^k_f\|
$$

###### Update the fixedFluxPressure BCs to ensure flux consistency

```cpp
    // Update the fixedFluxPressure BCs to ensure flux consistency
    {
        surfaceScalarField::Boundary phib(phi.boundaryField());
        phib = 0;
        forAll(phases, phasei)
        {
            phaseModel& phase = phases[phasei];
            phib +=
                alphafs[phasei].boundaryField()*phase.phi()().boundaryField();
        }

        setSnGrad<fixedFluxPressureFvPatchScalarField>
        (
            p_rgh.boundaryFieldRef(),
            (
                phiHbyA.boundaryField() - phib
            )/(mesh.magSf().boundaryField()*rAUf.boundaryField())
        );
    }
```

* define $phib$
* for every phase:
  * $$phib = \sum_{k=1}^N \alpha_f^k \phi^k$$
* set boundary of $p_{rgh}$
  * $$p_{rgh} = \frac{phiHbyA - phib}{\|\mathbf{S}_f\| rAUf} = \frac{\sum_{k=1}^N \alpha^k phiHbyA^k - \sum_{k=1}^N \alpha_f^k \phi^k}{\|\mathbf{S}_f\| \|\sum_{k=1}^N \alpha_f^k (\frac{\alpha}{A_p})^k_f\|}  $$

###### Compressible pressure equations

```cpp
    // Compressible pressure equations
    PtrList<fvScalarMatrix> pEqnComps(phases.size());
    PtrList<volScalarField> dmdts(fluid.dmdts());
```

* define pointer list of compressible pressure equations `pEqnComps` and `dmdts`
  
```cpp
    forAll(phases, phasei)
    {
        phaseModel& phase = phases[phasei];
        const volScalarField& alpha = phase;
        volScalarField& rho = phase.thermoRef().rho();

        pEqnComps.set(phasei, new fvScalarMatrix(p_rgh, dimVolume/dimTime));
        fvScalarMatrix& pEqnComp = pEqnComps[phasei];

        ...

    }
```

* for every phase:
  * define and get $\alpha$ and $\rho$
  * initialize `pEqnComps` and define current `pEqnComp`
    * $$pEqnComp = p_{rgh}$$

Density variation

```cpp
        // Density variation
        if (!phase.isochoric() || !phase.pure())
        {
            pEqnComp +=
                (
                    fvc::ddt(alpha, rho) + fvc::div(phase.alphaRhoPhi())
                  - fvc::Sp(fvc::ddt(alpha) + fvc::div(phase.alphaPhi()), rho)
                )/rho;
        }
```

`isochoric()` can be found in `applications\solvers\multiphase\multiphaseEulerFoam\phaseSystems\phaseModel\ThermoPhaseModel\ThermoPhaseModel.H`, is to return whether the phase is constant density

`pure` can be found in `applications\solvers\multiphase\multiphaseEulerFoam\phaseSystems\phaseModel\MultiComponentPhaseModel\MultiComponentPhaseModel.H`, is to return whether the phase is pure (i.e., not multi-component)

  * if current phase is not constant density or not pure, then:
    * $$pEqnComp = pEqnComp + \left[\frac{\partial \alpha \rho}{\partial t} + \nabla (\alpha \rho \phi) - Sp(\frac{\partial \alpha}{\partial t} + \nabla (\alpha \rho), \rho)\right] /\rho \\ = p_{rgh} + \left[\frac{\partial \alpha \rho}{\partial t} + \nabla (\alpha \rho \phi) - Sp(\frac{\partial \alpha}{\partial t} + \nabla (\alpha \rho), \rho)\right] /\rho $$

Compressibility

```cpp
        // Compressibility
        if (!phase.incompressible())
        {
            if (pimple.transonic())
            {
                const surfaceScalarField phid
                (
                    IOobject::groupName("phid", phase.name()),
                    fvc::interpolate(phase.thermo().psi())*phase.phi()
                );

                pEqnComp +=
                    correction
                    (
                        (alpha/rho)*
                        (
                            phase.thermo().psi()*fvm::ddt(p_rgh)
                          + fvm::div(phid, p_rgh)
                          - fvm::Sp(fvc::div(phid), p_rgh)
                        )
                    );

                pEqnComps[phasei].relax();
            }
            else
            {
                pEqnComp +=
                    (alpha*phase.thermo().psi()/rho)
                   *correction(fvm::ddt(p_rgh));
            }
        }
```

  * if current phase in not incompressible, in other words, it's compressible, then
    * if current phase is transonic, then
      * define $phid$ as
        * $$phid = \psi_f \phi$$
      * $$pEqnComp = pEqnComp + correction(\frac{\alpha}{\rho} (\psi \frac{\partial p_{rgh}}{\partial t} + \nabla \cdot (phid p_{rgh}) - Sp(\nabla phid, p_{rgh}) ) ) \\ = p_{rgh} + \left[\frac{\partial \alpha \rho}{\partial t} + \nabla (\alpha \rho \phi) - Sp(\frac{\partial \alpha}{\partial t} + \nabla (\alpha \rho), \rho)\right] /\rho + correction(\frac{\alpha}{\rho} (\psi \frac{\partial p_{rgh}}{\partial t} + \nabla \cdot (\psi_f \phi p_{rgh}) - Sp(\nabla \psi_f \phi, p_{rgh}) ) ) $$ 
    * if current phase in not transonic, then
      * $$pEqnComp = pEqnComp + \frac{\alpha \psi}{\rho} correction(\frac{\partial p_{rgh}}{\partial t}) \\ = p_{rgh} + \left[\frac{\partial \alpha \rho}{\partial t} + \nabla (\alpha \rho \phi) - Sp(\frac{\partial \alpha}{\partial t} + \nabla (\alpha \rho), \rho)\right] /\rho + \frac{\alpha \psi}{\rho} correction(\frac{\partial p_{rgh}}{\partial t})$$

Option sources and Mass transfer

```cpp
        // Option sources
        if (fvOptions.appliesToField(rho.name()))
        {
            pEqnComp -= (fvOptions(alpha, rho) & rho)/rho;
        }

        // Mass transfer
        if (dmdts.set(phasei))
        {
            pEqnComp -= dmdts[phasei]/rho;
        }
```

  * if there is Option sources, then
    * $$pEqnComp = pEqnComp - [fvOptions(\alpha, \rho) \cdot \rho] / \rho$$
  * if ther is Mass transfer, then
    * $$pEqnComp = pEqnComp - \frac{dmdt}{\rho}$$

`dmdts` can be found in `applications\solvers\multiphase\multiphaseEulerFoam\phaseSystems\phaseSystem\phaseSystem.H`, is to return the mass transfer rates for each phase

where $dmdt$ is the masss transfer for current phase

###### Cache p prior to solve for density update

```cpp
    // Cache p prior to solve for density update
    volScalarField p_rgh_0(p_rgh);
```

  * cache $p_{rgh}$ in $p_{rgh0}$

###### Iterate over the pressure equation to correct for non-orthogonality
    
```cpp
    // Iterate over the pressure equation to correct for non-orthogonality
    while (pimple.correctNonOrthogonal())
    {
        ...
    }
```

start the loop to solve pressure equation

Construct the transport part of the pressure equation
        
```cpp
        // Construct the transport part of the pressure equation
        fvScalarMatrix pEqnIncomp
        (
            fvc::div(phiHbyA)
          - fvm::laplacian(rAUf, p_rgh)
        );
```

* define $pEqnIncomp$
  * $$pEqnIncomp = \nabla phiHbyA - \nabla \cdot (rAUf \nabla p_{rgh})$$

Solve

```cpp
        // Solve
        {
            fvScalarMatrix pEqn(pEqnIncomp);

            forAll(phases, phasei)
            {
                pEqn += pEqnComps[phasei];
            }

            if (fluid.incompressible())
            {
                pEqn.setReference(pRefCell, getRefCellValue(p_rgh, pRefCell));
            }

            pEqn.solve();
        }
```

* define `pEqn` equal to `pEqnIncomp`
* for every phase:
  * pEqn = pEqn + pEqnComp
  * namely, $$pEqn = pEqn + \sum_{k=1}^N pEqnComp^k$$
* if it's incompressible, then
  * set reference for `pEqn`
* solve `pEqn` 

Correct fluxes and velocities on last non-orthogonal iteration

```cpp
        // Correct fluxes and velocities on last non-orthogonal iteration
        if (pimple.finalNonOrthogonalIter())
        {
            phi = phiHbyA + pEqnIncomp.flux();

            surfaceScalarField mSfGradp("mSfGradp", pEqnIncomp.flux()/rAUf);

            forAll(fluid.movingPhases(), movingPhasei)
            {
                phaseModel& phase = fluid.movingPhases()[movingPhasei];

                phase.phiRef() =
                    phiHbyAs[phase.index()]
                  + alpharAUfs[phase.index()]*mSfGradp;

                // Set the phase dilatation rate
                phase.divU(-pEqnComps[phase.index()] & p_rgh);
            }

            // Optionally relax pressure for velocity correction
            p_rgh.relax();

            mSfGradp = pEqnIncomp.flux()/rAUf;

            forAll(fluid.movingPhases(), movingPhasei)
            {
                phaseModel& phase = fluid.movingPhases()[movingPhasei];

                phase.URef() =
                    HbyAs[phase.index()]
                  + fvc::reconstruct
                    (
                        alpharAUfs[phase.index()]*mSfGradp
                      - phigFs[phase.index()]
                    );
            }

            if (partialElimination)
            {
                fluid.partialElimination(rAUs, KdUByAs, alphafs, phiKdPhis);
            }
            else
            {
                forAll(fluid.movingPhases(), movingPhasei)
                {
                    phaseModel& phase = fluid.movingPhases()[movingPhasei];
                    MRF.makeRelative(phase.phiRef());
                }
            }

            forAll(fluid.movingPhases(), movingPhasei)
            {
                phaseModel& phase = fluid.movingPhases()[movingPhasei];

                phase.URef().correctBoundaryConditions();
                fvOptions.correct(phase.URef());
            }
        }
```

* if there is final non orthogonal correction
  * $$\phi = phiHbyA + pEqnIncomp.flux()$$
  * define `nSfGradp`
    * $$mSfGradp = \frac{pEqnIncomp.flux()}{rAUf}$$
  * for every moving phase
    * get `phase`, or $\alpha$
    * set `phase.phiRef()`
      * $$phaze.phiRef() = phiHbyA^k + alpharAUf^k mSfGradp$$
    * set the phase dilatation rate
      * $$phase.divU = -pEqnComp^k \cdot p_{rgh}$$
  * relax $p_{rgh}$
  * recalculate `mSfGradp`
    * $$mSfGradp = \frac{pEqnIncomp.flux()}{rAUf}$$
  * for every moving phase:
    * get current moving phase `phase` or $\alpha$
    * set `phase.URef()`
      * $$phase.URef() = HbyA^k + reconstruct(alpharAUf^k mSfGradp - phigF^k)$$
  * if the drag system should be sovled for the velocities and fluxes, then
    * solve it
  * else:
    * for every moving phase:
      * get `phase`
      * get relative velocity in the MRF region
  * for every moving phase
    * get `phase`
    * correct boundary condition for `phase.URef()`
    * correct `fvOptions`of `phase.URef()`

`reconstruct()` can be found in `src\finiteVolume\finiteVolume\fvc\fvcReconstruct.H`, is to reconstruct volField from a face flux field

`MRF.makeRelative()` can be found in `src\finiteVolume\cfdTools\general\MRF\MRFZone.H`, is to make the given absolute velocity relative within the MRF region

###### Update and limit the static pressure

```cpp
    // Update and limit the static pressure
    p = max(p_rgh + rho*gh, pMin);
```

$$
p = \max(p_{rgh} + \rho \mathbf{g}h, pMin)
$$

###### Account for static pressure reference

```cpp
    // Account for static pressure reference
    if (p_rgh.needReference() && fluid.incompressible())
    {
        p += dimensionedScalar
        (
            "p",
            p.dimensions(),
            pRefValue - getRefCellValue(p, pRefCell)
        );
    }
```

* if $p_{rgh}$ needs reference and the fluid is incompressible, then
  * $$p = p + pRefValue - getRefCellValue(p, pRefCell)$$

###### Limit p_rgh

```cpp
    // Limit p_rgh
    p_rgh = p - rho*gh;
```

$$
p_{rgh} = p - \rho \mathbf{g} h
$$

###### Update densities from change in p_rgh

```cpp
    // Update densities from change in p_rgh
    forAll(phases, phasei)
    {
        phaseModel& phase = phases[phasei];
        phase.thermoRef().rho() += phase.thermo().psi()*(p_rgh - p_rgh_0);
    }
```

* for every phase:
  * get `phase`
  * $$\rho = rho + \psi (p_{rgh} - p_{rgh0})$$

###### Correct p_rgh for consistency with p and the updated densities

```cpp
    // Correct p_rgh for consistency with p and the updated densities
    rho = fluid.rho();
    p_rgh = p - rho*gh;
    p_rgh.correctBoundaryConditions();
```

* correct $\rho$, $p_{rgh}$ and the boundary conditions

##### clear rAUs

```cpp
if (!fluid.implicitPhasePressure())
{
    rAUs.clear();
}
```

* if the phase pressure is not treated implicitly, then
  * clear `rAUs`

`fluid.implicitPhasePressure()` can be found in `applications\solvers\multiphase\multiphaseEulerFoam\phaseSystems\phaseSystem\phaseSystem.H`, is to returns true if the phase pressure is treated implicitly in the phase fraction equation for any phase

### PUf

Not available now.

### multiphaseSystems

#### make

##### files

##### options

#### multiphaseSystems.C

##### basicMultiphaseSystem

```cpp
    typedef
        PhaseTransferPhaseSystem
        <
            OneResistanceHeatTransferPhaseSystem
            <
                MomentumTransferPhaseSystem<phaseSystem>
            >
        >
        basicMultiphaseSystem;

    addNamedToRunTimeSelectionTable
    (
        phaseSystem,
        basicMultiphaseSystem,
        dictionary,
        basicMultiphaseSystem
    );
```

* `PhaseTransferPhaseSystem` in `applications\solvers\multiphase\multiphaseEulerFoam\phaseSystems\PhaseSystems\PhaseTransferPhaseSystem\PhaseTransferPhaseSystem.H`, is the class which models non-thermally-coupled or weakly thermally coupled mass transfers.
* `OneResistanceHeatTransferPhaseSystem` in `applications\solvers\multiphase\multiphaseEulerFoam\phaseSystems\PhaseSystems\OneResistanceHeatTransferPhaseSystem\OneResistanceHeatTransferPhaseSystem.H`, is the class which models interfacial heat transfer between a number of phases. **A single heat transfer model is used for each interface.**
* `MomentumTransferPhaseSystem` in `applications\solvers\multiphase\multiphaseEulerFoam\phaseSystems\PhaseSystems\MomentumTransferPhaseSystem\MomentumTransferPhaseSystem.H`, is the class which models **interfacial momentum transfer** between a number of phases. **Drag, virtual mass, lift, wall lubrication and turbulent dispersion are all modelled.** The explicit contribution from the drag is omitted from the transfer matrices, as this forms part of the solution of the pressure equation.

define the type of `basicMultiphaseSystem` and add it to run time selection table.

Tutorials:

* `Laminar/bubbleColumn`

##### interfaceCompositionPhaseChangeMultiphaseSystem

```cpp
    typedef
        InterfaceCompositionPhaseChangePhaseSystem
        <
            PhaseTransferPhaseSystem
            <
                TwoResistanceHeatTransferPhaseSystem
                <
                    MomentumTransferPhaseSystem<phaseSystem>
                >
            >
        >
        interfaceCompositionPhaseChangeMultiphaseSystem;

    addNamedToRunTimeSelectionTable
    (
        phaseSystem,
        interfaceCompositionPhaseChangeMultiphaseSystem,
        dictionary,
        interfaceCompositionPhaseChangeMultiphaseSystem
    );
```

* `InterfaceCompositionPhaseChangePhaseSystem` in `applications\solvers\multiphase\multiphaseEulerFoam\phaseSystems\PhaseSystems\InterfaceCompositionPhaseChangePhaseSystem\InterfaceCompositionPhaseChangePhaseSystem.H`, is the class to provide **interfacial heat and mass transfer** between a number of phases according to a interface composition model. The interface temperature is calculated such that the net rate at which the heat is transferred to the interface is equal to the latent heat consumed by the mass transfer.
* `TwoResistanceHeatTransferPhaseSystem` in `applications\solvers\multiphase\multiphaseEulerFoam\phaseSystems\PhaseSystems\TwoResistanceHeatTransferPhaseSystem\TwoResistanceHeatTransferPhaseSystem.H`, is the class which models interfacial heat transfer between a number of phases. **Two heat transfer models are stored at each interface, one for each phase.** This permits definition of an interface temperature with which heat transfer occurs. It also allows derived systems to define other thermodynamic properties at the interface and therefore represent **phase changes**.

define the type of `interfaceCompositionPhaseChangeMultiphaseSystem` and add it to run time selection table. This can deal with phase change problems.

Tutorials:

* `Laminar/bubbleColumnEvaporating`
* `Laminar/bubbleColumnEvaporatingDissolving`
* `RAS/bubbleColumnEvaporatingReacting`

##### thermalPhaseChangeMultiphaseSystem

```cpp
    typedef
        ThermalPhaseChangePhaseSystem
        <
            PhaseTransferPhaseSystem
            <
                TwoResistanceHeatTransferPhaseSystem
                <
                    MomentumTransferPhaseSystem<phaseSystem>
                >
            >
        >
        thermalPhaseChangeMultiphaseSystem;

    addNamedToRunTimeSelectionTable
    (
        phaseSystem,
        thermalPhaseChangeMultiphaseSystem,
        dictionary,
        thermalPhaseChangeMultiphaseSystem
    );
```

* `ThermalPhaseChangePhaseSystem` in `applications\solvers\multiphase\multiphaseEulerFoam\phaseSystems\PhaseSystems\ThermalPhaseChangePhaseSystem\ThermalPhaseChangePhaseSystem.H`, is the class to provide **interfacial heat and mass transfer** between a number of phases according the interfacial temperature approximated by **the saturation temperature**.

Based on the implemnetaion described in:

    \verbatim
        Peltola, J., Pättikangas, T., Bainbridge, W., Lehnigk, R., Schlegel, F. (2019). On Development and validation of subcooled nucleate boiling models for OpenFOAM Foundation Release. NURETH-18 Conference Proceedings, Portland, Oregon, United States, 2019.
    \endverbatim

The present implementation includes simplified support for **non-volatile
components** in addition to a single volatile component in order to account compressibility effects when non-volatile gas bubbles of non-volatile gas filled pressure reservoirs are present.

The phase change mass transfer calculation is still only dependent on the interfacial temperature estimate and interfacial heat transfer models. The mass diffusion effects in presence of non-volatile components at the interface are neglected.

define the type of `thermalPhaseChangeMultiphaseSystem` and add it to run time selection table.

The difference between `thermalPhaseChangeMultiphaseSystem` and `interfaceCompositionPhaseChangeMultiphaseSystem` is that the former uses `ThermalPhaseChangePhaseSystem` while the latter uses `InterfaceCompositionPhaseChangePhaseSystem`.

Tutorials:

* `Laminar/steamInjection`
* `RAS/wallBoiling`
* `RAS/wallBoilingIATE`

##### populationBalanceMultiphaseSystem

```cpp
    typedef
        PopulationBalancePhaseSystem
        <
            PhaseTransferPhaseSystem
            <
                OneResistanceHeatTransferPhaseSystem
                <
                    MomentumTransferPhaseSystem<phaseSystem>
                >
            >
        >
        populationBalanceMultiphaseSystem;

    addNamedToRunTimeSelectionTable
    (
        phaseSystem,
        populationBalanceMultiphaseSystem,
        dictionary,
        populationBalanceMultiphaseSystem
    );
```

* `PopulationBalancePhaseSystem` in `applications\solvers\multiphase\multiphaseEulerFoam\phaseSystems\PhaseSystems\PopulationBalancePhaseSystem\PopulationBalancePhaseSystem.H`, is the class which provides population balance functionality.

define the type of `populationBalanceMultiphaseSystem` and add it to run time selection table.

Tutorials:

* `Laminar/bubbleColumnPolydisperse`
* `Laminar/titaniaSynthesis`
* `Laminar/titaniaSynthesisSurface`
* `RAS/bubbleColumnPolydisperse`

##### thermalPhaseChangePopulationBalanceMultiphaseSystem

```cpp
    typedef
        ThermalPhaseChangePhaseSystem
        <
            PopulationBalancePhaseSystem
            <
                PhaseTransferPhaseSystem
                <
                    TwoResistanceHeatTransferPhaseSystem
                    <
                        MomentumTransferPhaseSystem<phaseSystem>
                    >
                >
            >
        >
        thermalPhaseChangePopulationBalanceMultiphaseSystem;

    addNamedToRunTimeSelectionTable
    (
        phaseSystem,
        thermalPhaseChangePopulationBalanceMultiphaseSystem,
        dictionary,
        thermalPhaseChangePopulationBalanceMultiphaseSystem
    );
```

define the type of `thermalPhaseChangePopulationBalanceMultiphaseSystem` and add it to run time selection table.

Tutorials:

* `RAS/wallBoilingPolydisperse`
* `RAS/wallBoilingPolydisperseTwoGroups`

| phase type | PhaseTransferPhaseSystem | OneResistanceHeatTransferPhaseSystem | TwoResistanceHeatTransferPhaseSystem | MomentumTransferPhaseSystem | InterfaceCompositionPhaseChangePhaseSystem | ThermalPhaseChangePhaseSystem | PopulationBalancePhaseSystem | function |
| - | - | - | - | - | - | - | - | - |
| basicMultiphaseSystem | yes | yes |  | yes |  |  |  |  |
| interfaceCompositionPhaseChangeMultiphaseSystem | yes |  | yes | yes | yes |  |  |  |
| thermalPhaseChangeMultiphaseSystem | yes |  | yes | yes |  | yes |  |  |
| populationBalanceMultiphaseSystem | yes | yes |  | yes |  |  | yes |  |
| thermalPhaseChangePopulationBalanceMultiphaseSystem | yes |  | yes | yes |  | yes | yes |  |














