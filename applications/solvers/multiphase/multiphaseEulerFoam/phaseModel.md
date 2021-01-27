# phaseModel

- [phaseModel](#phasemodel)
  - [phaseModel.H](#phasemodelh)
    - [private data](#private-data)
    - [public](#public)
    - [public member functions](#public-member-functions)
      - [return protected data and correction](#return-protected-data-and-correction)
      - [density variation and compressibility](#density-variation-and-compressibility)
      - [thermo](#thermo)
      - [species](#species)
      - [momentum](#momentum)
      - [transport](#transport)
  - [phaseModel.C](#phasemodelc)
    - [incnlude](#incnlude)
    - [Static Data Members](#static-data-members)
    - [constructors](#constructors)
    - [clone](#clone)
    - [destructors](#destructors)
    - [member functions](#member-functions)
      - [return](#return)
      - [correct](#correct)
        - [correct](#correct-1)
        - [other correct](#other-correct)
        - [read](#read)
        - [correct flow](#correct-flow)
  - [phaseModelNew.C](#phasemodelnewc)
  - [phaseModels.C](#phasemodelsc)
    - [purePhaseModel](#purephasemodel)
    - [pureStationaryPhaseModel](#purestationaryphasemodel)
    - [pureIsothermalPhaseModel](#pureisothermalphasemodel)
    - [pureStationaryIsothermalPhaseModel](#purestationaryisothermalphasemodel)
    - [multiComponentPhaseModel](#multicomponentphasemodel)
    - [multiComponentIsothermalPhaseModel](#multicomponentisothermalphasemodel)
    - [reactingPhaseModel](#reactingphasemodel)

## phaseModel.H

### private data

```cpp
    // Private Data

        //- Reference to the phaseSystem to which this phase belongs
        const phaseSystem& fluid_;

        //- Name of phase
        word name_;

        //- Index of phase
        label index_;

        //- Return the residual phase-fraction for given phase
        //  Used to stabilize the phase momentum as the phase-fraction -> 0
        dimensionedScalar residualAlpha_;

        //- Optional maximum phase-fraction (e.g. packing limit)
        scalar alphaMax_;

        //- Diameter model
        autoPtr<diameterModel> diameterModel_;
```

### public

```cpp
public:

    //- Runtime type information
    ClassName("phaseModel");


    // Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            phaseModel,
            phaseSystem,
            (
                const phaseSystem& fluid,
                const word& phaseName,
                const bool referencePhase,
                const label index
            ),
            (fluid, phaseName, referencePhase, index)
        );


    // Constructors

        phaseModel
        (
            const phaseSystem& fluid,
            const word& phaseName,
            const bool referencePhase,
            const label index
        );

        //- Return clone
        autoPtr<phaseModel> clone() const;


    // Selectors

        static autoPtr<phaseModel> New
        (
            const phaseSystem& fluid,
            const word& phaseName,
            const bool referencePhase,
            const label index
        );

        //- Return a pointer to a new phase created on freestore
        //  from Istream
        class iNew
        {
            const phaseSystem& fluid_;
            const word& referencePhaseName_;
            mutable label indexCounter_;

        public:

            iNew
            (
                const phaseSystem& fluid,
                const word& referencePhaseName
            )
            :
                fluid_(fluid),
                referencePhaseName_(referencePhaseName),
                indexCounter_(-1)
            {}

            autoPtr<phaseModel> operator()(Istream& is) const
            {
                indexCounter_++;

                const word phaseName(is);

                return autoPtr<phaseModel>
                (
                    phaseModel::New
                    (
                        fluid_,
                        phaseName,
                        phaseName == referencePhaseName_,
                        indexCounter_
                    )
                );
            }
        };


    //- Destructor
    virtual ~phaseModel();
```

define name, constructors, selectors and destructors

### public member functions

#### return protected data and correction

```cpp
    // Member Functions

        //- Return the name of this phase
        const word& name() const;

        //- Return the name of the phase for use as the keyword in PtrDictionary
        const word& keyword() const;

        //- Return the index of the phase
        label index() const;

        //- Return the system to which this phase belongs
        const phaseSystem& fluid() const;

        //- Return the residual phase-fraction for given phase
        //  Used to stabilize the phase momentum as the phase-fraction -> 0
        const dimensionedScalar& residualAlpha() const;

        //- Return the maximum phase-fraction (e.g. packing limit)
        scalar alphaMax() const;

        //- Return the Sauter-mean diameter
        tmp<volScalarField> d() const;

        //- Return const-reference to diameterModel of the phase
        const autoPtr<diameterModel>& dPtr() const;

        //- Correct the phase properties
        virtual void correct();

        //- Correct the continuity error
        virtual void correctContinuityError(const volScalarField& source);

        //- Correct the kinematics
        virtual void correctKinematics();

        //- Correct the thermodynamics
        virtual void correctThermo();

        //- Correct the reactions
        virtual void correctReactions();

        //- Correct the species concentrations
        virtual void correctSpecies();

        //- Correct the turbulence
        virtual void correctTurbulence();

        //- Correct the energy transport
        virtual void correctEnergyTransport();

        //- Ensure that the flux at inflow/outflow BCs is preserved
        void correctInflowOutflow(surfaceScalarField& alphaPhi) const;

        //- Read phase properties dictionary
        virtual bool read();
```

define some correct functions

#### density variation and compressibility

```cpp
        // Density variation and compressibility

            //- Return true if the phase is incompressible otherwise false
            virtual bool incompressible() const = 0;

            //- Return true if the phase is constant density otherwise false
            virtual bool isochoric() const = 0;

            //- Return the phase dilatation rate (d(alpha)/dt + div(alpha*phi))
            virtual tmp<volScalarField> divU() const = 0;

            //- Set the phase dilatation rate (d(alpha)/dt + div(alpha*phi))
            virtual void divU(tmp<volScalarField> divU) = 0;
```

#### thermo

```cpp
        // Thermo

            //- Return the thermophysical model
            virtual const rhoThermo& thermo() const = 0;

            //- Access the thermophysical model
            virtual rhoThermo& thermoRef() = 0;

            //- Return the density field
            virtual tmp<volScalarField> rho() const = 0;

            //- Return whether the phase is isothermal
            virtual bool isothermal() const = 0;

            //- Return the enthalpy equation
            virtual tmp<fvScalarMatrix> heEqn() = 0;
```

#### species

```cpp
        // Species

            //- Return whether the phase is pure (i.e., not multi-component)
            virtual bool pure() const = 0;

            //- Return the species fraction equation
            virtual tmp<fvScalarMatrix> YiEqn(volScalarField& Yi) = 0;

            //- Return the species mass fractions
            virtual const PtrList<volScalarField>& Y() const = 0;

            //- Return a species mass fraction by name
            virtual const volScalarField& Y(const word& name) const = 0;

            //- Access the species mass fractions
            virtual PtrList<volScalarField>& YRef() = 0;

            //- Return the active species mass fractions
            virtual const UPtrList<volScalarField>& YActive() const = 0;

            //- Access the active species mass fractions
            virtual UPtrList<volScalarField>& YActiveRef() = 0;

            //- Return the fuel consumption rate matrix
            virtual tmp<fvScalarMatrix> R(volScalarField& Yi) const = 0;
```

#### momentum

```cpp
        // Momentum

            //- Return whether the phase is stationary
            virtual bool stationary() const = 0;

            //- Return the momentum equation
            virtual tmp<fvVectorMatrix> UEqn() = 0;

            //- Return the momentum equation for the face-based algorithm
            virtual tmp<fvVectorMatrix> UfEqn() = 0;

            //- Return the velocity
            virtual tmp<volVectorField> U() const = 0;

            //- Access the velocity
            virtual volVectorField& URef() = 0;

            //- Return the volumetric flux
            virtual tmp<surfaceScalarField> phi() const = 0;

            //- Access the volumetric flux
            virtual surfaceScalarField& phiRef() = 0;

            //- Return the volumetric flux of the phase
            virtual tmp<surfaceScalarField> alphaPhi() const = 0;

            //- Access the volumetric flux of the phase
            virtual surfaceScalarField& alphaPhiRef() = 0;

            //- Return the mass flux of the phase
            virtual tmp<surfaceScalarField> alphaRhoPhi() const = 0;

            //- Access the mass flux of the phase
            virtual surfaceScalarField& alphaRhoPhiRef() = 0;

            //- Return the substantive acceleration
            virtual tmp<volVectorField> DUDt() const = 0;

            //- Return the substantive acceleration on the faces
            virtual tmp<surfaceScalarField> DUDtf() const = 0;

            //- Return the continuity error
            virtual tmp<volScalarField> continuityError() const = 0;

            //- Return the phase kinetic energy
            virtual tmp<volScalarField> K() const = 0;
```

#### transport

```cpp
        // Transport

            //- Effective thermal turbulent diffusivity for temperature
            //  of mixture for patch [W/m/K]
            virtual tmp<scalarField> kappaEff(const label patchi) const = 0;

            //- Return the turbulent kinetic energy
            virtual tmp<volScalarField> k() const = 0;

            //- Return the phase-pressure'
            //  (derivative of phase-pressure w.r.t. phase-fraction)
            virtual tmp<volScalarField> pPrime() const = 0;
```

## phaseModel.C

### incnlude

```cpp
#include "phaseModel.H"
#include "phaseSystem.H"
#include "diameterModel.H"
```

### Static Data Members

```cpp
namespace Foam
{
    defineTypeNameAndDebug(phaseModel, 0);
    defineRunTimeSelectionTable(phaseModel, phaseSystem);
}
```

### constructors

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

it should be noted that there is only one definition in class `phaseSystem` owns the type of `volScalarField`, which is also the parent class.

### clone

```cpp
Foam::autoPtr<Foam::phaseModel> Foam::phaseModel::clone() const
{
    NotImplemented;
    return autoPtr<phaseModel>(nullptr);
}
```

### destructors

```cpp
Foam::phaseModel::~phaseModel()
{}
```

### member functions

#### return

```cpp
const Foam::word& Foam::phaseModel::name() const
{
    return name_;
}


const Foam::word& Foam::phaseModel::keyword() const
{
    return name_;
}


Foam::label Foam::phaseModel::index() const
{
    return index_;
}


const Foam::phaseSystem& Foam::phaseModel::fluid() const
{
    return fluid_;
}


const Foam::dimensionedScalar& Foam::phaseModel::residualAlpha() const
{
    return residualAlpha_;
}


Foam::scalar Foam::phaseModel::alphaMax() const
{
    return alphaMax_;
}


Foam::tmp<Foam::volScalarField> Foam::phaseModel::d() const
{
    return diameterModel_().d();
}


const Foam::autoPtr<Foam::diameterModel>& Foam::phaseModel::dPtr() const
{
    return diameterModel_;
}
```

return:

* name_
* index_
* fluid_
* residualAlpha_
* alphaMax_
* diameter of diameter model
* diameter model

#### correct

##### correct

```cpp
void Foam::phaseModel::correct()
{
    diameterModel_->correct();
}
```

correct by diameter model

##### other correct

```cpp
void Foam::phaseModel::correctContinuityError(const volScalarField& source)
{}


void Foam::phaseModel::correctKinematics()
{}


void Foam::phaseModel::correctThermo()
{}

void Foam::phaseModel::correctReactions()
{}

void Foam::phaseModel::correctSpecies()
{}

void Foam::phaseModel::correctTurbulence()
{}


void Foam::phaseModel::correctEnergyTransport()
{}
```

empty definition

##### read

```cpp
bool Foam::phaseModel::read()
{
    return diameterModel_->read(fluid_.subDict(name_));
}
```

read by diameter model

##### correct flow

```cpp
void Foam::phaseModel::correctInflowOutflow(surfaceScalarField& alphaPhi) const
{
    surfaceScalarField::Boundary& alphaPhiBf = alphaPhi.boundaryFieldRef();
    const volScalarField::Boundary& alphaBf = boundaryField();
    const surfaceScalarField::Boundary& phiBf = phi()().boundaryField();

    forAll(alphaPhiBf, patchi)
    {
        fvsPatchScalarField& alphaPhip = alphaPhiBf[patchi];

        if (!alphaPhip.coupled())
        {
            alphaPhip = phiBf[patchi]*alphaBf[patchi];
        }
    }
}
```

if `alphaPhiBf` is not coupled, then

$$
alphaPhiBf = \phi_{Bf} \alpha_{Bf}
$$

## phaseModelNew.C

## phaseModels.C

### purePhaseModel

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

### pureStationaryPhaseModel

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

### pureIsothermalPhaseModel

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

### pureStationaryIsothermalPhaseModel

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

### multiComponentPhaseModel

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

### multiComponentIsothermalPhaseModel

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

### reactingPhaseModel

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
