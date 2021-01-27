# MovingPhaseModel

- [MovingPhaseModel](#movingphasemodel)
  - [MovingPhaseModel.H](#movingphasemodelh)
    - [include](#include)
    - [template](#template)
    - [protected](#protected)
    - [private](#private)
    - [public](#public)
      - [constructors and destructors](#constructors-and-destructors)
      - [member functions](#member-functions)
        - [correct](#correct)
        - [momentum](#momentum)
        - [Compressibility (variable density)](#compressibility-variable-density)
        - [Momentum transport](#momentum-transport)
        - [Thermophysical transport](#thermophysical-transport)
  - [MovingPhaseModel.C](#movingphasemodelc)
    - [include](#include-1)
    - [static member function](#static-member-function)
    - [constructors](#constructors)
    - [destructor](#destructor)
    - [member functions](#member-functions-1)
      - [correctContinuityError()](#correctcontinuityerror)
  - [phaseCompressibleMomentumTransportModel.H](#phasecompressiblemomentumtransportmodelh)
    - [include](#include-2)
    - [typedef](#typedef)
  - [phaseCompressibleMomentumTransportModelFwd.H](#phasecompressiblemomentumtransportmodelfwdh)

## MovingPhaseModel.H

### include

```cpp
#ifndef MovingPhaseModel_H
#define MovingPhaseModel_H

#include "phaseModel.H"
#include "phaseCompressibleMomentumTransportModel.H"
#include "PhaseThermophysicalTransportModel.H"
```

### template

```cpp
template<class BasePhaseModel>
class MovingPhaseModel
:
    public BasePhaseModel
{
    ...
}
```

define a template

### protected

```cpp
protected:

    // Protected data

        //- Velocity field
        volVectorField U_;

        //- Flux
        surfaceScalarField phi_;

        //- Volumetric flux
        surfaceScalarField alphaPhi_;

        //- Mass flux
        surfaceScalarField alphaRhoPhi_;

        //- Lagrangian acceleration field (needed for virtual-mass)
        mutable tmp<volVectorField> DUDt_;

        //- Lagrangian acceleration field on the faces (needed for virtual-mass)
        mutable tmp<surfaceScalarField> DUDtf_;

        //- Dilatation rate
        tmp<volScalarField> divU_;

        //- Turbulence model
        autoPtr<phaseCompressibleMomentumTransportModel> turbulence_;

        //- Thermophysical transport model
        autoPtr
        <
            PhaseThermophysicalTransportModel
            <
                phaseCompressibleMomentumTransportModel,
                typename BasePhaseModel::thermoModel
            >
        > thermophysicalTransport_;

        //- Continuity error
        volScalarField continuityError_;

        //- Kinetic Energy
        mutable tmp<volScalarField> K_;
```

define protected data

### private

```cpp
private:

    // Private static member functions

        //- Calculate and return the flux field
        tmp<surfaceScalarField> phi(const volVectorField& U) const;
```

### public

#### constructors and destructors

```cpp
public:

    // Constructors

        MovingPhaseModel
        (
            const phaseSystem& fluid,
            const word& phaseName,
            const bool referencePhase,
            const label index
        );


    //- Destructor
    virtual ~MovingPhaseModel();
```

#### member functions

##### correct

```cpp
        //- Correct the phase properties other than the thermo and turbulence
        virtual void correct();

        //- Correct the continuity error
        virtual void correctContinuityError(const volScalarField& source);

        //- Correct the kinematics
        virtual void correctKinematics();

        //- Correct the turbulence
        virtual void correctTurbulence();

        //- Correct the energy transport e.g. alphat
        virtual void correctEnergyTransport();
```

##### momentum

```cpp
        // Momentum

            //- Return whether the phase is stationary
            virtual bool stationary() const;

            //- Return the momentum equation
            virtual tmp<fvVectorMatrix> UEqn();

            //- Return the momentum equation for the face-based algorithm
            virtual tmp<fvVectorMatrix> UfEqn();

            //- Return the velocity
            virtual tmp<volVectorField> U() const;

            //- Access the velocity
            virtual volVectorField& URef();

            //- Return the volumetric flux
            virtual tmp<surfaceScalarField> phi() const;

            //- Access the volumetric flux
            virtual surfaceScalarField& phiRef();

            //- Return the volumetric flux of the phase
            virtual tmp<surfaceScalarField> alphaPhi() const;

            //- Access the volumetric flux of the phase
            virtual surfaceScalarField& alphaPhiRef();

            //- Return the mass flux of the phase
            virtual tmp<surfaceScalarField> alphaRhoPhi() const;

            //- Access the mass flux of the phase
            virtual surfaceScalarField& alphaRhoPhiRef();

            //- Return the substantive acceleration
            virtual tmp<volVectorField> DUDt() const;

            //- Return the substantive acceleration on the faces
            virtual tmp<surfaceScalarField> DUDtf() const;

            //- Return the continuity error
            virtual tmp<volScalarField> continuityError() const;

            //- Return the phase kinetic energy
            virtual tmp<volScalarField> K() const;
```

##### Compressibility (variable density)

```cpp
        // Compressibility (variable density)

            //- Return the phase dilatation rate (d(alpha)/dt + div(alpha*phi))
            virtual tmp<volScalarField> divU() const;

            //- Set the phase dilatation rate (d(alpha)/dt + div(alpha*phi))
            virtual void divU(tmp<volScalarField> divU);
```

##### Momentum transport

```cpp
        // Momentum transport

            //- Return the effective thermal conductivity on a patch
            virtual tmp<scalarField> kappaEff(const label patchi) const;

            //- Return the turbulent kinetic energy
            virtual tmp<volScalarField> k() const;

            //- Return the phase-pressure'
            //  (derivative of phase-pressure w.r.t. phase-fraction)
            virtual tmp<volScalarField> pPrime() const;
```

##### Thermophysical transport

```cpp
        // Thermophysical transport

            //- Return the source term for the energy equation
            virtual tmp<fvScalarMatrix> divq(volScalarField& he) const;

            //- Return the source term for the given specie mass-fraction
            //  equation
            virtual tmp<fvScalarMatrix> divj(volScalarField& Yi) const;

```

## MovingPhaseModel.C

### include

```cpp
#include "MovingPhaseModel.H"
#include "phaseSystem.H"
#include "fixedValueFvPatchFields.H"
#include "slipFvPatchFields.H"
#include "partialSlipFvPatchFields.H"

#include "fvmDdt.H"
#include "fvmDiv.H"
#include "fvmSup.H"
#include "fvcDdt.H"
#include "fvcDiv.H"
#include "fvcFlux.H"
```

### static member function

```cpp
template<class BasePhaseModel>
Foam::tmp<Foam::surfaceScalarField>
Foam::MovingPhaseModel<BasePhaseModel>::phi(const volVectorField& U) const
{
    word phiName(IOobject::groupName("phi", this->name()));

    IOobject phiHeader
    (
        phiName,
        U.mesh().time().timeName(),
        U.mesh(),
        IOobject::NO_READ
    );

    if (phiHeader.typeHeaderOk<surfaceScalarField>(true))
    {
        Info<< "Reading face flux field " << phiName << endl;

        return tmp<surfaceScalarField>
        (
            new surfaceScalarField
            (
                IOobject
                (
                    phiName,
                    U.mesh().time().timeName(),
                    U.mesh(),
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                U.mesh()
            )
        );
    }
    else
    {
        Info<< "Calculating face flux field " << phiName << endl;

        wordList phiTypes
        (
            U.boundaryField().size(),
            calculatedFvPatchScalarField::typeName
        );

        forAll(U.boundaryField(), patchi)
        {
            if (!U.boundaryField()[patchi].assignable())
            {
                phiTypes[patchi] = fixedValueFvPatchScalarField::typeName;
            }
        }

        return tmp<surfaceScalarField>
        (
            new surfaceScalarField
            (
                IOobject
                (
                    phiName,
                    U.mesh().time().timeName(),
                    U.mesh(),
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                fvc::flux(U),
                phiTypes
            )
        );
    }
}
```

### constructors

```cpp
template<class BasePhaseModel>
Foam::MovingPhaseModel<BasePhaseModel>::MovingPhaseModel
(
    const phaseSystem& fluid,
    const word& phaseName,
    const bool referencePhase,
    const label index
)
:
    BasePhaseModel(fluid, phaseName, referencePhase, index),
    U_
    (
        IOobject
        (
            IOobject::groupName("U", this->name()),
            fluid.mesh().time().timeName(),
            fluid.mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        fluid.mesh()
    ),
    phi_(phi(U_)),
    alphaPhi_
    (
        IOobject
        (
            IOobject::groupName("alphaPhi", this->name()),
            fluid.mesh().time().timeName(),
            fluid.mesh()
        ),
        fluid.mesh(),
        dimensionedScalar(dimensionSet(0, 3, -1, 0, 0), 0)
    ),
    alphaRhoPhi_
    (
        IOobject
        (
            IOobject::groupName("alphaRhoPhi", this->name()),
            fluid.mesh().time().timeName(),
            fluid.mesh()
        ),
        fluid.mesh(),
        dimensionedScalar(dimensionSet(1, 0, -1, 0, 0), 0)
    ),
    DUDt_(nullptr),
    DUDtf_(nullptr),
    divU_(nullptr),
    turbulence_
    (
        phaseCompressibleMomentumTransportModel::New
        (
            *this,
            this->thermo().rho(),
            U_,
            alphaRhoPhi_,
            phi_,
            *this
        )
    ),
    thermophysicalTransport_
    (
        PhaseThermophysicalTransportModel
        <
            phaseCompressibleMomentumTransportModel,
            typename BasePhaseModel::thermoModel
        >::New(turbulence_, this->thermo_)
    ),
    continuityError_
    (
        IOobject
        (
            IOobject::groupName("continuityError", this->name()),
            fluid.mesh().time().timeName(),
            fluid.mesh()
        ),
        fluid.mesh(),
        dimensionedScalar(dimDensity/dimTime, 0)
    ),
    K_(nullptr)
{
    phi_.writeOpt() = IOobject::AUTO_WRITE;

    correctKinematics();
}
```

initialize

### destructor

```cpp
template<class BasePhaseModel>
Foam::MovingPhaseModel<BasePhaseModel>::~MovingPhaseModel()
{}
```

### member functions

#### correctContinuityError()

```cpp
template<class BasePhaseModel>
void Foam::MovingPhaseModel<BasePhaseModel>::correctContinuityError
(
    const volScalarField& source
)
{
    volScalarField& rho = this->thermoRef().rho();

    continuityError_ = fvc::ddt(*this, rho) + fvc::div(alphaRhoPhi_) - source;
}
```

$$
continuityError_  = \frac{\partial \alpha^k \rho^k}{\partial t} + \nabla \cdot (\alpha^k \rho^k \phi^k) - source
$$


##  phaseCompressibleMomentumTransportModel.H

Typedef for phaseCompressibleMomentumTransportModel

### include

```cpp
#ifndef phaseCompressibleMomentumTransportModel_H
#define phaseCompressibleMomentumTransportModel_H

#include "phaseCompressibleMomentumTransportModelFwd.H"
#include "PhaseCompressibleMomentumTransportModel.H"
#include "phaseModel.H"
```

### typedef

```cpp
namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

typedef PhaseCompressibleMomentumTransportModel<phaseModel>
    phaseCompressibleMomentumTransportModel;

// template specialisation for TransportTraits<phaseModel>
template<>
class TransportTraits<phaseModel>
{
public:

    typedef fluidThermo thermoModel;

    static const thermoModel& thermo(const phaseModel& pm)
    {
        return pm.thermo();
    }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam
```

## phaseCompressibleMomentumTransportModelFwd.H

Forward declaration of typedef for phaseCompressibleMomentumTransportModel

```cpp
namespace Foam
{
    class phaseModel;

    template<class TransportModel>
    class PhaseCompressibleMomentumTransportModel;

    typedef PhaseCompressibleMomentumTransportModel<phaseModel>
          phaseCompressibleMomentumTransportModel;
}
```

