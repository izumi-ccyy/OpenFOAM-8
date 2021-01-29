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
      - [correct()](#correct-1)
      - [correctKinematics()](#correctkinematics)
      - [correctTurbulence()](#correctturbulence)
      - [correctEnergyTransport()](#correctenergytransport)
      - [stationary()](#stationary)
      - [UEqn()](#ueqn)
      - [UfEqn()](#ufeqn)
      - [U(), URef(), phi(), phiRef(), alphaPhi(), alphaPhiRef(), alphaRhoPhi(), alphaRhoPhiRef()](#u-uref-phi-phiref-alphaphi-alphaphiref-alpharhophi-alpharhophiref)
      - [DUDt()](#dudt)
      - [DUDtf()](#dudtf)
      - [continuityError()](#continuityerror)
      - [K()](#k)
      - [divU()](#divu)
      - [divU(divU)](#divudivu)
      - [kappaEff(), k(), pPrime()](#kappaeff-k-pprime)
      - [divq(), divj()](#divq-divj)
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
continuityError_ = \frac{\partial \alpha^k \rho^k}{\partial t} + \nabla \cdot (\alpha^k \rho^k \phi^k) - source
$$

this is embeded in definition of `correctContinuityError()` as in `applications\solvers\multiphase\multiphaseEulerFoam\phaseSystems\phaseSystem\phaseSystem.C`, where `source` is calculated and passed to this function as parameter

```cpp
void Foam::phaseSystem::correctContinuityError()
{
    const PtrList<volScalarField> dmdts = this->dmdts();

    forAll(movingPhaseModels_, movingPhasei)
    {
        phaseModel& phase = movingPhaseModels_[movingPhasei];
        const volScalarField& alpha = phase;
        volScalarField& rho = phase.thermoRef().rho();

        volScalarField source
        (
            volScalarField::New
            (
                IOobject::groupName("source", phase.name()),
                mesh_,
                dimensionedScalar(dimDensity/dimTime, 0)
            )
        );

        if (fvOptions().appliesToField(rho.name()))
        {
            source += fvOptions()(alpha, rho)&rho;
        }

        if (dmdts.set(phase.index()))
        {
            source += dmdts[phase.index()];
        }

        phase.correctContinuityError(source);
    }
}
```

in `applications\solvers\multiphase\multiphaseEulerFoam\multiphaseEulerFoam\multiphaseEulerFoam.C`, `correctContinuityError()` is performed to get `continuityError_`

```cpp
        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            fluid.solve(rAUs, rAUfs);
            fluid.correct();
            fluid.correctContinuityError();
```

#### correct()

```cpp
template<class BasePhaseModel>
void Foam::MovingPhaseModel<BasePhaseModel>::correct()
{
    BasePhaseModel::correct();
    this->fluid().MRF().correctBoundaryVelocity(U_);
}
```

* correct with `BasePhaseModel` and `MRF`

#### correctKinematics()

```cpp
template<class BasePhaseModel>
void Foam::MovingPhaseModel<BasePhaseModel>::correctKinematics()
{
    BasePhaseModel::correctKinematics();

    if (DUDt_.valid())
    {
        DUDt_.clear();
        DUDt();
    }

    if (DUDtf_.valid())
    {
        DUDtf_.clear();
        DUDtf();
    }

    if (K_.valid())
    {
        K_.clear();
        K();
    }
}
```

* correct with `BasePhaseModel`
* clear existed `DUDt_`, `DUDtf` and `K_` and recalculate them

#### correctTurbulence()

```cpp
template<class BasePhaseModel>
void Foam::MovingPhaseModel<BasePhaseModel>::correctTurbulence()
{
    BasePhaseModel::correctTurbulence();

    turbulence_->correct();
}
```

* correct turbulence with `BasePhaseModel` and `turbulence_`

#### correctEnergyTransport()

```cpp
template<class BasePhaseModel>
void Foam::MovingPhaseModel<BasePhaseModel>::correctEnergyTransport()
{
    BasePhaseModel::correctEnergyTransport();
    thermophysicalTransport_->correct();
}
```

* correct turbulence with `BasePhaseModel` and `thermophysicalTransport_`

#### stationary() 

```cpp
template<class BasePhaseModel>
bool Foam::MovingPhaseModel<BasePhaseModel>::stationary() const
{
    return false;
}
```

* return `false`, since it's moving phase

#### UEqn()

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


The `UEqn` is defined as

$$
\frac{\partial \alpha^k \rho^k \mathbf{U}^k}{\partial t} + \nabla \cdot (\alpha^k \rho^k \phi \mathbf{U}^k) + SuSp(continuityError_, \mathbf{U}^k) + MRF(\alpha^k \rho^k \mathbf{U}^k) -\nabla \cdot \left[\alpha^k \rho^k \nu_{Eff}^k \left((\nabla \mathbf{U}^k+(\nabla \mathbf{U}^k)^T) - \frac{2}{3} (\nabla \cdot \mathbf{U}^k) \mathbf{I}\right)\right]
$$

or

$$
\frac{\partial \alpha^k \rho^k \mathbf{U}^k}{\partial t} + \nabla \cdot (\alpha^k \rho^k \phi \mathbf{U}^k) + SuSp(continuityError_, \mathbf{U}^k) + MRF(\alpha^k \rho^k \mathbf{U}^k) -\nabla \cdot \tau^k
$$

`continuityError_` is defined above as

$$
continuityError_  = \frac{\partial \alpha^k \rho^k}{\partial t} + \nabla \cdot (\alpha^k \rho^k \phi^k) - source
$$

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

#### UfEqn()

```cpp
template<class BasePhaseModel>
Foam::tmp<Foam::fvVectorMatrix>
Foam::MovingPhaseModel<BasePhaseModel>::UfEqn()
{
    // As the "normal" U-eqn but without the ddt terms

    const volScalarField& alpha = *this;
    const volScalarField& rho = this->thermo().rho();

    return
    (
        fvm::div(alphaRhoPhi_, U_)
      + fvm::SuSp(fvc::ddt(*this, rho) - this->continuityError(), U_)
      + this->fluid().MRF().DDt(alpha*rho, U_)
      + turbulence_->divDevTau(U_)
    );
}
```

#### U(), URef(), phi(), phiRef(), alphaPhi(), alphaPhiRef(), alphaRhoPhi(), alphaRhoPhiRef()

```cpp
template<class BasePhaseModel>
Foam::tmp<Foam::volVectorField>
Foam::MovingPhaseModel<BasePhaseModel>::U() const
{
    return U_;
}


template<class BasePhaseModel>
Foam::volVectorField&
Foam::MovingPhaseModel<BasePhaseModel>::URef()
{
    return U_;
}


template<class BasePhaseModel>
Foam::tmp<Foam::surfaceScalarField>
Foam::MovingPhaseModel<BasePhaseModel>::phi() const
{
    return phi_;
}


template<class BasePhaseModel>
Foam::surfaceScalarField&
Foam::MovingPhaseModel<BasePhaseModel>::phiRef()
{
    return phi_;
}


template<class BasePhaseModel>
Foam::tmp<Foam::surfaceScalarField>
Foam::MovingPhaseModel<BasePhaseModel>::alphaPhi() const
{
    return alphaPhi_;
}


template<class BasePhaseModel>
Foam::surfaceScalarField&
Foam::MovingPhaseModel<BasePhaseModel>::alphaPhiRef()
{
    return alphaPhi_;
}


template<class BasePhaseModel>
Foam::tmp<Foam::surfaceScalarField>
Foam::MovingPhaseModel<BasePhaseModel>::alphaRhoPhi() const
{
    return alphaRhoPhi_;
}


template<class BasePhaseModel>
Foam::surfaceScalarField&
Foam::MovingPhaseModel<BasePhaseModel>::alphaRhoPhiRef()
{
    return alphaRhoPhi_;
}
```

return the variables and return the references of the variables to modify them:

* U(), 
* URef(),
* phi(), 
* phiRef(),
* alphaPhi(), 
* alphaPhiRef(), 
* alphaRhoPhi(), 
* alphaRhoPhiRef()


#### DUDt()

```cpp
template<class BasePhaseModel>
Foam::tmp<Foam::volVectorField>
Foam::MovingPhaseModel<BasePhaseModel>::DUDt() const
{
    if (!DUDt_.valid())
    {
        DUDt_ = fvc::ddt(U_) + fvc::div(phi_, U_) - fvc::div(phi_)*U_;
    }

    return tmp<volVectorField>(DUDt_());
}
```

* return the substantive acceleration

* if `DUDt_` is not calculated, then
  * $$DUDt_ = \frac{\partial \mathbf{U}^k}{\partial t} + \nabla \cdot (\phi \mathbf{U}^k) - \nabla \phi \cdot \mathbf{U}^k = \frac{\partial \mathbf{U}^k}{\partial t} + \phi \nabla \mathbf{U}^k$$

material derivative or  substantial derivative:

$$
\frac{DF}{Dt} = \frac{\partial F}{\partial t} + \mathbf{U} \cdot \nabla F
$$

acceleration:

$$
\frac{D \mathbf{U}}{Dt} = \frac{\partial \mathbf{U}}{\partial t} + \mathbf{U} \cdot \nabla \mathbf{U}
$$

#### DUDtf()

```cpp
template<class BasePhaseModel>
Foam::tmp<Foam::surfaceScalarField>
Foam::MovingPhaseModel<BasePhaseModel>::DUDtf() const
{
    if (!DUDtf_.valid())
    {
        DUDtf_ = byDt(phi_ - phi_.oldTime());
    }

    return tmp<surfaceScalarField>(DUDtf_());
}
```

* Return the substantive acceleration on the faces
* if `DUDtf_` is not calculated, then
  * $$DUDtf_ = \frac{\phi - \phi_{old}}{\Delta t}$$

`byDt()` can be found in `applications\solvers\multiphase\multiphaseEulerFoam\phaseSystems\phaseSystem\phaseSystem.C`

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

$$
byDt = \frac{s_f}{\Delta t}
$$

            
#### continuityError()

```cpp
template<class BasePhaseModel>
Foam::tmp<Foam::volScalarField>
Foam::MovingPhaseModel<BasePhaseModel>::continuityError() const
{
    return continuityError_;
}
```

* return `continuityError_`

#### K()

```cpp
template<class BasePhaseModel>
Foam::tmp<Foam::volScalarField>
Foam::MovingPhaseModel<BasePhaseModel>::K() const
{
    if (!K_.valid())
    {
        K_ = volScalarField::New
        (
            IOobject::groupName("K", this->name()),
            0.5*magSqr(this->U())
        );
    }

    return tmp<volScalarField>(K_());
}
```

* Return the phase kinetic energy
  * $$K = \frac{\|\mathbf{U}^k\|^2}{2}$$

#### divU()

```cpp
template<class BasePhaseModel>
Foam::tmp<Foam::volScalarField>
Foam::MovingPhaseModel<BasePhaseModel>::divU() const
{
    return divU_.valid() ? tmp<volScalarField>(divU_()) : tmp<volScalarField>();
}
```

* if `divU_` exits, return it
* else return an empty

#### divU(divU)

```cpp
template<class BasePhaseModel>
void Foam::MovingPhaseModel<BasePhaseModel>::divU(tmp<volScalarField> divU)
{
    divU_ = divU;
}
```

* set `divU_`  as `divU`
* $$divU_ = divU$$

#### kappaEff(), k(), pPrime()

```cpp
template<class BasePhaseModel>
Foam::tmp<Foam::scalarField>
Foam::MovingPhaseModel<BasePhaseModel>::kappaEff(const label patchi) const
{
    return thermophysicalTransport_->kappaEff(patchi);
}


template<class BasePhaseModel>
Foam::tmp<Foam::volScalarField>
Foam::MovingPhaseModel<BasePhaseModel>::k() const
{
    return turbulence_->k();
}


template<class BasePhaseModel>
Foam::tmp<Foam::volScalarField>
Foam::MovingPhaseModel<BasePhaseModel>::pPrime() const
{
    return turbulence_->pPrime();
}
```

* return variables from thermal or turbulence model
  * kappaEff(), 
  * k(), 
  * pPrime()

#### divq(), divj()

```cpp
template<class BasePhaseModel>
Foam::tmp<Foam::fvScalarMatrix>
Foam::MovingPhaseModel<BasePhaseModel>::divq(volScalarField& he) const
{
    return thermophysicalTransport_->divq(he);
}


template<class BasePhaseModel>
Foam::tmp<Foam::fvScalarMatrix>
Foam::MovingPhaseModel<BasePhaseModel>::divj(volScalarField& Yi) const
{
    return thermophysicalTransport_->divj(Yi);
}
```

* return divq(), the source term for the energy equation, and divj(), the source term for the given specie mass-fraction equation
* they can be found in `src\ThermophysicalTransportModels\laminar\Fourier\Fourier.C` or `src\ThermophysicalTransportModels\turbulence\eddyDiffusivity\eddyDiffusivity.C`

in the former:

```cpp
template<class BasicThermophysicalTransportModel>
tmp<fvScalarMatrix>
Fourier<BasicThermophysicalTransportModel>::divq(volScalarField& he) const
{
    return -fvm::laplacian(this->alpha()*this->thermo().alpha(), he);
}


template<class BasicThermophysicalTransportModel>
tmp<fvScalarMatrix>
Fourier<BasicThermophysicalTransportModel>::divj(volScalarField& Yi) const
{
    return -fvm::laplacian(this->alpha()*this->thermo().alpha(), Yi);
}
```

in the later:

```cpp
template<class TurbulenceThermophysicalTransportModel>
tmp<fvScalarMatrix>
eddyDiffusivity<TurbulenceThermophysicalTransportModel>::divq
(
    volScalarField& he
) const
{
    return -fvm::laplacian(this->alpha()*this->alphaEff(), he);
}


template<class TurbulenceThermophysicalTransportModel>
tmp<fvScalarMatrix>
eddyDiffusivity<TurbulenceThermophysicalTransportModel>::divj
(
    volScalarField& Yi
) const
{
    return -fvm::laplacian(this->alpha()*this->DEff(Yi), Yi);
}
```

The later should be adopted, so

$$
divq = -\nabla \cdot \left[ \alpha^k \alpha_{Eff, thermo}^k \nabla he^k\right]
$$

$$
divj = -\nabla \cdot \left[ \alpha^k D_{Eff}^k \nabla Yi^k\right]
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

