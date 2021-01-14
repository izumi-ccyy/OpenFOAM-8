# multiphaseEulerFoam

- [multiphaseEulerFoam](#multiphaseeulerfoam)
  - [phaseSystems](#phasesystems)
    - [phaseSystem](#phasesystem)
      - [phaseSystem.H](#phasesystemh)
      - [phaseSystemI.H](#phasesystemih)
      - [phaseSystemNew.C](#phasesystemnewc)
      - [phaseSystem.C](#phasesystemc)
        - [calcPhi](#calcphi)
        - [generatePairs](#generatepairs)
        - [sumAlphaMoving](#sumalphamoving)
        - [setMixtureU](#setmixtureu)
        - [setMixturePhi](#setmixturephi)
        - [nHatfv](#nhatfv)
        - [nHatf](#nhatf)
        - [correctContactAngle](#correctcontactangle)
        - [K](#k)
        - [Constructors 1: phaseSystem](#constructors-1-phasesystem)
        - [Destructor ~phaseSystem](#destructor-phasesystem)
        - [rho](#rho)
        - [U](#u)
        - [E](#e)
        - [sigma1](#sigma1)
        - [sigma2](#sigma2)
        - [nearInterface](#nearinterface)
        - [dmdtf](#dmdtf)
        - [dmdts](#dmdts)
        - [incompressible](#incompressible)
        - [implicitPhasePressure](#implicitphasepressure)
        - [surfaceTension](#surfacetension)
        - [correct](#correct)
        - [correctContinuityError](#correctcontinuityerror)
        - [correctKinematics](#correctkinematics)
        - [correctThermo](#correctthermo)
        - [correctReactions](#correctreactions)
        - [correctSpecies](#correctspecies)
        - [correctTurbulelnce](#correctturbulelnce)
        - [correctEnergyTransport](#correctenergytransport)
        - [read](#read)
        - [byDt1](#bydt1)
        - [byDt2](#bydt2)
      - [phaseSystemSolve.C](#phasesystemsolvec)
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
    - [PUf](#puf)

## phaseSystems

### phaseSystem

#### phaseSystem.H

```cpp
#include "IOdictionary.H"

#include "phaseModel.H"
#include "phasePair.H"
#include "orderedPhasePair.H"
#include "HashPtrTable.H"
#include "PtrListDictionary.H"

#include "IOMRFZoneList.H"
#include "fvOptions.H"

#include "volFields.H"
#include "surfaceFields.H"
#include "fvMatricesFwd.H"
```

include `phaseModel.H` and `phasePair.H` 

```cpp
// Public Typedefs

    typedef HashPtrTable<fvVectorMatrix> momentumTransferTable;

    typedef HashPtrTable<fvScalarMatrix> heatTransferTable;

    typedef HashPtrTable<fvScalarMatrix> specieTransferTable;

    typedef PtrListDictionary<phaseModel> phaseModelList;

    typedef UPtrList<phaseModel> phaseModelPartialList;

    typedef
        HashTable<autoPtr<phasePair>, phasePairKey, phasePairKey::hash>
        phasePairTable;

    typedef
        HashPtrTable<volScalarField, phasePairKey, phasePairKey::hash>
        dmdtfTable;

    typedef
        HashPtrTable
        <
            HashPtrTable<volScalarField>,
            phasePairKey,
            phasePairKey::hash
        >
        dmidtfTable;
```

define list to strore properties

#### phaseSystemI.H

define some simple inline functions, such as:

```cpp
inline const Foam::fvMesh& Foam::phaseSystem::mesh() const
{
    return mesh_;
}
```

#### phaseSystemNew.C

create a new `phaseSystem`

#### phaseSystem.C

##### calcPhi

```cpp
Foam::tmp<Foam::surfaceScalarField> Foam::phaseSystem::calcPhi
(
    const phaseModelList& phaseModels
) const
{
    tmp<surfaceScalarField> tmpPhi
    (
        surfaceScalarField::New
        (
            "phi",
            fvc::interpolate(phaseModels[0])*phaseModels[0].phi()
        )
    );

    for (label phasei=1; phasei<phaseModels.size(); phasei++)
    {
        tmpPhi.ref() +=
            fvc::interpolate(phaseModels[phasei])*phaseModels[phasei].phi();
    }

    return tmpPhi;
}
```

**what is `phaseModels`?**

`tmpPhi` is the `phi` of `phase[0]`

##### generatePairs

```cpp
void Foam::phaseSystem::generatePairs
(
    const dictTable& modelDicts
)
{
    forAllConstIter(dictTable, modelDicts, iter)
    {
        const phasePairKey& key = iter.key();

        // pair already exists
        if (phasePairs_.found(key))
        {}

        // new ordered pair
        else if (key.ordered())
        {
            phasePairs_.insert
            (
                key,
                autoPtr<phasePair>
                (
                    new orderedPhasePair
                    (
                        phaseModels_[key.first()],
                        phaseModels_[key.second()]
                    )
                )
            );
        }

        // new unordered pair
        else
        {
            phasePairs_.insert
            (
                key,
                autoPtr<phasePair>
                (
                    new phasePair
                    (
                        phaseModels_[key.first()],
                        phaseModels_[key.second()]
                    )
                )
            );
        }
    }
}
```

**what is `phasePairs`?**

##### sumAlphaMoving

```cpp
Foam::tmp<Foam::volScalarField> Foam::phaseSystem::sumAlphaMoving() const
{
    tmp<volScalarField> sumAlphaMoving
    (
        volScalarField::New
        (
            "sumAlphaMoving",
            movingPhaseModels_[0],
            calculatedFvPatchScalarField::typeName
        )
    );

    for
    (
        label movingPhasei=1;
        movingPhasei<movingPhaseModels_.size();
        movingPhasei++
    )
    {
        sumAlphaMoving.ref() += movingPhaseModels_[movingPhasei];
    }

    return sumAlphaMoving;
}
```

Return the sum of the phase fractions of the moving phases

##### setMixtureU

```cpp
void Foam::phaseSystem::setMixtureU(const volVectorField& Um0)
{
    // Calculate the mean velocity difference with respect to Um0
    // from the current velocity of the moving phases
    volVectorField dUm(Um0);

    forAll(movingPhaseModels_, movingPhasei)
    {
        dUm -=
            movingPhaseModels_[movingPhasei]
           *movingPhaseModels_[movingPhasei].U();
    }

    forAll(movingPhaseModels_, movingPhasei)
    {
        movingPhaseModels_[movingPhasei].URef() += dUm;
    }
}
```

Re-normalise the velocity of the phases around the specified mixture mean

$$
dUm = dUm - \sum \mathbf{U}_{moving phase}
$$

##### setMixturePhi

```cpp
void Foam::phaseSystem::setMixturePhi
(
    const PtrList<surfaceScalarField>& alphafs,
    const surfaceScalarField& phim0
)
{
    // Calculate the mean flux difference with respect to phim0
    // from the current flux of the moving phases
    surfaceScalarField dphim(phim0);

    forAll(movingPhaseModels_, movingPhasei)
    {
        dphim -=
            alphafs[movingPhaseModels_[movingPhasei].index()]
           *movingPhaseModels_[movingPhasei].phi();
    }

    forAll(movingPhaseModels_, movingPhasei)
    {
        movingPhaseModels_[movingPhasei].phiRef() += dphim;
    }
}
```

Re-normalise the flux of the phases around the specified mixture mean

##### nHatfv

```cpp
Foam::tmp<Foam::surfaceVectorField> Foam::phaseSystem::nHatfv
(
    const volScalarField& alpha1,
    const volScalarField& alpha2
) const
{
    /*
    // Cell gradient of alpha
    volVectorField gradAlpha =
        alpha2*fvc::grad(alpha1) - alpha1*fvc::grad(alpha2);

    // Interpolated face-gradient of alpha
    surfaceVectorField gradAlphaf = fvc::interpolate(gradAlpha);
    */

    surfaceVectorField gradAlphaf
    (
        fvc::interpolate(alpha2)*fvc::interpolate(fvc::grad(alpha1))
      - fvc::interpolate(alpha1)*fvc::interpolate(fvc::grad(alpha2))
    );

    // Face unit interface normal
    return gradAlphaf/(mag(gradAlphaf) + deltaN_);
}
```

Normal to interface between two phases Used for interface compression

$$
gradAlpha = \alpha_2 \nabla \alpha_1 - \alpha_1 \nabla \alpha_2
$$

$$
nHatfv = \frac{gradAlpha}{|gradAlpha| + deltaN_}
$$

##### nHatf

```cpp
Foam::tmp<Foam::surfaceScalarField> Foam::phaseSystem::nHatf
(
    const volScalarField& alpha1,
    const volScalarField& alpha2
) const
{
    // Face unit interface normal flux
    return nHatfv(alpha1, alpha2) & mesh_.Sf();
}
```

$$
nHatf = nHatfv \cdot \mathbf{S}_f
$$

##### correctContactAngle

```cpp
void Foam::phaseSystem::correctContactAngle
(
    const phaseModel& phase1,
    const phaseModel& phase2,
    surfaceVectorField::Boundary& nHatb
) const
{
    const volScalarField::Boundary& gbf
        = phase1.boundaryField();

    const fvBoundaryMesh& boundary = mesh_.boundary();

    forAll(boundary, patchi)
    {
        if (isA<alphaContactAngleFvPatchScalarField>(gbf[patchi]))
        {
            const alphaContactAngleFvPatchScalarField& acap =
                refCast<const alphaContactAngleFvPatchScalarField>(gbf[patchi]);

            vectorField& nHatPatch = nHatb[patchi];

            vectorField AfHatPatch
            (
                mesh_.Sf().boundaryField()[patchi]
               /mesh_.magSf().boundaryField()[patchi]
            );

            alphaContactAngleFvPatchScalarField::thetaPropsTable::
                const_iterator tp =
                    acap.thetaProps()
                   .find(phasePairKey(phase1.name(), phase2.name()));

            if (tp == acap.thetaProps().end())
            {
                FatalErrorInFunction
                    << "Cannot find interface "
                    << phasePairKey(phase1.name(), phase2.name())
                    << "\n    in table of theta properties for patch "
                    << acap.patch().name()
                    << exit(FatalError);
            }

            bool matched = (tp.key().first() == phase1.name());

            scalar theta0 = degToRad(tp().theta0(matched));
            scalarField theta(boundary[patchi].size(), theta0);

            scalar uTheta = tp().uTheta();

            // Calculate the dynamic contact angle if required
            if (uTheta > small)
            {
                scalar thetaA = degToRad(tp().thetaA(matched));
                scalar thetaR = degToRad(tp().thetaR(matched));

                // Calculated the component of the velocity parallel to the wall
                vectorField Uwall
                (
                    phase1.U()().boundaryField()[patchi].patchInternalField()
                  - phase1.U()().boundaryField()[patchi]
                );
                Uwall -= (AfHatPatch & Uwall)*AfHatPatch;

                // Find the direction of the interface parallel to the wall
                vectorField nWall
                (
                    nHatPatch - (AfHatPatch & nHatPatch)*AfHatPatch
                );

                // Normalise nWall
                nWall /= (mag(nWall) + small);

                // Calculate Uwall resolved normal to the interface parallel to
                // the interface
                scalarField uwall(nWall & Uwall);

                theta += (thetaA - thetaR)*tanh(uwall/uTheta);
            }


            // Reset nHatPatch to correspond to the contact angle

            scalarField a12(nHatPatch & AfHatPatch);

            scalarField b1(cos(theta));

            scalarField b2(nHatPatch.size());

            forAll(b2, facei)
            {
                b2[facei] = cos(acos(a12[facei]) - theta[facei]);
            }

            scalarField det(1 - a12*a12);

            scalarField a((b1 - a12*b2)/det);
            scalarField b((b2 - a12*b1)/det);

            nHatPatch = a*AfHatPatch + b*nHatPatch;

            nHatPatch /= (mag(nHatPatch) + deltaN_.value());
        }
    }
}
```

##### K

```cpp
Foam::tmp<Foam::volScalarField> Foam::phaseSystem::K
(
    const phaseModel& phase1,
    const phaseModel& phase2
) const
{
    tmp<surfaceVectorField> tnHatfv = nHatfv(phase1, phase2);

    correctContactAngle(phase1, phase2, tnHatfv.ref().boundaryFieldRef());

    // Simple expression for curvature
    return -fvc::div(tnHatfv & mesh_.Sf());
}
```

$$
K = \nabla \cdot (nHatfv \cdot \mathbf{S}_f)
$$

##### Constructors 1: phaseSystem

```cpp
Foam::phaseSystem::phaseSystem
(
    const fvMesh& mesh
)
:
    IOdictionary
    (
        IOobject
        (
            "phaseProperties",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),

    mesh_(mesh),

    referencePhaseName_
    (
        // Temporary hack for backward compatibility with
        // reactingTwoPhaseEulerFoam
        lookup<word>("type").find("TwoPhase") != string::npos
      ? lookup<wordList>("phases")[1]
      : lookupOrDefault("referencePhase", word::null)
    ),

    phaseModels_
    (
        lookup("phases"),
        phaseModel::iNew(*this, referencePhaseName_)
    ),

    phi_(calcPhi(phaseModels_)),

    dpdt_
    (
        IOobject
        (
            "dpdt",
            mesh.time().timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar(dimPressure/dimTime, 0)
    ),

    MRF_(mesh_),

    cAlphas_(lookupOrDefault("interfaceCompression", cAlphaTable())),

    deltaN_
    (
        "deltaN",
        1e-8/pow(average(mesh_.V()), 1.0/3.0)
    )
```

+ Define `phaseProperties` ro read phase properties;
+ Define and get `mesh`;
+ Define `referencePhaseName_`;
+ Define and read phase models;
+ Define and calculate `phi`;
+ Define and initialize `dpdt_`;
+ Deifne `MRF_` models;
+ Define `cAlphas_` for interface compression;
+ Define `deltaN` as $\frac{10^{-8}}{V^{1/3}}$.

```cpp
{
    // Groupings
    label movingPhasei = 0;
    label stationaryPhasei = 0;
    label anisothermalPhasei = 0;
    label multiComponentPhasei = 0;
    forAll(phaseModels_, phasei)
    {
        phaseModel& phase = phaseModels_[phasei];
        movingPhasei += !phase.stationary();
        stationaryPhasei += phase.stationary();
        anisothermalPhasei += !phase.isothermal();
        multiComponentPhasei += !phase.pure();
    }
    movingPhaseModels_.resize(movingPhasei);
    stationaryPhaseModels_.resize(stationaryPhasei);
    anisothermalPhaseModels_.resize(anisothermalPhasei);
    multiComponentPhaseModels_.resize(multiComponentPhasei);

    movingPhasei = 0;
    stationaryPhasei = 0;
    anisothermalPhasei = 0;
    multiComponentPhasei = 0;
    forAll(phaseModels_, phasei)
    {
        phaseModel& phase = phaseModels_[phasei];
        if (!phase.stationary())
        {
            movingPhaseModels_.set(movingPhasei++, &phase);
        }
        if (phase.stationary())
        {
            stationaryPhaseModels_.set(stationaryPhasei++, &phase);
        }
        if (!phase.isothermal())
        {
            anisothermalPhaseModels_.set(anisothermalPhasei++, &phase);
        }
        if (!phase.pure())
        {
            multiComponentPhaseModels_.set(multiComponentPhasei++, &phase);
        }
    }

    // Write phi
    phi_.writeOpt() = IOobject::AUTO_WRITE;

    // Blending methods
    forAllConstIter(dictionary, subDict("blending"), iter)
    {
        blendingMethods_.insert
        (
            iter().keyword(),
            blendingMethod::New
            (
                iter().keyword(),
                iter().dict(),
                phaseModels_.toc()
            )
        );
    }

    // Sub-models
    generatePairsAndSubModels("surfaceTension", surfaceTensionModels_);
    generatePairsAndSubModels("aspectRatio", aspectRatioModels_);

    // Update motion fields
    correctKinematics();

    // Set the optional reference phase fraction from the other phases
    if (referencePhaseName_ != word::null)
    {
        phaseModel* referencePhasePtr = &phases()[referencePhaseName_];
        volScalarField& referenceAlpha = *referencePhasePtr;

        referenceAlpha = 1;

        forAll(phaseModels_, phasei)
        {
            if (&phaseModels_[phasei] != referencePhasePtr)
            {
                referenceAlpha -= phaseModels_[phasei];
            }
        }
    }

    forAll(phases(), phasei)
    {
        const volScalarField& alphai = phases()[phasei];
        mesh_.setFluxRequired(alphai.name());
    }
}
```

##### Destructor ~phaseSystem 

```cpp
Foam::phaseSystem::~phaseSystem()
{}
```

##### rho

```cpp
Foam::tmp<Foam::volScalarField> Foam::phaseSystem::rho() const
{
    tmp<volScalarField> rho(movingPhaseModels_[0]*movingPhaseModels_[0].rho());

    for
    (
        label movingPhasei=1;
        movingPhasei<movingPhaseModels_.size();
        movingPhasei++
    )
    {
        rho.ref() +=
            movingPhaseModels_[movingPhasei]
           *movingPhaseModels_[movingPhasei].rho();
    }

    if (stationaryPhaseModels_.empty())
    {
        return rho;
    }
    else
    {
        return rho/sumAlphaMoving();
    }
}
```

$$
\rho = 
$$

##### U

```cpp
Foam::tmp<Foam::volVectorField> Foam::phaseSystem::U() const
{
    tmp<volVectorField> U(movingPhaseModels_[0]*movingPhaseModels_[0].U());

    for
    (
        label movingPhasei=1;
        movingPhasei<movingPhaseModels_.size();
        movingPhasei++
    )
    {
        U.ref() +=
            movingPhaseModels_[movingPhasei]
           *movingPhaseModels_[movingPhasei].U();
    }

    if (stationaryPhaseModels_.empty())
    {
        return U;
    }
    else
    {
        return U/sumAlphaMoving();
    }
}
```

##### E

```cpp
Foam::tmp<Foam::volScalarField>
Foam::phaseSystem::E(const phasePairKey& key) const
{
    if (aspectRatioModels_.found(key))
    {
        return aspectRatioModels_[key]->E();
    }
    else
    {
        return volScalarField::New
        (
            aspectRatioModel::typeName + ":E",
            mesh_,
            dimensionedScalar(dimless, 1)
        );
    }
}
```

Return the aspect-ratio for a pair.

##### sigma1

```cpp
Foam::tmp<Foam::volScalarField>
Foam::phaseSystem::sigma(const phasePairKey& key) const
{
    if (surfaceTensionModels_.found(key))
    {
        return surfaceTensionModels_[key]->sigma();
    }
    else
    {
        return volScalarField::New
        (
            surfaceTensionModel::typeName + ":sigma",
            mesh_,
            dimensionedScalar(surfaceTensionModel::dimSigma, 0)
        );
    }
}
```

Return the surface tension coefficient for a pair

##### sigma2

```cpp
Foam::tmp<Foam::scalarField>
Foam::phaseSystem::sigma(const phasePairKey& key, label patchi) const
{
    if (surfaceTensionModels_.found(key))
    {
        return surfaceTensionModels_[key]->sigma(patchi);
    }
    else
    {
        return tmp<scalarField>
        (
            new scalarField(mesh_.boundary()[patchi].size(), 0)
        );
    }
}
```

Return the surface tension coefficient for a pair on a patch

##### nearInterface

```cpp
Foam::tmp<Foam::volScalarField>
Foam::phaseSystem::nearInterface() const
{
    tmp<volScalarField> tnearInt
    (
        volScalarField::New
        (
            "nearInterface",
            mesh_,
            dimensionedScalar(dimless, 0)
        )
    );

    forAll(phases(), phasei)
    {
        tnearInt.ref() = max
        (
            tnearInt(),
            pos0(phases()[phasei] - 0.01)*pos0(0.99 - phases()[phasei])
        );
    }

    return tnearInt;
}
```

`pos0(x)`: if $x$ greater-equal zero: 1.0 else 0.0

pos0(phases()[phasei] - 0.01)*pos0(0.99 - phases()[phasei])

* if phases()[phasei] > 0.01, then pos0(phases()[phasei] - 0.01) = 1;
* if phases()[phasei] < 0.99, then pos0(0.99 - phases()[phasei]) = 1;
* so if 0.01 < phases()[phasei] < 0.99, it equals to 1

So only if 0.01 < phases()[phasei] < 0.99, tnearInt = 1, otherwise tneatInt = 0

##### dmdtf

```cpp
Foam::tmp<Foam::volScalarField> Foam::phaseSystem::dmdtf
(
    const phasePairKey& key
) const
{
    const phasePair pair
    (
        phaseModels_[key.first()],
        phaseModels_[key.second()]
    );

    return volScalarField::New
    (
        IOobject::groupName("dmdtf", pair.name()),
        mesh(),
        dimensionedScalar(dimDensity/dimTime, 0)
    );
}
```

$$
dmdtf = 0
$$

##### dmdts

```cpp
Foam::PtrList<Foam::volScalarField> Foam::phaseSystem::dmdts() const
{
    return PtrList<volScalarField>(phaseModels_.size());
}
```

##### incompressible

```cpp
bool Foam::phaseSystem::incompressible() const
{
    forAll(phaseModels_, phasei)
    {
        if (!phaseModels_[phasei].incompressible())
        {
            return false;
        }
    }

    return true;
}
```

if incompressible return `True`, else return `False`

##### implicitPhasePressure

```cpp
bool Foam::phaseSystem::implicitPhasePressure(const phaseModel& phase) const
{
    return false;
}


bool Foam::phaseSystem::implicitPhasePressure() const
{
    return false;
}
```

##### surfaceTension

```cpp
Foam::tmp<Foam::surfaceScalarField> Foam::phaseSystem::surfaceTension
(
    const phaseModel& phase1
) const
{
    tmp<surfaceScalarField> tSurfaceTension
    (
        surfaceScalarField::New
        (
            "surfaceTension",
            mesh_,
            dimensionedScalar(dimensionSet(1, -2, -2, 0, 0), 0)
        )
    );

    forAll(phases(), phasej)
    {
        const phaseModel& phase2 = phases()[phasej];

        if (&phase2 != &phase1)
        {
            phasePairKey key12(phase1.name(), phase2.name());

            cAlphaTable::const_iterator cAlpha(cAlphas_.find(key12));

            if (cAlpha != cAlphas_.end())
            {
                tSurfaceTension.ref() +=
                    fvc::interpolate(sigma(key12)*K(phase1, phase2))
                   *(
                        fvc::interpolate(phase2)*fvc::snGrad(phase1)
                      - fvc::interpolate(phase1)*fvc::snGrad(phase2)
                    );
            }
        }
    }

    return tSurfaceTension;
}
```

$$
tSurfaceTension = \sum_{i=2}^N \sigma K (phase_i \nabla phase_1 - phase_1 \nabla phase_i)
$$

##### correct

```cpp
void Foam::phaseSystem::correct()
{
    forAll(phaseModels_, phasei)
    {
        phaseModels_[phasei].correct();
    }
}
```

##### correctContinuityError

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

##### correctKinematics

```cpp
void Foam::phaseSystem::correctKinematics()
{
    bool updateDpdt = false;

    forAll(phaseModels_, phasei)
    {
        phaseModels_[phasei].correctKinematics();

        updateDpdt = updateDpdt || phaseModels_[phasei].thermo().dpdt();
    }

    // Update the pressure time-derivative if required
    if (updateDpdt)
    {
        dpdt_ = fvc::ddt(phaseModels_.begin()().thermo().p());
    }
}
```

##### correctThermo

```cpp
void Foam::phaseSystem::correctThermo()
{
    forAll(phaseModels_, phasei)
    {
        phaseModels_[phasei].correctThermo();
    }
}
```

##### correctReactions

```cpp
void Foam::phaseSystem::correctReactions()
{
    forAll(phaseModels_, phasei)
    {
        phaseModels_[phasei].correctReactions();
    }
}
```

##### correctSpecies

```cpp
void Foam::phaseSystem::correctSpecies()
{
    forAll(phaseModels_, phasei)
    {
        phaseModels_[phasei].correctSpecies();
    }
}
```

##### correctTurbulelnce

```cpp
void Foam::phaseSystem::correctTurbulence()
{
    forAll(phaseModels_, phasei)
    {
        phaseModels_[phasei].correctTurbulence();
    }
}
```

##### correctEnergyTransport

```cpp
void Foam::phaseSystem::correctEnergyTransport()
{
    forAll(phaseModels_, phasei)
    {
        phaseModels_[phasei].correctEnergyTransport();
    }
}
```

##### read

```cpp
bool Foam::phaseSystem::read()
{
    if (regIOobject::read())
    {
        bool readOK = true;

        forAll(phaseModels_, phasei)
        {
            readOK &= phaseModels_[phasei].read();
        }

        // models ...

        return readOK;
    }
    else
    {
        return false;
    }
}
```

##### byDt1

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
```

$$
byDt = \Delta t v_f
$$

or

$$
byDt = \frac{v_f}{\Delta t}
$$

##### byDt2

```cpp
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
byDt = \Delta t_f s_f
$$

or

$$
byDt = \frac{s_f}{\Delta t_f}
$$

#### phaseSystemSolve.C

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



#### pEqn.H

### PUf











