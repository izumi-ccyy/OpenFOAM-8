# phaseSystem

- [phaseSystem](#phasesystem)
  - [phaseSystem.H](#phasesystemh)
    - [include](#include)
    - [public typedef](#public-typedef)
    - [protected typedefs](#protected-typedefs)
    - [protected data](#protected-data)
  - [phaseSystemI.H](#phasesystemih)
      - [phaseSystemNew.C](#phasesystemnewc)
      - [phaseSystemTemplates.C](#phasesystemtemplatesc)
        - [addfield 1](#addfield-1)
      - [phaseSystem.C](#phasesystemc)
        - [calcPhi](#calcphi)
        - [generatePairs](#generatepairs)
        - [sumAlphaMoving](#sumalphamoving)
        - [setMixtureU](#setmixtureu)
        - [setMixturePhi](#setmixturephi)
        - [nHatfv](#nhatfv)
        - [nHatf](#nhatf)
        - [correctContactAngle (uncompleted)](#correctcontactangle-uncompleted)
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
        - [Definitions to solve the equation](#definitions-to-solve-the-equation)
        - [The loop to solve the equation](#the-loop-to-solve-the-equation)
          - [Source](#source)
          - [Controls](#controls)
          - [Solving loop](#solving-loop)

## phaseSystem.H

### include

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

### public typedef

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

### protected typedefs

```cpp
protected:

    // Protected typedefs

        typedef
            HashTable<dictionary, phasePairKey, phasePairKey::hash>
            dictTable;

        typedef
            HashTable<autoPtr<blendingMethod>, word, word::hash>
            blendingMethodTable;

        typedef
            HashTable
            <
                autoPtr<surfaceTensionModel>,
                phasePairKey,
                phasePairKey::hash
            >
            surfaceTensionModelTable;

        typedef
            HashTable
            <
                autoPtr<aspectRatioModel>,
                phasePairKey,
                phasePairKey::hash
            >
            aspectRatioModelTable;

        typedef HashTable<scalar, phasePairKey, phasePairKey::hash>
            cAlphaTable;
```

### protected data



## phaseSystemI.H

define some simple inline functions, such as:

```cpp
inline const Foam::fvMesh& Foam::phaseSystem::mesh() const
{
    return mesh_;
}
```

The following functions are defined here:

* mesh()
* phases()
* movingPhases()
* stationaryPhases()
* anisothermalPhases()
* multiComponentPhases()
* phasePairs()
* otherPhase(const phaseModel& phase)
* phi()
* dpdt()
* MRF()
* fvOptions()

#### phaseSystemNew.C

```cpp
Foam::autoPtr<Foam::phaseSystem> Foam::phaseSystem::New
(
    const fvMesh& mesh
)
{
    const word phaseSystemType
    (
        IOdictionary
        (
            IOobject
            (
                propertiesName,
                mesh.time().constant(),
                mesh,
                IOobject::MUST_READ_IF_MODIFIED,
                IOobject::NO_WRITE,
                false
            )
        ).lookup("type")
    );

    Info<< "Selecting phaseSystem "
        << phaseSystemType << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(phaseSystemType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown phaseSystemType type "
            << phaseSystemType << endl << endl
            << "Valid phaseSystem types are : " << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return cstrIter()(mesh);
}
```

create a new `phaseSystem` containing one or more phases types

#### phaseSystemTemplates.C

define some protected member function:

* createSubModels()
* generatePairsAndSubModels()

define some member functions:

* fillFields()
* foundSubModel()
* lookupSubModel()
* foundBlendedSubModel()
* lookupBlendedSubModel()

define some inline function

* addField()

##### addfield 1

```cpp
template<class GeoField, class Group>
inline void addField
(
    const Group& group,
    const word& name,
    tmp<GeoField> field,
    PtrList<GeoField>& fieldList
)
{
    if (fieldList.set(group.index()))
    {
        fieldList[group.index()] += field;
    }
    else
    {
        fieldList.set
        (
            group.index(),
            new GeoField
            (
                IOobject::groupName(name, group.name()),
                field
            )
        );
    }
}
```

* if `group` already exists in `fieldList`, then
  * `fieldList[group.index()] += field`
* else, create `group` in `fieldList`
  * with `name`, and `field`

* the first parameter `group` is to provide index of the new field
* the second parameter `name` is to provide the name of the new field
* the third parameter `field` is to provide the field of the new field
* the fouth parameter `fieldList` is to provide the list where the new field is added to

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

calculate and return the mixture flux

* define `tmpPhi`
  * $$tmpPhi = \phi = \sum_{k = 1}^N (\alpha^k)_f \phi^k_f$$

where $N$ is the phase number

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

for the interfering phases

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

* define `sumAlphaMoving`
* $$sumAlphaMoving = \sum_{k = 1}^N \alpha^k$$

where $N$ is the number of moving phases

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
dUm = U_{m0} - \sum_{k=1}^N \alpha^k\mathbf{U}^k
$$

$$
\mathbf{U}^k = \mathbf{U}^k + dUm
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

* set `dphim = phim0`
* $$dphim = phim0 - \sum_{k=1}^N \alpha^k_f \phi^k$$
  * where $N$ is the moving phase number
* $$\phi^k = \phi_k + dphim$$

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
gradAlphaf = (\alpha_2)_f \nabla (\alpha_1)_f - (\alpha_1)_f \nabla (\alpha_2)_f
$$

$$
nHatfv = \frac{gradAlphaf}{\|gradAlphaf\| + deltaN\_}
$$

$$
nHatfv = \frac{(\alpha_2)_f \nabla (\alpha_1)_f - (\alpha_1)_f \nabla (\alpha_2)_f}{\|(\alpha_2)_f \nabla (\alpha_1)_f - (\alpha_1)_f \nabla (\alpha_2)_f\| + deltaN\_}
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
nHatf = nHatfv \cdot \mathbf{S}_f = \frac{(\alpha_2)_f \nabla (\alpha_1)_f - (\alpha_1)_f \nabla (\alpha_2)_f}{\|(\alpha_2)_f \nabla (\alpha_1)_f - (\alpha_1)_f \nabla (\alpha_2)_f\| + deltaN\_} \cdot \mathbf{S}_f
$$

##### correctContactAngle (uncompleted)

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

* get boundary of `phase1`, `gbf`
* get `boundary` as mesh boundary
* for every patch in boundaries:
  * if current patchi is alphaContackAngle...
    * get `nHatPatch` as `nHat` on patchi
    * get `AfHatPatch` as
      * $$AfHatPatch = \frac{\mathbf{S}_f}{\|\mathbf{S}_f\|}$$

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

Curvature of interface between two phases, used for interface compression

$$
K = \nabla \cdot (nHatfv \cdot \mathbf{S}_f) = \nabla \cdot \left(\frac{(\alpha_2)_f \nabla (\alpha_1)_f - (\alpha_1)_f \nabla (\alpha_2)_f}{\|(\alpha_2)_f \nabla (\alpha_1)_f - (\alpha_1)_f \nabla (\alpha_2)_f\| + deltaN\_} \cdot \mathbf{S}_f\right)
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

return the mixture density

if there is no stationary phases:

$$
\rho = \sum_{k = 1}^{N} \alpha^k \rho^k
$$

else

$$
\rho = \frac{\sum_{k = 1}^{N} \alpha^k \rho^k}{\sum_{l=1}^{M} \alpha^l}
$$

where $M$ is the number of moving phases

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

return the mixture velocity

if there is no stationary phases:

$$
\mathbf{U} = \sum_{k = 1}^{N} \alpha^k \mathbf{U}^k
$$

else

$$
\mathbf{U} = \frac{\sum_{k = 1}^{N} \alpha^k \mathbf{U}^k}{\sum_{l=1}^{M} \alpha^l}
$$

where $M$ is the number of moving phases

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

Return the aspect-ratio $E_{k,l}$ for a pair of $k$th phase and $l$th phase.

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

Return the surface tension coefficient $\sigma_{k,l}$ for a pair $k$th phase and $l$th phase

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

Return the surface tension coefficient $\sigma_{k,l}$ for a pair $k$th phase and $l$th phase on a patch

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
return the mass transfer rate for an interface

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

return a mass transfer rates for each phase

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

is the phase pressure treated implicit

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
tSurfaceTension = \sum_{i=2}^N \sigma_{1,i} K_{1,i} (phase_i \nabla phase_1 - phase_1 \nabla phase_i)
$$

or

$$
tSurfaceTension = \sum_{i=2}^N \sigma_{1,i} K_{1,i} ((\alpha^i)_f (\nabla \alpha^1)_f - (\alpha^1)_f (\nabla \alpha^i)_f)
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
byDt = \frac{v_f}{\Delta t}
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
byDt = \frac{s_f}{\Delta t_f}
$$

or

$$
byDt = \frac{s_f}{\Delta t_f}
$$

#### phaseSystemSolve.C

##### Definitions to solve the equation

```cpp
void Foam::phaseSystem::solve
(
    const PtrList<volScalarField>& rAUs,
    const PtrList<surfaceScalarField>& rAUfs
)
{
    ...
}
```

define the solve function

```cpp
    const dictionary& alphaControls = mesh_.solverDict("alpha");

    const label nAlphaSubCycles(alphaControls.lookup<label>("nAlphaSubCycles"));
    const label nAlphaCorr(alphaControls.lookup<label>("nAlphaCorr"));

    const bool LTS = fv::localEulerDdt::enabled(mesh_);

    // Optional reference phase which is not solved for
    // but obtained from the sum of the other phases
    phaseModel* referencePhasePtr = nullptr;

    // The phases which are solved
    // i.e. the moving phases less the optional reference phase
    phaseModelPartialList solvePhases;
```

define and get some control variables:

* get `alphaControls` from dictionary
* define and get `nAlphaSubCycles` and `nAlphaCorr`
* get `LTS`
* get reference phase pointer `referencePhasePtr`
  * Optional reference phase which is not solved for but obtained from the sum of the other phases
* get the list of phase that are to be solved `solvePhases`

```cpp
    if (referencePhaseName_ != word::null)
    {
        referencePhasePtr = &phases()[referencePhaseName_];

        solvePhases.setSize(movingPhases().size() - 1);
        label solvePhasesi = 0;
        forAll(movingPhases(), movingPhasei)
        {
            if (&movingPhases()[movingPhasei] != referencePhasePtr)
            {
                solvePhases.set(solvePhasesi++, &movingPhases()[movingPhasei]);
            }
        }
    }
    else
    {
        solvePhases = movingPhases();
    }
```

get `solvePhases` from `movingPhases`

* if there is a reference phase
  * reduce 1 from the size of solved phases (because reference phase is not solved)
  * for every moving phase:
    * if current moving phase is not the reference phase, then
      * set current moving phase as `slovePhases`
* else
  * set `solvePhases` as `movingPhases`

```cpp
    // The phases included in the flux sum limit
    // which is all moving phases if the number of solved phases is > 1
    // otherwise it is just the solved phases
    // as the flux sum limit is not needed in this case
    phaseModelPartialList fluxPhases;
    if (solvePhases.size() == 1)
    {
        fluxPhases = solvePhases;
    }
    else
    {
        fluxPhases = movingPhases();
    }
```

* The phases included in the flux sum limit which is all moving phases if the number of solved phases is > 1
* otherwise it is just the solved phases as the flux sum limit is not needed in this case

```cpp
    forAll(phases(), phasei)
    {
        phases()[phasei].correctBoundaryConditions();
    }
```

* for every phase:
  * correct boundary conditions

```cpp
    PtrList<surfaceScalarField> alphaPhiDbyA0s(phases().size());
    if (implicitPhasePressure() && (rAUs.size() || rAUfs.size()))
    {
        const PtrList<surfaceScalarField> DByAfs(this->DByAfs(rAUs, rAUfs));

        forAll(solvePhases, solvePhasei)
        {
            phaseModel& phase = solvePhases[solvePhasei];
            volScalarField& alpha = phase;

            alphaPhiDbyA0s.set
            (
                phase.index(),
                DByAfs[phase.index()]
               *fvc::snGrad(alpha, "bounded")*mesh_.magSf()
            );
        }
    }
```

define a pointer list of `surfaceScalarField` variables with the size of phase size, `alphaPhiDbyA0s`

* if phase pressure is treated implicitly and `rAUs` or `rAUfs` exists:
  * define a pointer list of `surfaceScalarField` variables `DByAfs`
  * for every phase to be solved (solvePhases):
    * define and get $\alpha$
    * set `alphaPhiDbyA0s` as
      * $$alphaPhiDbyA0^k = DByAfs^k (\nabla \alpha^k)_f \|\mathbf{S}_f\|$$

`DByAfs()` can be found in `applications\solvers\multiphase\multiphaseEulerFoam\phaseSystems\phaseSystem\phaseSystem.H`, is to return the phase diffusivity divided by the momentum central coefficient

```cpp
    // Calculate the void fraction
    volScalarField alphaVoid
    (
        IOobject
        (
            "alphaVoid",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedScalar(dimless, 1)
    );
    forAll(stationaryPhases(), stationaryPhasei)
    {
        alphaVoid -= stationaryPhases()[stationaryPhasei];
    }
```

* define `alphaVoid() = 1` as the void fraction
* for every stationary phases:
  * $$alphaVoid - alphaVoid - \alpha_{stationary, k}$$
  * namely
    * $$alphaVoid = 1 - \sum_{k = 1}^N \alpha_{stationary, k}$$

```cpp
    bool dilatation = false;
    forAll(fluxPhases, fluxPhasei)
    {
        if (fluxPhases[fluxPhasei].divU().valid())
        {
            dilatation = true;
            break;
        }
    }
```

* define `dilatation` as `false`
* for every flux phase:
  * if there is a valid $\nabla \cdot \mathbf{U}$ of flux phases, in other words, $\nabla \cdot \mathbf{U} \neq 0$, then
    * `dilatation` is `true` 

##### The loop to solve the equation

```cpp
    for (int acorr=0; acorr<nAlphaCorr; acorr++)
    {
        ...
    }
```

* for $acorr < nAlphaCorr$, to solve $\alpha$ equation

###### Source

```cpp
        PtrList<volScalarField::Internal> Sps(phases().size());
        PtrList<volScalarField::Internal> Sus(phases().size());

        forAll(fluxPhases, fluxPhasei)
        {
            phaseModel& phase = fluxPhases[fluxPhasei];
            volScalarField& alpha = phase;
            const label phasei = phase.index();

            Sps.set
            (
                phasei,
                new volScalarField::Internal
                (
                    IOobject
                    (
                        "Sp",
                        mesh_.time().timeName(),
                        mesh_
                    ),
                    mesh_,
                    dimensionedScalar(dimless/dimTime, 0)
                )
            );

            Sus.set
            (
                phasei,
                new volScalarField::Internal
                (
                    "Su",
                    min(alpha, scalar(1))
                    *fvc::div(fvc::absolute(phi_, phase.U()))
                )
            );

            if (dilatation)
            {
                // Construct the dilatation rate source term
                volScalarField::Internal dgdt
                (
                    volScalarField::Internal::New
                    (
                        "dgdt",
                        mesh_,
                        dimensionedScalar(dimless/dimTime, 0)
                    )
                );

                forAll(phases(), phasej)
                {
                    const phaseModel& phase2 = phases()[phasej];
                    const volScalarField& alpha2 = phase2;

                    if (&phase2 != &phase)
                    {
                        if (phase.divU().valid())
                        {
                            dgdt += alpha2()*phase.divU()()();
                        }

                        if (phase2.divU().valid())
                        {
                            dgdt -= alpha()*phase2.divU()()();
                        }
                    }
                }

                volScalarField::Internal& Sp = Sps[phasei];
                volScalarField::Internal& Su = Sus[phasei];

                forAll(dgdt, celli)
                {
                    if (dgdt[celli] > 0)
                    {
                        Sp[celli] -= dgdt[celli]/max(1 - alpha[celli], 1e-4);
                        Su[celli] += dgdt[celli]/max(1 - alpha[celli], 1e-4);
                    }
                    else if (dgdt[celli] < 0)
                    {
                        Sp[celli] += dgdt[celli]/max(alpha[celli], 1e-4);
                    }
                }
            }
        }
```

* define a pointer list of `volScalarField` variables internal fields, `Sps` and `Sus`
* for every flux phases:
  * define and get $\alpha$
  * get index of current phase `phasei`
  * initialize `Sps` and `Sus`:
    * $Sp^k = 0$
    * $Su^k = \min(\alpha^k, 1) \nabla \mathbf{U}^k_{absolute}$
    * if `dilatation` is `true`, namely for one phase, $\nabla \cdot \mathbf{U} \neq 0$
      * initialize `dgdt`:
        * $$dgdt = 0$$
      * for every phase:
        * define and get current phase (`phase2`) volume fraction as $\alpha_2$
        * if `phase2` is not `phase`
          * if $\nabla \cdot \mathbf{U} \neq 0$ for `phase`, then
            * $$dgdt = dgdt + \alpha_2 \nabla \cdot \mathbf{U}^k $$
          * if $\nabla \cdot \mathbf{U} \neq 0$ for `phase2`, then
            * $$dgdt = dgdt - \alpha \nabla \cdot \mathbf{U}^l$$
        * namely,
          * $$dgdt = \sum_{k = 1}^n \sum_{l = 1, l \neq k}^N \alpha_l \nabla \cdot \mathbf{U}^k - \sum_{k = 1}^N \sum_{l = 1, l \neq k}^n \alpha_k \nabla \cdot \mathbf{U}^l$$
          * where, $n$ is the number of phases whose $\nabla \cdot \mathbf{U} \neq 0$ while $N$ is the total number of flux phases
      * get $Sp$ and $Su$ for current phase
        * $Sp^k$ and $Su^k$
      * for $dgdt$ of every grid cell $i$:
        * if $dgdt[i] > 0$
          * $$Sp^k[i] = Sp^k[i] - \frac{dgdt[i]}{\max(1-\alpha^k[i], 10^{-4})}$$
          * $$Su^k[i] = Su^k[i] + \frac{dgdt[i]}{\max(1-\alpha^k[i], 10^{-4})}$$
        * else, if $dgdt[i] < 0$
          * $$Sp^k[i] = Sp^k[i] + \frac{dgdt[i]}{\max(\alpha^k[i], 10^{-4})}$$

###### Controls

```cpp
        tmp<volScalarField> trSubDeltaT;

        if (LTS && nAlphaSubCycles > 1)
        {
            trSubDeltaT =
                fv::localEulerDdt::localRSubDeltaT(mesh_, nAlphaSubCycles);
        }

        List<volScalarField*> alphaPtrs(phases().size());
        forAll(phases(), phasei)
        {
            alphaPtrs[phasei] = &phases()[phasei];
        }
```

* define `trSubDeltaT`
* if LTS and `nAlphaSubCycles` > 1, then
  * $$trSubDeltaT = \frac{nAlphaSubCycles}{\Delta t}$$
* define a list of pointer for `volScalarField` variables $\alpha$, `alphaPtrs`
* for every phase:
  * get the pointer for $\alpha$

`localRSubDeltaT` can be found in `src\finiteVolume\finiteVolume\ddtSchemes\localEulerDdtScheme\localEulerDdtScheme.H` and `src\finiteVolume\finiteVolume\ddtSchemes\localEulerDdtScheme\localEulerDdt.C`, is to calculate and return the reciprocal of the local sub-cycling time-step

```cpp
Foam::tmp<Foam::volScalarField> Foam::fv::localEulerDdt::localRSubDeltaT
(
    const fvMesh& mesh,
    const label nAlphaSubCycles
)
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            rSubDeltaTName,
            nAlphaSubCycles
           *mesh.objectRegistry::lookupObject<volScalarField>
            (
                rDeltaTName
            )
        )
    );
}
```

$$
localRSubDeltaT = \frac{nAlphaSubCycles}{\Delta t}
$$

###### Solving loop

```cpp
        for
        (
            subCycle<volScalarField, subCycleFields> alphaSubCycle
            (
                alphaPtrs,
                nAlphaSubCycles
            );
            !(++alphaSubCycle).end();
        )
        {
            ...
        }
```

solve the equations in the subcycle

Generate face-alphas

```cpp
            // Generate face-alphas
            PtrList<surfaceScalarField> alphafs(phases().size());
            if (solvePhases.size() > 1)
            {
                forAll(phases(), phasei)
                {
                    phaseModel& phase = phases()[phasei];
                    alphafs.set
                    (
                        phasei,
                        new surfaceScalarField
                        (
                            IOobject::groupName("alphaf", phase.name()),
                            upwind<scalar>(mesh_, phi_).interpolate(phase)
                        )
                    );
                }
            }
```

* define a pointer list of $\alpha$ on surface, $\alpha_f$
* if the number of phases to be solved is larger than 1, then:
  * for every phase:
    * get $k$th `phase`
    * set `alphafs` using interpolate with upwind schemes:
      * $\alpha_f^k$

Create correction fluxes

```cpp
            // Create correction fluxes
            PtrList<surfaceScalarField> alphaPhiCorrs(phases().size());

            if (solvePhases.size() > 1)
            {
                forAll(stationaryPhases(), stationaryPhasei)
                {
                    phaseModel& phase = stationaryPhases()[stationaryPhasei];

                    alphaPhiCorrs.set
                    (
                        phase.index(),
                        new surfaceScalarField
                        (
                            IOobject::groupName("alphaPhiCorr", phase.name()),
                          - upwind<scalar>(mesh_, phi_).flux(phase)
                        )
                    );
                }
            }
```

* defien correction fluxes `alphaPhiCorrs`
* if the number of phases to be solved is larger than 1, then:
  * for every stationary phases
    * get stationary phase
    * set `alphaPhiCorrs` with upwind schemes
      * $$alphaPhiCorrs = -\phi$$

```cpp
            forAll(fluxPhases, fluxPhasei)
            {
                phaseModel& phase = fluxPhases[fluxPhasei];
                volScalarField& alpha = phase;

                alphaPhiCorrs.set
                (
                    phase.index(),
                    new surfaceScalarField
                    (
                        IOobject::groupName("alphaPhiCorr", phase.name()),
                        fvc::flux(phi_, alpha, "div(phi," + alpha.name() + ')')
                    )
                );

                surfaceScalarField& alphaPhiCorr = alphaPhiCorrs[phase.index()];

                forAll(phases(), phasei)
                {
                    phaseModel& phase2 = phases()[phasei];
                    volScalarField& alpha2 = phase2;

                    if (&phase2 == &phase) continue;

                    surfaceScalarField phir(phase.phi() - phase2.phi());

                    cAlphaTable::const_iterator cAlpha
                    (
                        cAlphas_.find(phasePairKey(phase.name(), phase2.name()))
                    );

                    if (cAlpha != cAlphas_.end())
                    {
                        surfaceScalarField phic
                        (
                            (mag(phi_) + mag(phir))/mesh_.magSf()
                        );

                        phir +=
                            min(cAlpha()*phic, max(phic))
                           *nHatf(alpha, alpha2);
                    }

                    word phirScheme
                    (
                        "div(phir," + alpha2.name() + ',' + alpha.name() + ')'
                    );

                    alphaPhiCorr += fvc::flux
                    (
                       -fvc::flux(-phir, alpha2, phirScheme),
                        alpha,
                        phirScheme
                    );
                }

                if (alphaPhiDbyA0s.set(phase.index()))
                {
                    alphaPhiCorr +=
                        fvc::interpolate(max(alpha, scalar(0)))
                       *fvc::interpolate(max(1 - alpha, scalar(0)))
                       *alphaPhiDbyA0s[phase.index()];
                }

                phase.correctInflowOutflow(alphaPhiCorr);

                MULES::limit
                (
                    geometricOneField(),
                    alpha,
                    phi_,
                    alphaPhiCorr,
                    Sps[phase.index()],
                    Sus[phase.index()],
                    min(alphaVoid.primitiveField(), phase.alphaMax())(),
                    zeroField(),
                    true
                );
            }
```

* for every `fluxPhase`:
  * define and get $\alpha^k$
  * set flux correction `alphaPhiCorrs`
    * $$alphaPhiCorrs = \phi$$
  * get `alphaPhiCorrs` fro current phase
    * $$alphaPhiCorr = alphaPhiCorr^k = \phi^k$$
  * for every phase:
    * get `phase2`, indexed by `l`
    * get `alpha2`, $\alpha^l$
    * if $l = k$, namely `phase2` is `phase`, then
      * stop this step and enter next step of the loop
    * if $l \neq k$, namely, different phases
      * define `phir` as
        * $$phir = \phi^k - \phi^l$$
      * get `cAlpha` from the list with the name of `phase` and `phase2`
      * if `cAlpha` is not end, then
        * defien `phic`
          * $$phic = \frac{\|\phi\| + \|phir\|}{\|\mathbf{S}_f\|} = \frac{\|\phi\| + \|\phi^k - \phi^l\|}{\|\mathbf{S}_f\|}$$
        * calculate `phir`
          * $$phir = phir + \min(cAlpha phic, \max(phic)) nHatf$$
          * $$phir = \phi^k - \phi^l + \min\left(c_{\alpha} \frac{\|\phi\| + \|\phi^k - \phi^l\|}{\|\mathbf{S}_f\|}, \max(\frac{\|\phi\| + \|\phi^k - \phi^l\|}{\|\mathbf{S}_f\|})\right) \times \left[\frac{(\alpha_2)_f \nabla (\alpha_1)_f - (\alpha_1)_f \nabla (\alpha_2)_f}{\|(\alpha_2)_f \nabla (\alpha_1)_f - (\alpha_1)_f \nabla (\alpha_2)_f\| + deltaN\_} \cdot \mathbf{S}_f\right]$$
      * define word `phirScheme`
      * calculate `alphaPhiCorr` as
        * $$alphaPhiCorr = alphaPhiCorr + phir = \phi + \phi_r$$
  * if $alphaPhiDbyA0^k$ is set, then
    * $$alphaPhiCorr = alphaPhiCorr + (\max(\alpha^k, 0))_f \cdot (\max(1-\alpha^k, 0))_f \cdot alphaPhiDbyA0^k = \phi + \phi_r  + (\max(\alpha^k, 0))_f \cdot (\max(1-\alpha^k, 0))_f \cdot DByAfs^k (\nabla \alpha^k)_f \|\mathbf{S}_f\|$$
  * correct flow using `alphaPhiCorr`
  * define `MULES` as limit

`cAlpha`, $c_{\alpha}$ is the Compression coefficient

`end()` can be found in `applications\solvers\multiphase\multiphaseEulerFoam\phaseSystems\phasePair\phasePair\phasePairI.H`, is to return const_iterator set to beyond the end of the pair

`phi_` can be found in `applications\solvers\multiphase\multiphaseEulerFoam\phaseSystems\phaseSystem\phaseSystem.H`, is the total volumetric flux

`flux` can be found in `src\finiteVolume\finiteVolume\fvc\fvcFluxTemplates.C`, is to reutrn flux

```cpp
template<class Type>
tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
flux
(
    const surfaceScalarField& phi,
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    Istream& schemeData
)
{
    return fv::convectionScheme<Type>::New
    (
        vf.mesh(),
        phi,
        schemeData
    )().flux(phi, vf);
}
```

the second parameter is to provide mesh information

limit flux sums

```cpp
            if (solvePhases.size() > 1)
            {
                // Limit the flux sums, fixing those of the stationary phases
                labelHashSet fixedAlphaPhiCorrs;
                forAll(stationaryPhases(), stationaryPhasei)
                {
                    fixedAlphaPhiCorrs.insert
                    (
                        stationaryPhases()[stationaryPhasei].index()
                    );
                }
                MULES::limitSum(alphafs, alphaPhiCorrs, fixedAlphaPhiCorrs);
            }
```

* if the number of solved phases is larger than 1:
  * define `fixedAlphaPhiCorrs`
  * for every stationary phase:
    * insert stationary phases into `fixedAlphaPhiCorrs`
  * using `MULES` to limit

solve

```cpp
            forAll(solvePhases, solvePhasei)
            {
                phaseModel& phase = solvePhases[solvePhasei];
                volScalarField& alpha = phase;

                surfaceScalarField& alphaPhi = alphaPhiCorrs[phase.index()];
                alphaPhi += upwind<scalar>(mesh_, phi_).flux(phase);
                phase.correctInflowOutflow(alphaPhi);

                MULES::explicitSolve
                (
                    geometricOneField(),
                    alpha,
                    alphaPhi,
                    Sps[phase.index()],
                    Sus[phase.index()]
                );

                if (alphaSubCycle.index() == 1)
                {
                    phase.alphaPhiRef() = alphaPhi;
                }
                else
                {
                    phase.alphaPhiRef() += alphaPhi;
                }
            }
```

* for every phases to be solved
  * get $\alpah^k$
  * define $alphaPhi = alphaPhiCorr^k$, then
  * $$alphaPhi = alphaPhi + \phi$$
  *  correct flow with `alphaPhi`
  *  use `MULES` to solve
  *  if `alphaSubCycle.index()` = 1,
     * $$phase.alphaPhiRef() = alphaPhi$$
  * else
    * $$phase.alphaPhiRef() = phase.alphaPhiRef() + alphaPhi$$

```cpp
            if (implicitPhasePressure() && (rAUs.size() || rAUfs.size()))
            {
                const PtrList<surfaceScalarField> DByAfs
                (
                    this->DByAfs(rAUs, rAUfs)
                );

                forAll(solvePhases, solvePhasei)
                {
                    phaseModel& phase = solvePhases[solvePhasei];
                    volScalarField& alpha = phase;

                    const surfaceScalarField alphaDbyA
                    (
                        fvc::interpolate(max(alpha, scalar(0)))
                       *fvc::interpolate(max(1 - alpha, scalar(0)))
                       *DByAfs[phase.index()]
                    );

                    fvScalarMatrix alphaEqn
                    (
                        fvm::ddt(alpha) - fvc::ddt(alpha)
                      - fvm::laplacian(alphaDbyA, alpha, "bounded")
                    );

                    alphaEqn.solve();

                    phase.alphaPhiRef() += alphaEqn.flux();
                }
            }
```

* if phase pressure is treat implicitly and `rAUs` or `rAUfs` exists:
  * define and get `DByAfs`
  * for every phase to be solved:
    * get $\alpha^k$
    * define `alphaDbyA` as
      * $$alphaDbyA = (\max(\alpha^k, 0))_f \cdot (\max(1-\alpha^k, 0))_f \cdot DByAf^k$$
    * define `alphaEqn` as
      * $$\left(\frac{\partial \alpha^k}{\partial t}\right)_{implicit} - \left(\frac{\partial \alpha^k}{\partial t}\right)_{explicit} - \nabla \cdot (alphaDbyA \cdot \nabla \alpha)$$
    * solve `alphaEqn`
    * $$alphaPhiRef() = alphaPhiRef() + alphaEqn.flux()$$

```cpp
            // Report the phase fractions and the phase fraction sum
            forAll(solvePhases, solvePhasei)
            {
                phaseModel& phase = solvePhases[solvePhasei];

                Info<< phase.name() << " fraction, min, max = "
                    << phase.weightedAverage(mesh_.V()).value()
                    << ' ' << min(phase).value()
                    << ' ' << max(phase).value()
                    << endl;
            }
```

* print the phase fractions of phases to be solved

```cpp
            if (!referencePhasePtr)
            {
                volScalarField sumAlphaMoving
                (
                    IOobject
                    (
                        "sumAlphaMoving",
                        mesh_.time().timeName(),
                        mesh_
                    ),
                    mesh_,
                    dimensionedScalar(dimless, 0)
                );
                forAll(movingPhases(), movingPhasei)
                {
                    sumAlphaMoving += movingPhases()[movingPhasei];
                }

                Info<< "Phase-sum volume fraction, min, max = "
                    << (sumAlphaMoving + 1 - alphaVoid)()
                      .weightedAverage(mesh_.V()).value()
                    << ' ' << min(sumAlphaMoving + 1 - alphaVoid).value()
                    << ' ' << max(sumAlphaMoving + 1 - alphaVoid).value()
                    << endl;

                // Correct the sum of the phase fractions to avoid drift
                forAll(movingPhases(), movingPhasei)
                {
                    movingPhases()[movingPhasei] *= alphaVoid/sumAlphaMoving;
                }
            }
```