# dragModel

- [dragModel](#dragmodel)
  - [dragModel](#dragmodel-1)
    - [dragModel.H](#dragmodelh)
      - [include](#include)
      - [protected](#protected)
      - [public](#public)
        - [constructor and destructor](#constructor-and-destructor)
        - [member functions](#member-functions)
    - [dragModelNew.C](#dragmodelnewc)
    - [dragModel.C](#dragmodelc)
      - [include and static data members](#include-and-static-data-members)
      - [constructor and destructor](#constructor-and-destructor-1)
      - [member functions](#member-functions-1)
        - [Ki()](#ki)
        - [K()](#k)
        - [Kf()](#kf)
        - [writeDta()](#writedta)
  - [SchillerNaumann](#schillernaumann)
    - [SchillerNaumann.H](#schillernaumannh)
      - [CdRe()](#cdre)
    - [SchillerNaumann.C](#schillernaumannc)
      - [CdRe()](#cdre-1)
  - [IshiiZuber](#ishiizuber)
    - [IshiiZuber.H](#ishiizuberh)
      - [CdRe()](#cdre-2)
    - [IshiiZuber.C](#ishiizuberc)
      - [CdRe()](#cdre-3)

## dragModel

### dragModel.H

#### include

```cpp
#ifndef dragModel_H
#define dragModel_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "volFields.H"
#include "dictionary.H"
#include "runTimeSelectionTables.H"

namespace Foam
{

class phasePair;
class swarmCorrection;
```

#### protected

```cpp
protected:

    // Protected data

        //- Phase pair
        const phasePair& pair_;

        //- Swarm correction
        autoPtr<swarmCorrection> swarmCorrection_;
```

define pair_ and swarmCorrection_

#### public

##### constructor and destructor

```cpp
public:

    //- Runtime type information
    TypeName("dragModel");


    // Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            dragModel,
            dictionary,
            (
                const dictionary& dict,
                const phasePair& pair,
                const bool registerObject
            ),
            (dict, pair, registerObject)
        );


    // Static Data Members

        //- Coefficient dimensions
        static const dimensionSet dimK;


    // Constructors

        // Construct without residual constants
        dragModel
        (
            const phasePair& pair,
            const bool registerObject
        );

        // Construct with residual constants
        dragModel
        (
            const dictionary& dict,
            const phasePair& pair,
            const bool registerObject
        );


    //- Destructor
    virtual ~dragModel();


    // Selectors

        static autoPtr<dragModel> New
        (
            const dictionary& dict,
            const phasePair& pair
        );

```

define name, runtime selection table, static data member (coefficient dimensions), constructor, destructor, selector

##### member functions

```cpp
    // Member Functions

        //- Drag coefficient
        virtual tmp<volScalarField> CdRe() const = 0;

        //- Return the phase-intensive drag coefficient Ki
        //  used in the momentum equations
        //    ddt(alpha1*rho1*U1) + ... = ... alphad*K*(U1-U2)
        //    ddt(alpha2*rho2*U2) + ... = ... alphad*K*(U2-U1)
        virtual tmp<volScalarField> Ki() const;

        //- Return the drag coefficient K
        //  used in the momentum equations
        //    ddt(alpha1*rho1*U1) + ... = ... K*(U1-U2)
        //    ddt(alpha2*rho2*U2) + ... = ... K*(U2-U1)
        virtual tmp<volScalarField> K() const;

        //- Return the drag coefficient Kf
        //  used in the face-momentum equations
        virtual tmp<surfaceScalarField> Kf() const;

        //- Dummy write for regIOobject
        bool writeData(Ostream& os) const;
```

* drag coefficient $CdRe$
* phase-intensive drag coefficient Ki in momentum equation, $k_i$
* drag coefficient K in momentum equation, $K$
* drag coefficient Kf in the face-momentum equation, $K_f$

### dragModelNew.C

```cpp
#include "dragModel.H"
#include "phasePair.H"

// * * * * * * * * * * * * * * * * Selector  * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::dragModel> Foam::dragModel::New
(
    const dictionary& dict,
    const phasePair& pair
)
{
    word dragModelType(dict.lookup("type"));

    Info<< "Selecting dragModel for "
        << pair << ": " << dragModelType << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(dragModelType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown dragModelType type "
            << dragModelType << endl << endl
            << "Valid dragModel types are : " << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return cstrIter()(dict, pair, true);
}
```

read dict to construct drag model

### dragModel.C

#### include and static data members

```cpp
#include "dragModel.H"
#include "phasePair.H"
#include "noSwarm.H"
#include "surfaceInterpolate.H"
#include "BlendedInterfacialModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(dragModel, 0);
    defineBlendedInterfacialModelTypeNameAndDebug(dragModel, 0);
    defineRunTimeSelectionTable(dragModel, dictionary);
}

const Foam::dimensionSet Foam::dragModel::dimK(1, -3, -1, 0, 0);
```

define the dimension `dimK` as $kg \cdot m^{-3} \cdot s^{-1}$

#### constructor and destructor

```cpp
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dragModel::dragModel
(
    const phasePair& pair,
    const bool registerObject
)
:
    regIOobject
    (
        IOobject
        (
            IOobject::groupName(typeName, pair.name()),
            pair.phase1().mesh().time().timeName(),
            pair.phase1().mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            registerObject
        )
    ),
    pair_(pair)
{}


Foam::dragModel::dragModel
(
    const dictionary& dict,
    const phasePair& pair,
    const bool registerObject
)
:
    regIOobject
    (
        IOobject
        (
            IOobject::groupName(typeName, pair.name()),
            pair.phase1().mesh().time().timeName(),
            pair.phase1().mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            registerObject
        )
    ),
    pair_(pair),
    swarmCorrection_
    (
        dict.found("swarmCorrection")
      ? swarmCorrection::New
        (
            dict.subDict("swarmCorrection"),
            pair
        )
      : autoPtr<swarmCorrection>(new swarmCorrections::noSwarm(dict, pair))
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::dragModel::~dragModel()
{}
```

#### member functions

**CdRe() is not defined here**

##### Ki()

```cpp
Foam::tmp<Foam::volScalarField> Foam::dragModel::Ki() const
{
    return
        0.75
       *CdRe()
       *swarmCorrection_->Cs()
       *pair_.continuous().rho()
       *pair_.continuous().thermo().nu()
       /sqr(pair_.dispersed().d());
}
```

$$
K_i = \frac{0.75 C_{d, Re} C_s \rho_{continuous} \nu_{contimuous}}{d^2_{dispersed}}
$$

##### K()

```cpp
Foam::tmp<Foam::volScalarField> Foam::dragModel::K() const
{
    return max(pair_.dispersed(), pair_.dispersed().residualAlpha())*Ki();
}
```

$$
K = \max(\alpha_{dispersed}, residualAlpha_{dispersed}) K_i
$$

##### Kf()

```cpp
Foam::tmp<Foam::surfaceScalarField> Foam::dragModel::Kf() const
{
    return
        max
        (
            fvc::interpolate(pair_.dispersed()),
            pair_.dispersed().residualAlpha()
        )*fvc::interpolate(Ki());
}
```

$$
K_f = \max((\alpha_{dispersed})_f, residualAlpha_{dispersed}) (K_i)_f
$$

##### writeDta()

```cpp
bool Foam::dragModel::writeData(Ostream& os) const
{
    return os.good();
}
```

## SchillerNaumann

### SchillerNaumann.H

```cpp
    // Private Data

        //- Residual Reynolds Number
        const dimensionedScalar residualRe_;
```

define a private data `residualRe_`

#### CdRe()

```cpp
    // Member Functions

        //- Drag coefficient
        virtual tmp<volScalarField> CdRe() const;
```

define `CdRe()`

### SchillerNaumann.C

#### CdRe()

```cpp
Foam::tmp<Foam::volScalarField> Foam::dragModels::SchillerNaumann::CdRe() const
{
    volScalarField Re(pair_.Re());

    return
        neg(Re - 1000)*24*(1.0 + 0.15*pow(Re, 0.687))
      + pos0(Re - 1000)*0.44*max(Re, residualRe_);
}
```

$Re$ in `pair` is defined as

$$
Re = \frac{\|\mathbf{U}_1 - \mathbf{U}_2\| \cdot d_{dispersed}}{\nu_{continuous}}
$$

`neg()` and `pos0()` is defined as

```cpp
inline label neg(const label s)
{
    return (s < 0)? 1: 0;
}

inline label pos0(const label s)
{
    return (s >= 0)? 1: 0;
}
```

So,

$$
C_{d, Re} = 24 (1 + 0.15 Re^{0.687}), Re < 1000 \\
C_{d, Re} = 0.44 \max(Re, residualRe_), Re > 1000
$$

$$
C_{d, Re} = C_d \cdot Re
$$

## IshiiZuber

### IshiiZuber.H

#### CdRe()

```cpp
        //- Drag coefficient
        virtual tmp<volScalarField> CdRe() const;
```

### IshiiZuber.C

#### CdRe()

```cpp
Foam::tmp<Foam::volScalarField>
Foam::dragModels::IshiiZuber::CdRe() const
{
    const volScalarField Re(pair_.Re());
    const volScalarField Eo(pair_.Eo());

    const volScalarField mud(pair_.dispersed().thermo().mu());
    const volScalarField muc(pair_.continuous().thermo().mu());

    const volScalarField muStar((mud + 0.4*muc)/(mud + muc));

    const volScalarField muMix
    (
        muc*pow(max(1 - pair_.dispersed(), scalar(1e-3)), -2.5*muStar)
    );

    const volScalarField ReM(Re*muc/muMix);
    const volScalarField CdRe
    (
        pos0(1000 - ReM)*24*(1 + 0.1*pow(ReM, 0.75))
      + neg(1000 - ReM)*0.44*ReM
    );

    volScalarField F((muc/muMix)*sqrt(1 - pair_.dispersed()));
    F.max(1e-3);

    const volScalarField Ealpha((1 + 17.67*pow(F, 0.8571428))/(18.67*F));

    const volScalarField CdReEllipse(Ealpha*0.6666*sqrt(Eo)*Re);

    return
        pos0(CdReEllipse - CdRe)
       *min(CdReEllipse, Re*sqr(1 - pair_.dispersed())*2.66667)
      + neg(CdReEllipse - CdRe)*CdRe;
}
```

get $Re$ and $Eo$ as

$$
Re = \frac{\|\mathbf{U}_1 - \mathbf{U}_2\| \cdot d_{dispersed}}{\nu_{continuous}}
$$

$$
Eo = \frac{(\rho_{dispersed} - \rho_{continuous}) \|\mathbf{g}\| d_{dispersed}^2}{\sigma}
$$

$$
\mu_d = \mu_{dispersed}
$$

$$
\mu_c = \mu_{continuous}
$$

$$
\mu^* = \frac{\mu_d + 0.4 \mu_c}{\mu_d + \mu_c}
$$

$$
\alpha_d = \alpha_{dispersed}
$$

$$
\mu_{Mix} = \mu_c \left[\max(1 - \alpha_d, 10^{-3})\right]^{-2.5 \mu^*}
$$

$$
Re_{M} = \frac{Re \mu_c}{\mu_d}
$$

$$
C_{d, Re} = 24 (1 + 0.1 Re_{M}^{0.75}), Re_M > 1000
$$

$$
C_{d, Re} = 0.44 Re_{M}, Re_M < 1000
$$

$$
F = \frac{\mu_c}{\mu_{Mix}} \sqrt{1 - \alpha_d}
$$

$$
F = \max(F, 10^{-3})
$$

$$
E_\alpha = \frac{1 + 17.67 F^{0.8571428}}{18.67 F}
$$

$$
C_{d, Re, Ellipse} = 0.6666 E_\alpha  \sqrt{Eo} Re
$$

if $C_{d, Re, Ellipse} > C_{d, Re}$, return:

$$
\min(C_{d, Re, Ellipse}, 2.66667 Re (1 - \alpha_d)^2)
$$

else if $C_{d, Re, Ellipse} < C_{d, Re}$, return:

$$
C_{d, Re}
$$








