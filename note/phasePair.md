# phasePair

- [phasePair](#phasepair)
  - [phasePair](#phasepair-1)
    - [phasePair.H](#phasepairh)
      - [include](#include)
      - [inheritance](#inheritance)
      - [public data](#public-data)
      - [private data](#private-data)
      - [private member function](#private-member-function)
      - [public](#public)
        - [constructor and destructor](#constructor-and-destructor)
        - [member functions](#member-functions)
          - [calc](#calc)
          - [access](#access)
          - [const iterator](#const-iterator)
    - [phasePairI.H](#phasepairih)
      - [access](#access-1)
      - [interators](#interators)
    - [phasePair.C](#phasepairc)
      - [private member function](#private-member-function-1)
      - [constructor and destructor](#constructor-and-destructor-1)
      - [dispersed() and continuous()](#dispersed-and-continuous)
      - [name() and ortherName()](#name-and-orthername)
      - [rho(), MagUr() ...](#rho-magur-)
  - [orderedPhasePair](#orderedphasepair)
    - [orderedPhasePair.H](#orderedphasepairh)
      - [include](#include-1)
      - [inheritance](#inheritance-1)
      - [public](#public-1)
    - [orderedPhasePair.C](#orderedphasepairc)
      - [constructor and destructor](#constructor-and-destructor-2)
      - [dispersed(), continuous(), name(), othername(), E()](#dispersed-continuous-name-othername-e)
  - [phasePairKey](#phasepairkey)
    - [phasePairKey.H](#phasepairkeyh)
      - [include](#include-2)
      - [friend functions](#friend-functions)
      - [inheritance](#inheritance-2)
      - [data and member function](#data-and-member-function)
      - [constructor and destructor](#constructor-and-destructor-3)
      - [access and friend functions](#access-and-friend-functions)
    - [phasePairKey.C](#phasepairkeyc)
      - [constructor and destructor](#constructor-and-destructor-4)
      - [order()](#order)
      - [Member Operators](#member-operators)
      - [friend operators and istream operators](#friend-operators-and-istream-operators)

## phasePair

### phasePair.H

#### include

```cpp
#ifndef phasePair_H
#define phasePair_H

#include "phaseModel.H"
#include "phasePairKey.H"
#include "uniformDimensionedFields.H"
```

#### inheritance

```cpp
class phasePair
:
    public phasePairKey
{
    ...
}
```

inherited from `phasePairKey`

#### public data

```cpp
public:

    // Hash table types

        //- Dictionary hash table
        typedef HashTable<dictionary, phasePairKey, phasePairKey::hash>
            dictTable;

        //- Scalar hash table
        typedef HashTable<scalar, phasePairKey, phasePairKey::hash>
            scalarTable;
```

define `dictTable` and `scalarTable`

#### private data

```cpp
private:

    // Private Data

        //- Phase 1
        const phaseModel& phase1_;

        //- Phase 2
        const phaseModel& phase2_;

        //- Gravitational acceleration
        const uniformDimensionedVectorField& g_;
```

define 
* `pahse1_`
* `phase_2`
* `g_`

#### private member function

```cpp
    // Private Member Functions

        // Etvos number for given diameter
        tmp<volScalarField> EoH(const volScalarField& d) const;
```

define the `EoH`, the Eötvös number (Eo), also called the Bond number (Bo), is a dimensionless number measuring the importance of gravitational forces compared to surface tension forces and is used (together with Morton number) to characterize the shape of bubbles or drops moving in a surrounding fluid. 

$$
Eo = Bo = \frac{\Delta \rho g d^2}{\sigma}
$$

where: 
* $\Delta \rho$: difference in density of the two phases, 
* $g$: gravitational acceleration
* $d$: characteristic length, (SI units : m) (for example the radii of curvature for a drop)
* $\sigma$: surface tension

#### public

##### constructor and destructor

```cpp
    // Constructors

        //- Construct from two phases and gravity
        phasePair
        (
            const phaseModel& phase1,
            const phaseModel& phase2,
            const bool ordered = false
        );


    //- Destructor
    virtual ~phasePair();
```

##### member functions

###### calc

```cpp
        //- Dispersed phase
        virtual const phaseModel& dispersed() const;

        //- Continuous phase
        virtual const phaseModel& continuous() const;

        //- Pair name
        virtual word name() const;

        //- Other pair name
        virtual word otherName() const;

        //- Average density
        tmp<volScalarField> rho() const;

        //- Relative velocity magnitude
        tmp<volScalarField> magUr() const;

        //- Relative velocity
        tmp<volVectorField> Ur() const;

        //- Reynolds number
        tmp<volScalarField> Re() const;

        //- Prandtl number
        tmp<volScalarField> Pr() const;

        //- Eotvos number
        tmp<volScalarField> Eo() const;

        //- Eotvos number based on hydraulic diameter type 1
        tmp<volScalarField> EoH1() const;

        //- Eotvos number based on hydraulic diameter type 2
        tmp<volScalarField> EoH2() const;

        //- Surface tension coefficient
        tmp<volScalarField> sigma() const;

        //- Morton Number
        tmp<volScalarField> Mo() const;

        //- Takahashi Number
        tmp<volScalarField> Ta() const;

        //- Aspect ratio
        virtual tmp<volScalarField> E() const;
```

###### access

```cpp
            //- Return phase 1
            inline const phaseModel& phase1() const;

            //- Return phase 2
            inline const phaseModel& phase2() const;

            //- Return true if this phasePair contains the given phase
            inline bool contains(const phaseModel& phase) const;

            //- Return the other phase relative to the given phase
            //  Generates a FatalError if this phasePair does not contain
            //  the given phase
            inline const phaseModel& otherPhase(const phaseModel& phase) const;

            //- Return the index of the given phase. Generates a FatalError if
            //  this phasePair does not contain the given phase
            inline label index(const phaseModel& phase) const;

            //- Return gravitation acceleration
            inline const uniformDimensionedVectorField& g() const;

            //- Return the mesh
            inline const fvMesh& mesh() const;
```

return some variables

###### const iterator

```cpp
        //- STL const_iterator
        class const_iterator
        {
            // Private Data

                //- Reference to the pair for which this is an iterator
                const phasePair& pair_;

                //- Current index
                label index_;

                //- Construct an iterator with the given index
                inline const_iterator(const phasePair&, const label index);

        public:

            friend class phasePair;

            // Constructors

                //- Construct from pair, moving to its 'begin' position
                inline explicit const_iterator(const phasePair&);


            // Access

                //- Return the current index
                inline label index() const;


            // Member Operators

                inline bool operator==(const const_iterator&) const;

                inline bool operator!=(const const_iterator&) const;

                inline const phaseModel& operator*() const;
                inline const phaseModel& operator()() const;

                inline const phaseModel& otherPhase() const;

                inline const_iterator& operator++();
                inline const_iterator operator++(int);
        };


        //- const_iterator set to the beginning of the pair
        inline const_iterator cbegin() const;

        //- const_iterator set to beyond the end of the pair
        inline const_iterator cend() const;

        //- const_iterator set to the beginning of the pair
        inline const_iterator begin() const;

        //- const_iterator set to beyond the end of the pair
        inline const_iterator end() const;
```






### phasePairI.H

#### access

```cpp
inline const Foam::phaseModel& Foam::phasePair::phase1() const
{
    return phase1_;
}


inline const Foam::phaseModel& Foam::phasePair::phase2() const
{
    return phase2_;
}


inline bool Foam::phasePair::contains(const phaseModel& phase) const
{
    return &phase1_ == &phase || & phase2_ == &phase;
}


inline const Foam::phaseModel& Foam::phasePair::otherPhase
(
    const phaseModel& phase
) const
{
    if (&phase1_ == &phase)
    {
        return phase2_;
    }
    else if (&phase2_ == &phase)
    {
        return phase1_;
    }
    else
    {
        FatalErrorInFunction
            << "this phasePair does not contain phase " << phase.name()
            << exit(FatalError);

        return phase;
    }
}


inline Foam::label Foam::phasePair::index(const phaseModel& phase) const
{
    if (&phase1_ == &phase)
    {
        return 0;
    }
    else if (&phase2_ == &phase)
    {
        return 1;
    }
    else
    {
        FatalErrorInFunction
            << "this phasePair does not contain phase " << phase.name()
            << exit(FatalError);

        return -1;
    }
}


inline const Foam::uniformDimensionedVectorField& Foam::phasePair::g() const
{
    return g_;
}


inline const Foam::fvMesh& Foam::phasePair::mesh() const
{
    return phase1_.mesh();
}
```

return:

* phase1_
* phase2_
* if the phase in the pair
* the other phase
* index of the phase
  * 0 for phase1_
  * 1 for phase2)
* gravity g_
* mesh_

#### interators

```cpp
inline Foam::phasePair::const_iterator::const_iterator
(
    const phasePair& pair,
    const label index
)
:
    pair_(pair),
    index_(index)
{}


inline Foam::phasePair::const_iterator::const_iterator(const phasePair& pair)
:
    const_iterator(pair, 0)
{}


inline Foam::label Foam::phasePair::const_iterator::index() const
{
    return index_;
}


inline bool Foam::phasePair::const_iterator::operator==
(
    const const_iterator& iter
) const
{
    return (this->index_ == iter.index_);
}


inline bool Foam::phasePair::const_iterator::operator!=
(
    const const_iterator& iter
) const
{
    return !(this->operator==(iter));
}


inline const Foam::phaseModel&
Foam::phasePair::const_iterator::operator*() const
{
    if (index_ == 0)
    {
        return pair_.phase1_;
    }
    else
    {
        return pair_.phase2_;
    }
}


inline const Foam::phaseModel&
Foam::phasePair::const_iterator::operator()() const
{
    return operator*();
}


inline const Foam::phaseModel&
Foam::phasePair::const_iterator::otherPhase() const
{
    if (index_ == 0)
    {
        return pair_.phase2_;
    }
    else
    {
        return pair_.phase1_;
    }
}


inline Foam::phasePair::const_iterator&
Foam::phasePair::const_iterator::operator++()
{
    index_++;
    return *this;
}


inline Foam::phasePair::const_iterator
Foam::phasePair::const_iterator::operator++(int)
{
    const_iterator old = *this;
    this->operator++();
    return old;
}


inline Foam::phasePair::const_iterator Foam::phasePair::cbegin() const
{
    return const_iterator(*this);
}


inline Foam::phasePair::const_iterator Foam::phasePair::cend() const
{
    return const_iterator(*this, 2);
}


inline Foam::phasePair::const_iterator Foam::phasePair::begin() const
{
    return const_iterator(*this);
}


inline Foam::phasePair::const_iterator Foam::phasePair::end() const
{
    return const_iterator(*this, 2);
}
```

### phasePair.C

#### private member function

```cpp
Foam::tmp<Foam::volScalarField> Foam::phasePair::EoH
(
    const volScalarField& d
) const
{
    return
        mag(dispersed().rho() - continuous().rho())
       *mag(g())
       *sqr(d)
       /sigma();
}
```

calculate $Eo$ as definition:

$$
Eo = \frac{(\rho_{dispersed} - \rho_{continuous}) \|\mathbf{g}\| d^2}{\sigma}
$$

#### constructor and destructor

```cpp
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::phasePair::phasePair
(
    const phaseModel& phase1,
    const phaseModel& phase2,
    const bool ordered
)
:
    phasePairKey(phase1.name(), phase2.name(), ordered),
    phase1_(phase1),
    phase2_(phase2),
    g_(phase1.mesh().lookupObject<uniformDimensionedVectorField>("g"))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::phasePair::~phasePair()
{}
```

#### dispersed() and continuous()

```cpp
const Foam::phaseModel& Foam::phasePair::dispersed() const
{
    FatalErrorInFunction
        << "Requested dispersed phase from an unordered pair."
        << exit(FatalError);

    return phase1_;
}


const Foam::phaseModel& Foam::phasePair::continuous() const
{
    FatalErrorInFunction
        << "Requested continuous phase from an unordered pair."
        << exit(FatalError);

    return phase1_;
}
```

if dispersed() and continuous() defined here are performed, the pair was not ordered and error emerge. The definitions in `orderedPhasePair`
should be used

#### name() and ortherName()

```cpp
Foam::word Foam::phasePair::name() const
{
    word name2(second());
    name2[0] = toupper(name2[0]);
    return first() + "And" + name2;
}


Foam::word Foam::phasePair::otherName() const
{
    word name1(first());
    name1[0] = toupper(name1[0]);
    return second() + "And" + name1;
}
```

* output the name of the two phases in the pair
* output the name of other phase in the pair

#### rho(), MagUr() ...

```cpp
Foam::tmp<Foam::volScalarField> Foam::phasePair::rho() const
{
    return phase1()*phase1().rho() + phase2()*phase2().rho();
}


Foam::tmp<Foam::volScalarField> Foam::phasePair::magUr() const
{
    return mag(phase1().U() - phase2().U());
}


Foam::tmp<Foam::volVectorField> Foam::phasePair::Ur() const
{
    return dispersed().U() - continuous().U();
}


Foam::tmp<Foam::volScalarField> Foam::phasePair::Re() const
{
    return magUr()*dispersed().d()/continuous().thermo().nu();
}


Foam::tmp<Foam::volScalarField> Foam::phasePair::Pr() const
{
    return
         continuous().thermo().nu()
        *continuous().thermo().Cpv()
        *continuous().rho()
        /continuous().thermo().kappa();
}


Foam::tmp<Foam::volScalarField> Foam::phasePair::Eo() const
{
    return EoH(dispersed().d());
}


Foam::tmp<Foam::volScalarField> Foam::phasePair::EoH1() const
{
    return
        EoH
        (
            dispersed().d()
           *cbrt(1 + 0.163*pow(Eo(), 0.757))
        );
}


Foam::tmp<Foam::volScalarField> Foam::phasePair::EoH2() const
{
    return
        EoH
        (
            dispersed().d()
           /cbrt(E())
        );
}


Foam::tmp<Foam::volScalarField> Foam::phasePair::sigma() const
{
    return
        phase1().fluid().lookupSubModel<surfaceTensionModel>
        (
            phasePair(phase1(), phase2())
        ).sigma();
}


Foam::tmp<Foam::volScalarField> Foam::phasePair::Mo() const
{
    return
        mag(g())
       *continuous().thermo().nu()
       *pow3
        (
            continuous().thermo().nu()
           *continuous().rho()
           /sigma()
        );
}


Foam::tmp<Foam::volScalarField> Foam::phasePair::Ta() const
{
    return Re()*pow(Mo(), 0.23);
}


Foam::tmp<Foam::volScalarField> Foam::phasePair::E() const
{
    FatalErrorInFunction
        << "Requested aspect ratio of the dispersed phase in an unordered pair"
        << exit(FatalError);

    return phase1();
}
```

$$
\rho = \alpha_1 \rho_1 + \alpha_2 \rho_2
$$

$$
magUr = \|\mathbf{U}_r\| = \|\mathbf{U}_1 - \mathbf{U}_2\|
$$

$$
\mathbf{U}_r = \mathbf{U}_{dispersied} - \mathbf{U}_{conrinuous}
$$

$$
Re = \frac{\|\mathbf{U}_1 - \mathbf{U}_2\| \cdot d_{dispersed}}{\nu_{continuous}}
$$

$$
Pr = \frac{\alpha}{\nu} = \frac{c_p \mu}{k}
$$

* $\nu$: momentum diffusivity (kinematic viscosity), $\nu =\mu /\rho$, (SI units: m2/s)
* $\alpha$: thermal diffusivity, $\alpha =k/(\rho c_{p})$, (SI units: m2/s)
* $\mu$: dynamic viscosity, (SI units: Pa s = N s/m2)
* $k$: thermal conductivity, (SI units: W/(m·K))
* $c_{p}$: specific heat, (SI units: J/(kg·K))
* $\rho$: density, (SI units: kg/m3).

$$
Pr = \frac{\nu_{continuous} c_{p, continuous} \rho_{continuous}}{\kappa_{continuous}}
$$

$$
Eo = \frac{(\rho_{dispersed} - \rho_{continuous}) \|\mathbf{g}\| d_{dispersed}^2}{\sigma}
$$

$$
Eo_1 = \frac{(\rho_{dispersed} - \rho_{continuous}) \|\mathbf{g}\| [d_{dispersed}(1 + 0.163 Eo^{0.757})^{1/3}]^2}{\sigma}
$$

$$
Eo_2 = \frac{(\rho_{dispersed} - \rho_{continuous}) \|\mathbf{g}\| [d_{dispersed}E^{1/3}]^2}{\sigma}
$$

where $E$ is the aspect ratio

$$
sigma = \sigma
$$

the Morton number (Mo) is a dimensionless number used together with the Eötvös number or Bond number to characterize the shape of bubbles or drops moving in a surrounding fluid or continuous phase, c. It is named after Rose Morton, who described it with W. L. Haberman in 1953.

$$
Mo = \frac{g \mu_{continuous}^4 \Delta \rho}{\rho_{continuous}^2 \sigma^3}
$$

where $\sigma$ is the surface tension coefficient

For the case of a bubble with a negligible inner density the Morton number can be simplified to

$$
{Mo} = \frac{g\mu_{continuous}^4}{\rho_{continuous} \sigma^3}.
$$

$$
Mo = \|\mathbf{g}\| \nu_{continuous}\left(\frac{\nu_{continuous} \rho}{\sigma}\right)^3
$$

$$
Ta = Re Mo^{0.23}
$$

the aspect ratio $E$ is calc from `orderedPhasePair`



## orderedPhasePair

### orderedPhasePair.H

#### include

```cpp
#ifndef orderedPhasePair_H
#define orderedPhasePair_H

#include "phasePair.H"
```

#### inheritance

inherited from `phasePair`

#### public

```cpp
public:

    // Constructors

        //- Construct from two phases and gravity
        orderedPhasePair
        (
            const phaseModel& dispersed,
            const phaseModel& continuous
        );


    //- Destructor
    virtual ~orderedPhasePair();


    // Member Functions

        //- Dispersed phase
        virtual const phaseModel& dispersed() const;

        //- Continuous phase
        virtual const phaseModel& continuous() const;

        //- Pair name
        virtual word name() const;

        //- Other pair name
        virtual word otherName() const;

        //- Aspect ratio
        virtual tmp<volScalarField> E() const;
```

### orderedPhasePair.C

#### constructor and destructor

```cpp
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::orderedPhasePair::orderedPhasePair
(
    const phaseModel& dispersed,
    const phaseModel& continuous
)
:
    phasePair
    (
        dispersed,
        continuous,
        true
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::orderedPhasePair::~orderedPhasePair()
{}
```

#### dispersed(), continuous(), name(), othername(), E()

```cpp
const Foam::phaseModel& Foam::orderedPhasePair::dispersed() const
{
    return phase1();
}


const Foam::phaseModel& Foam::orderedPhasePair::continuous() const
{
    return phase2();
}

Foam::word Foam::orderedPhasePair::name() const
{
    word namec(second());
    namec[0] = toupper(namec[0]);
    return first() + "In" + namec;
}


Foam::word Foam::orderedPhasePair::otherName() const
{
    FatalErrorInFunction
        << "Requested other name phase from an ordered pair."
        << exit(FatalError);

    return word::null;
}


Foam::tmp<Foam::volScalarField> Foam::orderedPhasePair::E() const
{
    return phase1().fluid().E(*this);
}
```

in `orderedPhasePair`, `phase1_` is dispersed while `phase2_` is continuous

name of phase1_ "IN" name of phase2_, so the entries containing "in", the former is dispersed and the latter is continuous. In `phasePair`, it is "and" rather than "in"

aspect ratio is defined for dispersed phase, so it is calculated with `phase1_`

**ATTENTION: It should be noted that in `orderedPhasePair`, dispersed and continuous are defined but how to initialize or how to specify which phase is dispersed or continuous is not mentioned here!**

## phasePairKey

### phasePairKey.H

#### include

```cpp
#ifndef phasePairKey_H
#define phasePairKey_H

#include "Pair.H"
```

#### friend functions

```cpp
// Forward declaration of friend functions and operators

class phasePairKey;

bool operator==(const phasePairKey&, const phasePairKey&);
bool operator!=(const phasePairKey&, const phasePairKey&);

Istream& operator>>(Istream&, phasePairKey&);
Ostream& operator<<(Ostream&, const phasePairKey&);
```

#### inheritance

inherited from `Pair`

#### data and member function

```cpp
        class hash
        :
            public Hash<phasePairKey>
        {
        public:

            // Constructors

                // Construct null
                hash();


            // Member Operators

                // Generate a hash from a phase pair key
                label operator()(const phasePairKey& key) const;
        };
```

define a class

```cpp
private:

    // Private Data

        //- Flag to indicate whether ordering is important
        bool ordered_;
```

define a flag

#### constructor and destructor

```cpp
    // Constructors

        //- Construct null
        phasePairKey();

        //- Construct from names and the ordering flag
        phasePairKey
        (
            const word& name1,
            const word& name2,
            const bool ordered = false
        );


    // Destructor
    virtual ~phasePairKey();
```

#### access and friend functions

```cpp
    // Access

        //- Return the ordered flag
        bool ordered() const;


    // Friend Operators

        //- Test if keys are equal
        friend bool operator==(const phasePairKey& a, const phasePairKey& b);

        //- Test if keys are unequal
        friend bool operator!=(const phasePairKey& a, const phasePairKey& b);

        //- Read from stdin
        friend Istream& operator>>(Istream& is, phasePairKey& key);

        //- Write to stdout
        friend Ostream& operator<<(Ostream& os, const phasePairKey& key);
```

constructor and destructor

### phasePairKey.C

#### constructor and destructor

```cpp
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::phasePairKey::hash::hash()
{}


Foam::phasePairKey::phasePairKey()
{}


Foam::phasePairKey::phasePairKey
(
    const word& name1,
    const word& name2,
    const bool ordered
)
:
    Pair<word>(name1, name2),
    ordered_(ordered)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::phasePairKey::~phasePairKey()
{}
```

#### order()

```cpp
bool Foam::phasePairKey::ordered() const
{
    return ordered_;
}
```

return if ordered

#### Member Operators

```cpp
Foam::label Foam::phasePairKey::hash::operator()
(
    const phasePairKey& key
) const
{
    if (key.ordered_)
    {
        return
            word::hash()
            (
                key.first(),
                word::hash()(key.second())
            );
    }
    else
    {
        return
            word::hash()(key.first())
          + word::hash()(key.second());
    }
}
```

Generate a hash from a phase pair key, 
* if ordered:
  * generate one hash
* else:
  * generate two hash

#### friend operators and istream operators

```cpp
// * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * * //

bool Foam::operator==
(
    const phasePairKey& a,
    const phasePairKey& b
)
{
    const label c = Pair<word>::compare(a, b);

    return
        (a.ordered_ == b.ordered_)
     && (
            (a.ordered_ && (c == 1))
         || (!a.ordered_ && (c != 0))
        );
}


bool Foam::operator!=
(
    const phasePairKey& a,
    const phasePairKey& b
)
{
    return !(a == b);
}


// * * * * * * * * * * * * * * Istream Operator  * * * * * * * * * * * * * * //

Foam::Istream& Foam::operator>>(Istream& is, phasePairKey& key)
{
    const FixedList<word, 3> temp(is);

    key.first() = temp[0];

    if (temp[1] == "and")
    {
        key.ordered_ = false;
    }
    else if (temp[1] == "in")
    {
        key.ordered_ = true;
    }
    else
    {
        FatalErrorInFunction
            << "Phase pair type is not recognised. "
            << temp
            << "Use (phaseDispersed in phaseContinuous) for an ordered"
            << "pair, or (phase1 and pase2) for an unordered pair."
            << exit(FatalError);
    }

    key.second() = temp[2];

    return is;
}


// * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const phasePairKey& key)
{
    os  << token::BEGIN_LIST
        << key.first()
        << token::SPACE
        << (key.ordered_ ? "in" : "and")
        << token::SPACE
        << key.second()
        << token::END_LIST;

    return os;
}
```

* Test if keys are equal
* Test if keys are unequal
* Read from stdin
* Write to stdout

from `operator>>()` it can be noted that:

* "in" means ordered pair, the former is the dispersed phase
* "and" means unordered pair