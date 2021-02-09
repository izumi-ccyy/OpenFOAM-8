# diameterModel

- [diameterModel](#diametermodel)
  - [diameterModel](#diametermodel-1)
    - [diameterModel.H](#diametermodelh)
      - [include](#include)
      - [private data](#private-data)
      - [protected access](#protected-access)
      - [public name, constructor and destructor](#public-name-constructor-and-destructor)
      - [public member functions](#public-member-functions)
    - [diameterModelNew.C](#diametermodelnewc)
    - [diameterModel.C](#diametermodelc)
      - [static data members](#static-data-members)
      - [protected](#protected)
        - [dRef()](#dref)
        - [aRef()](#aref)
      - [constructor and destructor](#constructor-and-destructor)
      - [member functions](#member-functions)
        - [d()](#d)
        - [a()](#a)
        - [correct()](#correct)
        - [read()](#read)
  - [constantDiameter](#constantdiameter)
    - [constantDiameter.H](#constantdiameterh)
      - [inheritance](#inheritance)
      - [private data](#private-data-1)
      - [protected](#protected-1)
      - [public](#public)
    - [constantDiameter.C](#constantdiameterc)
      - [calcD()](#calcd)
      - [constructor and destructor](#constructor-and-destructor-1)
      - [read()](#read-1)
  - [isothermalDiameter](#isothermaldiameter)
    - [isothermalDiameter.H](#isothermaldiameterh)
      - [inheritance and private data](#inheritance-and-private-data)
      - [protected](#protected-2)
      - [public](#public-1)
    - [isothermalDiameter.C](#isothermaldiameterc)
      - [calcD()](#calcd-1)
      - [con- and destructor](#con--and-destructor)
      - [correct()](#correct-1)
      - [read()](#read-2)

## diameterModel

### diameterModel.H

#### include 

```cpp
#ifndef diameterModel_H
#define diameterModel_H

#include "dictionary.H"
#include "phaseModel.H"
#include "runTimeSelectionTables.H"
```

#### private data

```cpp
    // Private Data

        //- The phase diameter properties dictionary
        dictionary diameterProperties_;

        //- The phase that this model applies
        const phaseModel& phase_;

        //- Optionally stored diameter field
        autoPtr<volScalarField> dPtr_;

        //- Optionally stored surface area per unit volume field
        autoPtr<volScalarField> aPtr_;
```

define dict and phase model

#### protected access

```cpp
    // Access

        //- Get a reference to the stored diameter field
        volScalarField& dRef();

        //- Get a reference to the stored surface area per unit volume field
        volScalarField& aRef();

        //- Get the diameter field
        virtual tmp<volScalarField> calcD() const = 0;

        //- Get the surface area per unit volume field
        virtual tmp<volScalarField> calcA() const = 0;
```

get access to:

* a reference to the stored diameter field
* a reference to the stored surface area per unit volume field
* the diameter field
* the surface area per unit volume field

#### public name, constructor and destructor

```cpp
public:

    //- Runtime type information
    TypeName("diameterModel");


    // Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            diameterModel,
            dictionary,
            (
                const dictionary& diameterProperties,
                const phaseModel& phase
            ),
            (diameterProperties, phase)
        );


    // Constructors

        diameterModel
        (
            const dictionary& diameterProperties,
            const phaseModel& phase
        );


    //- Destructor
    virtual ~diameterModel();
    

    // Selectors

        static autoPtr<diameterModel> New
        (
            const dictionary& diameterProperties,
            const phaseModel& phase
        );
```

#### public member functions

```cpp
    // Member Functions

        //- Return the phase diameter properties dictionary
        const dictionary& diameterProperties() const
        {
            return diameterProperties_;
        }

        //- Return the phase
        const phaseModel& phase() const
        {
            return phase_;
        }

        //- Return the diameter
        tmp<volScalarField> d() const;

        //- Return the surface area per unit volume
        tmp<volScalarField> a() const;

        //- Correct the diameter field
        virtual void correct();

        //- Read phaseProperties dictionary
        virtual bool read(const dictionary& phaseProperties);
```

to:

* Return the phase diameter properties dictionary
* Return the phase
* Return the diameter
* Return the surface area per unit volume
* Correct the diameter field
* Read phaseProperties dictionary


### diameterModelNew.C

```cpp
Foam::autoPtr<Foam::diameterModel> Foam::diameterModel::New
(
    const dictionary& dict,
    const phaseModel& phase
)
{
    word diameterModelType
    (
        dict.lookup("diameterModel")
    );

    Info << "Selecting diameterModel for phase "
        << phase.name()
        << ": "
        << diameterModelType << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(diameterModelType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorInFunction
           << "Unknown diameterModelType type "
           << diameterModelType << endl << endl
           << "Valid diameterModel types are : " << endl
           << dictionaryConstructorTablePtr_->sortedToc()
           << exit(FatalError);
    }

    return cstrIter()
    (
        dict.optionalSubDict(diameterModelType + "Coeffs"),
        phase
    );
}
```

read diameter model from the dict

### diameterModel.C

#### static data members

```cpp
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(diameterModel, 0);
    defineRunTimeSelectionTable(diameterModel, dictionary);
}
```

#### protected

##### dRef()

```cpp
Foam::volScalarField& Foam::diameterModel::dRef()
{
    if (!dPtr_.valid())
    {
        dPtr_.reset
        (
            new volScalarField
            (
                IOobject
                (
                    IOobject::groupName("d", phase_.name()),
                    phase_.time().timeName(),
                    phase_.mesh()
                ),
                phase_.mesh(),
                dimensionedScalar(dimLength, 0)
            )
        );
    }

    return dPtr_();
}
```

* if dPtr_ does not exists:
  * create a new field of diameter $d  = 0$ 
* return the existed diameter field or the new created one

##### aRef()

```cpp
Foam::volScalarField& Foam::diameterModel::aRef()
{
    if (!aPtr_.valid())
    {
        aPtr_.reset
        (
            new volScalarField
            (
                IOobject
                (
                    IOobject::groupName("a", phase_.name()),
                    phase_.time().timeName(),
                    phase_.mesh()
                ),
                phase_.mesh(),
                dimensionedScalar(dimless/dimLength, 0)
            )
        );
    }

    return aPtr_();
}
```

Return the existed surface area field or return a new reated zero field

#### constructor and destructor

```cpp
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::diameterModel::diameterModel
(
    const dictionary& diameterProperties,
    const phaseModel& phase
)
:
    diameterProperties_(diameterProperties),
    phase_(phase),
    dPtr_(nullptr),
    aPtr_(nullptr)
{
    if (diameterProperties.lookupOrDefault("storeD", false))
    {
        dRef();
    }
    if (diameterProperties.lookupOrDefault("storeA", false))
    {
        aRef();
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::diameterModel::~diameterModel()
{}
```

#### member functions

##### d()

```cpp
Foam::tmp<Foam::volScalarField> Foam::diameterModel::d() const
{
    if (dPtr_.valid())
    {
        return dPtr_();
    }
    else
    {
        return calcD();
    }
}
```

* if dPtr_ exists
  * return dPtr_
* else 
  * return the result of calcD(), the diameter field

**calcD() is not defined here, but in specific diameter models**

##### a()

```cpp
Foam::tmp<Foam::volScalarField> Foam::diameterModel::a() const
{
    if (aPtr_.valid())
    {
        return aPtr_();
    }
    else
    {
        return calcA();
    }
}
```

* if aPtr_ exists
  * return aPtr_
* else 
  * return the result of calcA(), the diameter field

**calcA() is not defined here, but in specific diameter models**

##### correct()

```cpp
void Foam::diameterModel::correct()
{
    if (dPtr_.valid())
    {
        tmp<volScalarField> td = calcD();
        if (td.isTmp())
        {
            dPtr_() = td;
        }
    }
    if (aPtr_.valid())
    {
        tmp<volScalarField> tA = calcA();
        if (tA.isTmp())
        {
            aPtr_() = tA;
        }
    }
}
```

* if dPtr_ exists
  * set dPtr() = calcD()
  * namely, the diameter field obtained by diameter models
* * if aPtr_ exists
  * set aPtr() = calcA()
  * namely, the surface area field obtained by diameter models

##### read()

```cpp
bool Foam::diameterModel::read(const dictionary& phaseProperties)
{
    diameterProperties_ = phaseProperties.optionalSubDict(type() + "Coeffs");

    return true;
}
```

## constantDiameter

### constantDiameter.H

#### inheritance

```cpp
class constant
:
    public spherical
```

inherited from spherical

#### private data

```cpp
    // Private Data

        //- The constant diameter of the phase
        dimensionedScalar d_;
```

define the constant diameter `d_`

#### protected

```cpp
protected:

    // Protected Member Functions

        //- Get the diameter field
        virtual tmp<volScalarField> calcD() const;
```

#### public

```cpp
public:

    //- Runtime type information
    TypeName("constant");


    // Constructors

        //- Construct from components
        constant
        (
            const dictionary& diameterProperties,
            const phaseModel& phase
        );


    //- Destructor
    virtual ~constant();


    // Member Functions

        //- Read diameterProperties dictionary
        virtual bool read(const dictionary& diameterProperties);
```

### constantDiameter.C

#### calcD()

```cpp
Foam::tmp<Foam::volScalarField> Foam::diameterModels::constant::calcD() const
{
    return volScalarField::New
    (
        IOobject::groupName("d", phase().name()),
        phase().mesh(),
        d_
    );
}
```

define and set the diameter field as the constant `d_`

#### constructor and destructor

```cpp
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::diameterModels::constant::constant
(
    const dictionary& diameterProperties,
    const phaseModel& phase
)
:
    spherical(diameterProperties, phase),
    d_("d", dimLength, diameterProperties)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::diameterModels::constant::~constant()
{}
```

construc from `d_`

#### read()

```cpp
bool Foam::diameterModels::constant::read(const dictionary& phaseProperties)
{
    spherical::read(phaseProperties);

    diameterProperties().lookup("d") >> d_;

    return true;
}
```

read `d_` 

## isothermalDiameter

### isothermalDiameter.H

#### inheritance and private data

```cpp
class isothermal
:
    public spherical
{
    // Private Data

        //- Reference diameter for the isothermal expansion
        dimensionedScalar d0_;

        //- Reference pressure for the isothermal expansion
        dimensionedScalar p0_;

        //- Diameter field
        volScalarField& d_;
```

* also inherited from spherical
* define private data 
  * d0_, specified in dict
  * p0_, specified in dict
  * d_

#### protected

```cpp
    // Protected Member Functions

        //- Get the diameter field
        virtual tmp<volScalarField> calcD() const;
```

#### public

```cpp
public:

    //- Runtime type information
    TypeName("isothermal");


    // Constructors

        //- Construct from components
        isothermal
        (
            const dictionary& diameterProperties,
            const phaseModel& phase
        );


    //- Destructor
    virtual ~isothermal();


    // Member Functions

        //- Correct the diameter field
        virtual void correct();

        //- Read phaseProperties dictionary
        virtual bool read(const dictionary& phaseProperties);
```

compared with `constantDiamter`, the `correct()` is redefiend here. This is because in `constantDiameter`, the diamter will bot change, there is no need to update the diamter. But, in `isothermalDiameter`, the diameter will change, updates are required.

### isothermalDiameter.C

#### calcD()

```cpp
Foam::tmp<Foam::volScalarField> Foam::diameterModels::isothermal::calcD() const
{
    return d_;
}
```

return the calculated diameter field

#### con- and destructor

```cpp
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::diameterModels::isothermal::isothermal
(
    const dictionary& diameterProperties,
    const phaseModel& phase
)
:
    spherical(diameterProperties, phase),
    d0_("d0", dimLength, diameterProperties),
    p0_("p0", dimPressure, diameterProperties),
    d_(dRef())
{
    d_ = d0_;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::diameterModels::isothermal::~isothermal()
{}
```

construct from `d0_`

#### correct()

```cpp
void Foam::diameterModels::isothermal::correct()
{
    const volScalarField& p = phase().db().lookupObject<volScalarField>("p");

    d_ = d0_*pow(p0_/p, 1.0/3.0);
}
```

$$
d = d_0 \left(\frac{p_0}{p}\right)^{1/3}
$$

#### read()

```cpp
bool Foam::diameterModels::isothermal::read(const dictionary& phaseProperties)
{
    spherical::read(phaseProperties);

    diameterProperties().lookup("d0") >> d0_;
    diameterProperties().lookup("p0") >> p0_;

    return true;
}
```

read `d0_` and `p0_`