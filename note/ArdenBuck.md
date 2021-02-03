# ArdenBuck model

## ArdenBuck.H

### Description

ArdenBuck equation for the vapour pressure of moist air.

### include

```cpp
#ifndef ArdenBuck_H
#define ArdenBuck_H

#include "saturationModel.H"
```

### private

```cpp
tmp<volScalarField> xByTC(const volScalarField& TC) const;
```

Exponent divided by the temperature

### public

```cpp
    //- Runtime type information
    TypeName("ArdenBuck");

    // Constructors

        //- Construct from a dictionary
        ArdenBuck(const dictionary& dict, const phasePair& pair);


    //- Destructor
    virtual ~ArdenBuck();
```

constructors and destructors

### member functions

```cpp
    // Member Functions

        //- Saturation pressure
        virtual tmp<volScalarField> pSat(const volScalarField& T) const;

        //- Saturation pressure derivative w.r.t. temperature
        virtual tmp<volScalarField> pSatPrime(const volScalarField& T) const;

        //- Natural log of the saturation pressure
        virtual tmp<volScalarField> lnPSat(const volScalarField& T) const;

        //- Saturation temperature
        virtual tmp<volScalarField> Tsat(const volScalarField& p) const;
```

## ArdenBuck.C

### static members

```cpp
static const Foam::dimensionedScalar zeroC("", Foam::dimTemperature, 273.15);
static const Foam::dimensionedScalar A("", Foam::dimPressure, 611.21);
static const Foam::dimensionedScalar B("", Foam::dimless, 18.678);
static const Foam::dimensionedScalar C("", Foam::dimTemperature, 234.5);
static const Foam::dimensionedScalar D("", Foam::dimTemperature, 257.14);
```

$$
zeroC = 273.15, A = 611.21, B = 18.678, C = 234.5, D = 257.14
$$

### private member functions

```cpp
Foam::tmp<Foam::volScalarField>
Foam::saturationModels::ArdenBuck::xByTC
(
    const volScalarField& TC
) const
{
    return (B - TC/C)/(D + TC);
}
```

$$
xByTC = \frac{B - \frac{T_C}{C}}{D + T_C} = \frac{18.678 - \frac{T_C}{234.5}}{257.14 + T_C}
$$

### constructor and destructor

```cpp
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::saturationModels::ArdenBuck::ArdenBuck
(
    const dictionary& dict,
    const phasePair& pair
)
:
    saturationModel(pair)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::saturationModels::ArdenBuck::~ArdenBuck()
{}
```

### member fuctions

#### pSat()

```cpp
Foam::tmp<Foam::volScalarField>
Foam::saturationModels::ArdenBuck::pSat
(
    const volScalarField& T
) const
{
    volScalarField TC(T - zeroC);

    return A*exp(TC*xByTC(TC));
}
```

$$
T_C = T - zeroC = T - 273.15
$$

$$
p_{Sat} = A \exp(T_C  \frac{18.678 - \frac{T_C}{234.5}}{257.14 + T_C})= 611.21 \exp \left((18.678 - \frac{T_C}{234.5}) (\frac{T_C}{257.14 + T_C})\right)
$$

#### pSatPrime()

```cpp
Foam::tmp<Foam::volScalarField>
Foam::saturationModels::ArdenBuck::pSatPrime
(
    const volScalarField& T
) const
{
    volScalarField TC(T - zeroC);

    volScalarField x(xByTC(TC));

    return A*exp(TC*x)*(D*x - TC/C)/(D + TC);
}
```

$$
x = \frac{18.678 - \frac{T_C}{234.5}}{257.14 + T_C}
$$

$$
\frac{\partial p_{Sat}}{\partial T} = A*exp(TC*x)*(D*x - TC/C)/(D + TC)
$$

#### lnPSat()

```cpp
Foam::tmp<Foam::volScalarField>
Foam::saturationModels::ArdenBuck::lnPSat
(
    const volScalarField& T
) const
{
    volScalarField TC(T - zeroC);

    return log(A.value()) + TC*xByTC(TC);
}
```

$$
\ln{p_{Sat}} = \ln(611.21) + (18.678 - \frac{T_C}{234.5}) (\frac{T_C}{257.14 + T_C})
$$

#### Tsat()

```cpp
Foam::tmp<Foam::volScalarField>
Foam::saturationModels::ArdenBuck::Tsat
(
    const volScalarField& p
) const
{
    NotImplemented;

    return volScalarField::null();
}
```

no definition





