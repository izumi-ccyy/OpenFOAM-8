# BlendedInterfacialModel

- [BlendedInterfacialModel](#blendedinterfacialmodel)
  - [BlendedInterfacialModel.H](#blendedinterfacialmodelh)
    - [include](#include)
    - [private data](#private-data)
    - [private member functions](#private-member-functions)
    - [public](#public)
      - [constructor and destructor](#constructor-and-destructor)
      - [member functions and operator](#member-functions-and-operator)
  - [BlendedInterfacialModel.C](#blendedinterfacialmodelc)
    - [interpolate](#interpolate)
    - [private member functions](#private-member-functions-1)
      - [calculateBlendingCoeffs()](#calculateblendingcoeffs)
      - [correctFixedFluxBCs()](#correctfixedfluxbcs)
      - [evaluate() - tmp](#evaluate---tmp)
      - [evaluate() - HashPtrTable](#evaluate---hashptrtable)
    - [constructors and destructor](#constructors-and-destructor)
      - [constructor - from two phase](#constructor---from-two-phase)
      - [constructor - from model table](#constructor---from-model-table)
      - [Destructor](#destructor)
    - [member functions](#member-functions)
      - [K(), Kf(), F(), Ff(), D()](#k-kf-f-ff-d)
      - [mixture()](#mixture)
      - [dmdtf()](#dmdtf)
      - [species()](#species)
      - [dmidtf()](#dmidtf)
      - [writeData()](#writedata)
  - [blendingMethods](#blendingmethods)
    - [blendingMethod](#blendingmethod)
      - [blendingMethod.H](#blendingmethodh)
        - [include](#include-1)
        - [public](#public-1)
          - [type name](#type-name)
          - [constructor, selector and destructor](#constructor-selector-and-destructor)
          - [member functions](#member-functions-1)
      - [blendingMethodNew.C](#blendingmethodnewc)
      - [blendingMethod.C](#blendingmethodc)
        - [include and static data member](#include-and-static-data-member)
        - [constructor and destructor](#constructor-and-destructor-1)
    - [noBlending](#noblending)
      - [noBlending.H](#noblendingh)
        - [beginning](#beginning)
        - [constructor and destructor](#constructor-and-destructor-2)
        - [member functions](#member-functions-2)
      - [noBlending.C](#noblendingc)
        - [include and static data member](#include-and-static-data-member-1)
        - [constructor and destructor](#constructor-and-destructor-3)
        - [f1() and f2()](#f1-and-f2)
    - [linear](#linear)
      - [linear.H](#linearh)
        - [include](#include-2)
        - [private data](#private-data-1)
        - [name and constructor and destructor](#name-and-constructor-and-destructor)
        - [f1() and f2()](#f1-and-f2-1)
      - [linear.C](#linearc)
        - [include](#include-3)
        - [constructor and destructor](#constructor-and-destructor-4)
        - [f1() and f2()](#f1-and-f2-2)
    - [hyperbolic](#hyperbolic)
      - [hyperbolic.H](#hyperbolich)
        - [private data](#private-data-2)
      - [hyperbolic.C](#hyperbolicc)
        - [f1() and f2()](#f1-and-f2-3)
  - [Discussion](#discussion)

## BlendedInterfacialModel.H

### include

```cpp
#ifndef BlendedInterfacialModel_H
#define BlendedInterfacialModel_H

#include "blendingMethod.H"
#include "phasePair.H"
#include "orderedPhasePair.H"
#include "HashPtrTable.H"
#include "hashedWordList.H"
#include "geometricZeroField.H"
```

### private data

```cpp
        //- Reference to phase 1
        const phaseModel& phase1_;

        //- Reference to phase 2
        const phaseModel& phase2_;

        //- Blending model
        const blendingMethod& blending_;

        //- Model for region with no obvious dispersed phase
        autoPtr<ModelType> model_;

        //- Model for dispersed phase 1 in continuous phase 2
        autoPtr<ModelType> model1In2_;

        //- Model for dispersed phase 2 in continuous phase 1
        autoPtr<ModelType> model2In1_;

        //- If true set coefficients and forces to 0 at fixed-flux BCs
        bool correctFixedFluxBCs_;
```

### private member functions

```cpp
        //- Calculate the blending coefficients
        template<template<class> class PatchField, class GeoMesh>
        void calculateBlendingCoeffs
        (
            tmp<GeometricField<scalar, PatchField, GeoMesh>>& f1,
            tmp<GeometricField<scalar, PatchField, GeoMesh>>& f2,
            const bool subtract
        ) const;

        //- Correct coeff/value on fixed flux boundary conditions
        template<class Type, template<class> class PatchField, class GeoMesh>
        void correctFixedFluxBCs
        (
            GeometricField<Type, PatchField, GeoMesh>& field
        ) const;

        //- Return the blended coeff/value
        template
        <
            class Type,
            template<class> class PatchField,
            class GeoMesh,
            class ... Args
        >
        tmp<GeometricField<Type, PatchField, GeoMesh>> evaluate
        (
            tmp<GeometricField<Type, PatchField, GeoMesh>>
            (ModelType::*method)(Args ...) const,
            const word& name,
            const dimensionSet& dims,
            const bool subtract,
            Args ... args
        ) const;

        //- Return the blended coeff/value
        template
        <
            class Type,
            template<class> class PatchField,
            class GeoMesh,
            class ... Args
        >
        HashPtrTable<GeometricField<Type, PatchField, GeoMesh>> evaluate
        (
            HashPtrTable<GeometricField<Type, PatchField, GeoMesh>>
            (ModelType::*method)(Args ...) const,
            const word& name,
            const dimensionSet& dims,
            const bool subtract,
            Args ... args
        ) const;
```

### public

#### constructor and destructor

```cpp
    //- Runtime type information
    TypeName("BlendedInterfacialModel");


    // Constructors

        //- Construct from two phases, blending method and three models
        BlendedInterfacialModel
        (
            const phaseModel& phase1,
            const phaseModel& phase2,
            const blendingMethod& blending,
            autoPtr<ModelType> model,
            autoPtr<ModelType> model1In2,
            autoPtr<ModelType> model2In1,
            const bool correctFixedFluxBCs = true
        );


        //- Construct from the model table, dictionary and pairs
        BlendedInterfacialModel
        (
            const phasePair::dictTable& modelTable,
            const blendingMethod& blending,
            const phasePair& pair,
            const orderedPhasePair& pair1In2,
            const orderedPhasePair& pair2In1,
            const bool correctFixedFluxBCs = true
        );

        //- Disallow default bitwise copy construction
        BlendedInterfacialModel
        (
            const BlendedInterfacialModel<ModelType>&
        ) = delete;


    //- Destructor
    ~BlendedInterfacialModel();
```

#### member functions and operator

```cpp
    // Member Functions

        //- Return the blended force coefficient
        tmp<volScalarField> K() const;

        //- Return the blended force coefficient with a specified residual alpha
        tmp<volScalarField> K(const scalar residualAlpha) const;

        //- Return the face blended force coefficient
        tmp<surfaceScalarField> Kf() const;

        //- Return the blended force
        template<class Type>
        tmp<GeometricField<Type, fvPatchField, volMesh>> F() const;

        //- Return the face blended force
        tmp<surfaceScalarField> Ff() const;

        //- Return the blended diffusivity
        tmp<volScalarField> D() const;

        //- Return the list of individual species that are transferred
        bool mixture() const;

        //- Return the blended mass transfer rate
        tmp<volScalarField> dmdtf() const;

        //- Return the list of individual species that are transferred
        hashedWordList species() const;

        //- Return the blended mass transfer rates for individual species
        HashPtrTable<volScalarField> dmidtf() const;

        //- Dummy write for regIOobject
        bool writeData(Ostream& os) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const BlendedInterfacialModel<ModelType>&) = delete;

```

## BlendedInterfacialModel.C

### interpolate

```cpp
namespace Foam
{
namespace blendedInterfacialModel
{

template<class GeoField>
inline tmp<GeoField> interpolate(tmp<volScalarField> f);

template<>
inline tmp<Foam::volScalarField> interpolate(tmp<volScalarField> f)
{
    return f;
}

template<>
inline tmp<Foam::surfaceScalarField> interpolate(tmp<volScalarField> f)
{
    return fvc::interpolate(f);
}

} // End namespace blendedInterfacialModel
} // End namespace Foam
```

define some interpolate functions

### private member functions

#### calculateBlendingCoeffs()

```cpp
template<class ModelType>
template<template<class> class PatchField, class GeoMesh>
void Foam::BlendedInterfacialModel<ModelType>::calculateBlendingCoeffs
(
    tmp<GeometricField<scalar, PatchField, GeoMesh>>& f1,
    tmp<GeometricField<scalar, PatchField, GeoMesh>>& f2,
    const bool subtract
) const
{
    typedef GeometricField<scalar, PatchField, GeoMesh> scalarGeoField;

    if (model_.valid() && subtract)
    {
        FatalErrorInFunction
            << "Cannot treat an interfacial model with no distinction between "
            << "continuous and dispersed phases as signed"
            << exit(FatalError);
    }

    if (model_.valid() || model1In2_.valid())
    {
        f1 =
            blendedInterfacialModel::interpolate<scalarGeoField>
            (
                blending_.f1(phase1_, phase2_)
            );
    }

    if (model_.valid() || model2In1_.valid())
    {
        f2 =
            (subtract ? -1 : +1)
           *blendedInterfacialModel::interpolate<scalarGeoField>
            (
                blending_.f2(phase1_, phase2_)
            );
    }
}
```

Calculate the blending coefficients, three situations:

* if **model for region with no obvious dispersed phase** exists and `subtract` is true, then
  * ERROR!
* if **model for region with no obvious dispersed phase** exists and dispersed phase 1 in continuous phase 2, then
  * get f1 from the blending method using phase1 and phase2 and then interpolate it as defined in above to obtain $f_1$
* if **model for region with no obvious dispersed phase** exists and dispersed phase 2 in continuous phase 1, then
  * get f2 from the blending method using phase1 and phase2 and then interpolate it as defined in above to obtain $f_2$
  * if `subtract` is true, then
    * $f_2 = - f_2$
  * else, `subtract` is false, then
    * $f_2 = f_2$

**what is `subtract`?**

#### correctFixedFluxBCs()

```cpp
template<class ModelType>
template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::BlendedInterfacialModel<ModelType>::correctFixedFluxBCs
(
    GeometricField<Type, PatchField, GeoMesh>& field
) const
{
    typedef GeometricField<Type, PatchField, GeoMesh> typeGeoField;

    typename typeGeoField::Boundary& fieldBf = field.boundaryFieldRef();

    forAll(fieldBf, patchi)
    {
        if
        (
            (
                !phase1_.stationary()
             && isA<fixedValueFvsPatchScalarField>
                (
                    phase1_.phi()().boundaryField()[patchi]
                )
            )
         || (
                !phase2_.stationary()
             && isA<fixedValueFvsPatchScalarField>
                (
                    phase2_.phi()().boundaryField()[patchi]
                )
            )
        )
        {
            fieldBf[patchi] = Zero;
        }
    }
}
```

Correct coeff/value on fixed flux boundary conditions:

* get boundary of the field as `fieldBf`
* for every `patch` in `fieldBf`
  * if phase1_ is not a stationary phase and the BC for phase1_ on this patch is fixedValue 
  * Or phase2_ is not a stationary phase and the BC for phase2_ on this patch is fixedValue, then
    * fieldBf[patchi] = Zero

#### evaluate() - tmp

```cpp
template<class ModelType>
template
<
    class Type,
    template<class> class PatchField,
    class GeoMesh,
    class ... Args
>
Foam::tmp<Foam::GeometricField<Type, PatchField, GeoMesh>>
Foam::BlendedInterfacialModel<ModelType>::evaluate
(
    tmp<GeometricField<Type, PatchField, GeoMesh>>
    (ModelType::*method)(Args ...) const,
    const word& name,
    const dimensionSet& dims,
    const bool subtract,
    Args ... args
) const
{
    typedef GeometricField<scalar, PatchField, GeoMesh> scalarGeoField;
    typedef GeometricField<Type, PatchField, GeoMesh> typeGeoField;

    tmp<scalarGeoField> f1, f2;
    calculateBlendingCoeffs(f1, f2, subtract);

    tmp<typeGeoField> x =
        typeGeoField::New
        (
            ModelType::typeName + ":"
          + IOobject::groupName(name, phasePair(phase1_, phase2_).name()),
            phase1_.mesh(),
            dimensioned<Type>(dims, Zero)
        );

    if (model_.valid())
    {
        x.ref() += (scalar(1) - f1() - f2())*(model_().*method)(args ...);
    }

    if (model1In2_.valid())
    {
        x.ref() += f1*(model1In2_().*method)(args ...);
    }

    if (model2In1_.valid())
    {
        x.ref() += f2*(model2In1_().*method)(args ...);
    }

    if
    (
        correctFixedFluxBCs_
     && (model_.valid() || model1In2_.valid() || model2In1_.valid())
    )
    {
        correctFixedFluxBCs(x.ref());
    }

    return x;
}
```

* define type
* calculate $f_1$, $f_2$
* define x
* if model_ is not empty, region with no obvious dispersed phase exists, then
  * $x = x + (1 - f_1 - f_2) * model_().*method$
* if model1In2_ is not empty, dispersed phase 1 in continuous phase 2
  * $$x = x + f_1 * model1In2_.*method$$
* if model2In1_ is not empty, dispersed phase 2 in continuous phase 1
  * $x = x + f_2 * model22In1_.*method$
* if it is required to set coefficients and forces to 0 at fixed-flux BCs and model_ or model1In2_ or model2In1_ is not empty, then
  * correct the fixed flux BCs of every element in xs, set zero
* return `x`

#### evaluate() - HashPtrTable

```cpp
template<class ModelType>
template
<
    class Type,
    template<class> class PatchField,
    class GeoMesh,
    class ... Args
>
Foam::HashPtrTable<Foam::GeometricField<Type, PatchField, GeoMesh>>
Foam::BlendedInterfacialModel<ModelType>::evaluate
(
    HashPtrTable<GeometricField<Type, PatchField, GeoMesh>>
    (ModelType::*method)(Args ...) const,
    const word& name,
    const dimensionSet& dims,
    const bool subtract,
    Args ... args
) const
{
    typedef GeometricField<scalar, PatchField, GeoMesh> scalarGeoField;
    typedef GeometricField<Type, PatchField, GeoMesh> typeGeoField;

    tmp<scalarGeoField> f1, f2;
    calculateBlendingCoeffs(f1, f2, subtract);

    HashPtrTable<typeGeoField> xs;

    auto addToXs = [&]
    (
        const scalarGeoField& f,
        const HashPtrTable<typeGeoField>& dxs
    )
    {
        forAllConstIter(typename HashPtrTable<typeGeoField>, dxs, dxIter)
        {
            if (xs.found(dxIter.key()))
            {
                *xs[dxIter.key()] += f**dxIter();
            }
            else
            {
                xs.insert
                (
                    dxIter.key(),
                    typeGeoField::New
                    (
                        ModelType::typeName + ':'
                      + IOobject::groupName
                        (
                            IOobject::groupName(name, dxIter.key()),
                            phasePair(phase1_, phase2_).name()
                        ),
                        f**dxIter()
                    ).ptr()
                );
            }
        }
    };

    if (model_.valid())
    {
        addToXs(scalar(1) - f1() - f2(), (model_().*method)(args ...));
    }

    if (model1In2_.valid())
    {
        addToXs(f1, (model1In2_().*method)(args ...));
    }

    if (model2In1_.valid())
    {
        addToXs(f2, (model2In1_().*method)(args ...));
    }

    if
    (
        correctFixedFluxBCs_
     && (model_.valid() || model1In2_.valid() || model2In1_.valid())
    )
    {
        forAllIter(typename HashPtrTable<typeGeoField>, xs, xIter)
        {
            correctFixedFluxBCs(*xIter());
        }
    }

    return xs;
}
```

* define type
* calculate blending coefficient $f_1$ and $f_2$
* define hashPtrTable `xs`
* define `addToXs` to add element to `xs` using two parameters f and dxs
  * for every element `dxIter` in `dxs`
    * if dxIter.key can be found in xs, then
      * $$xs[dxIter.key] = xs[dxIter] + f * dxIter$$
    * else, dxIter.key is not in xs, then
      * insert element to xs as
        * xs[dxIter.key] = f * dxIter
* using `addToXs` to add elements from the three models
  * if model_ is not empty, region with no obvious dispersed phase exists, then
    * $f = 1 - f_1 - f_2$
  * if model1In2_ is not empty, dispersed phase 1 in continuous phase 2
    * $f = f_1$
  * if model2In1_ is not empty, dispersed phase 2 in continuous phase 1
    * $f = f_2$
* if it is required to set coefficients and forces to 0 at fixed-flux BCs and model_ or model1In2_ or model2In1_ is not empty, then
  * correct the fixed flux BCs of every element in xs, set zero
* return `xs`

ATTENTION: `f**dxIter()` the two star (*), the first is multiply, the second is a pointer

### constructors and destructor

#### constructor - from two phase

```cpp
template<class ModelType>
Foam::BlendedInterfacialModel<ModelType>::BlendedInterfacialModel
(
    const phaseModel& phase1,
    const phaseModel& phase2,
    const blendingMethod& blending,
    autoPtr<ModelType> model,
    autoPtr<ModelType> model1In2,
    autoPtr<ModelType> model2In1,
    const bool correctFixedFluxBCs
)
:
    regIOobject
    (
        IOobject
        (
            IOobject::groupName(typeName, phasePair(phase1, phase2).name()),
            phase1.mesh().time().timeName(),
            phase1.mesh()
        )
    ),
    phase1_(phase1),
    phase2_(phase2),
    blending_(blending),
    model_(model),
    model1In2_(model1In2),
    model2In1_(model2In1),
    correctFixedFluxBCs_(correctFixedFluxBCs)
{}
```

Construct from two phases, blending method and three models

#### constructor - from model table

```cpp
template<class ModelType>
Foam::BlendedInterfacialModel<ModelType>::BlendedInterfacialModel
(
    const phasePair::dictTable& modelTable,
    const blendingMethod& blending,
    const phasePair& pair,
    const orderedPhasePair& pair1In2,
    const orderedPhasePair& pair2In1,
    const bool correctFixedFluxBCs
)
:
    regIOobject
    (
        IOobject
        (
            IOobject::groupName(typeName, pair.name()),
            pair.phase1().mesh().time().timeName(),
            pair.phase1().mesh()
        )
    ),
    phase1_(pair.phase1()),
    phase2_(pair.phase2()),
    blending_(blending),
    correctFixedFluxBCs_(correctFixedFluxBCs)
{
    if (modelTable.found(pair))
    {
        model_.set
        (
            ModelType::New
            (
                modelTable[pair],
                pair
            ).ptr()
        );
    }

    if (modelTable.found(pair1In2))
    {
        model1In2_.set
        (
            ModelType::New
            (
                modelTable[pair1In2],
                pair1In2
            ).ptr()
        );
    }

    if (modelTable.found(pair2In1))
    {
        model2In1_.set
        (
            ModelType::New
            (
                modelTable[pair2In1],
                pair2In1
            ).ptr()
        );
    }
}
```

Construct from the model table, dictionary and pairs:

* since there is no pair, so first define a pair
* if found the model, then set the pair

#### Destructor

```cpp
template<class ModelType>
Foam::BlendedInterfacialModel<ModelType>::~BlendedInterfacialModel()
{}
```

### member functions

#### K(), Kf(), F(), Ff(), D()

```cpp
template<class ModelType>
Foam::tmp<Foam::volScalarField>
Foam::BlendedInterfacialModel<ModelType>::K() const
{
    tmp<volScalarField> (ModelType::*k)() const = &ModelType::K;

    return evaluate(k, "K", ModelType::dimK, false);
}


template<class ModelType>
Foam::tmp<Foam::volScalarField>
Foam::BlendedInterfacialModel<ModelType>::K(const scalar residualAlpha) const
{
    tmp<volScalarField> (ModelType::*k)(const scalar) const = &ModelType::K;

    return evaluate(k, "K", ModelType::dimK, false, residualAlpha);
}


template<class ModelType>
Foam::tmp<Foam::surfaceScalarField>
Foam::BlendedInterfacialModel<ModelType>::Kf() const
{
    return evaluate(&ModelType::Kf, "Kf", ModelType::dimK, false);
}


template<class ModelType>
template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvPatchField, Foam::volMesh>>
Foam::BlendedInterfacialModel<ModelType>::F() const
{
    return evaluate(&ModelType::F, "F", ModelType::dimF, true);
}


template<class ModelType>
Foam::tmp<Foam::surfaceScalarField>
Foam::BlendedInterfacialModel<ModelType>::Ff() const
{
    return evaluate(&ModelType::Ff, "Ff", ModelType::dimF*dimArea, true);
}


template<class ModelType>
Foam::tmp<Foam::volScalarField>
Foam::BlendedInterfacialModel<ModelType>::D() const
{
    return evaluate(&ModelType::D, "D", ModelType::dimD, false);
}
```

using evaluate() to calculate K(), Kf(), F(), Ff(), D()

#### mixture()

```cpp
template<class ModelType>
bool Foam::BlendedInterfacialModel<ModelType>::mixture() const
{
    return
        (model1In2_.valid() && model1In2_->mixture())
     || (model2In1_.valid() && model2In1_->mixture())
     || (model_.valid() && model_->mixture());
}
```

if there is mixture in ther three models: model1IN2_, model2In1_ or model_, return `true`

#### dmdtf()

```cpp
template<class ModelType>
Foam::tmp<Foam::volScalarField>
Foam::BlendedInterfacialModel<ModelType>::dmdtf() const
{
    return evaluate(&ModelType::dmdtf, "dmdtf", ModelType::dimDmdt, true);
}
```

using evaluate() to calculate dmdtf

#### species()

```cpp
template<class ModelType>
Foam::hashedWordList Foam::BlendedInterfacialModel<ModelType>::species() const
{
    wordList species;

    if (model1In2_.valid())
    {
        species.append(model1In2_->species());
    }
    if (model2In1_.valid())
    {
        species.append(model2In1_->species());
    }
    if (model_.valid())
    {
        species.append(model_->species());
    }

    return hashedWordList(move(species));
}
```

add the species of the three model to the list of `species`

#### dmidtf()

```cpp
template<class ModelType>
Foam::HashPtrTable<Foam::volScalarField>
Foam::BlendedInterfacialModel<ModelType>::dmidtf() const
{
    return evaluate(&ModelType::dmidtf, "dmidtf", ModelType::dimDmdt, true);
}
```

using evaluate() to calculate dmidtf

#### writeData()

```cpp
template<class ModelType>
bool Foam::BlendedInterfacialModel<ModelType>::writeData(Ostream& os) const
{
    return os.good();
}
```

## blendingMethods

### blendingMethod

#### blendingMethod.H

##### include

```cpp
#ifndef blendingMethod_H
#define blendingMethod_H

#include "dictionary.H"
#include "runTimeSelectionTables.H"
#include "phaseModel.H"
```

##### public

###### type name

```cpp
    //- Runtime type information
    TypeName("blendingMethod");


    // Declare runtime construction
    declareRunTimeSelectionTable
    (
        autoPtr,
        blendingMethod,
        dictionary,
        (
            const dictionary& dict,
            const wordList& phaseNames
        ),
        (dict, phaseNames)
    );
```

typename and runtime construction

###### constructor, selector and destructor

```cpp
    // Constructors

        //- Construct from a dictionary
        blendingMethod
        (
            const dictionary& dict
        );


    // Selector

        static autoPtr<blendingMethod> New
        (
            const word& modelName,
            const dictionary& dict,
            const wordList& phaseNames
        );


    //- Destructor
    virtual ~blendingMethod();
```

###### member functions

```cpp
    // Member Functions

        //- Factor for first phase
        virtual tmp<volScalarField> f1
        (
            const phaseModel& phase1,
            const phaseModel& phase2
        ) const = 0;

        //- Factor for second phase
        virtual tmp<volScalarField> f2
        (
            const phaseModel& phase1,
            const phaseModel& phase2
        ) const = 0;
```

f1() and f2()

#### blendingMethodNew.C

define selector

```cpp
#include "blendingMethod.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::blendingMethod> Foam::blendingMethod::New
(
    const word& modelName,
    const dictionary& dict,
    const wordList& phaseNames
)
{
    word blendingMethodType(dict.lookup("type"));

    Info<< "Selecting " << modelName << " blending method: "
        << blendingMethodType << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(blendingMethodType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown blendingMethodType type "
            << blendingMethodType << endl << endl
            << "Valid blendingMethod types are : " << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return cstrIter()(dict, phaseNames);
}
```

read dict to create a blending method

#### blendingMethod.C

##### include and static data member

```cpp
#include "blendingMethod.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(blendingMethod, 0);
    defineRunTimeSelectionTable(blendingMethod, dictionary);
}

```

##### constructor and destructor

```cpp
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::blendingMethod::blendingMethod
(
    const dictionary& dict
)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::blendingMethod::~blendingMethod()
{}
```

ATTENTION: f1() and f2() are defined in specific blending methods instead of here

### noBlending

#### noBlending.H

##### beginning

```cpp
#ifndef noBlending_H
#define noBlending_H

#include "blendingMethod.H"
```

```cpp
class noBlending
:
    public blendingMethod
{
    // Private Data

        //- Name of the continuous phase
        const word continuousPhase_;
```

inherited from `blendingMethod`

define `continuousPhase_` as the name of the continuous phase

##### constructor and destructor

```cpp
    // Constructors

        //- Construct from a dictionary and two phases
        noBlending
        (
            const dictionary& dict,
            const wordList& phaseNames
        );


    //- Destructor
    ~noBlending();
```

##### member functions

```cpp
    // Member Functions

        //- Factor for primary phase
        virtual tmp<volScalarField> f1
        (
            const phaseModel& phase1,
            const phaseModel& phase2
        ) const;

        //- Factor for secondary phase
        virtual tmp<volScalarField> f2
        (
            const phaseModel& phase1,
            const phaseModel& phase2
        ) const;
```

#### noBlending.C

##### include and static data member

```cpp
#include "noBlending.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace blendingMethods
{
    defineTypeNameAndDebug(noBlending, 0);

    addToRunTimeSelectionTable
    (
        blendingMethod,
        noBlending,
        dictionary
    );
}
}
```

##### constructor and destructor

```cpp
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::blendingMethods::noBlending::noBlending
(
    const dictionary& dict,
    const wordList& phaseNames
)
:
    blendingMethod(dict),
    continuousPhase_(dict.lookup("continuousPhase"))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::blendingMethods::noBlending::~noBlending()
{}
```

##### f1() and f2()

```cpp
Foam::tmp<Foam::volScalarField> Foam::blendingMethods::noBlending::f1
(
    const phaseModel& phase1,
    const phaseModel& phase2
) const
{
    const fvMesh& mesh(phase1.mesh());

    return volScalarField::New
    (
        "f",
        mesh,
        dimensionedScalar(dimless, phase2.name() == continuousPhase_)
    );
}


Foam::tmp<Foam::volScalarField> Foam::blendingMethods::noBlending::f2
(
    const phaseModel& phase1,
    const phaseModel& phase2
) const
{
    const fvMesh& mesh(phase1.mesh());

    return volScalarField::New
    (
        "f",
        mesh,
        dimensionedScalar(dimless, phase1.name() == continuousPhase_)
    );
}
```

* for f1():
  * get mesh
  * define and return $f_1$
    * if phase2_ is continuousPhase_, then 
      * $f_1 = 1$
    * else phase2_ is not continuousPhase_, then 
      * $f_1 = 0$
* for f2():
  * get mesh
  * define and return $f_2$
    * if phase1_ is continuousPhase_, then 
      * $f_2 = 1$
    * else phase1_ is not continuousPhase_, then 
      * $f_2 = 0$

* If phase1_ is continuousPhase_, phase2_ is not continuousPhase_, then 
  * $f_1 = 0$, $f_2 = 1$
* If phase1_ is not continuousPhase_, phase2_ is continuousPhase_, then 
  * $f_1 = 1$, $f_2 = 0$

ATTENTION: this should be paid additional attention

### linear

#### linear.H

##### include

```cpp
#ifndef linear_H
#define linear_H

#include "blendingMethod.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace blendingMethods
{
```

##### private data

```cpp
        //- Minimum fraction of phases which can be considered fully continuous
        HashTable<dimensionedScalar, word, word::hash>
            minFullyContinuousAlpha_;

        //- Minimum fraction of phases which can be considered partly continuous
        HashTable<dimensionedScalar, word, word::hash>
            minPartlyContinuousAlpha_;
```

define `minFullyContinuousAlpha_` and `minPartlyContinuousAlpha_`

##### name and constructor and destructor

```cpp
    //- Runtime type information
    TypeName("linear");


    // Constructors

        //- Construct from a dictionary and two phases
        linear
        (
            const dictionary& dict,
            const wordList& phaseNames
        );


    //- Destructor
    ~linear();
```

##### f1() and f2()

```cpp
    // Member Functions

        //- Factor for primary phase
        virtual tmp<volScalarField> f1
        (
            const phaseModel& phase1,
            const phaseModel& phase2
        ) const;

        //- Factor for secondary phase
        virtual tmp<volScalarField> f2
        (
            const phaseModel& phase1,
            const phaseModel& phase2
        ) const;
```

#### linear.C

##### include

```cpp
#include "linear.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace blendingMethods
{
    defineTypeNameAndDebug(linear, 0);

    addToRunTimeSelectionTable
    (
        blendingMethod,
        linear,
        dictionary
    );
}
}
```

##### constructor and destructor

```cpp
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::blendingMethods::linear::linear
(
    const dictionary& dict,
    const wordList& phaseNames
)
:
    blendingMethod(dict)
{
    forAllConstIter(wordList, phaseNames, iter)
    {
        const word nameFull
        (
            IOobject::groupName("minFullyContinuousAlpha", *iter)
        );

        minFullyContinuousAlpha_.insert
        (
            *iter,
            dimensionedScalar
            (
                nameFull,
                dimless,
                dict.lookup(nameFull)
            )
        );

        const word namePart
        (
            IOobject::groupName("minPartlyContinuousAlpha", *iter)
        );

        minPartlyContinuousAlpha_.insert
        (
            *iter,
            dimensionedScalar
            (
                namePart,
                dimless,
                dict.lookup(namePart)
            )
        );

        if
        (
            minFullyContinuousAlpha_[*iter]
          < minPartlyContinuousAlpha_[*iter]
        )
        {
            FatalErrorInFunction
                << "The supplied fully continuous volume fraction for "
                << *iter
                << " is less than the partly continuous value."
                << endl << exit(FatalError);
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::blendingMethods::linear::~linear()
{}
```

for constructor:

* initialized with dict
* for every phase names in phaseNames:
  * read values from dict and insert them to `minFullyContinuousAlpha_` or `minPartlyContinuousAlpha_`
  * if `minFullyContinuousAlpha_` < `minPartlyContinuousAlpha_`
    * ERROR, because fully should not be smaller than partly

##### f1() and f2()

```cpp
Foam::tmp<Foam::volScalarField> Foam::blendingMethods::linear::f1
(
    const phaseModel& phase1,
    const phaseModel& phase2
) const
{
    const dimensionedScalar
        minFullAlpha(minFullyContinuousAlpha_[phase2.name()]);
    const dimensionedScalar
        minPartAlpha(minPartlyContinuousAlpha_[phase2.name()]);

    return
        min
        (
            max
            (
                (phase2 - minPartAlpha)
               /(minFullAlpha - minPartAlpha + small),
                scalar(0)
            ),
            scalar(1)
        );
}


Foam::tmp<Foam::volScalarField> Foam::blendingMethods::linear::f2
(
    const phaseModel& phase1,
    const phaseModel& phase2
) const
{
    const dimensionedScalar
        minFullAlpha(minFullyContinuousAlpha_[phase1.name()]);
    const dimensionedScalar
        minPartAlpha(minPartlyContinuousAlpha_[phase1.name()]);

    return
        min
        (
            max
            (
                (phase1 - minPartAlpha)
               /(minFullAlpha - minPartAlpha + small),
                scalar(0)
            ),
            scalar(1)
        );
}
```

* for f1()
  * get `minFullyContinuousAlpha_` and `minPartlyContinuousAlpha_` of phase2_ as minFullAlpha and minPartAlpha
  * define and return $f_1$ as
  * $$f_1 = \min(\max(\frac{\alpha_2 - minPartAlpha}{minFullAlpha - minPartAlpha + small}, 0), 1)$$
* for f2()
  * get `minFullyContinuousAlpha_` and `minPartlyContinuousAlpha_` of phase1_ as minFullAlpha and minPartAlpha
  * define and return $f_1$ as
  * $$f_1 = \min(\max(\frac{\alpha_1 - minPartAlpha}{minFullAlpha - minPartAlpha + small}, 0), 1)$$

according to `minFullyContinuousAlpha_` and `minPartlyContinuousAlpha_` and $\alpha$ to decide whether a phase is continuous or dispersed

* take f1() as example 
  * if $\alpha_2 < minPartAlpha$, phase2_ is not even partly continuous
    * $f_1 = 0$
  * if $minPartAlpha < \alpha_2 < minFullAlpha$, phase2_ is between fully and partly continuous
    * $f_1 = \frac{\alpha_2 - minPartAlpha}{minFullAlpha - minPartAlpha + small}$, which is between [0, 1]
  * if $\alpha_2 > minFullAlpha$, phase2_ is fully continuous
    * $f_1 = 1$, which is the same with `noBlending`
  
### hyperbolic

#### hyperbolic.H

##### private data

```cpp
        //- Minimum fraction of phases which can be considered continuous
        HashTable<dimensionedScalar, word, word::hash> minContinuousAlpha_;

        //- Width of the transition
        const dimensionedScalar transitionAlphaScale_;
```

define `minContinuousAlpha_` and `transitionAlphaScale_`

#### hyperbolic.C

##### f1() and f2()

```cpp
Foam::tmp<Foam::volScalarField> Foam::blendingMethods::hyperbolic::f1
(
    const phaseModel& phase1,
    const phaseModel& phase2
) const
{
    return
        (
            1
          + tanh
            (
                (4/transitionAlphaScale_)
               *(phase2 - minContinuousAlpha_[phase2.name()])
            )
        )/2;
}


Foam::tmp<Foam::volScalarField> Foam::blendingMethods::hyperbolic::f2
(
    const phaseModel& phase1,
    const phaseModel& phase2
) const
{
    return
        (
            1
          + tanh
            (
                (4/transitionAlphaScale_)
               *(phase1 - minContinuousAlpha_[phase1.name()])
            )
        )/2;
}
```

* for f1()
  * $$f_1 = \frac{1}{2}\left[1 + \tanh(\frac{4}{transitionAlphaScale\_}(\alpha_2 - minContinuousAlpha\__{phase2\_}))\right]$$
* for f2()
  * $$f_1 = \frac{1}{2}\left[1 + \tanh(\frac{4}{transitionAlphaScale\_}(\alpha_1 - minContinuousAlpha\__{phase1\_}))\right]$$

## Discussion

take noBlending as example, the coefficients are calculated as

* If phase1_ is continuousPhase_, phase2_ is not continuousPhase_, then 
  * $f_1 = 0$, $f_2 = 1$
* If phase1_ is not continuousPhase_, phase2_ is continuousPhase_, then 
  * $f_1 = 1$, $f_2 = 0$

then the coefficients are used in blending for three kind of model

* no dispersed and continuous phase, say, `model_`, then
  * the model is multiplied with a coefficient of $1 - f_1 - f_2$
* if phase1_ is dispersed and phase2_ is continuous, then
  * the model is multiplied with a coefficient of $f_1$, for no blending, it's $1$
* if phase2_ is dispersed and phase1_ is continuous, then
  * the model is multiplied with a coefficient of $f_2$, for no blending, it's also $1$

in summary, blending is to using models according to the relationship between the two phases


