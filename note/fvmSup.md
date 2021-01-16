# fvmSup

fvmSup is to calculate the matrix for implicit and explicit sources. In general, there are three kinds of sources: 
* explicit source, represented by `Su`; 
* implicit source, represented by `Sp`;
* implicit or explicit source depending on sign of coefficient, represented by `SuSp`.

## fvmSup.H

### include

```cpp
#ifndef fvmSup_H
#define fvmSup_H

#include "volFieldsFwd.H"
#include "fvMatrix.H"
#include "zeroField.H"
```

### name space

`Foam::fvm`

### Explicit source

```cpp
// Explicit source

    template<class Type>
    tmp<fvMatrix<Type>> Su
    (
        const DimensionedField<Type, volMesh>&,
        const GeometricField<Type, fvPatchField, volMesh>&
    );

    template<class Type>
    tmp<fvMatrix<Type>> Su
    (
        const tmp<DimensionedField<Type, volMesh>>&,
        const GeometricField<Type, fvPatchField, volMesh>&
    );

    template<class Type>
    tmp<fvMatrix<Type>> Su
    (
        const tmp<GeometricField<Type, fvPatchField, volMesh>>&,
        const GeometricField<Type, fvPatchField, volMesh>&
    );

    template<class Type>
    zeroField Su
    (
        const zero&,
        const GeometricField<Type, fvPatchField, volMesh>&
    );
```

There are four different definitions of `Su`.

| No. | type of function | type of Parameter No. 1 | type of Parameter No. 2 |
| - | - | - | - |
| 1 | tmp<fvMatrix<Type>> | const DimensionedField<Type, volMesh>& | const GeometricField<Type, fvPatchField, volMesh>& |
| 2 | tmp<fvMatrix<Type>> | const tmp<DimensionedField<Type, volMesh>>& | const GeometricField<Type, fvPatchField, volMesh>& |
| 3 | tmp<fvMatrix<Type>> | const tmp<GeometricField<Type, fvPatchField, volMesh>>& | const GeometricField<Type, fvPatchField, volMesh>& |
| 4 | zeroField | const zero& | const GeometricField<Type, fvPatchField, volMesh>& |

It can be found that types of parameter 1 are different for the four `Su` while those of parameter 2 are the same.


### Implicit source

```cpp
// Implicit source

    template<class Type>
    tmp<fvMatrix<Type>> Sp
    (
        const volScalarField::Internal&,
        const GeometricField<Type, fvPatchField, volMesh>&
    );

    template<class Type>
    tmp<fvMatrix<Type>> Sp
    (
        const tmp<volScalarField::Internal>&,
        const GeometricField<Type, fvPatchField, volMesh>&
    );

    template<class Type>
    tmp<fvMatrix<Type>> Sp
    (
        const tmp<volScalarField>&,
        const GeometricField<Type, fvPatchField, volMesh>&
    );


    template<class Type>
    tmp<fvMatrix<Type>> Sp
    (
        const dimensionedScalar&,
        const GeometricField<Type, fvPatchField, volMesh>&
    );


    template<class Type>
    zeroField Sp
    (
        const zero&,
        const GeometricField<Type, fvPatchField, volMesh>&
    );
```

There are five different definitions of `Su`.

| No. | type of function | type of Parameter No. 1 | type of Parameter No. 2 |
| - | - | - | - |
| 1 | tmp<fvMatrix<Type>> | const volScalarField::Internal& | const GeometricField<Type, fvPatchField, volMesh>& |
| 2 | tmp<fvMatrix<Type>> | const tmp<volScalarField::Internal>& | const GeometricField<Type, fvPatchField, volMesh>& |
| 3 | tmp<fvMatrix<Type>> | const tmp<volScalarField>& | const GeometricField<Type, fvPatchField, volMesh>& |
| 4 | tmp<fvMatrix<Type>> | const dimensionedScalar& | const GeometricField<Type, fvPatchField, volMesh>& |
| 5 | zeroField | const zero& | const GeometricField<Type, fvPatchField, volMesh>& |

It can be found that types of parameter 1 are different for the five `Sp` while those of parameter 2 are the same.

### Implicit/Explicit source depending on sign of coefficient

```cpp
// Implicit/Explicit source depending on sign of coefficient

    template<class Type>
    tmp<fvMatrix<Type>> SuSp
    (
        const volScalarField::Internal&,
        const GeometricField<Type, fvPatchField, volMesh>&
    );

    template<class Type>
    tmp<fvMatrix<Type>> SuSp
    (
        const tmp<volScalarField::Internal>&,
        const GeometricField<Type, fvPatchField, volMesh>&
    );

    template<class Type>
    tmp<fvMatrix<Type>> SuSp
    (
        const tmp<volScalarField>&,
        const GeometricField<Type, fvPatchField, volMesh>&
    );

    template<class Type>
    zeroField SuSp
    (
        const zero&,
        const GeometricField<Type, fvPatchField, volMesh>&
    );
```

There are four different definitions of `Su`.

| No. | type of function | type of Parameter No. 1 | type of Parameter No. 2 |
| - | - | - | - |
| 1 | tmp<fvMatrix<Type>> | const volScalarField::Internal& | const GeometricField<Type, fvPatchField, volMesh>& |
| 2 | tmp<fvMatrix<Type>> | const tmp<volScalarField::Internal>& | const GeometricField<Type, fvPatchField, volMesh>& |
| 3 | tmp<fvMatrix<Type>> | const tmp<volScalarField>& | const GeometricField<Type, fvPatchField, volMesh>& |
| 4 | zeroField | const zero& | const GeometricField<Type, fvPatchField, volMesh>& |

It can be found that types of parameter 1 are different for the four `SuSp` while those of parameter 2 are the same.

### Questions

Why these three kinds of source have different type of parameter No. 1?

## fvmSup.C

### Explicit Source Su

#### definition of Su 1

```cpp
template<class Type>
Foam::tmp<Foam::fvMatrix<Type>>
Foam::fvm::Su
(
    const DimensionedField<Type, volMesh>& su,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    const fvMesh& mesh = vf.mesh();

    tmp<fvMatrix<Type>> tfvm
    (
        new fvMatrix<Type>
        (
            vf,
            dimVol*su.dimensions()
        )
    );
    fvMatrix<Type>& fvm = tfvm.ref();

    fvm.source() -= mesh.V()*su.field();

    return tfvm;
}
```

* get `mesh`
* define a temporary matrix `tfvm`, as
  * $$tfvm = vf$$, with the unit of volume times su (the first parameter)
* define a reference of `tfvm` as `fvm`
* set the source of `fvm` as:
  * $$fvm.source() = fvm.source() - V \cdot su$$
* return `tfvm`, whose source has been replaced with $$fvm.source() - V \cdot su$$

`tmp` is a class for managing temporary objects

since `tfvm` is a `tmp` variable, so `tfvm.ref()` can be found in `src\OpenFOAM\memory\tmpNrc\tmpNrcI.H` or `src\OpenFOAM\memory\tmpNrc\tmpNrc.H`, is to return non-const reference or generate a fatal error if the object is const.

#### definition of Su 2 and 3

```cpp
template<class Type>
Foam::tmp<Foam::fvMatrix<Type>>
Foam::fvm::Su
(
    const tmp<DimensionedField<Type, volMesh>>& tsu,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    tmp<fvMatrix<Type>> tfvm = fvm::Su(tsu(), vf);
    tsu.clear();
    return tfvm;
}


template<class Type>
Foam::tmp<Foam::fvMatrix<Type>>
Foam::fvm::Su
(
    const tmp<GeometricField<Type, fvPatchField, volMesh>>& tsu,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    tmp<fvMatrix<Type>> tfvm = fvm::Su(tsu(), vf);
    tsu.clear();
    return tfvm;
}
```

Su 2 and 3 are just using Su 1 with different parameters

#### definition of Su 4

```cpp
template<class Type>
Foam::zeroField
Foam::fvm::Su
(
    const zero&,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    return zeroField();
}
```

* return a zero field


### Implicit Source Sp

#### definition of Sp 1

```cpp
template<class Type>
Foam::tmp<Foam::fvMatrix<Type>>
Foam::fvm::Sp
(
    const volScalarField::Internal& sp,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    const fvMesh& mesh = vf.mesh();

    tmp<fvMatrix<Type>> tfvm
    (
        new fvMatrix<Type>
        (
            vf,
            dimVol*sp.dimensions()*vf.dimensions()
        )
    );
    fvMatrix<Type>& fvm = tfvm.ref();

    fvm.diag() += mesh.V()*sp.field();

    return tfvm;
}
```

* get `mesh`
* define a temporary matrix `tfvm`, as
  * $$tfvm = vf$$, with the unit of volume times sp (the first parameter) times vf (the second parameter)
* define a reference of `tfvm` as `fvm`
* set the diagonal elements of `fvm` as:
  * $$fvm.diag() = fvm.diag() + V \cdot sp$$
* return `tfvm`, whose diagonal has been replaced with $$fvm.diag() + V \cdot sp$$

Why the dimension is set to `dimVol*sp.dimensions()*vf.dimensions()`?

It should be noted that the type of the first parameter `sp` is `const volScalarField::Internal&`, which is a scalar. So, the dimension in `tfvm` is reasonable.

#### definition of Sp 2 and 3

```cpp
template<class Type>
Foam::tmp<Foam::fvMatrix<Type>>
Foam::fvm::Sp
(
    const tmp<volScalarField::Internal>& tsp,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    tmp<fvMatrix<Type>> tfvm = fvm::Sp(tsp(), vf);
    tsp.clear();
    return tfvm;
}


template<class Type>
Foam::tmp<Foam::fvMatrix<Type>>
Foam::fvm::Sp
(
    const tmp<volScalarField>& tsp,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    tmp<fvMatrix<Type>> tfvm = fvm::Sp(tsp(), vf);
    tsp.clear();
    return tfvm;
}
```

Sp 2 and 3 are just using Sp 1 with different parameters

#### definition of Sp 4

```cpp
template<class Type>
Foam::tmp<Foam::fvMatrix<Type>>
Foam::fvm::Sp
(
    const dimensionedScalar& sp,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    const fvMesh& mesh = vf.mesh();

    tmp<fvMatrix<Type>> tfvm
    (
        new fvMatrix<Type>
        (
            vf,
            dimVol*sp.dimensions()*vf.dimensions()
        )
    );
    fvMatrix<Type>& fvm = tfvm.ref();

    fvm.diag() += mesh.V()*sp.value();

    return tfvm;
}
```

almost the same with Sp 1, except that the type of the first paramter is different

#### definition of Sp 5

```cpp
template<class Type>
Foam::zeroField
Foam::fvm::Sp
(
    const zero&,
    const GeometricField<Type, fvPatchField, volMesh>&
)
{
    return zeroField();
}
```

* return a zero field

### Implicit/Explicit source depending on sign of coefficient SuSp

#### definition of SuSp 1

```cpp
template<class Type>
Foam::tmp<Foam::fvMatrix<Type>>
Foam::fvm::SuSp
(
    const volScalarField::Internal& susp,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    const fvMesh& mesh = vf.mesh();

    tmp<fvMatrix<Type>> tfvm
    (
        new fvMatrix<Type>
        (
            vf,
            dimVol*susp.dimensions()*vf.dimensions()
        )
    );
    fvMatrix<Type>& fvm = tfvm.ref();

    fvm.diag() += mesh.V()*max(susp.field(), scalar(0));

    fvm.source() -= mesh.V()*min(susp.field(), scalar(0))
        *vf.primitiveField();

    return tfvm;
}
```

* get `mesh`
* define a temporary matrix `tfvm`, as
  * $$tfvm = vf$$, with the unit of volume times susp (the first parameter) times vf (the second parameter)
* define a reference of `tfvm` as `fvm`
* set the diagonal elements of `fvm` as:
  * $$fvm.diag() = fvm.diag() + V \cdot \max(susp, 0)$$
  * namely, if $susp$ is positive, then it is added to the diagonal elements of `fvm`; if it's negative, then it is not added 
* set the source of `fvm` as:
  * $$fvm.source() = fvm.source() - V \cdot \min(susp, 0) \cdot vf$$
  * namely, if $susp$ is negative, then it is added to the source, else, it is add to the diagonal as above
* return `tfvm`

`vf.primitiveField()` can be found in `src\OpenFOAM\fields\GeometricFields\GeometricField\GeometricField.H`, is to return a const-reference to the  internal field

It should be noted that the type of the first parameter `susp` is `const volScalarField::Internal&`, which is a scalar. So, the dimension in `tfvm` is reasonable.

#### definition of SuSp 2 and 3

```cpp
template<class Type>
Foam::tmp<Foam::fvMatrix<Type>>
Foam::fvm::SuSp
(
    const tmp<volScalarField::Internal>& tsusp,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    tmp<fvMatrix<Type>> tfvm = fvm::SuSp(tsusp(), vf);
    tsusp.clear();
    return tfvm;
}


template<class Type>
Foam::tmp<Foam::fvMatrix<Type>>
Foam::fvm::SuSp
(
    const tmp<volScalarField>& tsusp,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    tmp<fvMatrix<Type>> tfvm = fvm::SuSp(tsusp(), vf);
    tsusp.clear();
    return tfvm;
}
```

SuSp 2 and 3 are just using SuSp 1 with different parameters

#### definition of SuSp 4

```cpp
template<class Type>
Foam::zeroField
Foam::fvm::SuSp
(
    const zero&,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    return zeroField();
}
```

* return a zero field