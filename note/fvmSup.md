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

There are four different definations of `Su`.

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

There are five different definations of `Su`.

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

There are four different definations of `Su`.

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

