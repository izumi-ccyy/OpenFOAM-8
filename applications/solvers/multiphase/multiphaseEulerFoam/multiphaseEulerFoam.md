# multiphaseEulerFoam

## pphaseSystems

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