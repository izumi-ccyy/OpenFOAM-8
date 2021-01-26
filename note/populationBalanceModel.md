# Population Balance Model

## PopulationBalancePhaseSystem

### PopulationBalancePhaseSystem.H

#### include

```cpp
#ifndef PopulationBalancePhaseSystem_H
#define PopulationBalancePhaseSystem_H

#include "phaseSystem.H"
#include "populationBalanceModel.H"
```

#### template

```cpp
template<class BasePhaseSystem>
class PopulationBalancePhaseSystem
:
    public BasePhaseSystem
{
    ...
}
```

define a class template

#### Private data

```cpp
    // Private data

        //- Population balances
        PtrList<diameterModels::populationBalanceModel> populationBalances_;

        //- Mass transfer rates
        phaseSystem::dmdtfTable dmdtfs_;
```

define private data

* `populationBalances_`: population balances
* `dmdtfs`: mass transfer rates

#### public

```cpp
public:

    // Constructors

        //- Construct from fvMesh
        PopulationBalancePhaseSystem(const fvMesh&);


    //- Destructor
    virtual ~PopulationBalancePhaseSystem();


    // Member Functions

        //- Return the mass transfer rate for an interface
        virtual tmp<volScalarField> dmdtf(const phasePairKey& key) const;

        //- Return the mass transfer rates for each phase
        virtual PtrList<volScalarField> dmdts() const;

        //- Return the momentum transfer matrices for the cell-based algorithm
        virtual autoPtr<phaseSystem::momentumTransferTable> momentumTransfer();

        //- Return the momentum transfer matrices for the face-based algorithm
        virtual autoPtr<phaseSystem::momentumTransferTable> momentumTransferf();

        //- Return the heat transfer matrices
        virtual autoPtr<phaseSystem::heatTransferTable> heatTransfer() const;

        //- Return the specie transfer matrices
        virtual autoPtr<phaseSystem::specieTransferTable>
            specieTransfer() const;

        //- Read base phaseProperties dictionary
        virtual bool read();

        //- Solve all population balance equations
        virtual void solve
        (
            const PtrList<volScalarField>& rAUs,
            const PtrList<surfaceScalarField>& rAUfs
        );
```

* define constructors and destructor
* member function: all defined as `virtual`
  * `dmdtf()`: Return the mass transfer rate for an interface
  * `dmdts()`: Return the mass transfer rates for each phase
  * `momentumTransfer()`: Return the momentum transfer matrices for the cell-based algorithm
  * `momentumTransferf()`: Return the momentum transfer matrices for the face-based algorithm
  * `heatTransfer()`: Return the heat transfer matrices
  * `speciesTransfer()`: Return the specie transfer matrices
  * `read()`: Read base phaseProperties dictionary
  * `solve()`: Solve all population balance equations 

### PopulationBalancePhaseSystem.C

#### constructors

#### destructor

#### member function

##### dmdtf()

```cpp
template<class BasePhaseSystem>
Foam::tmp<Foam::volScalarField>
Foam::PopulationBalancePhaseSystem<BasePhaseSystem>::dmdtf
(
    const phasePairKey& key
) const
{
    tmp<volScalarField> tDmdt = BasePhaseSystem::dmdtf(key);

    if (!dmdtfs_.found(key))
    {
        const label dmdtSign(Pair<word>::compare(this->phasePairs_[key], key));

        tDmdt.ref() += dmdtSign**dmdtfs_[key];
    }

    return tDmdt;
}
```

* define `tDmdt` as `dmdtf`, the mass transfer rate for an interface
* get the sign according to the phase
* set `tDmdt`
  * $$tDmdt = tDmdt + sign dmdtfs_[key]$$
  * which is the mass transfer rate for an interface plus mass transfer rate (sign is considered)

`compare()` can be found in `src\OpenFOAM\primitives\Pair\Pair.H`, is to compare pairs and return a sign

```cpp
        //- Compare Pairs
        //  Returning:
        //  -  0: different
        //  - +1: identical
        //  - -1: same pair, but reversed order
        static inline int compare(const Pair<Type>& a, const Pair<Type>& b)
        {
            if (a == b)
            {
                return 1;
            }
            else if (a == reverse(b))
            {
                return -1;
            }
            else
            {
                return 0;
            }
        }
```

##### dmdts()

```cpp
template<class BasePhaseSystem>
Foam::PtrList<Foam::volScalarField>
Foam::PopulationBalancePhaseSystem<BasePhaseSystem>::dmdts() const
{
    PtrList<volScalarField> dmdts(BasePhaseSystem::dmdts());

    forAllConstIter(phaseSystem::dmdtfTable, dmdtfs_, dmdtfIter)
    {
        const phasePair& pair = this->phasePairs_[dmdtfIter.key()];
        const volScalarField& pDmdt = *dmdtfIter();

        addField(pair.phase1(), "dmdt", pDmdt, dmdts);
        addField(pair.phase2(), "dmdt", - pDmdt, dmdts);
    }

    return dmdts;
}
```

* get `dmdts`
* for every `dmdtf` in `dmdtfs`
  * get `pair`
  * get `dmdtfIter()` as `pDmdt`, in fact `dmdtfIter()` is $dmdtf^k$
  * add `pDmdt`, $dmdtf^k$, to `dmdt` for `phase1` of the `pair`
  * add `-pDmdt`, $-dmdtf^k$, to `dmdt` for `phase2` of the `pair`
  * then for `phase1`
    * $$dmdt^k_1 = dmdtf^k_1 - dmdtf^k_2$$
    * $$dmdt^k_2 = dmdtf^k_2 - dmdtf^k_1$$
    * where $k$th pair

`forAllConstIter` can be found in `src\OpenFOAM\containers\Lists\UList\UList.H`, as

```cpp
#define forAllConstIter(Container,container,iter)                              \
    for                                                                        \
    (                                                                          \
        Container::const_iterator iter = (container).begin();                  \
        iter != (container).end();                                             \
        ++iter                                                                 \
    )
```

##### momentumTransfer()

```cpp
template<class BasePhaseSystem>
Foam::autoPtr<Foam::phaseSystem::momentumTransferTable>
Foam::PopulationBalancePhaseSystem<BasePhaseSystem>::momentumTransfer()
{
    autoPtr<phaseSystem::momentumTransferTable> eqnsPtr =
        BasePhaseSystem::momentumTransfer();

    phaseSystem::momentumTransferTable& eqns = eqnsPtr();

    this->addDmdtUfs(dmdtfs_, eqns);

    return eqnsPtr;
}
```

* define a pointer `eqnsPtr` as `momentumTransfer()`
* define a `momentumTransferTable` `eqns` as `eqnsPtr()`
* add momentum transfer term of `eqns` to `dmdtfs_`
* return `momentumTransfer()`

`addDmdtUfs` can be found in `applications\solvers\multiphase\multiphaseEulerFoam\phaseSystems\PhaseSystems\MomentumTransferPhaseSystem\MomentumTransferPhaseSystem.H`, is to dd momentum transfer terms which result from bulk mass transfers

```cpp
        //- Add momentum transfer terms which result from bulk mass transfers
        void addDmdtUfs
        (
            const phaseSystem::dmdtfTable& dmdtfs,
            phaseSystem::momentumTransferTable& eqns
        );
```

##### momentumTransferf()

```cpp
template<class BasePhaseSystem>
Foam::autoPtr<Foam::phaseSystem::momentumTransferTable>
Foam::PopulationBalancePhaseSystem<BasePhaseSystem>::momentumTransferf()
{
    autoPtr<phaseSystem::momentumTransferTable> eqnsPtr =
        BasePhaseSystem::momentumTransferf();

    phaseSystem::momentumTransferTable& eqns = eqnsPtr();

    this->addDmdtUfs(dmdtfs_, eqns);

    return eqnsPtr;
}
```

similarly, is to add momentum transfer term to `dmdtfs_`

##### heatTransfer()

```cpp
template<class BasePhaseSystem>
Foam::autoPtr<Foam::phaseSystem::heatTransferTable>
Foam::PopulationBalancePhaseSystem<BasePhaseSystem>::heatTransfer() const
{
    autoPtr<phaseSystem::heatTransferTable> eqnsPtr =
        BasePhaseSystem::heatTransfer();

    phaseSystem::heatTransferTable& eqns = eqnsPtr();

    this->addDmdtHefs(dmdtfs_, eqns);

    return eqnsPtr;
}
```

similarly, is to add heat transfer term to `dmdtfs_`

##### specieTransfer()

```cpp
template<class BasePhaseSystem>
Foam::autoPtr<Foam::phaseSystem::specieTransferTable>
Foam::PopulationBalancePhaseSystem<BasePhaseSystem>::specieTransfer() const
{
    autoPtr<phaseSystem::specieTransferTable> eqnsPtr =
        BasePhaseSystem::specieTransfer();

    phaseSystem::specieTransferTable& eqns = eqnsPtr();

    this->addDmdtYfs(dmdtfs_, eqns);

    return eqnsPtr;
}

```

similarly, is to add specie transfer term to `dmdtfs_`

##### read()

```cpp
template<class BasePhaseSystem>
bool Foam::PopulationBalancePhaseSystem<BasePhaseSystem>::read()
{
    if (BasePhaseSystem::read())
    {
        bool readOK = true;

        // Models ...

        return readOK;
    }
    else
    {
        return false;
    }
}
```

read base phaseProperties dictionary

##### solve()

```cpp
template<class BasePhaseSystem>
void Foam::PopulationBalancePhaseSystem<BasePhaseSystem>::solve
(
    const PtrList<volScalarField>& rAUs,
    const PtrList<surfaceScalarField>& rAUfs
)
{
    BasePhaseSystem::solve(rAUs, rAUfs);

    forAll(populationBalances_, i)
    {
        populationBalances_[i].solve();
    }
}
```

solve all population balance equations

