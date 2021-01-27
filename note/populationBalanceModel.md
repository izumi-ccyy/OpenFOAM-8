# Population Balance Model

- [Population Balance Model](#population-balance-model)
  - [PopulationBalancePhaseSystem](#populationbalancephasesystem)
    - [PopulationBalancePhaseSystem.H](#populationbalancephasesystemh)
      - [include](#include)
      - [template](#template)
      - [Private data](#private-data)
      - [public](#public)
    - [PopulationBalancePhaseSystem.C](#populationbalancephasesystemc)
      - [constructors](#constructors)
      - [destructor](#destructor)
      - [member function](#member-function)
        - [dmdtf()](#dmdtf)
        - [dmdts()](#dmdts)
        - [momentumTransfer()](#momentumtransfer)
        - [momentumTransferf()](#momentumtransferf)
        - [heatTransfer()](#heattransfer)
        - [specieTransfer()](#specietransfer)
        - [read()](#read)
        - [solve()](#solve)
  - [populationBalanceModel](#populationbalancemodel)
    - [populationBalanceModel](#populationbalancemodel-1)
      - [populationBalanceModel.H](#populationbalancemodelh)
        - [description](#description)
        - [include](#include-1)
        - [namespace](#namespace)
        - [class](#class)
        - [private data](#private-data-1)
        - [private member function](#private-member-function)
        - [public](#public-1)
        - [public member function](#public-member-function)
      - [populationBalanceModelI.H](#populationbalancemodelih)
      - [populationBalanceModel.C](#populationbalancemodelc)
        - [include](#include-2)
        - [Static Data Members](#static-data-members)
        - [Private Member Functions](#private-member-functions)
          - [registerVelocityGroups()](#registervelocitygroups)
          - [registerSizeGroups()](#registersizegroups)
          - [createPhasePairs()](#createphasepairs)
          - [correct()](#correct)
          - [birthByCoalescence()](#birthbycoalescence)
          - [deathByCoalescence()](#deathbycoalescence)
          - [](#)
          - [](#-1)
          - [](#-2)
          - [](#-3)
          - [](#-4)
          - [](#-5)
          - [](#-6)
          - [](#-7)
          - [](#-8)
          - [](#-9)
          - [](#-10)
          - [](#-11)
        - [public member functions](#public-member-functions)
          - [clone()](#clone)
          - [writeData()](#writedata)
          - [eta()](#eta)
          - [](#-12)
          - [](#-13)

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

```cpp
template<class BasePhaseSystem>
Foam::PopulationBalancePhaseSystem<BasePhaseSystem>::
PopulationBalancePhaseSystem
(
    const fvMesh& mesh
)
:
    BasePhaseSystem(mesh),

    populationBalances_
    (
        this->lookup("populationBalances"),
        diameterModels::populationBalanceModel::iNew(*this, dmdtfs_)
    )
{
    forAll(populationBalances_, i)
    {
        const Foam::diameterModels::populationBalanceModel& popBal =
            populationBalances_[i];

        forAllConstIter(phaseSystem::phasePairTable, popBal.phasePairs(), iter)
        {
            const phasePairKey& key = iter.key();

            if (!this->phasePairs_.found(key))
            {
                this->phasePairs_.insert
                (
                    key,
                    autoPtr<phasePair>
                    (
                        new phasePair
                        (
                            this->phaseModels_[key.first()],
                            this->phaseModels_[key.second()]
                        )
                    )
                );
            }

            dmdtfs_.insert
            (
                key,
                new volScalarField
                (
                    IOobject
                    (
                        IOobject::groupName
                        (
                            "populationBalance:dmdtf",
                            this->phasePairs_[key]->name()
                        ),
                        this->mesh().time().timeName(),
                        this->mesh(),
                        IOobject::READ_IF_PRESENT,
                        IOobject::AUTO_WRITE
                    ),
                    this->mesh(),
                    dimensionedScalar(dimDensity/dimTime, 0)
                )
            );
        }
    }
}
```

* initialize with `mesh` and `populationBalances_` found in dictionary
* create `phasePairs`
* initialize and insert `dmdtfs`

#### destructor

```cpp
template<class BasePhaseSystem>
Foam::PopulationBalancePhaseSystem<BasePhaseSystem>::
~PopulationBalancePhaseSystem()
{}
```

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

## populationBalanceModel

### populationBalanceModel

#### populationBalanceModel.H

##### description

Class that solves the univariate population balance equation by means of a **class method** (also called sectional or discrete method). The internal coordinate is set to the particle volume, so the equation is based on a transport equation of the volume-based number density function. The discretization is done using the fixed pivot technique of Kumar and Ramkrishna (1996). The source terms are written in a way that particle  number and mass are preserved. Coalescence (aggregation), breakup, drift (growth and surface loss) as well as nucleation are supported. For the discrete breakup term two recipes are available, depending on the model choice. For models which state a total breakup rate and a separate daughter size distribution function, the formulation of Kumar and Ramkrishna (1996) is applied which is applicable for binary and multiple breakup events. The second formulation is given by Liao et al. (2018). It is useful for binary breakup models which give the breakup rate between a sizeGroup pair directly, without an explicit expression for the daughter size distribution. The drift term is implemented using a finite difference upwind scheme. Although it is diffusive, it ensures a stable and number-conservative solution.

##### include

```cpp
#ifndef populationBalanceModel_H
#define populationBalanceModel_H

#include "sizeGroup.H"
#include "phasePair.H"
#include "pimpleControl.H"
#include "phaseCompressibleMomentumTransportModelFwd.H"
#include "HashPtrTable.H"
```

* include `sizeGroup.H` etc.

##### namespace

```cpp
namespace Foam
{

class phaseSystem;

namespace diameterModels
{

class coalescenceModel;
class breakupModel;
class binaryBreakupModel;
class driftModel;
class nucleationModel;
```

define `namespace`

##### class

```cpp
class populationBalanceModel
:
    public regIOobject
{
    ...
}
```

define the class inherit from `regIOobject`

##### private data

```cpp
private:

    // Private Typedefs

        typedef
            HashTable<autoPtr<phasePair>, phasePairKey, phasePairKey::hash>
            phasePairTable;

        typedef
            HashPtrTable<volScalarField, phasePairKey, phasePairKey::hash>
            pDmdtTable;


    // Private Data

        //- Reference to the phaseSystem
        const phaseSystem& fluid_;

        //- Interfacial mass transfer rates
        pDmdtTable& pDmdt_;

        //- Reference to the mesh
        const fvMesh& mesh_;

        //- Name of the populationBalance
        word name_;

        //- Dictionary
        dictionary dict_;

        //- Reference to pimpleControl
        const pimpleControl& pimple_;

        //- Continuous phase
        const phaseModel& continuousPhase_;

        //- velocityGroups belonging to this populationBalance
        UPtrList<velocityGroup> velocityGroups_;

        //- sizeGroups belonging to this populationBalance
        UPtrList<sizeGroup> sizeGroups_;

        //- List of unordered phasePairs in this populationBalance
        phasePairTable phasePairs_;

        //- sizeGroup boundaries
        PtrList<dimensionedScalar> v_;

        //- Section width required for binary breakup formulation
        PtrList<PtrList<dimensionedScalar>> delta_;

        //- Explicitly treated sources
        PtrList<volScalarField> Su_;

        //- Sources treated implicitly or explicitly depending on sign
        PtrList<volScalarField> SuSp_;

        //- Field for caching sources
        volScalarField Sui_;

        //- Coalescence models
        PtrList<coalescenceModel> coalescence_;

        //- Coalescence rate
        autoPtr<volScalarField> coalescenceRate_;

        //- BreakupModels
        PtrList<breakupModel> breakup_;

        //- Breakup rate
        autoPtr<volScalarField> breakupRate_;

        //- Binary breakup models
        PtrList<binaryBreakupModel> binaryBreakup_;

        //- Binary breakup rate
        autoPtr<volScalarField> binaryBreakupRate_;

        //- Drift models
        PtrList<driftModel> drift_;

        //- Drift rate
        autoPtr<volScalarField> driftRate_;

        //- Ratio between successive representative volumes
        autoPtr<volScalarField> rx_;

        //- Ratio between successive class widths
        autoPtr<volScalarField> rdx_;

        //- Zeroeth order models
        PtrList<nucleationModel> nucleation_;

        //- Zeroeth order rate
        autoPtr<volScalarField> nucleationRate_;

        //- Total void fraction
        autoPtr<volScalarField> alphas_;

        //- Mean Sauter diameter
        autoPtr<volScalarField> dsm_;

        //- Average velocity
        autoPtr<volVectorField> U_;

        //- Counter for interval between source term updates
        label sourceUpdateCounter_;
```

define typedef and private data

##### private member function

```cpp
    // Private Member Functions

        void registerVelocityGroups();

        void registerSizeGroups(sizeGroup& group);

        void createPhasePairs();

        void correct();

        void birthByCoalescence(const label j, const label k);

        void deathByCoalescence(const label i, const label j);

        void birthByBreakup(const label k, const label model);

        void deathByBreakup(const label i);

        void calcDeltas();

        void birthByBinaryBreakup(const label i, const label j);

        void deathByBinaryBreakup(const label j, const label i);

        void drift(const label i, driftModel& model);

        void nucleation(const label i, nucleationModel& model);

        void sources();

        void dmdt();

        void calcAlphas();

        tmp<volScalarField> calcDsm();

        void calcVelocity();

        //- Return whether the sources should be updated on this iteration
        bool updateSources();

        //- Return the number of corrections
        inline label nCorr() const;

        //- Return the interval at which the sources are updated
        inline label sourceUpdateInterval() const;
```

##### public

```cpp
public:

    //- Runtime type information
    TypeName("populationBalanceModel");


    // Constructor

        populationBalanceModel
        (
            const phaseSystem& fluid,
            const word& name,
            HashPtrTable
            <
                volScalarField,
                phasePairKey,
                phasePairKey::hash
            >& pDmdt
        );

        //- Return clone
        autoPtr<populationBalanceModel> clone() const;

        //- Return a pointer to a new populationBalanceModel object created on
        //  freestore from Istream
        class iNew
        {
            const phaseSystem& fluid_;

            HashPtrTable<volScalarField, phasePairKey, phasePairKey::hash>&
                pDmdt_;

        public:

            iNew
            (
                const phaseSystem& fluid,
                HashPtrTable<volScalarField, phasePairKey, phasePairKey::hash>&
                    pDmdt
            )
            :
                fluid_(fluid),
                pDmdt_(pDmdt)
            {}

            autoPtr<populationBalanceModel> operator()(Istream& is) const
            {
                return autoPtr<populationBalanceModel>
                (
                    new populationBalanceModel(fluid_, word(is), pDmdt_)
                );
            }
        };


    //- Destructor
    virtual ~populationBalanceModel();
```

define 

* type name
* constructors
  * clone: Return clone
  * iNew: Return a pointer to a new populationBalanceModel object created on freestore from Istream
* destructor

##### public member function

```cpp
    // Member Functions

        //- Dummy write for regIOobject
        bool writeData(Ostream&) const;

        //- Return reference to the phaseSystem
        inline const phaseSystem& fluid() const;

        //- Return reference to the mesh
        inline const fvMesh& mesh() const;

        //- Return populationBalanceCoeffs dictionary
        inline const dictionary& dict() const;

        //- Return continuous phase
        inline const phaseModel& continuousPhase() const;

        //- Return the velocityGroups belonging to this populationBalance
        inline const UPtrList<velocityGroup>& velocityGroups() const;

        //- Return the sizeGroups belonging to this populationBalance
        inline const UPtrList<sizeGroup>& sizeGroups() const;

        //- Return semi-implicit source terms
        inline const volScalarField& SuSp(const label i) const;

        //- Return list of unordered phasePairs in this populationBalance
        inline const phasePairTable& phasePairs() const;

        //- Returns true if both phases are velocity groups and
        //  belong to this populationBalance
        inline bool isVelocityGroupPair(const phasePair& pair) const;

        //- Return the sizeGroup boundaries
        inline const PtrList<dimensionedScalar>& v() const;

        //- Return total void of phases belonging to this populationBalance
        inline const volScalarField& alphas() const;

        //- Return average velocity
        inline const volVectorField& U() const;

        //- Return allocation coefficient
        const dimensionedScalar eta
        (
            const label i,
            const dimensionedScalar& v
        ) const;

        //- Return the surface tension coefficient between a given dispersed
        //  and the continuous phase
        const tmp<volScalarField> sigmaWithContinuousPhase
        (
            const phaseModel& dispersedPhase
        ) const;

        //- Return reference to turbulence model of the continuous phase
        const phaseCompressibleMomentumTransportModel&
            continuousTurbulence() const;

        //- Solve the population balance equation
        void solve();
```

define inline functions etc.

#### populationBalanceModelI.H

```cpp
inline Foam::label Foam::diameterModels::populationBalanceModel::nCorr() const
{
    return mesh_.solverDict(name_).lookup<label>("nCorr");
}


inline Foam::label
Foam::diameterModels::populationBalanceModel::sourceUpdateInterval() const
{
    return
        mesh_.solverDict(name_)
       .lookupOrDefault<label>("sourceUpdateInterval", 1);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline const Foam::phaseSystem&
Foam::diameterModels::populationBalanceModel::fluid() const
{
    return fluid_;
}


inline const Foam::fvMesh&
Foam::diameterModels::populationBalanceModel::mesh() const
{
    return mesh_;
}


inline const Foam::dictionary&
Foam::diameterModels::populationBalanceModel::dict() const
{
    return dict_;
}


inline const Foam::phaseModel&
Foam::diameterModels::populationBalanceModel::continuousPhase() const
{
    return continuousPhase_;
}


inline const Foam::UPtrList<Foam::diameterModels::velocityGroup>&
Foam::diameterModels::populationBalanceModel::velocityGroups() const
{
    return velocityGroups_;
}


inline const Foam::UPtrList<Foam::diameterModels::sizeGroup>&
Foam::diameterModels::populationBalanceModel::sizeGroups() const
{
    return sizeGroups_;
}


inline const Foam::volScalarField&
Foam::diameterModels::populationBalanceModel::SuSp(const label i) const
{
    return SuSp_[i];
}


inline const Foam::diameterModels::populationBalanceModel::phasePairTable&
Foam::diameterModels::populationBalanceModel::phasePairs() const
{
    return phasePairs_;
}


inline bool Foam::diameterModels::populationBalanceModel::isVelocityGroupPair
(
    const phasePair& pair
) const
{
    if
    (
        isA<velocityGroup>(pair.phase1().dPtr()())
        &&
        isA<velocityGroup>(pair.phase2().dPtr()())
    )
    {
        const velocityGroup& velGroup1 =
            refCast<const velocityGroup>(pair.phase1().dPtr()());

        const velocityGroup& velGroup2 =
            refCast<const velocityGroup>(pair.phase2().dPtr()());

        return velGroup1.popBalName() == velGroup2.popBalName();
    }

    return false;
}


inline const Foam::PtrList<Foam::dimensionedScalar>&
Foam::diameterModels::populationBalanceModel::v() const
{
    return v_;
}


inline const Foam::volScalarField&
Foam::diameterModels::populationBalanceModel::alphas() const
{
    if (velocityGroups_.size() > 1)
    {
        return alphas_();
    }
    else
    {
        return velocityGroups_.first().phase();
    }
}


inline const Foam::volVectorField&
Foam::diameterModels::populationBalanceModel::U() const
{
    if (velocityGroups_.size() > 1)
    {
        return U_();
    }
    else
    {
        return velocityGroups_.first().phase().U();
    }
}
```

define private and public member functions to return variables

#### populationBalanceModel.C

##### include

```cpp
#include "populationBalanceModel.H"
#include "coalescenceModel.H"
#include "breakupModel.H"
#include "binaryBreakupModel.H"
#include "driftModel.H"
#include "nucleationModel.H"
#include "phaseSystem.H"
#include "surfaceTensionModel.H"
#include "fvmDdt.H"
#include "fvcDdt.H"
#include "fvmSup.H"
#include "fvcSup.H"
#include "fvcDiv.H"
#include "phaseCompressibleMomentumTransportModel.H"
#include "shapeModel.H"
```

##### Static Data Members

```cpp
namespace Foam
{
namespace diameterModels
{
    defineTypeNameAndDebug(populationBalanceModel, 0);
}
}
```

##### Private Member Functions

###### registerVelocityGroups()

```cpp
void Foam::diameterModels::populationBalanceModel::registerVelocityGroups()
{
    forAll(fluid_.phases(), phasei)
    {
        if (isA<velocityGroup>(fluid_.phases()[phasei].dPtr()()))
        {
            const velocityGroup& velGroup =
                refCast<const velocityGroup>(fluid_.phases()[phasei].dPtr()());

            if (velGroup.popBalName() == this->name())
            {
                velocityGroups_.resize(velocityGroups_.size() + 1);

                velocityGroups_.set
                (
                    velocityGroups_.size() - 1,
                    &const_cast<velocityGroup&>(velGroup)
                );

                forAll(velGroup.sizeGroups(), i)
                {
                    this->registerSizeGroups
                    (
                        const_cast<sizeGroup&>(velGroup.sizeGroups()[i])
                    );
                }
            }
        }
    }
}
```

register velocity groups for every phase

###### registerSizeGroups()

```cpp

```

###### createPhasePairs()

```cpp

```

###### correct()

```cpp

```

###### birthByCoalescence()

```cpp

```

###### deathByCoalescence()

```cpp

```

###### 

```cpp

```

###### 

```cpp

```

###### 

```cpp

```

###### 

```cpp

```

###### 

```cpp

```

###### 

```cpp

```

###### 

```cpp

```

###### 

```cpp

```

###### 

```cpp

```

###### 

```cpp

```

###### 

```cpp

```

###### 

```cpp

```

##### public member functions

###### clone()

```cpp
Foam::autoPtr<Foam::diameterModels::populationBalanceModel>
Foam::diameterModels::populationBalanceModel::clone() const
{
    notImplemented("populationBalance::clone() const");
    return autoPtr<populationBalanceModel>(nullptr);
}
```

return a null pointer

###### writeData()

```cpp
bool Foam::diameterModels::populationBalanceModel::writeData(Ostream& os) const
{
    return os.good();
}
```

###### eta()

```cpp
const Foam::dimensionedScalar
Foam::diameterModels::populationBalanceModel::eta
(
    const label i,
    const dimensionedScalar& v
) const
{
    const dimensionedScalar& x0 = sizeGroups_[0].x();
    const dimensionedScalar& xi = sizeGroups_[i].x();
    const dimensionedScalar& xm = sizeGroups_.last().x();
    dimensionedScalar lowerBoundary(x0);
    dimensionedScalar upperBoundary(xm);

    if (i != 0) lowerBoundary = sizeGroups_[i-1].x();

    if (i != sizeGroups_.size() - 1) upperBoundary = sizeGroups_[i+1].x();

    if ((i == 0 && v < x0) || (i == sizeGroups_.size() - 1 && v > xm))
    {
        return v/xi;
    }
    else if (v < lowerBoundary || v > upperBoundary)
    {
        return 0;
    }
    else if (v.value() == xi.value())
    {
        return 1.0;
    }
    else if (v > xi)
    {
        return (upperBoundary - v)/(upperBoundary - xi);
    }
    else
    {
        return (v - lowerBoundary)/(xi - lowerBoundary);
    }
}
```

Return allocation coefficient



`sizeGroup.x()` can be found in `applications\solvers\multiphase\multiphaseEulerFoam\phaseSystems\diameterModels\velocityGroup\sizeGroup\sizeGroup.H`, is the representative volume of the sizeGroup

###### 

```cpp

```

###### 

```cpp

```





