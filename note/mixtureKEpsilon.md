# mixtureKEpsilon

## mixtureKEpsilon.H

### Description

The model is based on Behzadi et al.[^1], but but an effective density for the gas is used in the averaging and an alternative model for bubble-generated turbulence from Lahey et al.[^2]

[^1]: Behzadi, A., Issa, R. I., & Rusche, H. (2004). Modelling of dispersed bubble and droplet flow at high phase fractions. Chemical Engineering Science, 59(4), 759-770.
[^2]: Lahey Jr, R. T. (2005). The simulation of multidimensional multiphase flows. Nuclear Engineering and Design, 235(10), 1043-1060.

The default coefficients are:

| Variable | Value | Remarks |
| --|-- | -- |
| $C_\mu$ | 0.09 | |
| $C_1$ | 1.44 | |
| $C_2$ | 1.92 | |
| $C_3$ | $C_2$ | |
| $C_p$ | 0.25 | Bubble generated turbulence |
| $\alpha_p$ | 0.3 | Gas phase fraction below which bubble generated turbulence is included |
| $\sigma_k$ | 1.0 | |
| $\sigma_\epsilon$ | 1.3 |

### Private

```cpp
    // Private Data

        mutable mixtureKEpsilon<BasicMomentumTransportModel>
            *liquidTurbulencePtr_;


    // Private Member Functions

        //- Return the turbulence model for the other phase
        mixtureKEpsilon<BasicMomentumTransportModel>& liquidTurbulence() const;
```

### Protected

```cpp
        // Model coefficients

            dimensionedScalar Cmu_;
            dimensionedScalar C1_;
            dimensionedScalar C2_;
            dimensionedScalar C3_;
            dimensionedScalar Cp_;
            dimensionedScalar alphap_;
            dimensionedScalar sigmak_;
            dimensionedScalar sigmaEps_;

        // Fields

            volScalarField k_;
            volScalarField epsilon_;

        // Mixture fields

            autoPtr<volScalarField> Ct2_;
            autoPtr<volScalarField> rhom_;
            autoPtr<volScalarField> km_;
            autoPtr<volScalarField> epsilonm_;
```

define model coefficients, fields and mixture fields

### protected member function

```cpp
    // Protected Member Functions

        wordList epsilonBoundaryTypes(const volScalarField& epsilon) const;

        void correctInletOutlet
        (
            volScalarField& vsf,
            const volScalarField& refVsf
        ) const;

        void initMixtureFields();

        virtual void correctNut();

        tmp<volScalarField> Ct2() const;

        tmp<volScalarField> rholEff() const;
        tmp<volScalarField> rhogEff() const;
        tmp<volScalarField> rhom() const;

        tmp<volScalarField> mix
        (
            const volScalarField& fc,
            const volScalarField& fd
        ) const;

        tmp<volScalarField> mixU
        (
            const volScalarField& fc,
            const volScalarField& fd
        ) const;

        tmp<surfaceScalarField> mixFlux
        (
            const surfaceScalarField& fc,
            const surfaceScalarField& fd
        ) const;

        tmp<volScalarField> bubbleG() const;
        virtual tmp<fvScalarMatrix> kSource() const;
        virtual tmp<fvScalarMatrix> epsilonSource() const;

        //- Return the effective diffusivity for k
        tmp<volScalarField> DkEff(const volScalarField& nutm) const
        {
            return volScalarField::New
            (
                "DkEff",
                nutm/sigmak_
            );
        }

        //- Return the effective diffusivity for epsilon
        tmp<volScalarField> DepsilonEff(const volScalarField& nutm) const
        {
            return volScalarField::New
            (
                "DepsilonEff",
                nutm/sigmaEps_
            );
        }
```

$$
D_{k, Eff} = \frac{\nu_{t, m}}{\sigma_k}
$$

$$
D_{\epsilon, Eff} = \frac{\nu_{t, m}}{\sigma_\epsilon}
$$

### public

```cpp
    typedef typename BasicMomentumTransportModel::alphaField alphaField;
    typedef typename BasicMomentumTransportModel::rhoField rhoField;
    typedef typename BasicMomentumTransportModel::transportModel transportModel;


    //- Runtime type information
    TypeName("mixtureKEpsilon");

    // Constructors

        //- Construct from components
        mixtureKEpsilon
        (
            const alphaField& alpha,
            const rhoField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const transportModel& transport,
            const word& type = typeName
        );

        //- Disallow default bitwise copy construction
        mixtureKEpsilon(const mixtureKEpsilon&) = delete;


    //- Destructor
    virtual ~mixtureKEpsilon()
    {}
```

define type and constructors and destructors

### public member functions

```cpp
    // Member Functions

        //- Re-read model coefficients if they have changed
        virtual bool read();

        //- Return the turbulence kinetic energy
        virtual tmp<volScalarField> k() const
        {
            return k_;
        }

        //- Return the turbulence kinetic energy dissipation rate
        virtual tmp<volScalarField> epsilon() const
        {
            return epsilon_;
        }

        //- Solve the turbulence equations and correct the turbulence viscosity
        virtual void correct();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const mixtureKEpsilon&) = delete;
```

read, return $k$ and $\epsilon$, and correct

## mixtureKEpsilon.C

### constructors

```cpp
...
```

initialize turbulence model with default coefficients or read values, define fields of $k$ and $\epsilon$, bound $k$ and $\epsilon$ and print coefficients

### epsilonBoundaryTypes()

```cpp
template<class BasicMomentumTransportModel>
wordList mixtureKEpsilon<BasicMomentumTransportModel>::epsilonBoundaryTypes
(
    const volScalarField& epsilon
) const
{
    const volScalarField::Boundary& ebf = epsilon.boundaryField();

    wordList ebt = ebf.types();

    forAll(ebf, patchi)
    {
        if (isA<fixedValueFvPatchScalarField>(ebf[patchi]))
        {
            ebt[patchi] = fixedValueFvPatchScalarField::typeName;
        }
    }

    return ebt;
}
```

define boundary type of $\epsilon$

### correctInletOutlet()

```cpp
template<class BasicMomentumTransportModel>
void mixtureKEpsilon<BasicMomentumTransportModel>::correctInletOutlet
(
    volScalarField& vsf,
    const volScalarField& refVsf
) const
{
    volScalarField::Boundary& bf = vsf.boundaryFieldRef();
    const volScalarField::Boundary& refBf =
        refVsf.boundaryField();

    forAll(bf, patchi)
    {
        if
        (
            isA<inletOutletFvPatchScalarField>(bf[patchi])
         && isA<inletOutletFvPatchScalarField>(refBf[patchi])
        )
        {
            refCast<inletOutletFvPatchScalarField>
            (bf[patchi]).refValue() =
            refCast<const inletOutletFvPatchScalarField>
            (refBf[patchi]).refValue();
        }
    }
}
```

set boundary of vsf as refVsf, it is used to set boundary for mixture fields

### initMixtureFields()

```cpp
template<class BasicMomentumTransportModel>
void mixtureKEpsilon<BasicMomentumTransportModel>::initMixtureFields()
{
    if (rhom_.valid()) return;

    // Local references to gas-phase properties
    const volScalarField& kg = this->k_;
    const volScalarField& epsilong = this->epsilon_;

    // Local references to liquid-phase properties
    mixtureKEpsilon<BasicMomentumTransportModel>& turbc =
        this->liquidTurbulence();
    const volScalarField& kl = turbc.k_;
    const volScalarField& epsilonl = turbc.epsilon_;

    word startTimeName
    (
        this->runTime_.timeName(this->runTime_.startTime().value())
    );

    Ct2_.set
    (
        new volScalarField
        (
            IOobject
            (
                "Ct2",
                startTimeName,
                this->mesh_,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            Ct2()
        )
    );

    rhom_.set
    (
        new volScalarField
        (
            IOobject
            (
                "rhom",
                startTimeName,
                this->mesh_,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            rhom()
        )
    );

    km_.set
    (
        new volScalarField
        (
            IOobject
            (
                "km",
                startTimeName,
                this->mesh_,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            mix(kl, kg),
            kl.boundaryField().types()
        )
    );
    correctInletOutlet(km_(), kl);

    epsilonm_.set
    (
        new volScalarField
        (
            IOobject
            (
                "epsilonm",
                startTimeName,
                this->mesh_,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            mix(epsilonl, epsilong),
            epsilonBoundaryTypes(epsilonl)
        )
    );
    correctInletOutlet(epsilonm_(), epsilonl);
}
```

* get references of $k$ and $\epsilon$ for gas phase as $k_g$ and $\epsilon_g$
* get reference of turbulence model, $k$ and $\epsilon$ for the other phase (liquid) as:
  * turbc
  * $k_l$
  * $\epsilon_l$
* define and set fields of:
  * $Ct_2$
  * $\rho_m$
  * $k_m$
  * $\epsilon_m$
* correct inlet boundary of $\epsilon_m$ to $\epsilon_l$

### read()

```cpp
template<class BasicMomentumTransportModel>
bool mixtureKEpsilon<BasicMomentumTransportModel>::read()
{
    if (eddyViscosity<RASModel<BasicMomentumTransportModel>>::read())
    {
        Cmu_.readIfPresent(this->coeffDict());
        C1_.readIfPresent(this->coeffDict());
        C2_.readIfPresent(this->coeffDict());
        C3_.readIfPresent(this->coeffDict());
        Cp_.readIfPresent(this->coeffDict());
        sigmak_.readIfPresent(this->coeffDict());
        sigmaEps_.readIfPresent(this->coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}
```

read model coefficients

### correctNut()

```cpp
template<class BasicMomentumTransportModel>
void mixtureKEpsilon<BasicMomentumTransportModel>::correctNut()
{
    this->nut_ = Cmu_*sqr(k_)/epsilon_;
    this->nut_.correctBoundaryConditions();
    fv::options::New(this->mesh_).correct(this->nut_);
}
```

$$
\nu_t = C_\mu \frac{k^2}{\epsilon}
$$

correct boundary conditions for $\nu_t$

### liquidTurbulence()

```cpp
template<class BasicMomentumTransportModel>
mixtureKEpsilon<BasicMomentumTransportModel>&
mixtureKEpsilon<BasicMomentumTransportModel>::liquidTurbulence() const
{
    if (!liquidTurbulencePtr_)
    {
        const volVectorField& U = this->U_;

        const transportModel& gas = this->transport();
        const phaseSystem& fluid = gas.fluid();
        const transportModel& liquid = fluid.otherPhase(gas);

        liquidTurbulencePtr_ =
           &const_cast<mixtureKEpsilon<BasicMomentumTransportModel>&>
            (
                U.db().lookupObject
                <
                    mixtureKEpsilon<BasicMomentumTransportModel>
                >
                (
                    IOobject::groupName
                    (
                        momentumTransportModel::typeName,
                        liquid.name()
                    )
                )
            );
    }

    return *liquidTurbulencePtr_;
}
```

get turbulence model of the other pahse

### Ct2()

```cpp
template<class BasicMomentumTransportModel>
tmp<volScalarField> mixtureKEpsilon<BasicMomentumTransportModel>::Ct2() const
{
    const mixtureKEpsilon<BasicMomentumTransportModel>& liquidTurbulence =
        this->liquidTurbulence();

    const transportModel& gas = this->transport();
    const phaseSystem& fluid = gas.fluid();
    const transportModel& liquid = fluid.otherPhase(gas);

    const dragModel& drag = fluid.lookupSubModel<dragModel>(gas, liquid);

    const volScalarField& alphag = this->alpha_;

    volScalarField magUr(mag(liquidTurbulence.U() - this->U()));

    volScalarField beta
    (
        (6*this->Cmu_/(4*sqrt(3.0/2.0)))
       *drag.K()/liquid.rho()
       *(liquidTurbulence.k_/liquidTurbulence.epsilon_)
    );
    volScalarField Ct0((3 + beta)/(1 + beta + 2*gas.rho()/liquid.rho()));
    volScalarField fAlphad((180 + (-4.71e3 + 4.26e4*alphag)*alphag)*alphag);

    return sqr(1 + (Ct0 - 1)*exp(-fAlphad));
}
```

* get liquid phase (the other phase)
* get drag model
* get volume fraction $\alpha_g$
* get relative velocity magnitude $\|\mathbf{U}_l - \mathbf{U}_g\|$
* get $\beta$ as
  * $$\beta = \frac{6 C_\mu}{4 \sqrt{1.5}} \cdot \frac{K}{\rho_l} \cdot  \frac{k_l}{\epsilon_l}$$ 
  * where $K$ is drag coefficient K in momentum equation
* define Ct_0 as:
  * $$Ct_0 = \frac{3 + \beta}{1 + \beta + 2 \rho_g/\rho_l}$$
* define $f_{\alpha, d}$ as
  * $$f_{\alpha, d} = (180 + (-4710 + 42600\alpha_g) \alpha_g) \alpha_g = 42600 \alpha_g^3 - 4710 \alpha_g^2 + 180 \alpha_g$$
* return $Ct_2$ as
  * $$Ct_2 = (1 + (Ct_0 -1)\exp(-f_{\alpha, d}) )^2 $$

### rholEff()

```cpp
template<class BasicMomentumTransportModel>
tmp<volScalarField>
mixtureKEpsilon<BasicMomentumTransportModel>::rholEff() const
{
    const transportModel& gas = this->transport();
    const phaseSystem& fluid = gas.fluid();
    return fluid.otherPhase(gas).rho();
}
```

get $\rho_{l, Eff} = \rho_l$

### rhoEff()

```cpp
template<class BasicMomentumTransportModel>
tmp<volScalarField>
mixtureKEpsilon<BasicMomentumTransportModel>::rhogEff() const
{
    const transportModel& gas = this->transport();
    const phaseSystem& fluid = gas.fluid();
    const virtualMassModel& virtualMass =
        fluid.lookupSubModel<virtualMassModel>(gas, fluid.otherPhase(gas));
    return gas.rho() + virtualMass.Cvm()*fluid.otherPhase(gas).rho();
}
```

* get gas and liquid phase, and virtural mass model
* $$\rho_{g, Eff} = \rho_g + C_{vm} \rho_l$$

### rhom()

```cpp
template<class BasicMomentumTransportModel>
tmp<volScalarField> mixtureKEpsilon<BasicMomentumTransportModel>::rhom() const
{
    const volScalarField& alphag = this->alpha_;
    const volScalarField& alphal = this->liquidTurbulence().alpha_;

    return alphal*rholEff() + alphag*rhogEff();
}
```

$$
rho_m = \alpha_g \rho_{g, Eff} + \alpha_l \rho_{l, Eff}
$$

### mix()

```cpp
template<class BasicMomentumTransportModel>
tmp<volScalarField> mixtureKEpsilon<BasicMomentumTransportModel>::mix
(
    const volScalarField& fc,
    const volScalarField& fd
) const
{
    const volScalarField& alphag = this->alpha_;
    const volScalarField& alphal = this->liquidTurbulence().alpha_;

    return (alphal*rholEff()*fc + alphag*rhogEff()*fd)/rhom_();
}
```

$$
mix = \frac{\alpha_l \rho_{l, Eff} f_c + \alpha_g \rho_{g, Eff} f_d}{\rho_m}
$$

### mixU()

```cpp
template<class BasicMomentumTransportModel>
tmp<volScalarField> mixtureKEpsilon<BasicMomentumTransportModel>::mixU
(
    const volScalarField& fc,
    const volScalarField& fd
) const
{
    const volScalarField& alphag = this->alpha_;
    const volScalarField& alphal = this->liquidTurbulence().alpha_;

    return
        (alphal*rholEff()*fc + alphag*rhogEff()*Ct2_()*fd)
       /(alphal*rholEff() + alphag*rhogEff()*Ct2_());
}
```

$$
mixU = \frac{\alpha_l \rho_{l, Eff} f_c + \alpha_g \rho_{g, Eff} Ct_2 f_d}{\alpha_l \rho_{l, Eff} + \alpha_g \rho_{g, Eff} Ct_2}
$$

### mixFlux()

```cpp
template<class BasicMomentumTransportModel>
tmp<surfaceScalarField> mixtureKEpsilon<BasicMomentumTransportModel>::mixFlux
(
    const surfaceScalarField& fc,
    const surfaceScalarField& fd
) const
{
    const volScalarField& alphag = this->alpha_;
    const volScalarField& alphal = this->liquidTurbulence().alpha_;

    surfaceScalarField alphalf(fvc::interpolate(alphal));
    surfaceScalarField alphagf(fvc::interpolate(alphag));

    surfaceScalarField rholEfff(fvc::interpolate(rholEff()));
    surfaceScalarField rhogEfff(fvc::interpolate(rhogEff()));

    return
       (alphalf*rholEfff*fc + alphagf*rhogEfff*fvc::interpolate(Ct2_())*fd)
      /(alphalf*rholEfff + alphagf*rhogEfff*fvc::interpolate(Ct2_()));
}
```

get surface value by interpolating, then

$$
mixFlux = \frac{(\alpha_l)_f (\rho_{l, Eff})_f f_c + (\alpha_g)_f (\rho_{g, Eff})_f (Ct_2)_f f_d}{(\alpha_l)_f (\rho_{l, Eff})_f + (\alpha_g)_f (\rho_{g, Eff})_f (Ct_2)_f}
$$

### bubbleG()

```cpp
template<class BasicMomentumTransportModel>
tmp<volScalarField>
mixtureKEpsilon<BasicMomentumTransportModel>::bubbleG() const
{
    const mixtureKEpsilon<BasicMomentumTransportModel>& liquidTurbulence =
        this->liquidTurbulence();

    const transportModel& gas = this->transport();
    const phaseSystem& fluid = gas.fluid();
    const transportModel& liquid = fluid.otherPhase(gas);

    const dragModel& drag = fluid.lookupSubModel<dragModel>(gas, liquid);

    volScalarField magUr(mag(liquidTurbulence.U() - this->U()));

    // Lahey model
    tmp<volScalarField> bubbleG
    (
        Cp_
       *pos(alphap_ - gas)*liquid*liquid.rho()
       *(
            pow3(magUr)
          + pow(drag.CdRe()*liquid.thermo().nu()/gas.d(), 4.0/3.0)
           *pow(magUr, 5.0/3.0)
        )
       *gas
       /gas.d()
    );

    // Simple model
    // tmp<volScalarField> bubbleG
    // (
    //     Cp_*liquid*drag.K()*sqr(magUr)
    // );

    return bubbleG;
}
```

if $\alpha_p > \alpha_g$, then

$$
bubbleG = C_p \alpha_l \rho_l \left( \|\mathbf{U}_l - \mathbf{U}_g\|^3 + (\frac{C_d Re \nu_l}{d_g})^{4/3} \|\mathbf{U}_l - \mathbf{U}_g\|^{5/3} \right) \frac{\alpha_g}{d_g}
$$

else

$$
bubbleG = 0
$$
