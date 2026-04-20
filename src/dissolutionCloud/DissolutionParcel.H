#ifndef DissolutionParcel_H
#define DissolutionParcel_H

#include "MomentumParcel.H"
#include "autoPtr.H"
#include "IOdictionary.H"
#include "mathematicalConstants.H"

#include <cmath>

namespace Foam
{

class DissolutionParcel
:
    public MomentumParcel<particle>
{
    typedef MomentumParcel<particle> parcelType;

public:

    TypeName("DissolutionParcel");


    class trackingData
    :
        public parcelType::trackingData
    {
        autoPtr<interpolation<scalar>> CInterp_;
        scalar Cc_;

    public:

        template<class TrackCloudType>
        inline trackingData(const TrackCloudType& cloud);

        inline const interpolation<scalar>& CInterp() const;

        inline scalar Cc() const;
        inline scalar& Cc();
    };


    DissolutionParcel
    (
        const polyMesh& mesh,
        const barycentric& coordinates,
        const label celli,
        const label tetFacei,
        const label tetPti
    );

    DissolutionParcel
    (
        const polyMesh& mesh,
        const vector& position,
        const label celli
    );

    DissolutionParcel
    (
        const polyMesh& mesh,
        Istream& is,
        bool readFields = true
    );

    DissolutionParcel(const DissolutionParcel& p);
    DissolutionParcel(const DissolutionParcel& p, const polyMesh& mesh);

    virtual autoPtr<particle> clone() const
    {
        return autoPtr<particle>(new DissolutionParcel(*this));
    }

    virtual autoPtr<particle> clone(const polyMesh& mesh) const
    {
        return autoPtr<particle>(new DissolutionParcel(*this, mesh));
    }

    class iNew
    {
        const polyMesh& mesh_;

    public:

        explicit iNew(const polyMesh& mesh)
        :
            mesh_(mesh)
        {}

        autoPtr<DissolutionParcel> operator()(Istream& is) const
        {
            return autoPtr<DissolutionParcel>
            (
                new DissolutionParcel(mesh_, is, true)
            );
        }
    };


    template<class TrackCloudType>
    void setCellValues(TrackCloudType& cloud, trackingData& td);

    template<class TrackCloudType>
    void calc(TrackCloudType& cloud, trackingData& td, const scalar dt);


    friend Ostream& operator<<(Ostream& os, const DissolutionParcel& p);
};


template<class TrackCloudType>
inline DissolutionParcel::trackingData::trackingData
(
    const TrackCloudType& cloud
)
:
    parcelType::trackingData(cloud),
    CInterp_
    (
        interpolation<scalar>::New
        (
            cloud.solution().interpolationSchemes(),
            cloud.mesh().template lookupObject<volScalarField>("C")
        )
    ),
    Cc_(0)
{}


inline const interpolation<scalar>& DissolutionParcel::trackingData::CInterp()
const
{
    return CInterp_();
}


inline scalar DissolutionParcel::trackingData::Cc() const
{
    return Cc_;
}


inline scalar& DissolutionParcel::trackingData::Cc()
{
    return Cc_;
}


template<class TrackCloudType>
void DissolutionParcel::setCellValues(TrackCloudType& cloud, trackingData& td)
{
    parcelType::setCellValues(cloud, td);

    const tetIndices tetIs = this->currentTetIndices();

    td.Cc() = td.CInterp().interpolate(this->coordinates(), tetIs);

    if (!std::isfinite(td.Cc()))
    {
        td.Cc() = cloud.mesh().template lookupObject<volScalarField>("C")
            [this->cell()];
    }

    td.Cc() = max(td.Cc(), scalar(0));
}


template<class TrackCloudType>
void DissolutionParcel::calc
(
    TrackCloudType& cloud,
    trackingData& td,
    const scalar dt
)
{
    const scalar np0 = this->nParticle();
    const scalar mass0 = this->mass();
    const scalar d0 = this->d();

    if
    (
        dt <= rootVSmall
     || np0 <= 0
     || mass0 <= cloud.constProps().minParcelMass()
     || d0 <= rootVSmall
    )
    {
        return;
    }

    const IOdictionary& transportProperties =
        cloud.mesh().template lookupObject<IOdictionary>("transportProperties");

    const dimensionedScalar D("D", dimViscosity, transportProperties);
    const dimensionedScalar Cs("Cs", dimDensity, transportProperties);

    const scalar DValue = max(D.value(), rootVSmall);
    const scalar nu = td.muc()/max(td.rhoc(), rootVSmall);
    const scalar Re = max(this->Re(td), scalar(0));
    const scalar Sc = max(nu/DValue, scalar(0));
    const scalar Sh = 2.0 + 0.6*sqrt(Re)*cbrt(Sc);
    const scalar kM = max(Sh*DValue/max(d0, rootVSmall), scalar(0));
    const scalar dC = max(Cs.value() - td.Cc(), scalar(0));
    const scalar dmDtSingle =
        -min(kM*this->areaS(d0)*dC, mass0/max(dt, rootVSmall));

    vector dUTrans = Zero;
    scalar Spu = 0.0;
    this->U() =
        this->calcVelocity(cloud, td, dt, Re, td.muc(), mass0, Zero, dUTrans, Spu);

    if (cloud.solution().coupled())
    {
        cloud.UTransRef()[this->cell()] += np0*dUTrans;
        cloud.UCoeffRef()[this->cell()] += np0*Spu;
    }

    const scalar newMass = max(mass0 + dmDtSingle*dt, scalar(0));
    const scalar lostMassParcel = np0*(mass0 - newMass);

    if (lostMassParcel > rootVSmall)
    {
        volScalarField& S_C =
            const_cast<volScalarField&>
            (
                cloud.mesh().template lookupObject<volScalarField>("S_C")
            );

        S_C[this->cell()] += lostMassParcel/(cloud.mesh().V()[this->cell()]*dt);
    }

    if (newMass <= cloud.constProps().minParcelMass())
    {
        td.keepParticle = false;
        return;
    }

    const scalar rho = max(this->rho(), rootVSmall);
    const scalar newD =
        cbrt(6.0*newMass/(constant::mathematical::pi*rho));

    if (!std::isfinite(newD) || newD <= rootVSmall)
    {
        td.keepParticle = false;
        return;
    }

    this->d() = newD;
    this->dTarget() = newD;
}

} // End namespace Foam

#endif
