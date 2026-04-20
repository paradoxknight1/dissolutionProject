#ifndef DissolutionCloud_H
#define DissolutionCloud_H

#include "DissolutionParcel.H"
#include "ParcelCloud.H"
#include "ParcelCloudBase.H"
#include "MomentumCloud.H"
#include "volFields.H"

namespace Foam
{

typedef ParcelCloud
<
    MomentumCloud
    <
        ParcelCloudBase
        <
            DissolutionParcel
        >
    >
> DissolutionCloudBase;


class DissolutionCloud
:
    public DissolutionCloudBase
{
    volScalarField& S_C_;

public:

    DissolutionCloud
    (
        const word& cloudName,
        const volScalarField& rho,
        const volVectorField& U,
        const volScalarField& mu,
        const dimensionedVector& g,
        volScalarField& S_C
    );

    DissolutionCloud(DissolutionCloud& c, const word& name);
    DissolutionCloud
    (
        const fvMesh& mesh,
        const word& name,
        const DissolutionCloud& c
    );

    virtual autoPtr<Cloud<parcelType>> clone(const word& name)
    {
        return autoPtr<Cloud<parcelType>>(new DissolutionCloud(*this, name));
    }

    virtual autoPtr<Cloud<parcelType>> cloneBare(const word& name) const
    {
        return autoPtr<Cloud<parcelType>>
        (
            new DissolutionCloud(this->mesh(), name, *this)
        );
    }

    void evolve();

    scalar totalParticleMass() const;
};

} // End namespace Foam

#endif
