#include "DissolutionCloud.H"

#include "PstreamReduceOps.H"

namespace Foam
{

DissolutionCloud::DissolutionCloud
(
    const word& cloudName,
    const volScalarField& rho,
    const volVectorField& U,
    const volScalarField& mu,
    const dimensionedVector& g,
    volScalarField& S_C
)
:
    DissolutionCloudBase(cloudName, rho, U, mu, g, true),
    S_C_(S_C)
{}


DissolutionCloud::DissolutionCloud
(
    DissolutionCloud& c,
    const word& name
)
:
    DissolutionCloudBase(c, name),
    S_C_(c.S_C_)
{}


DissolutionCloud::DissolutionCloud
(
    const fvMesh& mesh,
    const word& name,
    const DissolutionCloud& c
)
:
    DissolutionCloudBase(mesh, name, c),
    S_C_(const_cast<volScalarField&>(c.S_C_))
{}


void DissolutionCloud::evolve()
{
    S_C_ = dimensionedScalar("zero", S_C_.dimensions(), 0);

    if (this->solution().canEvolve())
    {
        parcelType::trackingData td(*this);
        this->solve(*this, td);
    }

    S_C_.correctBoundaryConditions();
}


scalar DissolutionCloud::totalParticleMass() const
{
    scalar total = 0;

    forAllConstIter(DissolutionCloud, *this, iter)
    {
        total += iter().nParticle()*iter().mass();
    }

    return returnReduce(total, sumOp<scalar>());
}

} // End namespace Foam
