#include "DissolutionParcel.H"

namespace Foam
{

defineTypeNameAndDebug(DissolutionParcel, 0);


DissolutionParcel::DissolutionParcel
(
    const polyMesh& mesh,
    const barycentric& coordinates,
    const label celli,
    const label tetFacei,
    const label tetPti
)
:
    parcelType(mesh, coordinates, celli, tetFacei, tetPti)
{}


DissolutionParcel::DissolutionParcel
(
    const polyMesh& mesh,
    const vector& position,
    const label celli
)
:
    parcelType(mesh, position, celli)
{}


DissolutionParcel::DissolutionParcel
(
    const polyMesh& mesh,
    Istream& is,
    bool readFields
)
:
    parcelType(mesh, is, readFields)
{}


DissolutionParcel::DissolutionParcel(const DissolutionParcel& p)
:
    parcelType(p)
{}


DissolutionParcel::DissolutionParcel
(
    const DissolutionParcel& p,
    const polyMesh& mesh
)
:
    parcelType(p, mesh)
{}


Ostream& operator<<(Ostream& os, const DissolutionParcel& p)
{
    os << static_cast<const DissolutionParcel::parcelType&>(p);
    os.check("Ostream& operator<<(Ostream&, const DissolutionParcel&)");
    return os;
}

} // End namespace Foam
