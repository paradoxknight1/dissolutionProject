#include "DissolutionCloud.H"

#include "makeParcelCloudFunctionObjects.H"
#include "makeParcelForces.H"
#include "makeParcelDispersionModels.H"
#include "makeParcelInjectionModels.H"
#include "makeParcelPatchInteractionModels.H"
#include "makeParcelStochasticCollisionModels.H"
#include "makeParcelSurfaceFilmModels.H"

namespace Foam
{
    typedef DissolutionCloud dissolutionCloud;
}

makeParcelCloudFunctionObjects(dissolutionCloud);
makeParcelForces(dissolutionCloud);
makeParcelDispersionModels(dissolutionCloud);
makeParcelInjectionModels(dissolutionCloud);
makeParcelPatchInteractionModels(dissolutionCloud);
makeParcelStochasticCollisionModels(dissolutionCloud);
makeParcelSurfaceFilmModels(dissolutionCloud);
