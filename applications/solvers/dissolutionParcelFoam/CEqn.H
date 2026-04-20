// The DPM does not resolve the local surface flux J(x_s). Instead it uses the
// reduced-order chain: local diffusion -> area-averaged flux -> k_m -> Sh
// closure -> parcel dm/dt -> Eulerian source S_C. Particle mass loss is added
// to C with multiplicity scaling through nParticle.

fvScalarMatrix CEqn
(
    fvm::ddt(C)
  + fvm::div(phi, C)
  - fvm::laplacian(molecularDiffusivity, C)
 ==
    S_C + fvModels.source(C)
);

CEqn.relax();
fvConstraints.constrain(CEqn);
CEqn.solve();
fvConstraints.constrain(C);
