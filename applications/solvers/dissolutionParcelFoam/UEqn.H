// Solve the momentum equation. The MRF hooks are kept to preserve
// compatibility with a later rotating-frame extension.

MRF.correctBoundaryVelocity(U);

fvVectorMatrix UEqn
(
    fvm::ddt(U)
  + fvm::div(phi, U)
  + MRF.DDt(U)
  + turbulence->divDevSigma(U)
 ==
    g + fvModels.source(U)
);

UEqn.relax();
fvConstraints.constrain(UEqn);

if (piso.momentumPredictor())
{
    solve(UEqn == -fvc::grad(p));
    fvConstraints.constrain(U);
}
