Info<< "Calculating polymeriation field\n" << endl;
//i - number of a layer

if (i == 1)
{
    volScalarField exposure = laser.E()*BeerLambert(Z - laser.height(i) - small, Dp);
    totalSqrtExposure += sqrt(exposure);
    polymerization = neg0(Z - laser.height(i) - small)*pMax*
                (1 - pow(1 - pGel/pMax, totalSqrtExposure/sqrt(Ec)))
                + pos(Z - laser.height(i) - small)*pGel/2;
    Info<< "Calculating Young's modulus\n" << endl;

    E = pos0(polymerization - pGel)*Emax*((1 - e0)*(polymerization - pGel)/(pMax - pGel)
      + e0) + neg(polymerization - pGel)*Emax*e0;
    E.write();

    if (gMin(E) < SMALL)
        {
            FatalError
            << "Somewhere the Young's modulus is equal to zero."
            << exit(FatalError);
        }
}

Info<< "Calculating Lame's coefficients\n" << endl;
volScalarField mu(E/(2.0*(1.0 + nu)));
volScalarField lambda(nu*E/((1.0 + nu)*(1.0 - 2.0*nu)));
volScalarField threeK(E/(1.0 - 2.0*nu));

volTensorField gradD(fvc::grad(D));

if (i > 1)
{
    epsilonResidual -= epsilon*pos(Z - laser.height(i - 1));
}

volTensorField testEpsilon("testEpsilon", gradD - epsilonResidual);
testEpsilon.write();
//epsilonResidual.write();

sigma = mu*twoSymm(gradD) + lambda*I*tr(gradD - epsilonResidual) - 2*mu*epsilonResidual
      - threeK*epsilonChemicalMax/pMax*polymerization*I;
sigma.write();

//sigma = neg0(Z - laser.height(i))*(mu*twoSymm(gradD) + lambda*I*tr(gradD)
//      - threeK*epsilonChemicalMax/pMax*polymerization*I);

Info<< "Calculating explicit part of div(sigma) divSigmaExp\n" << endl;
volVectorField divSigmaExp
(
    "divSigmaExp",
    fvc::div(sigma) - fvc::laplacian(2*mu + lambda, D, "laplacian(DD,D)")
);
