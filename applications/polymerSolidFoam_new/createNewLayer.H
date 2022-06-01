Info<< "Calculating polymerization field\n" << endl;
//i - number of a layer
if (i == 1){
    exposure = laser.E()*BeerLambert(Z - laser.height(i) - small, Dp);
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
    Info<< "Calculating Lame's coefficients\n" << endl;
    mu = E/(2.0*(1.0 + nu));
    lambda = nu*E/((1.0 + nu)*(1.0 - 2.0*nu));
    threeK = E/(1.0 - 2.0*nu);
    muf = fvc::interpolate(mu);
    lambdaf = fvc::interpolate(lambda);
}

gradD = fvc::grad(D);

sigma = mu*twoSymm(gradD) + lambda*I*tr(gradD)
      - threeK*epsilonChemicalMax/pMax*polymerization*I;
sigma.write();

Info<< "Calculating explicit part of div(sigma) divSigmaExp\n" << endl;
volVectorField divSigmaExp
(
    "divSigmaExp",
    fvc::div(sigma) - fvc::laplacian(2*mu + lambda, D, "laplacian(DD,D)")
);

/*if (i > 1)
{
    divSigmaExp -=
        fvc::div((lambdaf + 2*muf)*epsilonf*mesh.magSf()*pos0(Zf - laser.height(i - 1)));
}*/

if (i > 1)
{
    divSigmaExp +=
        fvc::laplacian(pos0(Zf - laser.height(i - 1))*(lambdaf+2*muf),
        D_old, "laplacian(DD,D)");
}
