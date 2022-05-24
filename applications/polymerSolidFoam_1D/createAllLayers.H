Info<< "Calculating polymeriation field\n" << endl;
//i - number of a layer
 for (label i = 1; i <= numberOfLayers.value(); i++)
    {
        volScalarField exposure = laser.E()*BeerLambert(Z - laser.height(i) - small, Dp);
        totalSqrtExposure += sqrt(exposure);
        polymerization = neg0(Z - laser.height(i) - small)*pMax*
                (1 - pow(1 - pGel/pMax, totalSqrtExposure/sqrt(Ec)))
                + pos(Z - laser.height(i) - small)*pGel/2;
    }
Info<< "Calculating Young's modulus\n" << endl;

E = pos0(polymerization - pGel)*Emax*((1 - e0)*(polymerization - pGel)/(pMax - pGel)
  + e0) + neg(polymerization - pGel)*Emax*e0;
//Info<< "hi3" << endl;
E.write();

if (gMin(E) < SMALL)
        {
            FatalError
            << "Somewhere the Young's modulus is equal to zero."
            << exit(FatalError);
        }

Info<< "Calculating Lame's coefficients\n" << endl;
volScalarField mu(E/(2.0*(1.0 + nu)));
volScalarField lambda(nu*E/((1.0 + nu)*(1.0 - 2.0*nu)));
volScalarField threeK(E/(1.0 - 2.0*nu));
volTensorField gradD(fvc::grad(D));

sigma = mu*twoSymm(gradD) + lambda*I*tr(gradD) + threeK*epsilonChemicalMax*polymerization*I;

Info<< "Calculating explicit part of div(sigma) divSigmaExp\n" << endl;
volVectorField divSigmaExp
(
    "divSigmaExp",
    fvc::div(sigma) - fvc::laplacian(2*mu + lambda, D, "laplacian(DD,D)")
);

while (simple.loop())
    {
        Info<< "Iteration: " << runTime.value() << nl << endl;
        fvVectorMatrix DEqn
        (
            fvm::laplacian(2*mu + lambda, D, "laplacian(DD,D)")
            + divSigmaExp
        );
        DEqn.solve();
        gradD = fvc::grad(D);
        epsilon = symm(gradD);
        /*sigma = mu*twoSymm(gradD) + lambda*I*tr(gradD)
              - threeK*(epsilonChemicalMax/pMax*polymerization*I - epsilon);*/
        sigma = mu*twoSymm(gradD) + lambda*I*tr(gradD)
              - threeK*(epsilonChemicalMax/pMax*polymerization*I);
        divSigmaExp = fvc::div(sigma - (2*mu + lambda)*gradD, "div(sigma)");

        int mas[6] = {500, 501, 502, 503, 504, 505};
        int size = sizeof(mas)/sizeof(mas[0]);

        for (int i = 0; i < size; i++){
            if(runTime.value() == mas[i]){
                runTime.writeNow();
            }
        }

        if (runTime.writeTime())
        {
                epsilon = symm(gradD);
                sigmaEq = sqrt(1.5*magSqr(dev(sigma)));
                Info<< "Max sigmaEq = " << max(sigmaEq).value() << endl;
                runTime.write();
            }

        runTime.printExecutionTime(Info);
    }
