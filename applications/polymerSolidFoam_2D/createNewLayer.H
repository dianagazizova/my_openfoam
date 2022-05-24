Info<< "Calculating polymerization field\n" << endl;
//i - number of a layer

if (i == 1)
{
    polymerization = neg0(Z - laser.height(i) - small)*pMax;
}

Info<< "Calculating Lame's coefficients\n" << endl;
volScalarField mu(E/(2.0*(1.0 + nu)));
volScalarField lambda(nu*E/((1.0 + nu)*(1.0 - 2.0*nu)));
volScalarField threeK(E/(1.0 - 2.0*nu));

volTensorField gradD(fvc::grad(D));

sigma = mu*twoSymm(gradD) + lambda*I*tr(gradD)
      - threeK*epsilonChemicalMax/pMax*polymerization*I;
sigma.write();

Info<< "Calculating explicit part of div(sigma) divSigmaExp\n" << endl;
volVectorField divSigmaExp
(
    "divSigmaExp",
    fvc::div(sigma) - fvc::laplacian(2*mu + lambda, D, "laplacian(DD,D)")
);
if (i > 1)
{
    divSigmaExp +=
        (lambda + 2*mu)*fvc::laplacian(pos0(Zf - laser.height(i - 1)), D, "laplacian(DD,D)");
}
divSigmaExp.write();

volVectorField divSigmaExp2
(
    "divSigmaExp2",
    fvc::div(sigma) - fvc::laplacian(2*mu + lambda, D, "laplacian(DD,D)")
);

if (i > 1)
{
    divSigmaExp2 -=
        (lambda + 2*mu)*fvc::div(epsilonf*mesh.magSf()*pos0(Zf - laser.height(i - 1)));
}
divSigmaExp2.write();

