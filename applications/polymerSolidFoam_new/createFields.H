//Stereolithography
volScalarField exposure = volScalarField::New
(
    "exposure",
    mesh,
    dimensionedScalar(dimEnergy/dimArea, 0)
);


volScalarField totalExposure = volScalarField::New
(
    "totalExposure",
    mesh,
    dimensionedScalar(dimEnergy/dimArea)
);

volScalarField totalSqrtExposure = volScalarField::New
(
    "totalSqrtExposure",
     mesh,
     dimensionedScalar(sqrt(dimEnergy/dimArea))
);

volScalarField polymerization = volScalarField::New
(
    "polymerization",
    mesh,
    dimensionedScalar()
);

volScalarField nLayer = volScalarField::New
(
    "nLayer",
    mesh,
    dimensionedScalar(dimless, 1)
);

Info<< "Reading polymerProperties" << endl;
IOdictionary polymerProperties
(
    IOobject
    (
        "polymerProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

const dimensionedScalar Ec("criticalExposure", totalExposure.dimensions(), polymerProperties);
const dimensionedScalar Dp("penetrationDepth", dimLength, polymerProperties);

const laserScanner laser(mesh);


//polymerSolidFoam
Info<< "Reading field D\n" << endl;

volVectorField D
(
    IOobject
    (
        "D",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);
volTensorField gradD(fvc::grad(D));
volVectorField D_old = D;

const dimensionedScalar Emax("Emax", dimPressure, polymerProperties);
const dimensionedScalar e0("stiffnessFactor", dimless, polymerProperties);
const dimensionedScalar epsilonChemicalMax("epsilonChemicalMax", dimless, polymerProperties);

Info<< "Calculating mechanical properties\n" << endl;

volScalarField E = volScalarField::New
(
    "YoungsModulus",
    mesh,
    dimensionedScalar(dimPressure)
);

volScalarField nu = volScalarField::New
(
    "nu",
    mesh,
    dimensionedScalar("nu", dimless, polymerProperties)
);

volScalarField mu(E/(2.0*(1.0 + nu)));
volScalarField lambda(nu*E/((1.0 + nu)*(1.0 - 2.0*nu)));
volScalarField threeK(E/(1.0 - 2.0*nu));

surfaceScalarField muf = fvc::interpolate(mu);
surfaceScalarField lambdaf = fvc::interpolate(lambda);

volSymmTensorField sigma
(
    IOobject
    (
        "sigma",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    mesh,
    E.dimensions()
);
sigma.write();

/*volSymmTensorField epsilonResidual = volSymmTensorField::New
(
    "epsilonResidual",
    mesh,
    dimless
);
epsilonResidual =  symm(fvc::grad(D));
*/

volSymmTensorField epsilon("epsilon", symm(fvc::grad(D)));
epsilon.write();
surfaceVectorField epsilonf("epsilonf", fvc::snGrad(D));
symmTensor Ixy = symmTensor::I; Ixy.zz() = 0;

mesh.setFluxRequired(D.name());

// Flag to write auxiliary fields
bool writeProperties(runTime.controlDict().get<bool>("writeProperties"));

if (writeProperties)
{
    polymerization.writeOpt() = IOobject::AUTO_WRITE;
    nLayer.writeOpt() = IOobject::AUTO_WRITE;
    E.writeOpt() = IOobject::AUTO_WRITE;
    sigma.writeOpt() = IOobject::AUTO_WRITE;
    epsilon.writeOpt() = IOobject::AUTO_WRITE;
}
