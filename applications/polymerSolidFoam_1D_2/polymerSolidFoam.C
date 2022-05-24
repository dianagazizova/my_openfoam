/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author(s)
     \\/     M anipulation  |
-------------------------------------------------------------------------------
      solidDisplacementFoam | Copyright (C) 2011-2016 OpenFOAM Foundation
           polymerSolidFoam | Copyright (C) 2020 Oleg Rogozin
-------------------------------------------------------------------------------
License
    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    polymerSolidFoam

Description
    Segregated finite-volume solver of elastic small-strain deformation of a
    solid body with chemical stresses.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "simpleControl.H"
#include "laserScanner.H"

using constant::mathematical::pi;

tmp<volScalarField> BeerLambert(const volScalarField& Z, const dimensionedScalar Dp)
{
    return neg0(Z)*exp(Z/Dp);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "postProcess.H"

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "createFields.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    //Stereolithography
    const volVectorField coord = mesh.C();
    const surfaceVectorField coordf = mesh.Cf();
    const volScalarField Z = coord.component(2);
    const surfaceScalarField Zf = coordf.component(2);

    const boundBox& bounds = mesh.bounds();
    const dimensionedScalar ymin("ymin", dimLength, bounds.min().y());
    const dimensionedScalar ymax("ymax", dimLength, bounds.max().y());
    const dimensionedScalar zmin("zmin", dimLength, bounds.min().z());
    const dimensionedScalar zmax("zmax", dimLength, bounds.max().z());
    const dimensionedScalar L = ymax - ymin;
    const dimensionedScalar small("small", dimLength, SMALL);

    dimensionedScalar endOfTime = runTime.endTime();

    if (mag(zmin) > small)
    {
        FatalError
            << "Coordinate zmin is not equal to zero."
            << exit(FatalError);
    }

    for (label i = 1; laser.height(i) < zmax + small; i++)
    {
        nLayer += pos0(Z - laser.height(i) - small);
    }
    const dimensionedScalar numberOfLayers = gMax(nLayer);
    Info<< "NUMBER OF LAYERS: " << numberOfLayers << endl;
    //const scalar pGel(polymerProperties.get<scalar>("pGel"));
    const scalar pMax(polymerProperties.get<scalar>("pMax"));

    for (label i = 1; i <= numberOfLayers.value(); i++)
    {
        Info<< endl << "LAYER NUMBER: " << i << endl;
        #include "createNewLayer.H"
        #include "createControl.H"

        Info<< "\nCalculating displacement field\n" << endl;
        while (simple.loop())
        {
            Info<< "Iteration: " << runTime.value() << nl << endl;
            Info<< "creating matrix for the layer " << i << endl;
            fvVectorMatrix DEqn
            (
                fvm::laplacian(2*mu + lambda, D, "laplacian(DD,D)")
                + divSigmaExp
            );
            //#include "printMatrix.H"
            DEqn.solve();
            gradD = fvc::grad(D); // D on the current layer
            sigma = mu*twoSymm(gradD) + lambda*I*tr(gradD)
                  - threeK*epsilonChemicalMax/pMax*polymerization*I;

            divSigmaExp = fvc::div(sigma - (2*mu + lambda)*gradD, "div(sigma)");
            /*if (i > 1)
            {
                divSigmaExp +=
            (lambda + 2*mu)*fvc::laplacian(pos0(Zf - laser.height(i - 1)), D_old, "laplacian(DD,D)");
            }
            */
            if (i > 1)
            {
                divSigmaExp -=
        (lambda + 2*mu)*fvc::div(epsilonf*mesh.magSf()*pos0(Zf - laser.height(i - 1)));
            }
            runTime.printExecutionTime(Info);
        }
        runTime.stopAt(Foam::Time::stopAtControls::saEndTime);
        runTime.setEndTime(endOfTime);
        epsilon = symm(gradD);
        epsilon.write();

        epsilonf = -fvc::snGrad(D);
        epsilonf.write();

        //Info<< mesh.Sf() << endl;
        //Info<< mesh.magSf() << endl;
        //Info<< epsilonf << endl;
        //epsilonff = epsilonf & mesh.Sf();
        //epsilonff.write();
        //epsilonff = epsilonf & mesh.Sf();
        //epsilonff.write();
        //Info<< mesh.magSf().dimensions() << endl;;

        epsilonResidual.write();
        runTime.writeNow();
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
