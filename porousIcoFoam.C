/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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
    porousIcoFoam

Description
    Transient solver for incompressible, laminar flow of Newtonian fluids with 
    porous zones.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "pisoControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"

    pisoControl piso(mesh);

    #include "createFields.H"
    #include "initContinuityErrs.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        #include "CourantNo.H"

        // Momentum predictor
        fvVectorMatrix UEqn
        (
            (1.0/pty)*fvm::ddt(U)                   // Local acceleration term
          + (1.0/sqr(pty))*fvm::div(phi, U)         // Convective acceleration term
          //+ (1.0/pty)*fvm::div(phi/fvc::interpolate(pty), U)       // Convective acceleration term
          - (1.0/pty)*fvm::laplacian(nu, U)         // Viscous term (with Brinkman modification)
          + fvm::SuSp(nu/K_, U)                     // Darcy term
          + fvm::SuSp((cf/sqrt(K_))*mag(U), U)      // Forchheimer term
        );

        // STEP 1: Solve the momentum predictor [ M*U = -grad(p) ]
        if (piso.momentumPredictor())
        {
            solve(UEqn == -fvc::grad(p)); // Solve for updated velocity field
            // Note: Computed velocity field does not satisfy the continuity equation
        }

        // PISO loop
        while (piso.correct())
        {
            // STEP 2: Extract the A matrix and compute the H matrix [ H = A*U - M*U ]
            volScalarField rAU(1.0/UEqn.A()); // Extract the inverse of the diagonal components of the M matrix [ rAU = A^-1 ]
            volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p)); // Calculate A^-1*H
            
            surfaceScalarField phiHbyA
            (
                "phiHbyA",
                fvc::flux(HbyA)
              + fvc::interpolate(rAU)*fvc::ddtCorr(U, phi)
            );

            adjustPhi(phiHbyA, U, p); // Adjust the balance of fluxes to obey continuity

            constrainPressure(p, U, phiHbyA, rAU); // Update the pressure BCs to ensure flux consistency

            // Non-orthogonal pressure corrector loop
            while (piso.correctNonOrthogonal())
            {
                // Pressure corrector
                fvScalarMatrix pEqn
                (
                    fvm::laplacian(rAU, p) == fvc::div(phiHbyA) // Pressure equation (a form of the Laplace equation)
                );

                pEqn.setReference(pRefCell, pRefValue);

                // STEP 3: Solve the pressure equation [ div(A^-1*grad(p)) = div(A^-1*H) ]
                pEqn.solve(mesh.solver(p.select(piso.finalInnerIter()))); // Solve for the updated pressure field - OF2106
                //pEqn.solve(); // Solve for the updated pressure field - OF9

                if (piso.finalNonOrthogonalIter())
                {
                    phi = phiHbyA - pEqn.flux();
                }
            }

            #include "continuityErrs.H"

            // STEP 4: Use updated pressure field to correct the velocity field
            U = HbyA - rAU*fvc::grad(p); // Explicit pressure-corrector stage [ U = A^-1*H - A^-1*grad(p) ]
            U.correctBoundaryConditions();
            // Note: Now that the velocity field has been updated the pressure field is no longer valid since the H matrix depends on the velocity field

            // STEP 5: Use the flux corrected U to update H directly by performing 'inner loops' of the Pressure-Velocity Coupling Algorithm
        }

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
