/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2018 OpenFOAM Foundation
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

Class
    Foam::solidificationSystem

Description
    A two-phase incompressible transport model for interacting phases
    having constant thermomechanical properties

SourceFiles
    solidificationSystem.C

\*---------------------------------------------------------------------------*/

#ifndef solidificationSystem_H
#define solidificationSystem_H

#include "twoPhaseMixture.H"
#include "IOdictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
         Class solidificationSystem Declaration
\*---------------------------------------------------------------------------*/

class solidificationSystem
:
    public IOdictionary,
    public twoPhaseMixture
{
protected:

    // Protected data

        dimensionedScalar rho1_;
        dimensionedScalar rho2_;

        dimensionedScalar Cp1_;
        dimensionedScalar Cp2_;

        dimensionedScalar kappa1_;
        dimensionedScalar kappa2_;

        dimensionedScalar mu1_;
        dimensionedScalar mu2_;

        dimensionedScalar D1_;
        dimensionedScalar D2_;

        dimensionedScalar DAS_;

        dimensionedScalar betaT_;
        dimensionedScalar betaC_;

        dimensionedScalar TRef_;
        dimensionedScalar CRef_;

        const volVectorField& U_;


public:

    TypeName("solidificationSystem");


    // Constructors

        //- Construct from components
        solidificationSystem
        (
            const volVectorField& U
        );


    //- Destructor
    virtual ~solidificationSystem()
    {}


    // Member Functions

        //- Return const-access to the dispersed-phase density
        const dimensionedScalar& rho1() const
        {
            return rho1_;
        }

        //- Return const-access to continuous-phase density
        const dimensionedScalar& rho2() const
        {
            return rho2_;
        };

        const dimensionedScalar& Cp1() const
        {
            return Cp1_;
        }

        //- Return const-access to continuous-phase density
        const dimensionedScalar& Cp2() const
        {
            return Cp2_;
        };

        const dimensionedScalar& kappa1() const
        {
            return kappa1_;
        }

        //- Return const-access to continuous-phase density
        const dimensionedScalar& kappa2() const
        {
            return kappa2_;
        };

        //- Return const-access to the dispersed-phase dynamic viscosity
        const dimensionedScalar& mu1() const
        {
            return mu1_;
        }

        //- Return const-access to continuous-phase dynamic viscosity
        const dimensionedScalar& mu2() const
        {
            return mu2_;
        };

        //- Return const-access to the dispersed-phase mass diffusivity
        const dimensionedScalar& D1() const
        {
            return D1_;
        }

        //- Return const-access to continuous-phase mass diffusivity
        const dimensionedScalar& D2() const
        {
            return D2_;
        };

        //- Return const-access to dnedrite arm spacing
        const dimensionedScalar& DAS() const
        {
            return DAS_;
        };

        //- Return const-access to liquid thermal expansion
        const dimensionedScalar& betaT() const
        {
            return betaT_;
        };

        //- Return const-access to liquid solutal expansion
        const dimensionedScalar& betaC() const
        {
            return betaC_;
        };

        //- Return const-access to reference temperature
        const dimensionedScalar& TRef() const
        {
            return TRef_;
        };

        //- Return const-access to reference composition
        const dimensionedScalar& CRef() const
        {
            return CRef_;
        };

        //- Return const-access to the mixture velocity
        const volVectorField& U() const
        {
            return U_;
        }

        //- Return the mixture density
        virtual tmp<volScalarField> rho() const
        {
            return alpha1_*rho1_ + alpha2_*rho2_;
        }

        //- Return the mixture density
        virtual tmp<volScalarField> Cp() const
        {
            return alpha1_*Cp1_ + alpha2_*Cp2_;
        }

        //- Return the mixture density
        virtual tmp<volScalarField> kappa() const
        {
            return alpha1_*kappa1_ + alpha2_*kappa2_;
        }

        //- Read base transportProperties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
