/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2024 OpenFOAM Foundation
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
    Foam::decompositionMethods::kahip

Description
    ParHIP (parallel KaHIP) redistribution in parallel

    Note: like zoltan/parMetis, this method does not support serial
    operation - it operates on the graph as already split across
    processors and is only usable as a 'distributor', not a 'decomposer'.
    Use a built-in method (e.g. 'simple') for the initial decomposer.

    Parameters
    - method : The KaHIP preconfiguration
      - ultrafast, fast, eco : for mesh-like graphs
      - ultrafastSocial, fastSocial, ecoSocial : for social-network-like
        graphs
      Default: fast

    - seed : random seed. Default: 0

    The imbalance tolerance is fixed at 0.02 (0 for a single domain), as in
    the parMetis distributor.

Usage
    To select the KaHIP (ParHIP) distributor add the following entries to
    decomposeParDict:

        distributor     kahip;
        libs            ("libkahipDecomp_aalto.so");

    The method/seed can be changed via the optional kahip sub-dictionary, e.g.:

        kahip
        {
            method      fast;
        }

SourceFiles
    kahip.C

\*---------------------------------------------------------------------------*/

#ifndef kahip_H
#define kahip_H

#include "decompositionMethod.H"

namespace Foam
{
namespace decompositionMethods
{

/*---------------------------------------------------------------------------*\
                           Class kahip Declaration
\*---------------------------------------------------------------------------*/

class kahip
:
    public decompositionMethod
{
    dictionary methodDict_;

    // Private Member Functions

        //- Call ParHIP with options from dictionary. Operates on the
        //  distributed graph, already split across processors.
        label decompose
        (
            const labelList& xadj,
            const labelList& adjncy,
            const pointField& cellCentres,
            const scalarField& cellWeights,
            labelList& finalDecomp
        );


public:

    //- Runtime type information
    TypeName("kahip");


    // Constructors

        //- Construct given the decomposition dictionary
        kahip
        (
            const dictionary& decompositionDict,
            const dictionary& methodDict
        );

        //- Disallow default bitwise copy construction
        kahip(const kahip&) = delete;


    //- Destructor
    virtual ~kahip()
    {}


    // Member Functions

        //- Inherit decompose from decompositionMethod
        using decompositionMethod::decompose;

        //- Return for every coordinate the wanted processor number. Use the
        //  mesh connectivity (if needed)
        //  Weights get normalised so the minimum value is 1 before truncation
        //  to an integer so the weights should be multiples of the minimum
        //  value. The overall sum of weights might otherwise overflow.
        virtual labelList decompose
        (
            const polyMesh& mesh,
            const pointField& points,
            const scalarField& pointWeights
        );

        //- Return for every coordinate the wanted processor number. Gets
        //  passed agglomeration map (from fine to coarse cells) and coarse cell
        //  location. Can be overridden by decomposers that provide this
        //  functionality natively.
        //  See note on weights above.
        virtual labelList decompose
        (
            const polyMesh& mesh,
            const labelList& cellToRegion,
            const pointField& regionPoints,
            const scalarField& regionWeights
        );

        //- Return for every coordinate the wanted processor number. Explicitly
        //  provided mesh connectivity.
        //  The connectivity is equal to mesh.cellCells() except for
        //  - in parallel the cell numbers are global cell numbers (starting
        //    from 0 at processor0 and then incrementing all through the
        //    processors)
        //  - the connections are across coupled patches
        //  See note on weights above.
        virtual labelList decompose
        (
            const labelListList& globalCellCells,
            const pointField& cellCentres,
            const scalarField& cellWeights
        );


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const kahip&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace decompositionMethods
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
