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

\*---------------------------------------------------------------------------*/

// Skip the deprecated MPI C++ bindings; must precede any mpi.h include
#ifndef MPICH_SKIP_MPICXX
#define MPICH_SKIP_MPICXX
#endif
#ifndef OMPI_SKIP_MPICXX
#define OMPI_SKIP_MPICXX
#endif

#include "kahip.H"
#include "Switch.H"
#include "globalIndex.H"
#include "PstreamGlobals.H"
#include "addToRunTimeSelectionTable.H"

#include "sigFpe.H"
#ifdef LINUX_GNUC
    #ifndef __USE_GNU
        #define __USE_GNU
    #endif
    #include <fenv.h>
#endif

extern "C"
{
    #include "parhip_interface.h"
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace decompositionMethods
{
    defineTypeNameAndDebug(kahip, 0);

    addToRunTimeSelectionTable
    (
        decompositionMethod,
        kahip,
        distributor
    );
}
}

const Foam::NamedEnum<Foam::decompositionMethods::kahip::kahipMethod, 6>
Foam::decompositionMethods::kahip::kahipMethodNames_
{
    "ultrafast",
    "fast",
    "eco",
    "ultrafastSocial",
    "fastSocial",
    "ecoSocial"
};

namespace
{
    // The ParHIP preconfiguration constant for each kahipMethod, in
    // kahipMethodNames_ order
    const int kahipModes[] =
    {
        ULTRAFASTMESH,
        FASTMESH,
        ECOMESH,
        ULTRAFASTSOCIAL,
        FASTSOCIAL,
        ECOSOCIAL
    };
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::label Foam::decompositionMethods::kahip::decompose
(
    const labelList& xadj,
    const labelList& adjncy,
    const label nCells,
    const scalarField& cellWeights,
    labelList& decomp
)
{
    // A distributor redistributes onto the running processors only
    if (nProcessors_ != Pstream::nProcs())
    {
        FatalErrorInFunction
            << "The kahip distributor redistributes onto the running "
            << Pstream::nProcs() << " processors, but numberOfSubdomains is "
            << nProcessors_ << "." << nl
            << "Set numberOfSubdomains equal to the number of processors."
            << exit(FatalError);
    }

    // Imbalance tolerance (fraction, as ParHIP expects)
    double imbalance = 0.02;

    // If only one processor there is no imbalance
    if (nProcessors_ == 1)
    {
        imbalance = 0;
    }

    // Are vertex weights in use on any processor? (collective)
    const bool useWeights =
        returnReduce(cellWeights.size(), sumOp<label>()) > 0;

    // Cell weights on the graph vertices; multiple weights per cell are
    // combined by summation as ParHIP supports a single weight only
    labelList intWeights;
    if (useWeights)
    {
        const label nWeights = nCells ? cellWeights.size()/nCells : 0;

        if (nCells && nWeights*nCells != cellWeights.size())
        {
            FatalErrorInFunction
                << "Number of cell weights " << cellWeights.size()
                << " is not a multiple of the number of cells " << nCells
                << exit(FatalError);
        }

        const label nGlobalWeights = returnReduce(nWeights, maxOp<label>());

        if (nCells && nWeights != nGlobalWeights)
        {
            FatalErrorInFunction
                << "Number of weights per cell " << nWeights
                << " differs between processors (max " << nGlobalWeights
                << ")" << exit(FatalError);
        }

        if (nGlobalWeights > 1)
        {
            Info<< "kahip : Combining " << nGlobalWeights
                << " weights per cell into one" << endl;
        }

        scalarField combinedWeights;
        if (nWeights > 1)
        {
            combinedWeights.setSize(nCells, 0);
            forAll(combinedWeights, i)
            {
                for (label wi = 0; wi < nWeights; wi++)
                {
                    combinedWeights[i] += cellWeights[nWeights*i + wi];
                }
            }
        }

        label nScaleWeights = 1;
        intWeights = scaleWeights
        (
            nWeights > 1 ? combinedWeights : cellWeights,
            nScaleWeights,
            true
        );
    }

    // Global cell offset per processor
    globalIndex globalMap(nCells);
    const labelList& cellOffsets = globalMap.offsets();

    // Build the vertex distribution over the cell-owning processors only
    labelList validProcs(Pstream::nProcs());
    List<idxtype> vtxdist(Pstream::nProcs() + 1);
    label nValidProcs = 0;
    for (label proci = 0; proci < Pstream::nProcs(); proci++)
    {
        if (cellOffsets[proci + 1] - cellOffsets[proci] > 0)
        {
            validProcs[nValidProcs] = proci;
            vtxdist[nValidProcs] = idxtype(cellOffsets[proci]);
            nValidProcs++;
        }
    }
    validProcs.setSize(nValidProcs);
    vtxdist[nValidProcs] = idxtype(cellOffsets[Pstream::nProcs()]);
    vtxdist.setSize(nValidProcs + 1);

    // Fresh communicator over the cell-owning processors only
    const label commi =
        Pstream::allocateCommunicator(UPstream::worldComm, validProcs);

    // Output: cell -> processor addressing
    List<idxtype> part(nCells, idxtype(0));

    // Output: number of cut edges
    int edgeCut = 0;

    if (nCells)
    {
        // Convert the CSR graph to ParHIP's 64-bit idxtype
        List<idxtype> xadjIndex(xadj);
        List<idxtype> adjncyIndex(adjncy);

        // Vertex weights; uniform where none were supplied
        List<idxtype> vwgt;
        if (useWeights)
        {
            vwgt.setSize(nCells, idxtype(1));
            forAll(intWeights, i)
            {
                vwgt[i] = idxtype(intWeights[i]);
            }
        }

        int nParts = nProcessors_;
        MPI_Comm comm = PstreamGlobals::MPICommunicators_[commi];

        // Switch off FPU error trapping to work around a divide-by-zero in
        // ParHIP's recursive initial partitioning (same issue as Zoltan's)
        #ifdef FE_NOMASK_ENV
        int oldExcepts = fedisableexcept
        (
            FE_DIVBYZERO
          | FE_INVALID
          | FE_OVERFLOW
        );
        #endif

        ParHIPPartitionKWay
        (
            vtxdist.begin(),
            xadjIndex.begin(),
            adjncyIndex.begin(),
            useWeights ? vwgt.begin() : nullptr,   // vwgt
            nullptr,                               // adjwgt: edge weights
            &nParts,
            &imbalance,
            suppressOutput_,
            seed_,
            kahipModes[unsigned(method_)],
            &edgeCut,
            part.begin(),
            &comm
        );

        #ifdef FE_NOMASK_ENV
        feenableexcept(oldExcepts);
        #endif
    }

    Pstream::freeCommunicator(commi);

    decomp.setSize(part.size());
    forAll(part, i)
    {
        decomp[i] = label(part[i]);
    }

    // Sum the number of cells allocated to each processor
    labelList nProcCells(Pstream::nProcs(), 0);

    forAll(decomp, i)
    {
        nProcCells[decomp[i]]++;
    }

    reduce(nProcCells, ListOp<sumOp<label>>());

    // If there are no cells allocated to this processor keep the first one
    // to ensure that all processors have at least one cell
    if (decomp.size() && nProcCells[Pstream::myProcNo()] == 0)
    {
        Pout<< "    No cells allocated to this processor"
               ", keeping first cell"
            << endl;
        decomp[0] = Pstream::myProcNo();
    }

    return label(edgeCut);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::decompositionMethods::kahip::kahip
(
    const dictionary& decompositionDict,
    const dictionary& methodDict
)
:
    decompositionMethod(decompositionDict),
    method_
    (
        kahipMethodNames_.lookupOrDefault
        (
            "method",
            methodDict,
            kahipMethod::fast
        )
    ),
    seed_(methodDict.lookupOrDefault<label>("seed", 0)),
    suppressOutput_(!methodDict.lookupOrDefault<Switch>("verbose", false))
{
    // ParHIP has no equivalent of parMetis's processorWeights
    if (methodDict.found("processorWeights"))
    {
        FatalIOErrorInFunction(methodDict)
            << "The kahip distributor does not support processorWeights"
            << exit(FatalIOError);
    }

    Info<< indent << "Using ParHIP method " << kahipMethodNames_[method_]
        << ", random seed " << seed_ << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::labelList Foam::decompositionMethods::kahip::decompose
(
    const polyMesh& mesh,
    const pointField& points,
    const scalarField& pointWeights
)
{
    if (points.size() != mesh.nCells())
    {
        FatalErrorInFunction
            << "Can use this decomposition method only for the whole mesh"
            << endl
            << "and supply one coordinate (cellCentre) for every cell." << endl
            << "The number of coordinates " << points.size() << endl
            << "The number of cells in the mesh " << mesh.nCells()
            << exit(FatalError);
    }

    // Make ParHIP distributed CSR (Compressed Storage Format) storage
    //   adjncy      : contains neighbours (= edges in graph), global indices
    //   xadj(celli) : start of information in adjncy for celli
    CompactListList<label> cellCells;
    calcCellCells
    (
        mesh,
        identityMap(mesh.nCells()),
        mesh.nCells(),
        true,
        cellCells
    );

    labelList decomp;
    decompose
    (
        cellCells.offsets(),
        cellCells.m(),
        points.size(),
        pointWeights,
        decomp
    );

    return decomp;
}


Foam::labelList Foam::decompositionMethods::kahip::decompose
(
    const polyMesh& mesh,
    const labelList& cellToRegion,
    const pointField& regionPoints,
    const scalarField& regionWeights
)
{
    if (cellToRegion.size() != mesh.nCells())
    {
        FatalErrorInFunction
            << "Size of cell-to-coarse map " << cellToRegion.size()
            << " differs from number of cells in mesh " << mesh.nCells()
            << exit(FatalError);
    }

    // Make ParHIP distributed CSR (Compressed Storage Format) storage
    //   adjncy      : contains neighbours (= edges in graph), global indices
    //   xadj(celli) : start of information in adjncy for celli
    CompactListList<label> cellCells;
    calcCellCells(mesh, cellToRegion, regionPoints.size(), true, cellCells);

    labelList decomp;
    decompose
    (
        cellCells.offsets(),
        cellCells.m(),
        regionPoints.size(),
        regionWeights,
        decomp
    );

    // Rework back into decomposition for original mesh
    labelList fineDistribution(cellToRegion.size());

    forAll(fineDistribution, i)
    {
        fineDistribution[i] = decomp[cellToRegion[i]];
    }

    return fineDistribution;
}


Foam::labelList Foam::decompositionMethods::kahip::decompose
(
    const labelListList& globalCellCells,
    const pointField& cellCentres,
    const scalarField& cellWeights
)
{
    if (cellCentres.size() != globalCellCells.size())
    {
        FatalErrorInFunction
            << "Inconsistent number of cells (" << globalCellCells.size()
            << ") and number of cell centres (" << cellCentres.size()
            << ")." << exit(FatalError);
    }

    // Make ParHIP distributed CSR (Compressed Storage Format) storage
    //   adjncy      : contains neighbours (= edges in graph), global indices
    //   xadj(celli) : start of information in adjncy for celli
    CompactListList<label> cellCells(globalCellCells);

    labelList decomp;
    decompose
    (
        cellCells.offsets(),
        cellCells.m(),
        cellCentres.size(),
        cellWeights,
        decomp
    );

    return decomp;
}


// ************************************************************************* //
