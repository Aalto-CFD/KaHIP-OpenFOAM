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

// Both PstreamGlobals.H and parhip_interface.h pull in mpi.h. Prevent it from
// also dragging in the deprecated MPI C++ bindings, which fail to compile on
// some MPI stacks. The OpenFOAM MPI build rules normally define these already;
// setting them here too keeps the file safe to compile in isolation. Must
// precede any include that reaches mpi.h.
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


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::label Foam::decompositionMethods::kahip::decompose
(
    const labelList& xadj,
    const labelList& adjncy,
    const pointField& cellCentres,
    const scalarField& cellWeights,
    labelList& decomp
)
{
    // A distributor redistributes the existing mesh across the running
    // processors, so the number of subdomains must equal the number of MPI
    // ranks. Guarding here also keeps every partition index ParHIP returns
    // within [0, nProcs), which the per-processor addressing below indexes.
    if (nProcessors_ != Pstream::nProcs())
    {
        FatalErrorInFunction
            << "The kahip distributor redistributes onto the running "
            << Pstream::nProcs() << " processors, but numberOfSubdomains is "
            << nProcessors_ << "." << nl
            << "Set numberOfSubdomains equal to the number of processors."
            << exit(FatalError);
    }

    // Cell weights, if any, must provide exactly one weight per cell
    if (cellWeights.size() && cellWeights.size() != cellCentres.size())
    {
        FatalErrorInFunction
            << "Number of cell weights " << cellWeights.size()
            << " does not equal number of cells " << cellCentres.size()
            << "." << nl
            << "The kahip distributor supports a single weight per cell only."
            << exit(FatalError);
    }

    // ParHIP preconfiguration method
    word method("fast");
    methodDict_.readIfPresent("method", method);

    int kahipMode = FASTMESH;
    if (method == "ultrafast")
    {
        kahipMode = ULTRAFASTMESH;
    }
    else if (method == "fast")
    {
        kahipMode = FASTMESH;
    }
    else if (method == "eco")
    {
        kahipMode = ECOMESH;
    }
    else if (method == "ultrafastSocial")
    {
        kahipMode = ULTRAFASTSOCIAL;
    }
    else if (method == "fastSocial")
    {
        kahipMode = FASTSOCIAL;
    }
    else if (method == "ecoSocial")
    {
        kahipMode = ECOSOCIAL;
    }
    else
    {
        FatalIOErrorInFunction(methodDict_)
            << "Unknown KaHIP method " << method << nl
            << "Valid methods are: ultrafast, fast, eco, "
            << "ultrafastSocial, fastSocial, ecoSocial"
            << exit(FatalIOError);
    }

    // Allowed imbalance between the resulting subdomains (fraction)
    double imbalance = methodDict_.lookupOrDefault<scalar>("imbalance", 0.03);

    // Random seed
    int seed = methodDict_.lookupOrDefault<label>("seed", 0);

    // Verbosity of the ParHIP library itself
    const Switch verbose(methodDict_.lookupOrDefault<Switch>("verbose", false));
    const bool suppressOutput = !verbose;

    Info<< "kahip : Using ParHIP method     " << method << nl
        << "        Allowed imbalance       " << imbalance << nl
        << "        Random seed             " << seed << nl << endl;

    // Cell weights on the graph vertices. scaleWeights performs a global
    // (collective) reduction to scale by the total weight, so it must run on
    // every processor, before any communicator sub-setting below.
    label nWeights = 1;
    const labelList intWeights(scaleWeights(cellWeights, nWeights, true));

    // Are vertex weights in use on any processor? (collective)
    const bool useWeights =
        returnReduce(cellWeights.size(), sumOp<label>()) > 0;

    // Distributed-CSR vertex distribution: the global cell offset per
    // processor. globalIndex construction is an all-gather, so on all ranks.
    globalIndex globalMap(cellCentres.size());
    const labelList& cellOffsets = globalMap.offsets();

    // Restrict ParHIP to the processors that actually own cells. ParHIP runs
    // MPI collectives internally (e.g. the vwgt all-reduce) and computes its
    // vertex range with unsigned arithmetic, so a processor with zero cells
    // would both diverge from the collectives the populated ranks make (a
    // hang) and underflow that range. Build vtxdist over the valid processors
    // only; empty processors contribute no cells and so do not shift offsets.
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

    // Communicator over the cell-owning processors only. Allocating a fresh
    // communicator (even when every processor owns cells) also isolates
    // ParHIP's raw MPI_ANY_SOURCE ghost-node exchange from OpenFOAM's own
    // Pstream traffic, which would otherwise silently corrupt the partition.
    const label commi =
        Pstream::allocateCommunicator(UPstream::worldComm, validProcs);

    // Output: cell -> processor addressing
    List<idxtype> part(cellCentres.size(), idxtype(0));

    // Output: number of cut edges
    int edgeCut = 0;

    if (cellCentres.size())
    {
        // Convert the CSR graph to ParHIP's fixed 64-bit idxtype, which
        // differs from the 32-bit OpenFOAM label and so cannot be aliased.
        List<idxtype> xadjIdx(xadj.size());
        forAll(xadj, i)
        {
            xadjIdx[i] = idxtype(xadj[i]);
        }

        List<idxtype> adjncyIdx(adjncy.size());
        forAll(adjncy, i)
        {
            adjncyIdx[i] = idxtype(adjncy[i]);
        }

        // Vertex weights. Pass a value on every participating processor so
        // ParHIP takes the vwgt!=NULL branch (and its collective) on all of
        // them; fall back to uniform weights where none were supplied.
        List<idxtype> vwgt;
        if (useWeights)
        {
            vwgt.setSize(cellCentres.size(), idxtype(1));
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
            xadjIdx.begin(),
            adjncyIdx.begin(),
            useWeights ? vwgt.begin() : nullptr,   // vwgt
            nullptr,                               // adjwgt: edge weights
            &nParts,
            &imbalance,
            suppressOutput,
            seed,
            kahipMode,
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
    methodDict_(methodDict)
{}


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
        points,
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
        regionPoints,
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
        cellCentres,
        cellWeights,
        decomp
    );

    return decomp;
}


// ************************************************************************* //
