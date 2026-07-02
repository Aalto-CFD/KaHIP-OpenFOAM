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
    // ParHIP preconfiguration mode
    word mode("fastMesh");
    methodDict_.readIfPresent("mode", mode);

    int kahipMode = FASTMESH;
    if (mode == "ultrafastMesh")
    {
        kahipMode = ULTRAFASTMESH;
    }
    else if (mode == "fastMesh")
    {
        kahipMode = FASTMESH;
    }
    else if (mode == "ecoMesh")
    {
        kahipMode = ECOMESH;
    }
    else if (mode == "ultrafastSocial")
    {
        kahipMode = ULTRAFASTSOCIAL;
    }
    else if (mode == "fastSocial")
    {
        kahipMode = FASTSOCIAL;
    }
    else if (mode == "ecoSocial")
    {
        kahipMode = ECOSOCIAL;
    }
    else
    {
        FatalIOErrorInFunction(methodDict_)
            << "Unknown KaHIP mode " << mode << nl
            << "Valid modes are: ultrafastMesh, fastMesh, ecoMesh, "
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

    if (!methodDict_.empty())
    {
        Info<< "kahip : Using ParHIP mode       " << mode << nl
            << "        Allowed imbalance       " << imbalance << nl
            << "        Random seed             " << seed << nl << endl;
    }

    // vtxdist: global cell offset per processor (size nProcs+1), the same
    // distributed-CSR convention as ParMETIS
    globalIndex globalMap(cellCentres.size());
    const labelList& offsets = globalMap.offsets();

    List<idxtype> vtxdist(offsets.size());
    forAll(offsets, i)
    {
        vtxdist[i] = idxtype(offsets[i]);
    }

    // ParHIP's graph arrays use a fixed 64-bit idxtype, distinct from the
    // 32-bit OpenFOAM label, so they must be converted rather than aliased.
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

    // Cell weights (so on the vertices of the graph). The graph is
    // distributed across processors, so the weights must be scaled using
    // the global (not per-processor local) sum.
    label nWeights = 1;
    const labelList intWeights(scaleWeights(cellWeights, nWeights, true));

    List<idxtype> vwgt;
    if (intWeights.size())
    {
        vwgt.setSize(intWeights.size());
        forAll(intWeights, i)
        {
            vwgt[i] = idxtype(intWeights[i]);
        }
    }

    int nParts = nProcessors_;

    // Output: number of cut edges
    int edgeCut = 0;

    // Output: cell -> processor addressing
    List<idxtype> part(cellCentres.size(), idxtype(0));

    // Give ParHIP its own communicator, duplicated from OpenFOAM's. ParHIP
    // uses raw MPI_ANY_SOURCE point-to-point calls internally for its
    // ghost-node exchange; sharing OpenFOAM's own communicator risks its
    // messages being intercepted by unrelated pending Pstream traffic on
    // the same tags, which silently corrupts the result rather than
    // crashing.
    MPI_Comm comm;
    MPI_Comm_dup(PstreamGlobals::MPI_COMM_FOAM, &comm);

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
        vwgt.size() ? vwgt.begin() : nullptr,  // vwgt
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

    MPI_Comm_free(&comm);

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
