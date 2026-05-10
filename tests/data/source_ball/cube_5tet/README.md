# cube_5tet reference mesh fixture

Unit cube split into 5 tetrahedra (Liu & Joe 1995, "Quality local
refinement of tetrahedral meshes based on bisection", J. Sci. Comput.
16(6)). Eight nodes, five cells. Used by
`Unit.SourceBall.DMPlexConstructionAndDistribution` to verify that
DMPlex partitioning and ghost-layer setup behave correctly at MPI=1
and MPI=2 (cell counts sum, no overlap of unique-owned cells).

This fixture is hand-written and does not require TetGen to
regenerate.
