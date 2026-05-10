# single_tet reference mesh fixture

Smallest valid TetGen `.node` / `.ele` pair the runtime parser
exercises. One tetrahedron, four nodes at the unit corner. Used by
`Unit.SourceBall.MeshLoadFromTetGen` to verify the file parser and
DMPlex construction at the smallest possible cell count.

This fixture is hand-written and does not require TetGen to
regenerate.
