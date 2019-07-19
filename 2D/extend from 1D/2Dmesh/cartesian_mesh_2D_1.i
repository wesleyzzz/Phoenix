[MeshGenerators]
  [./cmg]
    type = CartesianMeshGenerator
    dim = 2
    dx = '1 1 1'
    dy = '1 1'
    ix = '2 1 1'
    iy = '2 3'
    subdomain_id = '0 1 1 2 2 2'
  [../]
[]

[Mesh]
  type = MeshGeneratorMesh
[]
