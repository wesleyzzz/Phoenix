[MeshGenerators]
  [./cmg]
    type = CartesianMeshGenerator
    dim = 2
    dx = '1 1 1'
    dy = '1 1'
    ix = '20 20 20'
    iy = '20 20'
    subdomain_id = '0 1 2 0 0 2'
  [../]
[]

[Mesh]
  type = MeshGeneratorMesh
[]

