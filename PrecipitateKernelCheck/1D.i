#cannot run 

inactive = 'Functions'
[Mesh]
  type = GeneratedMesh
  nx = 10
  dim = 1
  xmax = 2
[]

[Variables]
  [Solute]
    order = FIRST
    family = LAGRANGE
  []
  [Precipitation]
    order = FIRST
    family = LAGRANGE
    inactive = 'InitialCondition'
    [InitialCondition]
      type = ConstantIC
      value = 0
      variable = Precipitation
    []
  []
[]

[Kernels]
  [Solute_diffusion]
    type = Diffusion
    variable = Solute
  []
  [Solute]
    type = ADPrecipitation
    unit_scalor = '1'
    scale_factor = 'scale_factor'
    variable = 'Solute'
    precipitate_variable = 'Precipitation'
    solubility = 'Solute_solubility'
  []
  [Precipitation]
    type = ADPrecipitation
    unit_scalor = '1'
    scale_factor = 'scale_factor'
    variable = 'Precipitation'
    dissolve_variable = 'Solute'
    solubility = 'Solute_solubility'
  []
[]

[Materials]
  [Solute_solubility]
    type = GenericConstantMaterial
    prop_names = 'Solute_solubility'
    prop_values = '1'
    outputs = 'exodus'
  []
  [scale_factor]
    type = GenericConstantMaterial
    prop_names = 'scale_factor'
    prop_values = '0.0795774729'
  []
[]

[BCs]
  [Solute_left]
    type = DirichletBC
    variable = Solute
    boundary = 'left'
    value = 0
  []
  [Solute_right]
    type = DirichletBC
    variable = Solute
    boundary = 'right'
    value = 2
  []
[]

[Executioner]
  type = Steady
  solve_type = PJFNK
[]

[Outputs]
  execute_on = 'timestep_end'
  exodus = true
[]

[VectorPostprocessors]
  [Solute]
    type = LineValueSampler
    num_points = 10
    end_point = '2 0 0'
    variable = 'Solute'
    sort_by = x
    start_point = '0 0 0'
  []
  [Precipitation]
    type = LineValueSampler
    num_points = 10
    end_point = '2 0 0'
    variable = 'Precipitation'
    sort_by = x
    start_point = '0 0 0'
  []
[]
