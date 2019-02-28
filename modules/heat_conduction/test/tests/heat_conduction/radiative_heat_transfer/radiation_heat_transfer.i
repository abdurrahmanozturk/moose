#radiation heat transfer calculations for 2 parallel blocks
[Mesh]
  type = FileMesh
  file = parallel_blocks.e
  # type = GeneratedMesh
  # xmax = 1
  # xmin = 0
  # ymax = 1
  # ymin = 0
  # zmax = 1
  # zmin = 0
  # dim = 3
  # nx = 1
  # ny = 1
  # nz = 1
  # elem_type = HEX8
[]

[Variables]
  [./temp]
    initial_condition = 500
  [../]
[]

[Kernels]
  [./HeatConduction]
    type = HeatConduction
    variable = temp
    diffusion_coefficient = thermal_conductivity
  [../]
  [./TimeDerivative]
    type = TimeDerivative
    variable = temp
  []
[]
[NodalNormals]
[]
[BCs]
  [./master]
    type = DirichletBC
    value = 1000 #K
    variable = temp
    boundary = 1
  [../]
  # [./slave]
  #   type = DirichletBC
  #   value = 900 #K
  #   variable = temp
  #   boundary = 7
  # [../]
  [./RadiationHeatTransfer]
    type = RadiativeHeatFluxBC
    variable = temp
    boundary = '2 7'
    emissivity = '1 1'
    viewfactor_userobject = ViewFactor
  [../]
  # [./RadiativeBC]
  #   type = RadiativeBC
  #   variable = temp
  #   viewfactor_method = MONTECARLO
  #   boundary = '2 7'
  #   emissivity = '1 1'
  #   sampling_number = 10
  #   source_number = 10
  # [../]
[]
[Materials]
  [./uo2_thermal]
    type = UO2
    temp = temp
  [../]
[]
[Executioner]
  type = Transient
  solve_type = PJFNK
  start_time = 0
  end_time = 1
  # dt = 1e-3
  dtmin = 1e-6
  nl_abs_tol = 1e-15
[]
[UserObjects]
  [./ViewFactor]
    type = ViewFactor
    boundary = '2 7'
    method = MONTECARLO
    sampling_number = 10
    source_number = 10
    print_screen = true
    execute_on = INITIAL
  [../]
[]
[Postprocessors]
  [./boundarytemp_1]
    type = SideAverageValue
    boundary = '1'
    variable = temp
  [../]
  [./boundarytemp_2]
    type = SideAverageValue
    boundary = '2'
    variable = temp
  [../]
  [./boundarytemp_7]
    type = SideAverageValue
    boundary = '7'
    variable = temp
  [../]
  [./boundarytemp_9]
    type = SideAverageValue
    boundary = '9'
    variable = temp
  [../]
[]

[Outputs]
  exodus = true
  file_base = radiation_heat_transfer
  console = true
[]
