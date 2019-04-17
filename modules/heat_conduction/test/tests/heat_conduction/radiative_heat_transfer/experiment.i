#radiation heat transfer calculations for 2 parallel blocks
[Mesh]
  type = FileMesh
  file = cylinder.e
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
    initial_condition = 400
  [../]
[]

[Kernels]
  [./HeatConduction]
    type = HeatConduction
    variable = temp
    diffusion_coefficient = thermal_conductivity
  [../]
  # [./TimeDerivative]
  #   type = TimeDerivative
  #   variable = temp
  # []
  # [./HeatSource]
  #   type = HeatSource
  #   variable = temp
  #   value = 1440
  #   block = 'innercylinder'
  # [../]
[]
[BCs]
  # [./master]
  #   type = DirichletBC
  #   value = 400 #K
  #   variable = temp
  #   boundary = 5
  # [../]
  # [./slave]
  #   type = DirichletBC
  #   value = 900 #K
  #   variable = temp
  #   boundary = 7
  # [../]
  # [./RadiationHeatTransfer]
  #   type = RadiativeHeatFluxBC
  #   variable = temp
  #   boundary = '5 2'
  #   emissivity = '1 1'
  #   viewfactor_userobject = ViewFactor
  # [../]
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
  type = Steady
  solve_type = PJFNK
  # start_time = 0
  # end_time = 10
  # # dt = 1e-3
  # dtmin = 1e-6
  # nl_abs_tol = 1e-15
[]
[UserObjects]
  [./ViewFactor]
    type = ViewFactor
    boundary = '5 2'
    method = MONTECARLO
    sampling_number = 10
    source_number = 10
    print_screen = true
    debug_mode = false
    execute_on = INITIAL
  [../]
[]
[Postprocessors]
  [./pellet_top]
    type = SideAverageValue
    boundary = '1'
    variable = temp
  [../]
  [./pellet_bottom]
    type = SideAverageValue
    boundary = '2'
    variable = temp
  [../]
  [./wall_bottom]
    type = SideAverageValue
    boundary = '5'
    variable = temp
  [../]
  [./wall_top]
    type = SideAverageValue
    boundary = '7'
    variable = temp
  [../]
[]

[Outputs]
  exodus = true
  file_base = experiment_out
  console = true
[]
