[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
[]

[Mesh]
  type = GeneratedMesh
  dim = 3
  elem_type = HEX8
[]

[AuxVariables]
  [./uncracked_pk2]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./uncracked_fp_zz]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./uncracked_e_zz]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./uncracked_gss]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./uncracked_slip_increment]
   order = CONSTANT
   family = MONOMIAL
  [../]
[]

[Physics/SolidMechanics/QuasiStatic/all]
  strain = FINITE
  add_variables = true
  generate_output = stress_zz
  base_name = uncracked
[]

[AuxKernels]
  [./uncracked_pk2]
   type = RankTwoAux
   variable = uncracked_pk2
   rank_two_tensor = second_piola_kirchhoff_stress
   index_j = 2
   index_i = 2
   execute_on = timestep_end
  [../]
  [./uncracked_fp_zz]
    type = RankTwoAux
    variable = uncracked_fp_zz
    rank_two_tensor = plastic_deformation_gradient
    index_j = 2
    index_i = 2
    execute_on = timestep_end
  [../]
  [./uncracked_e_zz]
    type = RankTwoAux
    variable = uncracked_e_zz
    rank_two_tensor = total_lagrangian_strain
    index_j = 2
    index_i = 2
    execute_on = timestep_end
  [../]
  [./uncracked_gss]
   type = MaterialStdVectorAux
   variable = uncracked_gss
   property = uncracked_slip_resistance
   index = 0
   execute_on = timestep_end
  [../]
  [./uncracked_slip_inc]
   type = MaterialStdVectorAux
   variable = uncracked_slip_increment
   property = uncracked_slip_increment
   index = 0
   execute_on = timestep_end
  [../]
[]

[BCs]
  [./symmy]
    type = DirichletBC
    variable = disp_y
    boundary = bottom
    value = 0
  [../]
  [./symmx]
    type = DirichletBC
    variable = disp_x
    boundary = left
    value = 0
  [../]
  [./symmz]
    type = DirichletBC
    variable = disp_z
    boundary = back
    value = 0
  [../]
  [./tdisp]
    type = FunctionDirichletBC
    variable = disp_z
    boundary = front
    function = '0.01*t'
  [../]
[]

[Materials]
  [./elasticity_tensor]
    type = ComputeElasticityTensorCP
    C_ijkl = '1.684e5 1.214e5 1.214e5 1.684e5 1.214e5 1.684e5 0.754e5 0.754e5 0.754e5'
    fill_method = symmetric9
    base_name = uncracked
  [../]
  [./trial_xtalpl]
    type = CrystalPlasticityKalidindiUpdate
    number_slip_systems = 12
    slip_sys_file_name = input_slip_sys.txt
    base_name = uncracked
  [../]
  [./stress]
    type = ComputeMultipleCrystalPlasticityStress
    crystal_plasticity_models = 'trial_xtalpl'
    tan_mod_type = exact
    base_name = uncracked
  [../]
  [./phi_pos]
    type = ComputeNeoHookeanTensileStrainEnergy
    dimension = 3
    nH1 = 1.0
    nH2 = 1.0
    base_name = uncracked
  [../]
[]

[Postprocessors]
  [./stress_zz]
    type = ElementAverageValue
    variable = uncracked_stress_zz
  [../]
  [./pk2]
   type = ElementAverageValue
   variable = uncracked_pk2
  [../]
  [./fp_zz]
    type = ElementAverageValue
    variable = uncracked_fp_zz
  [../]
  [./e_zz]
    type = ElementAverageValue
    variable = uncracked_e_zz
  [../]
  [./gss]
    type = ElementAverageValue
    variable = uncracked_gss
  [../]
  [./slip_increment]
   type = ElementAverageValue
   variable = uncracked_slip_increment
  [../]
[]

[Preconditioning]
  [./smp]
    type = SMP
    full = true
  [../]
[]

[Executioner]
  type = Transient
  solve_type = 'PJFNK'

  petsc_options_iname = '-pc_type -pc_asm_overlap -sub_pc_type -ksp_type -ksp_gmres_restart'
  petsc_options_value = ' asm      2              lu            gmres     200'
  nl_abs_tol = 1e-10
  nl_rel_tol = 1e-10
  nl_abs_step_tol = 1e-10

  dt = 0.05
  dtmin = 0.01
  dtmax = 10.0
  num_steps = 10
[]

[Outputs]
  exodus = true
[]
