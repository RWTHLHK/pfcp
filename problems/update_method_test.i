[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
[]

[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 3
    elem_type = HEX8
    nx = 4
    ny = 4
    nz = 4
  []
[]

[Variables]
  [./d]
   order = FIRST
   family = LAGRANGE
  []
 []

[AuxVariables]
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
  [./tau0]
    family = MONOMIAL
    order = CONSTANT
  [../]
  [./phi]
    family = MONOMIAL
    order = CONSTANT
  [../]
  [./normal_stress]
    family = MONOMIAL
    order = CONSTANT
  [../]
  [./c0]
    family = MONOMIAL
    order = CONSTANT
  [../]
  [./c1]
    family = MONOMIAL
    order = CONSTANT
  [../]
[]

[Physics/SolidMechanics/QuasiStatic/all]
  strain = FINITE
  add_variables = true
  generate_output = stress_zz
  base_name = uncracked
[]

[Kernels]
  [./dot_dm]
    type = TimeDerivative
    variable = d
  [../]
  [./ACbulkm]
    type = AllenCahn
    variable = d
    f_name = F
  [../]
  [./ACInterfacem]
    type = ACInterface
    variable = d
    kappa_name = kappa_op
    mob_name = L
  [../]
[]
[AuxKernels]
  [./uncracked_pk2]
   type = RankTwoAux
   variable = c0
   rank_two_tensor = second_piola_kirchhoff_stress
   index_j = 2
   index_i = 2
   execute_on = timestep_end
  [../]
  [./cracked_pk2]
    type = RankTwoAux
    variable = c1
    rank_two_tensor = cracked_pk2
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
  [./shear_stress0]
    type = MaterialStdVectorAux
    index = 0
    variable = tau0
    property = uncracked_applied_shear_stress
    execute_on = timestep_end
  [../]
  [./normal_stress0]
    type = MaterialStdVectorAux
    index = 0
    variable = normal_stress
    property = resolved_normal_stress
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
  [./cracked_stress]
    type = ComputeGeneralizedOrowanCrackedStress
    uncracked_base_name = uncracked
    number_slip_systems = 12
    dot_m0 = 0.2
    alpha = 1.0
    cm = 0.3
    tau_d = 71.0
    pm = 3
    d = d
    d_duc = 0.8

    dot_o0 = 0.5
    beta = 1.0
    co = 0.5
    sigma_d = 53.0
    po = 5
    slip_sys_file_name = input_slip_sys.txt
  [../]

  [./pfbulkmat]
    type = GenericConstantMaterial
    prop_names = 'gc_prop l visco d_duc'
    prop_values = '1e-3 0.05 1e-6 0.8'
  [../]

  [./define_mobility]
    type = ParsedMaterial
    material_property_names = 'gc_prop visco'
    property_name = L
    expression = '1.0/(gc_prop * visco)'
  [../]

  [./define_kappa]
    type = ParsedMaterial
    material_property_names = 'gc_prop l'
    property_name = kappa_op
    expression = 'gc_prop * l'
  [../]
  [./crack_formation_driving_energy]
    type = DerivativeParsedMaterial
    property_name = crack_formation_driving_energy
    material_property_names = 'd_duc mf'
    coupled_variables = 'd'
    constant_names = 'cm'
    constant_expressions = '0.3'
    expression = '1.0/d_duc^2 * (d_duc - d)^2 * cm * mf'
    derivative_order = 2
  [../]
  [./local_fracture_energy]
    type = DerivativeParsedMaterial
    property_name = local_fracture_energy
    coupled_variables = 'd'
    material_property_names = 'gc_prop l'
    expression = 'd^2 * gc_prop / 2 / l'
    derivative_order = 2
  [../]
  [./crack_open_driving_energy]
    type = DerivativeParsedMaterial
    property_name = crack_open_driving_energy
    material_property_names = 'ro'
    coupled_variables = 'd'
    constant_names = 'co'
    constant_expressions = '0.5'
    expression = '(1 - d)^2 * co * ro'
    derivative_order = 2
  [../]
  [./fracture_driving_energy]
    type = DerivativeSumMaterial
    coupled_variables = 'd'
    sum_materials = 'crack_formation_driving_energy crack_open_driving_energy local_fracture_energy'
    derivative_order = 2
    property_name = F
  [../]
[]

[Postprocessors]
  [./stress_zz]
    type = ElementAverageValue
    variable = uncracked_stress_zz
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
  [./tau]
    type = ElementAverageValue
    variable = tau0
  [../]
  [./mf]
    type = ElementAverageMaterialProperty
    mat_prop = mf
    execute_on = timestep_end
  [../]
  [./ro]
    type = ElementAverageMaterialProperty
    mat_prop = ro
    execute_on = timestep_end
  [../]
  [./d]
    type = ElementAverageValue
    variable = d
  [../]
  [./norm]
    type = ElementAverageValue
    variable = normal_stress
  [../]
  [./phi_pos]
    type = ElementAverageMaterialProperty
    mat_prop = cp_phi_pos
  [../]
  [./uncracked_pk2]
    type = ElementAverageValue
    variable = c0
  [../]
  [./cracked_pk2]
    type = ElementAverageValue
    variable = c1
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

  dt = 0.025
  dtmin = 0.001
  dtmax = 10.0
  num_steps = 50
  [./Predictor]
    type = SimplePredictor
    scale = 1.0
  [../]
[]

[Outputs]
  exodus = true
  csv = true
[]
