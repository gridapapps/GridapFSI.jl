function execute(problem::FSIProblem{:elasticFlagAggFem}; kwargs...)
 
   # Problem setting (Default FSI-2)
   println("Setting Elastic flag problem parameters")
   Um = _get_kwarg(:Um,kwargs,1.0)
   H = _get_kwarg(:H,kwargs,0.41)
   ⌀ = _get_kwarg(:D,kwargs,0.1)
   # Solid properties
   E_s = _get_kwarg(:E_s,kwargs,1.4e6)
   ν_s = _get_kwarg(:nu_s,kwargs,0.4)
   ρ_s = _get_kwarg(:rho_s,kwargs,1.0e4)
   # Fluid properties
   ρ_f = _get_kwarg(:rho_f,kwargs,1.0e3)
   Re = _get_kwarg(:Re,kwargs, 100.0)
   μ_f = ρ_f * Um * ⌀ / Re
   γ_f = _get_kwarg(:gamma_f,kwargs,1.0)
   # Mesh properties
   E_m = _get_kwarg(:E_m,kwargs,1.0)
   ν_m = _get_kwarg(:nu_m,kwargs,-0.1)
   α_m = _get_kwarg(:alpha_m,kwargs,1.0e-5)
   weight_strategy = _get_kwarg(:alpha_m_weight,kwargs,"constant")
   # Time stepping
   t0 = _get_kwarg(:t0,kwargs,0.0)
   tf = _get_kwarg(:tf,kwargs,0.1)
   dt = _get_kwarg(:dt,kwargs,0.1)
   θ  = _get_kwarg(:theta,kwargs,0.5)
   # Post-process
   is_vtk = _get_kwarg(:is_vtk,kwargs,false)

  ## Build geometry
  # Cilinder
  R = 0.05
  C = Point(0.2,0.2)
  cylinder = disk(R,x0=C)
  # Flag
  C1 = Point(0.2,0.19)
  d1 = VectorValue(0.4,0.0)
  d2 = VectorValue(0.0,0.02)
  flag = quadrilateral(x0=C1,d1=d1,d2=d2)
  # Fluid and solid domains
  solid_geo = setdiff(flag,cylinder,name="solid")
  not_fluid = union(cylinder,flag)
  fluid_geo = !(not_fluid,name="fluid")
  
  # Background model
  pmin = Point(0.0,0.0)
  pmax = Point(2.2,0.41)
  partition = (80,40)
  bgmodel = CartesianDiscreteModel(pmin,pmax,partition)
  dp = pmax - pmin
  h = dp[1]/80

  # Fluid and solid models
  cutgeo = cut(bgmodel,union(fluid_geo,solid_geo))
  model_fluid = DiscreteModel(cutgeo,"fluid")
  model_solid = DiscreteModel(cutgeo,"solid")
  writevtk(model_fluid,"tmp_fluid")
  writevtk(model_solid,"tmp_solid")
  
  # Boundary labels
  labels = get_face_labeling(bgmodel)
  add_tag_from_tags!(labels,"wall",[2,4,5,6])
  add_tag_from_tags!(labels,"inlet",[1,3,7])
  add_tag_from_tags!(labels,"outlet",[8])
  
  # Define BC functions
  println("Defining Boundary conditions")
  u_in(x, t) = VectorValue(1.5 * Um * x[2] * (H - x[2]) / ((H / 2)^2), 0.0)
  u_noSlip(x, t) = VectorValue(0.0, 0.0)
  u_in(t::Real) = x -> u_in(x, t)
  u_noSlip(t::Real) = x -> u_noSlip(x, t)
  
  # Triangulations
  println("Defining triangulations")
  Ω = Triangulation(bgmodel)
  Ωc = Triangulation(cutgeo)
  Ωs = Triangulation(model_solid)
  Ωf = Triangulation(model_fluid)
  #trian_Γi = InterfaceTriangulation(model_fluid,model_solid)
  ΓDi = EmbeddedBoundary(cutgeo,fluid_geo,solid_geo)
  ΓDf = EmbeddedBoundary(cutgeo,fluid_geo,cylinder)
  ΓDs = EmbeddedBoundary(cutgeo,solid_geo,cylinder)
  writevtk(ΓDi,"tmp_interface")
  writevtk(ΓDf,"tmp_fluid_boundary")
  writevtk(ΓDs,"tmp_solid_boundary")
  writevtk(Ωf,"tmp_fluid")
  writevtk(Ωs,"tmp_solid")
  writevtk(Ωc,"tmp_cut")

  # Quadratures
  println("Defining quadratures")
  order = _get_kwarg(:order,kwargs,2)
  degree = 2*order
  bdegree = 2*order
  dΩ = Measure(Ω,degree)
  dΩf = Measure(Ωf,degree)
  dΩs = Measure(Ωs,degree)
  dΓDi = Measure(ΓDi,bdegree)
  dΓDf = Measure(ΓDf,bdegree)
  dΓDs = Measure(ΓDs,bdegree)

  # Cell aggregation
  strategy = AggregateAllCutCells()
  aggregates = aggregate(strategy,cutgeo)
  
  ## FE Spaces
  println("Defining FE spaces")
  FEᵥ_std = FiniteElements(
    PhysicalDomain(),
    model_fluid,
    lagrangian,
    VectorValue{2,Float64},
    order
  )
  FEᵥ_ser = FiniteElements(
    PhysicalDomain(),
    model_fluid,
    lagrangian,
    VectorValue{2,Float64},
    order,
    space=:S
  )
  FEₚ_std = FiniteElements(
    PhysicalDomain(),
    model_fluid,
    lagrangian,
    Float64,
    order-1,
    space=:P
  )
  
  Vᵥ_std_ST = FESpace(model_fluid,FEᵥ_std,dirichlet_tags=["inlet","wall"])
  Vᵥ_ser_ST = FESpace(model_fluid,FEᵥ_ser,conformity=:L2)
  Q_std = FESpace(model_fluid,FEₚ_std,conformity=:L2)

  Vᵥ_ST = AgFEMSpace(Vᵥ_std_ST,aggregates,Vᵥ_ser_ST)
  Q = AgFEMSpace(Q_std,aggregates)

  Uᵥ_ST = TrialFESpace(Vᵥ_ST,[u_in(0.0),u_noSlip(0.0)])
  P = TrialFESpace(Q)

  Y_ST = MultiFieldFESpace([Vᵥ_ST,Q])
  X_ST = MultiFieldFESpace([Uᵥ_ST,P])

  ## Weak Forms
  # Stokes operator
  println("Defining Stokes operator")

  # Normals
  n_ΓDi = get_normal_vector(ΓDi)
  n_ΓDf = get_normal_vector(ΓDf)

  # Facet length
  h_ΓDi = CellField(get_array(∫(1)dΓDi),ΓDi)
  h_ΓDf = CellField(get_array(∫(1)dΓDf),ΓDf)
  # dim = num_cell_dims(ΓDi)
  # hΓᵢ = CellField( lazy_map(h->(h.^(-dim)),h_Γfs), Tₕ[:Γi])

  γ = 1.0
  f = VectorValue(0.0,0.0)
  function a_stokes(x,y)
    WeakForms.a_ST(x,y,μ_f,dΩf) + 
    WeakForms.a_ST_Γd(x,y,n_ΓDf,μ_f,γ,h_ΓDf,dΓDf) + 
    WeakForms.a_ST_Γd(x,y,n_ΓDi,μ_f,γ,h_ΓDi,dΓDi)
  end
  function l_stokes(y)
    WeakForms.l_ST(y,f,dΩf) + 
    WeakForms.l_ST_Γd(y,n_ΓDf,μ_f,γ,h_ΓDf,u_noSlip(0.0),dΓDf) + 
    WeakForms.l_ST_Γd(y,n_ΓDi,μ_f,γ,h_ΓDi,u_noSlip(0.0),dΓDi)
  end
  op_ST = AffineFEOperator(a_stokes,l_stokes,X_ST,Y_ST)

  # Solve Stokes problem
  xh = solve(op_ST)
  writevtk(trian_fluid,"stokes.vtu",cellfields=["vh"=>restrict(xh[1],trian_fluid),"ph"=>restrict(xh[2],trian_fluid)])
end
