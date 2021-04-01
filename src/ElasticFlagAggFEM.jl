function execute(problem::Problem{:elasticFlagAggFem}; kwargs...)
 
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
  C2 = Point(0.6,0.21)
  n1 = VectorValue(1.0,0.0)
  n2 = VectorValue(0.0,1.0)
  p1 = plane(x0=C1,v=-n1)
  p2 = plane(x0=C1,v=-n2)
  p3 = plane(x0=C2,v=n1)
  p4 = plane(x0=C2,v=n2)
  flag = intersect(intersect(p1,p2),intersect(p3,p4))
  #flag = !(quadrilateral(x0=C1,d1=d1,d2=d2))
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
  #cutgeo = cut(bgmodel,union(fluid_geo,solid_geo))
  cutgeo = cut(bgmodel,fluid_geo)
  model_fluid = DiscreteModel(cutgeo,"fluid")
  #model_solid = DiscreteModel(cutgeo,"solid")
  
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
  trian = Triangulation(bgmodel)
  trian_Ω = Triangulation(cutgeo)
  #trian_solid = Triangulation(model_solid)
  trian_fluid = Triangulation(model_fluid)
  #trian_Γi = InterfaceTriangulation(model_fluid,model_solid)
  trian_Γi = EmbeddedBoundary(cutgeo)
  #trian_boundary_Γi = get_left_boundary(trian_Γi)
  n_Γi = get_normal_vector(trian_Γi)
  hΓᵢ = reindex(cell_measure(trian_Γi,trian),trian_Γi)
  
  # Quadratures
  println("Defining quadratures")
  order = _get_kwarg(:order,kwargs,2)
  degree = 2*order
  bdegree = 2*order
  quad_Ω = CellQuadrature(trian_Ω,2*order)
  #quad_solid = CellQuadrature(trian_solid,degree)
  quad_fluid = CellQuadrature(trian_fluid,degree)
  quad_Γi = CellQuadrature(trian_Γi,bdegree)

  # Cell aggregation
  strategy = AggregateAllCutCells()
  aggregates = aggregate(strategy,cutgeo)
  
  # Test FE Spaces
  println("Defining FE spaces")
  Vf = TestFESpace(
  model=model_fluid,
  valuetype=VectorValue{2,Float64},
  reffe=:QLagrangian,
  order=order,
  conformity =:H1,
  dirichlet_tags=["inlet", "wall"],
  dof_space=:physical
  )
  # Vs = TestFESpace(
  # model=model_solid,
  # valuetype=VectorValue{2,Float64},
  # reffe=:Lagrangian,
  # order=order,
  # conformity =:H1,
  # dirichlet_tags=["fixed"]
  # )
  Vfser = TestFESpace(
  reffe=:SLagrangian,
  conformity=:L2, # we don't neet to impose continuity since we only use the cell dof basis / shapefuns
  valuetype=VectorValue{2,Float64},
  model=model_fluid,
  order=order,
  dof_space=:physical
  )
  Qf = TestFESpace(
  reffe=:PLagrangian,
  conformity=:L2,
  valuetype=Float64,
  model=model_fluid,
  order=order-1,
  dof_space=:physical
  )
  V = AgFEMSpace(Vf,aggregates,Vfser)
  Q = AgFEMSpace(Qf,aggregates)
  U = TrialFESpace(V,[u_in(0.0),u_noSlip(0.0)])
  P = TrialFESpace(Q)
  Y = MultiFieldFESpace([V,Q])
  X = MultiFieldFESpace([U,P])

  # Stokes operator
  println("Defining Stokes operator")
  function A_Ω(x,y)
    u, p = x
    v, q = y
    ∇(v)⊙∇(u) - q*(∇⋅u) - (∇⋅v)*p
  end
  function B_Ω(y)
    v, q = y
    0.0
  end
  γ = order*(order+1)
  function A_Γd(x,y)
    u, p = x
    v, q = y
    (γ/h)*v⋅u - v⋅(n_Γi⋅∇(u)) - (n_Γi⋅∇(v))⋅u + (p*n_Γi)⋅v + (q*n_Γi)⋅u
  end
  function B_Γd(y)
    v, q = y
    (γ/h)*v⋅u_noSlip(0.0) - (n_Γi⋅∇(v))⋅u_noSlip(0.0) + (q*n_Γi)⋅u_noSlip(0.0)
  end
  t_Ω = LinearFETerm(A_Ω,trian_Ω,quad_Ω)
  t_Γd = AffineFETerm(A_Γd,B_Γd,trian_Γi,quad_Γi)
  # res_ST(x,y) = WeakForms.stokes_residual(x,y,μ_f,VectorValue(0.0,0.0))
  # res_ST_Γd(x,y) = WeakForms.stokes_residual_Γd(x,y,n_Γi,μ_f,γ_f,hΓᵢ,u_noSlip(0.0))
  # jac_ST_Γd(x,dx,y) = WeakForms.stokes_residual_Γd(dx,y,n_Γi,μ_f,γ_f,hΓᵢ,u_noSlip(0.0))
  # t_ST_Ωf = FETerm(res_ST,trian_fluid,quad_fluid)
  # t_ST_Γd = LinearFETerm(res_ST_Γd,trian_Γi,quad_Γi)
  op_ST = FEOperator(X,Y,t_Ω,t_Γd)#,t_ST_Γd)

  # Solve Stokes problem
  xh = solve(op_ST)
  writevtk(trian_fluid,"stokes.vtu",cellfields=["vh"=>restrict(xh[1],trian_fluid),"ph"=>restrict(xh[2],trian_fluid)])
end
