"""
Executes a transient FSI driver with a constant function in space, linear in time:
  u(x,t) = [1.0, -1.0]^T * t
  v(x,t) = [1.0, -1.0]^T
  p(x,t) = 0.0
"""
function execute(problem::Problem{:analytical};kwargs...)

  # Problem setting
  println("Setting analytical fluid problem parameters")
  # Solid properties
  E_s = _get_kwarg(:E_s,kwargs,1.0)
  ν_s = _get_kwarg(:nu_s,kwargs,0.4)
  ρ_s = _get_kwarg(:rho_s,kwargs,1.0)
  # Fluid properties
  ρ_f = _get_kwarg(:rho_f,kwargs,1.0)
  μ_f = _get_kwarg(:mu_f,kwargs,1.0)
  # Mesh properties
  n_m = _get_kwarg(:n_m,kwargs,10)
  E_m = _get_kwarg(:E_m,kwargs,1.0)
  ν_m = _get_kwarg(:nu_m,kwargs,-0.1)
  α_m = _get_kwarg(:alpha_m,kwargs,1.0e-5)
  # Time stepping
  t0 = _get_kwarg(:t0,kwargs,0.0)
  tf = _get_kwarg(:tf,kwargs,0.5)
  dt = _get_kwarg(:dt,kwargs,0.1)
  # Post-process
  is_vtk = _get_kwarg(:is_vtk,kwargs,false)

  # Mesh strategy
  strategyName = _get_kwarg(:strategy,kwargs,"laplacian")
  strategy = WeakForms.MeshStrategy{Symbol(strategyName)}()

  # Define BC functions
  println("Defining Boundary conditions")
  #u(x, t) = 1.0e-2 * VectorValue( x[1]^2*x[2] , -x[1]*x[2]^2) * t
  #v(x, t) = 1.0e-2 * VectorValue( x[1]^2*x[2] , -x[1]*x[2]^2)
  u(x, t) = 1.0e-2 * VectorValue( 1.0 , -1.0) * t
  v(x, t) = 1.0e-2 * VectorValue( 1.0 , -1.0)
  u(t::Real) = x -> u(x,t)
  v(t::Real) = x -> v(x,t)
  ∂tu(t) = x -> VectorValue(ForwardDiff.derivative(t -> get_array(u(x,t)),t))
  ∂tv(t) = x -> VectorValue(ForwardDiff.derivative(t -> get_array(v(x,t)),t))
  ∂tu(x,t) = ∂tu(t)(x)
  ∂tv(x,t) = ∂tv(t)(x)
  T_u = typeof(u)
  T_v = typeof(v)
  @eval ∂t(::$T_u) = $∂tu
  @eval ∂t(::$T_v) = $∂tv
  p(x,t) = 0.0
  p(t::Real) = x -> p(x,t)
  bconds = get_boundary_conditions(problem,strategy,u,v)

  # Define Forcing terms
  I = TensorValue( 1.0, 0.0, 0.0, 1.0 )
  F(t) = x -> ∇(u(t))(x) + I
  J(t) = x -> det(F(t))(x)
  E(t) = x -> 0.5 * ((F(t)(x)')⋅F(t)(x) - I)
  (λ_s, μ_s) = WeakForms.lame_parameters(E_s,ν_s)
  (λ_m, μ_m) = WeakForms.lame_parameters(E_m,ν_m)
  S_SV(t) = x -> 2*μ_s*E(t)(x) + λ_s*tr(E(t)(x))*I
  fv_ST_Ωf(t) = x -> - μ_f*Δ(v(t))(x) + ∇(p(t))(x)
  function fu_closure(t,strategy)
    fu = if typeof(strategy) == WeakForms.MeshStrategy{:laplacian}
      x -> - α_m * Δ(u(t))(x)
    elseif typeof(strategy) == WeakForms.MeshStrategy{:linearElasticity}
      x -> - μ_m * Δ(u(t))(x)
    elseif typeof(strategy) == WeakForms.MeshStrategy{:neoHookean}
      x -> - μ_m * Δ(u(t))(x)
    else
      @notimplemented("The soruce term for $strategy strategy is not implemented")
    end
    return fu
  end
  fu_Ωf(t) = fu_closure(t,strategy)
  #fu_Ωf(t) = VectorValue(0.0, 0.0)
  # if ( typeof(strategy) == WeakForms.MeshStrategy{:laplacian} )
  #   fu_Ωf(t) = x -> - α_m * Δ(u(t))(x)
  # elseif ( typeof(strategy) == WeakForms.MeshStrategy{:linearElasticity} )
  #   fu_Ωf(t) = x -> - μ_m * Δ(u(t))(x)
  #
  fv_Ωf(t) = x -> ρ_f * ∂t(v)(t)(x) - μ_f * Δ(v(t))(x) + ∇(p(t))(x) + ρ_f*( (∇(v(t))(x)')⋅(v(t)(x) - ∂t(u)(t)(x)) )
  fp_Ωf(t) = x -> (∇⋅v(t))(x)
  fu_Ωs(t) = x -> ∂t(u)(t)(x) - v(t)(x)
  fv_Ωs(t) = x -> ρ_s * ∂t(v)(t)(x) #- (∇⋅(F(t)⋅S_SV(t)))(x)  # Divergence of a a doted function not working yet...

  # Discrete model
  println("Defining discrete model")
  domain = (-1,1,-1,1)
  partition = (n_m,n_m)
  model = CartesianDiscreteModel(domain,partition)
  trian = Triangulation(model)
  R = 0.5
  xc = 0.0
  yc = 0.0
  function is_in(coords)
    n = length(coords)
    x = (1/n)*sum(coords)
    d = (x[1]-xc)^2 + (x[2]-yc)^2 - R^2
    d < 1.0e-8
  end
  oldcell_to_coods = get_cell_coordinates(trian)
  oldcell_to_is_in = collect1d(apply(is_in,oldcell_to_coods))
  incell_to_cell = findall(oldcell_to_is_in)
  outcell_to_cell = findall(collect(Bool, .! oldcell_to_is_in))
  model_solid = DiscreteModel(model,incell_to_cell)
  model_fluid = DiscreteModel(model,outcell_to_cell)

  # Build fluid-solid interface labelling
  println("Defining Fluid-Solid interface")
  labeling = get_face_labeling(model_fluid)
  new_entity = num_entities(labeling) + 1
  topo = get_grid_topology(model_fluid)
  D = num_cell_dims(model_fluid)
  for d in 0:D-1
    fluid_boundary_mask = collect(Bool,get_isboundary_face(topo,d))
    fluid_outer_boundary_mask = get_face_mask(labeling,"boundary",d)
    fluid_interface_mask = collect(Bool,fluid_boundary_mask .!= fluid_outer_boundary_mask)
    dface_list = findall(fluid_interface_mask)
    for face in dface_list
      labeling.d_to_dface_to_entity[d+1][face] = new_entity
    end
  end
  add_tag!(labeling,"interface",[new_entity])

  # Triangulations
  println("Defining triangulations")
  trian_solid = Triangulation(model_solid)
  trian_fluid = Triangulation(model_fluid)
  trian_Γi = BoundaryTriangulation(model_fluid,labeling,["interface"])
  n_Γi = get_normal_vector(trian_Γi)

  # Compute cell area (auxiliar quantity for mesh motion eq.)
  volf = cell_measure(trian_fluid,trian)
  vols = cell_measure(trian_solid,trian)
  vol_fluid = reindex(volf,trian_fluid)
  vol_solid = reindex(vols,trian_solid)
  vol_Γi = reindex(volf,trian_Γi)

  # Quadratures
  println("Defining quadratures")
  order = _get_kwarg(:order,kwargs,2)
  degree = 2*order
  bdegree = 2*order
  quad = CellQuadrature(trian,degree)
  quad_solid = CellQuadrature(trian_solid,degree)
  quad_fluid = CellQuadrature(trian_fluid,degree)
  quad_Γi = CellQuadrature(trian_Γi,bdegree)

  # Test FE Spaces
  println("Defining FE spaces")
  Y_ST, X_ST, Y_FSI, X_FSI = get_FE_spaces(problem,strategy,model,model_fluid,order,bconds)

  # Stokes problem for initial solution
  println("Defining Stokes operator")
  res_ST(x,y) = WeakForms.stokes_residual(strategy,x,y,μ_f,fv_ST_Ωf(0.0))
  jac_ST(x,dx,y) = WeakForms.stokes_jacobian(strategy,dx,y,μ_f)
  t_ST_Ωf = FETerm(res_ST, jac_ST, trian_fluid, quad_fluid)
  op_ST = FEOperator(X_ST,Y_ST,t_ST_Ωf)

  # Setup equation parameters
  fsi_f_params = Dict("μ"=>μ_f,
  "ρ"=>ρ_f,
  "E"=>E_m,
  "ν"=>ν_m,
  "α"=>α_m,
  "fu"=>fu_Ωf,
  "fv"=>fv_Ωf,
  "vol"=>vol_fluid)
  fsi_s_params = Dict("ρ"=>ρ_s,
  "E"=>E_s,
  "ν"=>ν_s,
  "fu"=>fu_Ωs,
  "fv"=>fv_Ωs,
  "vol"=>vol_solid)
  fsi_Γi_params = Dict("n"=>n_Γi,
  "E"=>E_m,
  "ν"=>ν_m,
  "α"=>α_m,
  "vol"=>vol_Γi)

  # FSI problem
  println("Defining FSI operator")
  res_FSI_Ωf(t,x,xt,y) = WeakForms.fsi_residual_Ωf(strategy,t,x,xt,y,fsi_f_params)
  jac_FSI_Ωf(t,x,xt,dx,y) = WeakForms.fsi_jacobian_Ωf(strategy,x,xt,dx,y,fsi_f_params)
  jac_t_FSI_Ωf(t,x,xt,dxt,y) = WeakForms.fsi_jacobian_t_Ωf(strategy,x,xt,dxt,y,fsi_f_params)
  res_FSI_Ωs(t,x,xt,y) = WeakForms.fsi_residual_Ωs(strategy,t,x,xt,y,fsi_s_params)
  jac_FSI_Ωs(t,x,xt,dx,y) = WeakForms.fsi_jacobian_Ωs(strategy,x,xt,dx,y,fsi_s_params)
  jac_t_FSI_Ωs(t,x,xt,dxt,y) = WeakForms.fsi_jacobian_t_Ωs(strategy,x,xt,dxt,y,fsi_s_params)
  res_FSI_Γi(x,y) = WeakForms.fsi_residual_Γi(strategy,x,y,fsi_Γi_params)
  jac_FSI_Γi(x,dx,y) = WeakForms.fsi_jacobian_Γi(strategy,x,dx,y,fsi_Γi_params)
  t_FSI_Ωf = FETerm(res_FSI_Ωf, jac_FSI_Ωf, jac_t_FSI_Ωf, trian_fluid, quad_fluid)
  t_FSI_Ωs = FETerm(res_FSI_Ωs, jac_FSI_Ωs, jac_t_FSI_Ωs, trian_solid, quad_solid)
  t_FSI_Γi = FETerm(res_FSI_Γi,jac_FSI_Γi,trian_Γi,quad_Γi)
  op_FSI = TransientFEOperator(X_FSI,Y_FSI,t_FSI_Ωf,t_FSI_Ωs,t_FSI_Γi)

  # Setup output files
  folderName = "fsi-results"
  fileName = "fields"
  if !isdir(folderName)
    mkdir(folderName)
  end
  filePath = join([folderName, fileName], "/")

  # Solve Stokes problem
  @timeit "ST problem" begin
  println("Solving Stokes problem")
  xh = solve(op_ST)
  if(is_vtk)
    writePVD(filePath, trian_fluid, [(xh, 0.0)])
  end
end

# Compute Stokes solution L2-norm
l2(w) = w⋅w
eu_ST = u(0.0) - restrict(xh[1],trian_fluid)
ev_ST = v(0.0) - restrict(xh[2],trian_fluid)
ep_ST = p(0.0) - restrict(xh[3],trian_fluid)
eul2_ST = sqrt(sum( integrate(l2(eu_ST),trian_fluid,quad_fluid) ))
evl2_ST = sqrt(sum( integrate(l2(ev_ST),trian_fluid,quad_fluid) ))
epl2_ST = sqrt(sum( integrate(l2(ep_ST),trian_fluid,quad_fluid) ))
println("Stokes L2-norm u: ", eul2_ST)
println("Stokes L2-norm v: ", evl2_ST)
println("Stokes L2-norm p: ", epl2_ST)
@test eul2_ST < 1.0e-10
@test evl2_ST < 1.0e-10
@test epl2_ST < 1.0e-10

# Solve FSI problem
@timeit "FSI problem" begin
println("Solving FSI problem")
xh0 = interpolate_everywhere([u(0.0),v(0.0),p(0.0)],X_FSI(0.0))
nls = NLSolver(
show_trace = true,
method = :newton,
linesearch = BackTracking(),
ftol = 1.0e-10,
iterations = 50
)
odes =  ThetaMethod(nls, dt, 0.5)
solver = TransientFESolver(odes)
sol_FSI = solve(solver, op_FSI, xh0, t0, tf)
end

# Compute outputs
out_params = Dict("trian"=>trian,
"quad"=>quad,
"sol"=>sol_FSI,
"u"=>u,
"v"=>v,
"p"=>p,
"filePath"=>filePath,
"is_vtk"=>is_vtk)
output = computeOutputs(problem,strategy;params=out_params)

end

function get_boundary_conditions(problem::Problem{:analytical},strategy::WeakForms.MeshStrategy,u,v)
  boundary_conditions = (
  # Tags
  FSI_Vu_tags = ["boundary"],
  FSI_Vv_tags = ["boundary"],
  ST_Vu_tags = ["boundary","interface"],
  ST_Vv_tags = ["boundary","interface"],
  # Values,
  FSI_Vu_values = [u],
  FSI_Vv_values = [v],
  ST_Vu_values = [u(0.0),u(0.0)],
  ST_Vv_values = [v(0.0),v(0.0)],
  )
end

function get_FE_spaces(problem::Problem{:analytical},strategy::WeakForms.MeshStrategy,model,model_fluid,order,bconds)
  Vu_FSI = TestFESpace(
      model=model,
      valuetype=VectorValue{2,Float64},
      reffe=:Lagrangian,
      order=order,
      conformity =:H1,
      dirichlet_tags=bconds[:FSI_Vu_tags])
  Vv_FSI = TestFESpace(
      model=model,
      valuetype=VectorValue{2,Float64},
      reffe=:Lagrangian,
      order=order,
      conformity =:H1,
      dirichlet_tags=bconds[:FSI_Vv_tags])
  Vu_ST = TestFESpace(
      model=model_fluid,
      #model=model,
      valuetype=VectorValue{2,Float64},
      reffe=:Lagrangian,
      order=order,
      conformity =:H1,
      dirichlet_tags=bconds[:ST_Vu_tags])
  Vv_ST = TestFESpace(
      model=model_fluid,
      #model=model,
      valuetype=VectorValue{2,Float64},
      reffe=:Lagrangian,
      order=order,
      conformity =:H1,
      dirichlet_tags=bconds[:ST_Vv_tags])
  Q = TestFESpace(
      model=model_fluid,
      #model=model,
      valuetype=Float64,
      order=order-1,
      reffe=:Lagrangian,
      constraint=:zeromean,
      conformity=:C0)

  # Trial FE Spaces
  Uu_ST = TrialFESpace(Vu_ST,bconds[:ST_Vu_values])
  Uv_ST = TrialFESpace(Vv_ST,bconds[:ST_Vv_values])
  Uu_FSI = TransientTrialFESpace(Vu_FSI,bconds[:FSI_Vu_values])
  Uv_FSI = TransientTrialFESpace(Vv_FSI,bconds[:FSI_Vv_values])
  P = TrialFESpace(Q)

  # Multifield FE Spaces
  fe_spaces = (
      Y_ST = MultiFieldFESpace([Vu_ST,Vv_ST,Q]),
      X_ST = MultiFieldFESpace([Uu_ST,Uv_ST,P]),
      Y_FSI = MultiFieldFESpace([Vu_FSI,Vv_FSI,Q]),
      X_FSI = TransientMultiFieldFESpace([Uu_FSI,Uv_FSI,P])
  )
end

function computeOutputs(problem::Problem{:analytical},strategy::WeakForms.MeshStrategy;params=Dict())

  # Unpack parameters
  trian = params["trian"]
  quad = params["quad"]
  sol = params["sol"]
  u = params["u"]
  v = params["v"]
  p = params["p"]
  filePath = params["filePath"]
  is_vtk = params["is_vtk"]
  l2(w) = w⋅w


  ## Initialize arrays
  tpl = Real[]
  eupl = Real[]

  ## Loop over steps
  outfiles = paraview_collection(filePath, append=true) do pvd
    for (i,(xh, t)) in enumerate(sol)
      println("STEP: $i, TIME: $t")
      println("============================")

      # Compute errors
      eu = u(t) - restrict(xh[1],trian)
      ev = v(t) - restrict(xh[2],trian)
      ep = p(t) - restrict(xh[3],trian)
      eul2 = sqrt(sum( integrate(l2(eu),trian,quad )))
      evl2 = sqrt(sum( integrate(l2(ev),trian,quad )))
      epl2 = sqrt(sum( integrate(l2(ep),trian,quad )))

      # Write to PVD
      if(is_vtk)
        uh = restrict(xh[1],trian)
        vh = restrict(xh[2],trian)
        ph = restrict(xh[3],trian)
        pvd[t] = createvtk(
        trian,
        filePath * "_$t.vtu",
        cellfields = ["uh" => uh, "vh" => vh, "ph" => ph, "euh" => eu]
        )
      end

      # Test checks
      println("FSI L2-norm u: ", eul2)
      println("FSI L2-norm v: ", evl2)
      println("FSI L2-norm p: ", epl2)
      @test eul2 < 1.0e-10
      @test evl2 < 1.0e-10
      @test epl2 < 1.0e-10

      push!(tpl, t)
      push!(eupl, eul2)
    end
  end
  return (tpl, eupl)
end
