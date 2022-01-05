"""
Executes a transient FSI driver with a constant function in space, linear in time:
  u(x,t) = [1.0, -1.0]^T * t
  v(x,t) = [1.0, -1.0]^T
  p(x,t) = 0.0
"""
function execute(problem::FSIProblem{:analytical};kwargs...)

  # Problem setting
  println("Setting analytical fluid problem parameters")
  # Solid properties
  E_s = _get_kwarg(:E_s,kwargs,1.0)
  ν_s = _get_kwarg(:nu_s,kwargs,0.4)
  ρ_s = _get_kwarg(:rho_s,kwargs,1.0)
  # Fluid properties
  ρ_f = _get_kwarg(:rho_f,kwargs,1.0)
  μ_f = _get_kwarg(:mu_f,kwargs,1.0)
  γ_f = _get_kwarg(:gamma_f,kwargs,1.0)
  # Mesh properties
  n_m = _get_kwarg(:n_m,kwargs,10)
  E_m = _get_kwarg(:E_m,kwargs,1.0)
  ν_m = _get_kwarg(:nu_m,kwargs,-0.1)
  α_m = _get_kwarg(:alpha_m,kwargs,1.0e-5)
  weight_strategy = _get_kwarg(:alpha_m_weight,kwargs,"constant")
  # Time stepping
  t0 = _get_kwarg(:t0,kwargs,0.0)
  tf = _get_kwarg(:tf,kwargs,0.5)
  dt = _get_kwarg(:dt,kwargs,0.1)
  θ  = _get_kwarg(:theta,kwargs,0.5)
  # Post-process
  is_vtk = _get_kwarg(:is_vtk,kwargs,false)

  # Mesh strategy
  strategyName = _get_kwarg(:strategy,kwargs,"laplacian")
  strategy = WeakForms.MeshStrategy{Symbol(strategyName)}()

  # Fluid-Solid coupling type
  couplingName = _get_kwarg(:coupling,kwargs,"strong")
  coupling = WeakForms.Coupling{Symbol(couplingName)}()

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
  (λ_s, μ_s) = lame_parameters(E_s,ν_s)
  (λ_m, μ_m) = lame_parameters(E_m,ν_m)
  S_SV(t) = x -> 2*μ_s*E(t)(x) + λ_s*tr(E(t)(x))*I
  fv_ST_Ωf(t) = x -> - μ_f*Δ(v(t))(x) + ∇(p(t))(x)
  function fu_closure(t,strategy)
    fu = if typeof(strategy) == MeshStrategy{:laplacian}
      x -> - α_m * Δ(u(t))(x)
    elseif typeof(strategy) == MeshStrategy{:linearElasticity}
      x -> - μ_m * Δ(u(t))(x)
    elseif typeof(strategy) == MeshStrategy{:neoHookean}
      x -> - μ_m * Δ(u(t))(x)
    else
      @notimplemented("The soruce term for $strategy strategy is not implemented")
    end
    return fu
  end
  fu_Ωf(t) = fu_closure(t,strategy)
  fv_Ωf(t) = x -> ρ_f * ∂t(v)(t)(x) - μ_f * Δ(v(t))(x) + ∇(p(t))(x) + ρ_f*( (∇(v(t))(x)')⋅(v(t)(x) - ∂t(u)(t)(x)) )
  fp_Ωf(t) = x -> (∇⋅v(t))(x)
  fu_Ωs(t) = x -> ∂t(u)(t)(x) - v(t)(x)
  fv_Ωs(t) = x -> ρ_s * ∂t(v)(t)(x) #- (∇⋅(F(t)⋅S_SV(t)))(x)  # Divergence of a a doted function not working yet...

  # Discrete model
  println("Defining discrete model")
  domain = (-1,1,-1,1)
  partition = (n_m,n_m)
  𝒯 = CartesianDiscreteModel(domain,partition)
  Ω = Interior(𝒯)

  # Solid domain
  R = 0.5
  xc = 0.0
  yc = 0.0
  function is_in(coords)
    n = length(coords)
    x = (1/n)*sum(coords)
    d = (x[1]-xc)^2 + (x[2]-yc)^2 - R^2
    d < 1.0e-8
  end
  oldcell_to_coods = get_cell_coordinates(Ω)
  oldcell_to_is_in = collect1d(lazy_map(is_in,oldcell_to_coods))
  incell_to_cell = findall(oldcell_to_is_in)
  outcell_to_cell = findall(collect(Bool, .! oldcell_to_is_in))

  # Triangulations
  println("Defining triangulations")
  Ωs = Interior(Ω,incell_to_cell)
  Ωf = Interior(Ω,outcell_to_cell)
  Γi = InterfaceTriangulation(Ωf,Ωs)
  Γi_mask = Int64.(Γi.⁺.glue.face_to_bgface)
  if typeof(coupling) == Coupling{:strong}
    Γi = Γi.⁺
  end
  Tₕ = Dict(:Ω => Ω, :Ωf => Ωf, :Ωs => Ωs, :Γi=>Γi)

  # Add interface tag
  labels = get_face_labeling(𝒯)
  new_entity = num_entities(labels) + 1
  D = num_cell_dims(𝒯)
  for face in Γi_mask
    labels.d_to_dface_to_entity[D][face] = new_entity
  end
  add_tag!(labels,"interface",[new_entity])

  # Quadratures
  println("Defining quadratures")
  order = _get_kwarg(:order,kwargs,2)
  dTₕ = get_FSI_measures(Tₕ,order)

  # Test FE Spaces
  println("Defining FE spaces")
  Y_ST, X_ST, Y_FSI, X_FSI = get_FE_spaces(strategy,coupling,Tₕ,order,bconds,constraint=:zeromean)

  # Stokes problem for initial solution
  println("Defining Stokes operator")
  op_ST = get_Stokes_operator(X_ST,Y_ST,strategy,dTₕ[:Ωf],μ_f,fv_ST_Ωf(0.0))

  # Setup equation parameters
  mesh_params = Dict{Symbol,Any}(
    :w_strategy=>weight_strategy,
    :α=>α_m,
    :E=>E_m,
    :ν=>ν_m,
  )
  fluid_params = Dict{Symbol,Any}(
    :μ=>μ_f,
    :ρ=>ρ_f,
    :fu=>fu_Ωf,
    :fv=>fv_Ωf,
  )
  solid_params = Dict{Symbol,Any}(
    :ρ=>ρ_s,
    :E=>E_s,
    :ν=>ν_s,
    :fu=>fu_Ωs,
    :fv=>fv_Ωs,
  )
  Γi_params = Dict{Symbol,Any}(
    :μ=>μ_f,
    :γ=>γ_f,
    :dt=>dt
  )
  params = (mesh_params, fluid_params, solid_params, Γi_params)

  # FSI problem
  println("Defining FSI operator")
  op_FSI = get_FSI_operator(X_FSI,Y_FSI,coupling,strategy,Tₕ,dTₕ,params)

  # Setup output files
  folderName = "fsi-results"
  fileName = "fields"
  if !isdir(folderName)
    mkdir(folderName)
  end
  filePath = join([folderName, fileName], "/")

  # Solve Stokes problem
  @timeit "ST problem" begin
    println("Defining Stokes solver")
    xh = solve(op_ST)
    if(is_vtk)
      writePVD(filePath, Tₕ[:Ωf], [(xh, 0.0)])
    end
  end

  # Compute Stokes solution L2-norm
  l2(w) = w⋅w
  eu_ST = u(0.0) - xh[1]
  ev_ST = v(0.0) - xh[2]
  ep_ST = p(0.0) - xh[3]
  eul2_ST = sqrt(∑( ∫(l2(eu_ST))dTₕ[:Ωf] ))
  evl2_ST = sqrt(∑( ∫(l2(ev_ST))dTₕ[:Ωf] ))
  epl2_ST = sqrt(∑( ∫(l2(ep_ST))dTₕ[:Ωf] ))
  println("Stokes L2-norm u: ", eul2_ST)
  println("Stokes L2-norm v: ", evl2_ST)
  println("Stokes L2-norm p: ", epl2_ST)
  @test eul2_ST < 1.0e-10
  @test evl2_ST < 1.0e-10
  @test epl2_ST < 1.0e-10

  # Solve FSI problem
  @timeit "FSI problem" begin
  println("Defining FSI solver")
  xh0 = interpolate_everywhere([u(0.0),v(0.0),p(0.0)],X_FSI(0.0))
  nls = NLSolver(
  show_trace = true,
  method = :newton,
  linesearch = BackTracking(),
  ftol = 1.0e-10,
  iterations = 50
  )
  ode_solver =  ThetaMethod(nls, dt, 0.5)
  xht = solve(ode_solver, op_FSI, xh0, t0, tf)
  end

  # Compute outputs
  out_params = Dict(
  :u=>u,
  :v=>v,
  :p=>p,
  :filePath=>filePath,
  :is_vtk=>is_vtk
  )
  output = computeOutputs(xht,Tₕ,dTₕ,strategy,out_params)

end

function get_boundary_conditions(problem::FSIProblem{:analytical},strategy::WeakForms.MeshStrategy,u,v)
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

function get_FE_spaces(problem::FSIProblem{:analytical},strategy::WeakForms.MeshStrategy,model,model_fluid,order,bconds)
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

function computeOutputs(xht,Tₕ,dTₕ,strategy,params)

  # Unpack parameters
  u = params[:u]
  v = params[:v]
  p = params[:p]
  filePath = params[:filePath]
  is_vtk = params[:is_vtk]
  l2(w) = w⋅w

  ## Initialize arrays
  tpl = Real[]
  eupl = Real[]

  ## Loop over steps
  outfiles = paraview_collection(filePath, append=true) do pvd
    for (i,(xh, t)) in enumerate(xht)
      println("STEP: $i, TIME: $t")
      println("============================")

      # Compute errors
      eu = u(t) - xh[1]
      ev = v(t) - xh[2]
      ep = p(t) - xh[3]
      eul2 = sqrt(∑( ∫(l2(eu))dTₕ[:Ω] ))
      evl2 = sqrt(∑( ∫(l2(ev))dTₕ[:Ω] ))
      epl2 = sqrt(∑( ∫(l2(ep))dTₕ[:Ω] ))

      # Write to PVD
      if(is_vtk)
        uh = xh[1]
        vh = xh[2]
        ph = xh[3]
        pvd[t] = createvtk(
        Tₕ[:Ω],
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
