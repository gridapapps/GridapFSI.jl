function execute(problem::Problem{:oscillator}; kwargs...)

  # Define cylinder motion
  u_y, du_y = get_cylinder_motion()

  # Problem setting (Default FSI-2)
  println("Setting Forced Oscillator cylinder problem parameters")
  Um = _get_kwarg(:Um,kwargs,1.0)
  ⌀ = _get_kwarg(:D,kwargs,0.1)
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

  # Mesh strategy
  strategyName = _get_kwarg(:strategy,kwargs,"biharmonic")
  strategy = WeakForms.MeshStrategy{Symbol(strategyName)}()

  # Coupling (boundary conditions) type
  couplingName = _get_kwarg(:coupling,kwargs,"strong")
  coupling = WeakForms.Coupling{Symbol(couplingName)}()

  # Define BC functions
  println("Defining Boundary conditions")
  u_in(x, t) = VectorValue( Um , 0.0)
  u_sym(x, t) = 0.0
  u_cylinder(x,t) = VectorValue(0.0,u_y(t))
  u_in(t::Real) = x -> u_in(x, t)
  u_sym(t::Real) = x -> u_sym(x, t)
  u_cylinder(t::Real) = x -> u_cylinder(x, t)
  ∂tu_in(t) = x -> VectorValue(0.0, 0.0)
  ∂tu_sym(t) = x -> 0.0
  ∂tu_cylinder(t) = x -> VectorValue(0.0, du_y(t))
  ∂tu_in(x, t) = ∂tu_in(t)(x)
  ∂tu_sym(x, t) = ∂tu_sym(t)(x)
  ∂tu_cylinder(x, t) = ∂tu_cylinder(t)(x)
  T_in = typeof(u_in)
  T_sym = typeof(u_sym)
  T_cylinder = typeof(u_cylinder)
  @eval ∂t(::$T_in) = $∂tu_in
  @eval ∂t(::$T_sym) = $∂tu_sym
  @eval ∂t(::$T_cylinder) = $∂tu_cylinder
  bconds = get_boundary_conditions(problem,strategy,coupling,u_in,u_sym,u_cylinder,∂tu_cylinder)

  # Forcing terms
  f(t) = x -> VectorValue(0.0,0.0)

  # Discrete model
  println("Defining discrete model")
  modelName = _get_kwarg(:model,kwargs,"../models/cylinder_NSI.json")
  model = DiscreteModelFromFile(modelName)
  writevtk(model,"cylinder")

  # Triangulations
  println("Defining triangulations")
  trian = Triangulation(model)
  trian_Γc = BoundaryTriangulation(model,"cylinder")
  n_Γc = get_normal_vector(trian_Γc)

  # Compute cell area (auxiliar quantity for mesh motion eq.)
  if( weight_strategy == "volume")
    vol = cell_measure(trian,trian)
    α_Ω = α_m * vol
    α_Γc = α_m * reindex(vol,trian_Γc)
  else
    α_Ω = α_m; α_Γc = α_m
  end

  # Compute interface element size (if weak coupling)
  if ( typeof(coupling) == WeakForms.Coupling{:weak} )
    hΓ = reindex(cell_measure(trian_Γc,trian),trian_Γc)
  else
    hΓ = 0.0
  end

  # Quadratures
  println("Defining quadratures")
  order = _get_kwarg(:order,kwargs,2)
  degree = 2*order
  bdegree = 2*order
  quad = CellQuadrature(trian,degree)
  quad_Γc = CellQuadrature(trian_Γc,bdegree)

  # Test FE Spaces
  println("Defining FE spaces")
  Y_ST, X_ST, Y_NSI, X_NSI = get_FE_spaces(problem,strategy,coupling,model,order,bconds)

  # Stokes problem for initial solution
  println("Defining Stokes operator")
  res_ST(x,y) = WeakForms.stokes_residual(strategy,x,y,μ_f,f(0.0))
  jac_ST(x,dx,y) = WeakForms.stokes_jacobian(strategy,dx,y,μ_f)
  t_ST_Ω = FETerm(res_ST, jac_ST, trian, quad)
  op_ST = FEOperator(X_ST,Y_ST,t_ST_Ω)

  # Setup equation parameters
  nsi_f_params = Dict(
  :μ=>μ_f,
  :ρ=>ρ_f,
  "E"=>E_m,
  "ν"=>ν_m,
  "α"=>α_Ω,
  :fu=>f,
  :fv=>f,
  )
  nsi_Γc_params = Dict(
  "n"=>n_Γc,
  "E"=>E_m,
  "ν"=>ν_m,
  :μ=>μ_f,
  "α"=>α_Γc,
  :γ=>γ_f,
  :h=>hΓ,
  :vD=>du_y
  )

   # FSI problem
  println("Defining FSI operator")
  res_NSI_Ω(t,x,xt,y) = WeakForms.fluid_residual_Ω(strategy,t,x,xt,y,nsi_f_params)
  jac_NSI_Ω(t,x,xt,dx,y) = WeakForms.fluid_jacobian_Ω(strategy,x,xt,dx,y,nsi_f_params)
  jac_t_NSI_Ω(t,x,xt,dxt,y) = WeakForms.fluid_jacobian_t_Ω(strategy,x,xt,dxt,y,nsi_f_params)
  t_NSI_Ω = FETerm(res_NSI_Ω, jac_NSI_Ω, jac_t_NSI_Ω, trian, quad)
  if(typeof(coupling)==WeakForms.Coupling{:weak})
    res_NSI_Γc(t,x,xt,y) = WeakForms.fluid_residual_Γ(strategy,t,x,xt,y,nsi_Γc_params)
    jac_NSI_Γc(t,x,xt,dx,y) = WeakForms.fluid_jacobian_Γ(strategy,x,xt,dx,y,nsi_Γc_params)
    jac_t_NSI_Γc(t,x,xt,dxt,y) = WeakForms.fluid_jacobian_t_Γ(strategy,dxt,y,nsi_Γc_params)
    t_NSI_Γc = FETerm(res_NSI_Γc,jac_NSI_Γc,jac_t_NSI_Γc,trian_Γc,quad_Γc)
    op_FSI = TransientFEOperator(X_NSI,Y_NSI,t_NSI_Ω,t_NSI_Γc)
  else
    op_FSI = TransientFEOperator(X_NSI,Y_NSI,t_NSI_Ω)
  end

  # Setup output files
  folderName = "oscillator-results"
  fileName = "fields"
  if !isdir(folderName)
    mkdir(folderName)
  end
  filePath = join([folderName, fileName], "/")

  # Solve Stokes problem
  @timeit "ST problem" begin
    println("Solving Stokes problem")
    xh = Gridap.solve(op_ST)
    if(is_vtk)
      writePVD(filePath, trian, [(xh, 0.0)])
    end
  end

  # Solve FSI problem
  @timeit "FSI problem" begin
    println("Solving FSI problem")
    xh0  = interpolate(xh,X_NSI(0.0))
    nls = NLSolver(
    show_trace = true,
    method = :newton,
    linesearch = BackTracking(),
    ftol = 1.0e-6,
    iterations = 50
    )
    odes =  ThetaMethod(nls, dt, θ)
    solver = TransientFESolver(odes)
    sol_NSI = Gridap.solve(solver, op_FSI, xh0, t0, tf)
  end

  # Compute outputs
  out_params = Dict(
  :μ=>μ_f,
  "Um"=>Um,
  "⌀"=>⌀,
  :ρ=>ρ_f,
  "θ"=>θ,
  "model"=>model,
  "bdegree"=>bdegree,
  "trian"=>trian,
  "trian_Γc"=>trian_Γc,
  "quad_Γc"=>quad_Γc,
  "n_Γc"=>n_Γc,
  "xh0"=>xh0,
  "sol"=>sol_NSI,
  "filePath"=>filePath,
  "is_vtk"=>is_vtk,
  "coupling"=>coupling,
  )
  output = computeOutputs(problem,strategy;params=out_params)

end

function get_cylinder_motion()
  # Parameters
  A₀ = 1.5
  A₁ = 8.5
  A₂ = -11.1
  A₃ = 2.6
  B₀ = 4.2
  B₁ = 11.3
  B₂ = -68.7
  B₃ = 50.5

  # Physics Properties
  k = 1.0
  m = 1.0
  D = 0.1
  L = 1.0
  ρ = 1.0e3
  b = 1.0
  V = 1.0

  # Parameters
  q₁ = 2
  St = 0.2
  Cx0 = 2.0
  Cy0 = 0.3
  A = 12
  ε = 0.3

  # Auxiliar variables
  mₐ = π*ρ*D^2*L/4
  mʳ = m/mₐ
  ωₙ = sqrt(k/(m+mₐ))
  ωₛ = St*V/(2*π*D)
  Ωₙ = ωₙ/ωₛ
  ζ = b/(2*sqrt((m+mₐ)*k))

  function oscillator(x,p,t)
    y, q, dy, dq = x
    y_dot = dy
    q_dot = dq
    Cvy= (-2^π*St*dy*Cx0 + Cy0/q₁*q) * sqrt(1+4*π^2*St^2*dy^2)
    dy_dot = ρ*D^2*L/(m+mₐ)*1.0/(8*π^2*St^2)*Cvy - 2*ζ*Ωₙ*dy - Ωₙ^2*y
    dq_dot = A*dy_dot - ε*(q^2-1)*dq - q
    dx = [y_dot; q_dot; dy_dot; dq_dot]
  end

  prob = ODEProblem(oscillator,[0.0;2.0;0.0;0.0],(0.0,100.0))
  sol=DifferentialEquations.solve(prob)
  y(t) = sol(t)[1]
  dy(t) = sol(t)[3]

  return y,dy
end

function get_boundary_conditions(
  problem::Problem{:oscillator},
  strategy::WeakForms.MeshStrategy,
  coupling::WeakForms.Coupling{:weak},
  u_in,
  u_sym,
  u_cylinder,
  v_cylinder
  )

  u0 = VectorValue(0.0, 0.0)
  u0t(x,t) = u0
  u0t(t::Real) = x -> u0t(x,t)

  boundary_conditions = (
  # Tags
  NSI_Vw_tags = ["inlet", "noslip", "outlet"],
  NSI_Vu_tags = ["inlet", "noslip", "outlet"],
  NSI_Vv_tags = ["inlet", "noslip"],
  NSI_Vv_masks = [(true,true), (false,true)],
  ST_Vu_tags = ["inlet", "noslip", "cylinder","outlet"],
  ST_Vv_tags = ["inlet", "noslip", "cylinder"],
  ST_Vv_masks = [(true,true), (false,true), (true,true)],
  # Values,
  NSI_Vw_values = [u0, u0, u0],
  NSI_Vu_values = [u0t, u0t, u0t],
  NSI_Vv_values = [u_in, u_sym],
  ST_Vu_values = [u0, u0, u0, u0],
  ST_Vv_values = [u_in(0.0), u_sym(0.0), u_cylinder(0.0)],
  )
end

function get_boundary_conditions(
  problem::Problem{:oscillator},
  strategy::WeakForms.MeshStrategy,
  coupling::WeakForms.Coupling{:strong},
  u_in,
  u_sym,
  u_cylinder,
  v_cylinder
  )

  u0 = VectorValue(0.0, 0.0)
  u0t(x,t) = u0
  u0t(t::Real) = x -> u0t(x,t)

  boundary_conditions = (
  # Tags
  NSI_Vw_tags = ["inlet", "noslip", "cylinder", "outlet"],
  NSI_Vu_tags = ["inlet", "noslip", "cylinder", "outlet"],
  NSI_Vv_tags = ["inlet", "noslip", "cylinder"],
  NSI_Vv_masks = [(true,true), (false,true), (true,true)],
  ST_Vu_tags = ["inlet", "noslip", "cylinder","outlet"],
  ST_Vv_tags = ["inlet", "noslip", "cylinder"],
  ST_Vv_masks = [(true,true), (false,true), (true,true)],
  # Values,
  NSI_Vw_values = [u0, u0, u0, u0],
  NSI_Vu_values = [u0t, u0t, u_cylinder, u0t],
  NSI_Vv_values = [u_in, u_sym, v_cylinder],
  ST_Vu_values = [u0, u0, u_cylinder(0.0), u0],
  ST_Vv_values = [u_in(0.0), u_sym(0.0), v_cylinder(0.0)],
  )
end

function get_FE_spaces(
  problem::Problem{:oscillator},
  strategy::WeakForms.MeshStrategy{:biharmonic},
  coupling::WeakForms.Coupling,
  model,
  order,
  bconds)

  Vw_NSI = TestFESpace(
    model=model,
    valuetype=VectorValue{2,Float64},
    reffe=:Lagrangian,
    order=order,
    conformity =:H1,
    dirichlet_tags=bconds[:NSI_Vw_tags]
    )
  Vu_NSI = TestFESpace(
    model=model,
    valuetype=VectorValue{2,Float64},
    reffe=:Lagrangian,
    order=order,
    conformity =:H1,
    dirichlet_tags=bconds[:NSI_Vu_tags]
    )
  Vv_NSI = TestFESpace(
    model=model,
    valuetype=VectorValue{2,Float64},
    reffe=:Lagrangian,
    order=order,
    conformity =:H1,
    dirichlet_tags=bconds[:NSI_Vv_tags],
    dirichlet_masks=bconds[:NSI_Vv_masks]
    )
  Vw_ST = TestFESpace(
    model=model,
    valuetype=VectorValue{2,Float64},
    reffe=:Lagrangian,
    order=order,
    conformity =:H1,
    dirichlet_tags=bconds[:ST_Vu_tags]
    )
  Vu_ST = TestFESpace(
    model=model,
    valuetype=VectorValue{2,Float64},
    reffe=:Lagrangian,
    order=order,
    conformity =:H1,
    dirichlet_tags=bconds[:ST_Vu_tags]
    )
  Vv_ST = TestFESpace(
    model=model,
    valuetype=VectorValue{2,Float64},
    reffe=:Lagrangian,
    order=order,
    conformity =:H1,
    dirichlet_tags=bconds[:ST_Vv_tags],
    dirichlet_masks=bconds[:ST_Vv_masks]
    )
  Q = TestFESpace(
    model=model,
    valuetype=Float64,
    order=order-1,
    reffe=:Lagrangian,
    conformity=:C0
    )

  # Trial FE Spaces
  Uw_ST = TrialFESpace(Vu_ST,bconds[:ST_Vu_values])
  Uu_ST = TrialFESpace(Vu_ST,bconds[:ST_Vu_values])
  Uv_ST = TrialFESpace(Vv_ST,bconds[:ST_Vv_values])
  Uw_NSI = TrialFESpace(Vw_NSI,bconds[:NSI_Vw_values])
  Uu_NSI = TransientTrialFESpace(Vu_NSI,bconds[:NSI_Vu_values])
  Uv_NSI = TransientTrialFESpace(Vv_NSI,bconds[:NSI_Vv_values])
  P = TrialFESpace(Q)

  # Multifield FE Spaces
  fe_spaces = (
    Y_ST = MultiFieldFESpace([Vw_ST,Vu_ST,Vv_ST,Q]),
    X_ST = MultiFieldFESpace([Uw_ST,Uu_ST,Uv_ST,P]),
    Y_NSI = MultiFieldFESpace([Vw_NSI,Vu_NSI,Vv_NSI,Q]),
    X_NSI = TransientMultiFieldFESpace([Uw_NSI,Uu_NSI,Uv_NSI,P])
  )
end
function computeOutputs(problem::Problem{:oscillator},strategy::WeakForms.MeshStrategy;params=Dict())#, sol, xh0)

  # Unpack parameters
  model = params["model"]
  bdegree = params["bdegree"]
  xh0 = params["xh0"]
  sol = params["sol"]
  μ = params[:μ]
  trian = params["trian"]
  trian_Γc = params["trian_Γc"]
  quad_Γc = params["quad_Γc"]
  n_Γc = params["n_Γc"]
  Um = params["Um"]
  ⌀ = params["⌀"]
  ρ = params[:ρ]
  θ = params["θ"]
  filePath = params["filePath"]
  is_vtk = params["is_vtk"]
  coupling = params["coupling"]
  if( typeof(strategy) == WeakForms.MeshStrategy{:biharmonic} )
    uvpindex = [2,3,4]
  else
    uvpindex = [1,2,3]
  end

  #=
  # Aux function
  traction_boundary(n,u,v,p) = n ⋅ WeakForms.Pf_dev(μ,u,v)  + WeakForms.Pf_vol(u,p) * n
  function traction_closure(coupling,n,u,v,p)
    if typeof(coupling) == WeakForms.Coupling{:weak}
      traction_boundary(n,u,v,p).⁺
    else
      traction_boundary(n,u,v,p)
    end
  end
  traction_interface(n,u,v,p) = traction_closure(coupling,n,u,v,p)

  ## Initialize arrays
  tpl = Real[]
  FDpl = Real[]
  FLpl = Real[]
  uhn = xh0[uvpindex[1]]
  vhn = xh0[uvpindex[2]]
  phn = xh0[uvpindex[3]]
=#
  ## Loop over steps
  outfiles = paraview_collection(filePath, append=true
  ) do pvd
  for (i,(xh, t)) in enumerate(sol)
    println("STEP: $i, TIME: $t")
    println("============================")

    #=
    ## Get the solution at n+θ (where velocity and pressure are balanced)
    uh = xh[uvpindex[1]]
    vh = xh[uvpindex[2]]
    ph = xh[uvpindex[3]]
    uh_Γc = restrict(uh, trian_Γc)
    vh_Γc = restrict(vh, trian_Γc)
    ph_Γc = restrict(ph, trian_Γc)
    uhn_Γc = restrict(uhn, trian_Γc)
    vhn_Γc = restrict(vhn, trian_Γc)
    phn_Γc = restrict(phn, trian_Γc)
    uhθ_Γc = θ*uh_Γc + (1.0-θ)*uhn_Γc
    vhθ_Γc = θ*vh_Γc + (1.0-θ)*vhn_Γc
    phθ_Γc = θ*ph_Γc + (1.0-θ)*phn_Γc
    uh_Γi = restrict(uh, trian_Γi)
    vh_Γi = restrict(vh, trian_Γi)
    ph_Γi = restrict(ph, trian_Γi)
    uhn_Γi = restrict(uhn, trian_Γi)
    vhn_Γi = restrict(vhn, trian_Γi)
    phn_Γi = restrict(phn, trian_Γi)
    uhθ_Γi = θ*uh_Γi + (1.0-θ)*uhn_Γi
    vhθ_Γi = θ*vh_Γi + (1.0-θ)*vhn_Γi
    phθ_Γi = θ*ph_Γi + (1.0-θ)*phn_Γi

    # Integrate on the cylinder
    FDc, FLc = sum(integrate(traction_boundary(n_Γc,uhθ_Γc,vhθ_Γc,phθ_Γc), trian_Γc, quad_Γc))

    # Integrate on the interface
    FDi, FLi = sum(integrate(traction_interface(n_Γi,uhθ_Γi,vhθ_Γi,phθ_Γi), trian_Γi, quad_Γi))

    FD = FDc + FDi
    FL = FLc + FLi

    ## Drag and lift Forces
    push!(tpl, t)
    push!(FDpl, -FD)
    push!(FLpl, -FL)

    ## store step n
    uhn = uh
    vhn = vh
    phn = ph
=#
    # Write to PVD
    if(is_vtk)
      uh = xh[uvpindex[1]]
      vh = xh[uvpindex[2]]
      ph = xh[uvpindex[3]]
      pvd[t] = createvtk(
      trian,
      filePath * "_$t.vtu",
      cellfields = ["uh" => uh, "vh" => vh, "ph" => ph]
      )
    end
  end
end

#return (tpl, FDpl, FLpl)
end
