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
  couplingName = _get_kwarg(:coupling,kwargs,"weak")
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
  bconds = get_boundary_conditions(problem,strategy,coupling,u_in,u_sym,u_cylinder)

  # Forcing terms
  f(t) = x -> VectorValue(0.0,0.0)

  #= # Discrete model
  println("Defining discrete model")
  modelName = _get_kwarg(:model,kwargs,"../models/elasticFlag.json")
  model = DiscreteModelFromFile(modelName)
  model_solid = DiscreteModel(model,"solid")
  model_fluid = DiscreteModel(model,"fluid")

  # Triangulations
  println("Defining triangulations")
  trian = Triangulation(model)
  trian_solid = Triangulation(model_solid)
  trian_fluid = Triangulation(model_fluid)
  function Γi_triangulation(coupling)
    if typeof(coupling) == WeakForms.Coupling{:weak}
      InterfaceTriangulation(model_fluid,model_solid)
    else
      BoundaryTriangulation(model_fluid,"interface")
    end
  end
  trian_Γi = Γi_triangulation(coupling)
  n_Γi = get_normal_vector(trian_Γi)

  # Compute cell area (auxiliar quantity for mesh motion eq.)
  if( weight_strategy == "volume")
    volf = cell_measure(trian_fluid,trian)
    vols = cell_measure(trian_solid,trian)
    α_fluid = α_m * reindex(volf,trian_fluid)
    α_solid = α_m * reindex(vols,trian_solid)
    if ( typeof(coupling) == WeakForms.Coupling{:weak} )
      α_Γi = α_m * reindex(volf,get_left_boundary(trian_Γi))
    else
      α_Γi = α_m * reindex(volf,trian_Γi)
    end
  else
    α_fluid = α_m; α_solid = α_m; α_Γi = α_m
  end

  # Compute interface element size (if weak coupling)
  if ( typeof(coupling) == WeakForms.Coupling{:weak} )
    trian_boundary_Γi = get_left_boundary(trian_Γi)
    hΓᵢ = reindex(cell_measure(trian_boundary_Γi,trian),trian_boundary_Γi)
  else
    hΓᵢ = 0.0
  end

  # Quadratures
  println("Defining quadratures")
  order = _get_kwarg(:order,kwargs,2)
  degree = 2*order
  bdegree = 2*order
  quad_solid = CellQuadrature(trian_solid,degree)
  quad_fluid = CellQuadrature(trian_fluid,degree)
  quad_Γi = CellQuadrature(trian_Γi,bdegree)

  # Test FE Spaces
  println("Defining FE spaces")
  Y_ST, X_ST, Y_FSI, X_FSI = get_FE_spaces(problem,strategy,coupling,model,
  model_fluid,model_solid,order,bconds)

  # Stokes problem for initial solution
  println("Defining Stokes operator")
  res_ST(x,y) = WeakForms.stokes_residual(strategy,x,y,μ_f,f(0.0))
  jac_ST(x,dx,y) = WeakForms.stokes_jacobian(strategy,dx,y,μ_f)
  t_ST_Ωf = FETerm(res_ST, jac_ST, trian_fluid, quad_fluid)
  op_ST = FEOperator(X_ST,Y_ST,t_ST_Ωf)

  # Setup equation parameters
  fsi_f_params = Dict(
  "μ"=>μ_f,
  "ρ"=>ρ_f,
  "E"=>E_m,
  "ν"=>ν_m,
  "α"=>α_fluid,
  "fu"=>f,
  "fv"=>f,
  )
  fsi_s_params = Dict(
  "ρ"=>ρ_s,
  "E"=>E_s,
  "ν"=>ν_s,
  "fu"=>f,
  "fv"=>f,
  "α"=>α_solid,
  )
  fsi_Γi_params = Dict(
  "n"=>n_Γi,
  "E"=>E_m,
  "ν"=>ν_m,
  "μ"=>μ_f,
  "α"=>α_Γi,
  "γ"=>γ_f,
  "h"=>hΓᵢ,
  "dt"=>dt
  )

  # FSI problem
  println("Defining FSI operator")
  res_FSI_Ωf(t,x,xt,y) = WeakForms.fsi_residual_Ωf(strategy,coupling,t,x,xt,y,fsi_f_params)
  jac_FSI_Ωf(t,x,xt,dx,y) = WeakForms.fsi_jacobian_Ωf(strategy,coupling,x,xt,dx,y,fsi_f_params)
  jac_t_FSI_Ωf(t,x,xt,dxt,y) = WeakForms.fsi_jacobian_t_Ωf(strategy,coupling,x,xt,dxt,y,fsi_f_params)
  res_FSI_Ωs(t,x,xt,y) = WeakForms.fsi_residual_Ωs(strategy,coupling,t,x,xt,y,fsi_s_params)
  jac_FSI_Ωs(t,x,xt,dx,y) = WeakForms.fsi_jacobian_Ωs(strategy,coupling,x,xt,dx,y,fsi_s_params)
  jac_t_FSI_Ωs(t,x,xt,dxt,y) = WeakForms.fsi_jacobian_t_Ωs(strategy,coupling,x,xt,dxt,y,fsi_s_params)
  res_FSI_Γi(x,y) = WeakForms.fsi_residual_Γi(strategy,coupling,x,y,fsi_Γi_params)
  jac_FSI_Γi(x,dx,y) = WeakForms.fsi_jacobian_Γi(strategy,coupling,x,dx,y,fsi_Γi_params)
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

  # Solve FSI problem
  @timeit "FSI problem" begin
    println("Solving FSI problem")
    xh0  = interpolate(xh,X_FSI(0.0))
    nls = NLSolver(
    show_trace = true,
    method = :newton,
    linesearch = BackTracking(),
    ftol = 1.0e-6,
    iterations = 50
    )
    odes =  ThetaMethod(nls, dt, θ)
    solver = TransientFESolver(odes)
    sol_FSI = solve(solver, op_FSI, xh0, t0, tf)
  end

  # Compute outputs
  out_params = Dict(
  "μ"=>μ_f,
  "Um"=>Um,
  "⌀"=>⌀,
  "ρ"=>ρ_f,
  "θ"=>θ,
  "model"=>model,
  "bdegree"=>bdegree,
  "trian"=>trian,
  "trian_Γi"=>trian_Γi,
  "quad_Γi"=>quad_Γi,
  "n_Γi"=>n_Γi,
  "xh0"=>xh0,
  "sol"=>sol_FSI,
  "filePath"=>filePath,
  "is_vtk"=>is_vtk,
  "coupling"=>coupling,
  )
  output = computeOutputs(problem,strategy;params=out_params)
   =#
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
  u_cylinder
  )

  u0 = VectorValue(0.0, 0.0)
  u0t(t) = u0

  boundary_conditions = (
  # Tags
  NSI_Vw_f_tags = ["inlet", "noSlip", "outlet"],
  NSI_Vu_f_tags = ["inlet", "noSlip", "outlet"],
  NSI_Vv_f_tags = ["inlet", "noSlip"],
  ST_Vu_tags = ["inlet", "noSlip", "cylinder","outlet"],
  ST_Vv_tags = ["inlet", "noSlip", "cylinder",],
  # Values,
  NSI_Vw_f_values = [u0, u0, u0],
  NSI_Vu_f_values = [u0t, u0t, u0t],
  NSI_Vv_f_values = [u_in, u_sym],
  ST_Vu_values = [u0, u0, u0, u0],
  ST_Vv_values = [u_in(0.0), u_sym(0.0), u_cylinder(0.0)],
  )
end
