"""
Executes a elastic flag benchmark proposed by Turek et al in "Proposal for
numerical benchmarking of fluid-structure interaction between an elastic
object and laminar incompressible flow"

The default parameters correspond to FSI2 test case.
FS1 test can be run modifying the default values of the following parameters:
Um=0.2,
Re=20,
rho_s=1.0e-3
FS3 test can be run modifying the default values of the following parameters:
Um=2.0,
Re=200,
rho_s=1.0e-3
E_s=5.6e6

It supports the following mesh motion strategies:
 - Laplacian
 - Linear elasticity
 - Noehookean elasticity (not tested)
 - Biharmonic (Laplacian squared)

The fluid and structure can be strongly coupled, i.e. same variational space for
fluid and solid displacements and velocities, or weakly coupled with different
variational spaces and coupling enforced through Nitche's method.
"""
function execute(problem::Problem{:elasticFlag}; kwargs...)

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

  # Mesh strategy
  strategyName = _get_kwarg(:strategy,kwargs,"laplacian")
  strategy = WeakForms.MeshStrategy{Symbol(strategyName)}()

  # Fluid-Solid coupling type
  couplingName = _get_kwarg(:coupling,kwargs,"strong")
  coupling = WeakForms.Coupling{Symbol(couplingName)}()

  # Define BC functions
  println("Defining Boundary conditions")
  u_in(x, t) = VectorValue(1.5 * Um * x[2] * (H - x[2]) / ((H / 2)^2), 0.0)
  u_noSlip(x, t) = VectorValue(0.0, 0.0)
  u_in(t::Real) = x -> u_in(x, t)
  u_noSlip(t::Real) = x -> u_noSlip(x, t)
  ∂tu_in(t) = x -> VectorValue(0.0, 0.0)
  ∂tu_in(x, t) = ∂tu_in(t)(x)
  T_in = typeof(u_in)
  T_noSlip = typeof(u_noSlip)
  @eval ∂t(::$T_in) = $∂tu_in
  @eval ∂t(::$T_noSlip) = $∂tu_in
  bconds = get_boundary_conditions(problem,strategy,coupling,u_in,u_noSlip)

  # Forcing terms
  f(t) = x -> VectorValue(0.0,0.0)

  # Discrete model
  println("Defining discrete model")
  modelName = _get_kwarg(:model,kwargs,"../models/elasticFlag.json")
  model = DiscreteModelFromFile(modelName)
  model_solid = DiscreteModel(model,"solid")
  model_fluid = DiscreteModel(model,"fluid")
  models = Dict(:Ω => model, :Ωf => model_fluid, :Ωs => model_solid)

  # Triangulations
  println("Defining triangulations")
  Tₕ = get_FSI_triangulations(models,coupling)

  # Quadratures
  println("Defining quadratures")
  order = _get_kwarg(:order,kwargs,2)
  quads = get_FSI_quadratures(Tₕ,order)

  # Test FE Spaces
  println("Defining FE spaces")
  Y_ST, X_ST, Y_FSI, X_FSI = get_FE_spaces(strategy,coupling,models,order,bconds)

  # Stokes problem for initial solution
  println("Defining Stokes operator")
  op_ST = get_Stokes_operator([X_ST,Y_ST],strategy,Tₕ[:Ωf],quads[:Ωf],μ_f,f(0.0))

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
    :fu=>f,
    :fv=>f,
  )
  solid_params = Dict{Symbol,Any}(
    :ρ=>ρ_s,
    :E=>E_s,
    :ν=>ν_s,
    :fu=>f,
    :fv=>f,
  )
  Γi_params = Dict{Symbol,Any}(
    :μ=>μ_f,
    :γ=>γ_f,
    :dt=>dt
  )
  params = (mesh_params, fluid_params, solid_params, Γi_params)

  # FSI problem
  println("Defining FSI operator")
  op_FSI = get_FSI_operator([X_FSI,Y_FSI],coupling,strategy,Tₕ,quads,params)

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
  xh = Gridap.solve(op_ST)
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
xht = Gridap.solve(solver, op_FSI, xh0, t0, tf)
end

# Compute outputs
out_params = Dict{Symbol,Any}(
  :μ=>μ_f,
  :Um=>Um,
  :⌀=>⌀,
  :ρ=>ρ_f,
  :θ=>θ,
  :bdegree=>2*order,
  :filePath=>filePath,
  :is_vtk=>is_vtk,
  )
output = computeOutputs(xh0,xht,coupling,strategy,models,Tₕ,quads,out_params)

end

function get_boundary_conditions(
  problem::Problem{:elasticFlag},
  strategy::WeakForms.MeshStrategy,
  coupling::WeakForms.Coupling{:strong},
  u_in,
  u_noSlip
  )
  boundary_conditions = (
    # Tags
    FSI_Vu_tags = ["inlet", "noSlip", "cylinder","fixed","outlet"],
    FSI_Vv_tags = ["inlet", "noSlip", "cylinder","fixed"],
    ST_Vu_tags = ["inlet", "noSlip", "cylinder","interface","outlet"],
    ST_Vv_tags = ["inlet", "noSlip", "cylinder","interface"],
    # Values,
    FSI_Vu_values = [u_noSlip, u_noSlip, u_noSlip, u_noSlip, u_noSlip],
    FSI_Vv_values = [u_in, u_noSlip, u_noSlip, u_noSlip],
    ST_Vu_values = [u_noSlip(0.0), u_noSlip(0.0), u_noSlip(0.0), u_noSlip(0.0), u_noSlip(0.0)],
    ST_Vv_values = [u_in(0.0), u_noSlip(0.0), u_noSlip(0.0), u_noSlip(0.0)],
  )
end

function get_boundary_conditions(
  problem::Problem{:elasticFlag},
  strategy::WeakForms.MeshStrategy,
  coupling::WeakForms.Coupling{:weak},
  u_in,
  u_noSlip
  )
  boundary_conditions = (
    # Tags
    FSI_Vw_f_tags = ["inlet", "noSlip", "cylinder","outlet"],
    FSI_Vu_f_tags = ["inlet", "noSlip", "cylinder","outlet"],
    FSI_Vv_f_tags = ["inlet", "noSlip", "cylinder"],
    FSI_Vu_s_tags = ["fixed"],
    FSI_Vv_s_tags = ["fixed"],
    ST_Vu_tags = ["inlet", "noSlip", "cylinder","interface","outlet"],
    ST_Vv_tags = ["inlet", "noSlip", "cylinder","interface"],
    # Values,
    FSI_Vw_f_values = [u_noSlip(0.0), u_noSlip(0.0), u_noSlip(0.0), u_noSlip(0.0)],
    FSI_Vu_f_values = [u_noSlip, u_noSlip, u_noSlip, u_noSlip],
    FSI_Vv_f_values = [u_in, u_noSlip, u_noSlip],
    FSI_Vu_s_values = [u_noSlip],
    FSI_Vv_s_values = [u_noSlip],
    ST_Vu_values = [u_noSlip(0.0), u_noSlip(0.0), u_noSlip(0.0), u_noSlip(0.0), u_noSlip(0.0)],
    ST_Vv_values = [u_in(0.0), u_noSlip(0.0), u_noSlip(0.0), u_noSlip(0.0)],
  )
end

function computeOutputs(xh0,xht,coupling::WeakForms.Coupling,strategy::WeakForms.MeshStrategy,models,Tₕ,quads,params)

  # Unpack parameters
  bdegree = params[:bdegree]
  μ = params[:μ]
  Um = params[:Um]
  ⌀ = params[:⌀]
  ρ = params[:ρ]
  θ = params[:θ]
  filePath = params[:filePath]
  is_vtk = params[:is_vtk]
  if( typeof(strategy) == WeakForms.MeshStrategy{:biharmonic} )
    uvpindex = [2,3,4]
  else
    uvpindex = [1,2,3]
  end

  ## Surface triangulation
  trian_Γc = BoundaryTriangulation(models[:Ω], "cylinder")
  quad_Γc = CellQuadrature(trian_Γc, bdegree)
  n_Γc = get_normal_vector(trian_Γc)

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

  # Get interface normal Vector
  n_Γi = get_normal_vector(Tₕ[:Γi])

  ## Initialize arrays
  tpl = Real[]
  FDpl = Real[]
  FLpl = Real[]
  uhn = xh0[uvpindex[1]]
  vhn = xh0[uvpindex[2]]
  phn = xh0[uvpindex[3]]

  ## Loop over steps
  outfiles = paraview_collection(filePath, append=true
  ) do pvd
  for (i,(xh, t)) in enumerate(xht)
    println("STEP: $i, TIME: $t")
    println("============================")

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
    uh_Γi = restrict(uh, Tₕ[:Γi])
    vh_Γi = restrict(vh, Tₕ[:Γi])
    ph_Γi = restrict(ph, Tₕ[:Γi])
    uhn_Γi = restrict(uhn, Tₕ[:Γi])
    vhn_Γi = restrict(vhn, Tₕ[:Γi])
    phn_Γi = restrict(phn, Tₕ[:Γi])
    uhθ_Γi = θ*uh_Γi + (1.0-θ)*uhn_Γi
    vhθ_Γi = θ*vh_Γi + (1.0-θ)*vhn_Γi
    phθ_Γi = θ*ph_Γi + (1.0-θ)*phn_Γi

    # Integrate on the cylinder
    FDc, FLc = sum(integrate(traction_boundary(n_Γc,uhθ_Γc,vhθ_Γc,phθ_Γc), trian_Γc, quad_Γc))

    # Integrate on the interface
    FDi, FLi = sum(integrate(traction_interface(n_Γi,uhθ_Γi,vhθ_Γi,phθ_Γi), Tₕ[:Γi], quads[:Γi]))

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

return (tpl, FDpl, FLpl)
end
