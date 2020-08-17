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

    # Mesh strategy
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
    bconds = get_boundary_conditions(problem,strategy,u_in,u_noSlip)

    # Forcing terms
    f(t) = x -> VectorValue(0.0,0.0)

    # Discrete model
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
    trian_Γi = BoundaryTriangulation(model_fluid,"interface")
    n_Γi = get_normal_vector(trian_Γi)

    # Compute cell area (auxiliar quantity for mesh motion eq.)
    if( weight_strategy == "volume")
      volf = cell_measure(trian_fluid,trian)
      vols = cell_measure(trian_solid,trian)
      α_fluid = α_m * reindex(volf,trian_fluid)
      α_solid = α_m * reindex(vols,trian_solid)
      α_Γi = α_m * reindex(volf,trian_Γi)
    else
      α_fluid = α_m; α_solid = α_m; α_Γi = α_m
    end

    # Compute interface element size (if weak coupling)
    if ( typeof(coupling) == WeakForms.Coupling{:weak} )
      h_Γi = cell_measure(trian_Γi,trian_Γi)
    else
      h_Γi = 0.0
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
    Y_ST, X_ST, Y_FSI, X_FSI = get_FE_spaces(problem,strategy,model,model_fluid,order,bconds)

    # Stokes problem for initial solution
    println("Defining Stokes operator")
    res_ST(x,y) = WeakForms.stokes_residual(strategy,x,y,μ_f,f(0.0))
    jac_ST(x,dx,y) = WeakForms.stokes_jacobian(strategy,dx,y,μ_f)
    t_ST_Ωf = FETerm(res_ST, jac_ST, trian_fluid, quad_fluid)
    op_ST = FEOperator(X_ST,Y_ST,t_ST_Ωf)

    # Setup equation parameters
    fsi_f_params = Dict("μ"=>μ_f,
                        "ρ"=>ρ_f,
                        "E"=>E_m,
                        "ν"=>ν_m,
                        "α"=>α_fluid,
                        "fu"=>f,
                        "fv"=>f)
    fsi_s_params = Dict("ρ"=>ρ_s,
                        "E"=>E_s,
                        "ν"=>ν_s,
                        "fu"=>f,
                        "fv"=>f,
                        "α"=>α_solid)
    fsi_Γi_params = Dict("n"=>n_Γi,
                         "E"=>E_m,
                         "ν"=>ν_m,
                         "α"=>α_Γi,
                         "gamma"=>γ_f,
                         "h"=>h_Γi)

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

    # Solve FSI problem
    @timeit "FSI problem" begin
        println("Solving FSI problem")
		    xh0  = interpolate(X_FSI(0.0),xh)
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
    out_params = Dict("μ"=>μ_f,
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
                      "is_vtk"=>is_vtk)
    output = computeOutputs(problem,strategy;params=out_params)

end

function get_boundary_conditions(problem::Problem{:elasticFlag},strategy::WeakForms.MeshStrategy,u_in,u_noSlip)
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

function get_FE_spaces(problem::Problem,strategy::WeakForms.MeshStrategy,model,model_fluid,order,bconds)
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

function get_FE_spaces(problem::Problem,strategy::WeakForms.MeshStrategy{:biharmonic},model,model_fluid,order,bconds)
    Vw_FSI = TestFESpace(
        model=model,
        valuetype=VectorValue{2,Float64},
        reffe=:Lagrangian,
        order=order,
        conformity =:H1,
        dirichlet_tags=bconds[:FSI_Vu_tags])
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
    Vw_ST = TestFESpace(
        model=model_fluid,
		    #model=model,
        valuetype=VectorValue{2,Float64},
        reffe=:Lagrangian,
        order=order,
        conformity =:H1,
        dirichlet_tags=bconds[:ST_Vu_tags])
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
		    conformity=:C0)

    # Trial FE Spaces
    Uw_ST = TrialFESpace(Vu_ST,bconds[:ST_Vu_values])
    Uu_ST = TrialFESpace(Vu_ST,bconds[:ST_Vu_values])
    Uv_ST = TrialFESpace(Vv_ST,bconds[:ST_Vv_values])
    Uw_FSI = TrialFESpace(Vu_FSI,bconds[:ST_Vu_values])
    Uu_FSI = TransientTrialFESpace(Vu_FSI,bconds[:FSI_Vu_values])
    Uv_FSI = TransientTrialFESpace(Vv_FSI,bconds[:FSI_Vv_values])
    P = TrialFESpace(Q)

    # Multifield FE Spaces
    fe_spaces = (
        Y_ST = MultiFieldFESpace([Vw_ST,Vu_ST,Vv_ST,Q]),
        X_ST = MultiFieldFESpace([Uw_ST,Uu_ST,Uv_ST,P]),
        Y_FSI = MultiFieldFESpace([Vw_FSI,Vu_FSI,Vv_FSI,Q]),
        X_FSI = TransientMultiFieldFESpace([Uw_FSI,Uu_FSI,Uv_FSI,P])
    )
end

function computeOutputs(problem::Problem{:elasticFlag},strategy::WeakForms.MeshStrategy;params=Dict())#, sol, xh0)

    # Unpack parameters
    model = params["model"]
    bdegree = params["bdegree"]
    xh0 = params["xh0"]
    sol = params["sol"]
    μ = params["μ"]
    trian = params["trian"]
    trian_Γi = params["trian_Γi"]
    quad_Γi = params["quad_Γi"]
    n_Γi = params["n_Γi"]
    Um = params["Um"]
    ⌀ = params["⌀"]
    ρ = params["ρ"]
    θ = params["θ"]
    filePath = params["filePath"]
    is_vtk = params["is_vtk"]
    if( typeof(strategy) == WeakForms.MeshStrategy{:biharmonic} )
      uvpindex = [2,3,4]
    else
      uvpindex = [1,2,3]
    end
    ## Surface triangulation
    trian_Γc = BoundaryTriangulation(model, "cylinder")
    quad_Γc = CellQuadrature(trian_Γc, bdegree)
    n_Γc = get_normal_vector(trian_Γc)

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
        for (i,(xh, t)) in enumerate(sol)
            println("STEP: $i, TIME: $t")
            println("============================")

            ## Get the solution at n+θ (where velocity and pressure are balanced)
            uh = xh.blocks[uvpindex[1]]
            vh = xh.blocks[uvpindex[2]]
            ph = xh.blocks[uvpindex[3]]
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
            FDc, FLc = sum(integrate(
              (n_Γc ⋅ WeakForms.Pf_dev(μ,uhθ_Γc,vhθ_Γc)  + WeakForms.Pf_vol(uhθ_Γc,phθ_Γc) * n_Γc),
                trian_Γc,
                quad_Γc,
            ))
            # Integrate on the interface
            FDi, FLi = sum(integrate(
                (n_Γi ⋅ WeakForms.Pf_dev(μ,uhθ_Γi,vhθ_Γi)  + WeakForms.Pf_vol(uhθ_Γi,phθ_Γi) * n_Γi),
                trian_Γi,
                quad_Γi,
            ))
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
              uh = xh.blocks[uvpindex[1]]
              vh = xh.blocks[uvpindex[2]]
		  		    ph = xh.blocks[uvpindex[3]]
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
