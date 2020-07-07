function execute(problem::Problem{:elasticFlag}; kwargs...)


    # Problem setting (Default FSI-2)
    println("Setting Elastic flag problem parameters")
    Um = _get_kwarg(:Um,kwargs,1.0)
    H = _get_kwarg(:H,kwargs,0.41)
    ⌀ = _get_kwarg(:D,kwargs,0.1)
    Re = _get_kwarg(:Re,kwargs, 100.0)
    # Solid properties
    E_s = _get_kwarg(:E_s,kwargs,1.4e6)
    ν_s = _get_kwarg(:nu_s,kwargs,0.4)
    ρ_s = _get_kwarg(:rho_s,kwargs,1.0e4)
    # Fluid properties
    ρ_f = _get_kwarg(:rho_f,kwargs,1.0e3)
    μ_f = ρ_f * Um * ⌀ / Re
    # Mesh properties
    E_m = _get_kwarg(:E_m,kwargs,1.0)
    ν_m = _get_kwarg(:nu_m,kwargs,-0.1)
    # Time stepping
    t0 = _get_kwarg(:t0,kwargs,0.0)
    tf = _get_kwarg(:tf,kwargs,0.1)
    dt = _get_kwarg(:dt,kwargs,0.1)

    # Mesh strategy
    strategyName = _get_kwarg(:strategy,kwargs,"linearElasticity")
    strategy = WeakForms.MeshStrategy{Symbol(strategyName)}()

    # Define BC functions
    println("Defining Boundary conditions")

    #u_in(x, t) = VectorValue(1.5 * Um * x[2] * (H - x[2]) / ((H / 2)^2), 0.0)
    #u_noSlip(x, t) = VectorValue(0.0, 0.0)
    #ut_in(t::Real) = x -> u_in(x, t)
    #ut_noSlip(t::Real) = x -> u_noSlip(x, t)
    #∂tu_in(t) = x -> VectorValue(0.0, 0.0)
    #∂tu_in(x, t) = ∂tu_in(t)(x)
    #∂t(::typeof(ut_in)) = ∂tu_in
    #∂t(::typeof(ut_noSlip)) = ∂tu_in
    bconds = get_boundary_conditions(strategy,ut_in,ut_noSlip)

    # Discrete model
    println("Defining discrete model")
    modelName = _get_kwarg(:model,kwargs,"../models/elasticFlag_coarse.json")
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
    Y_ST, X_ST, Y_FSI, X_FSI = get_FE_spaces(strategy,model,model_fluid,order,bconds)

    # Stokes problem for initial solution
    println("Defining Stokes operator")
    res_ST(x,y) = WeakForms.stokes_uvp_residual(strategy,x,y,μ_f)
    jac_ST(x,dx,y) = WeakForms.stokes_uvp_residual(strategy,dx,y,μ_f)
    t_ST_Ωf = FETerm(res_ST, jac_ST, trian_fluid, quad_fluid)
    op_ST = FEOperator(X_ST,Y_ST,t_ST_Ωf)

    # FSI problem
    println("Defining FSI operator")
    res_FSI_Ωf(t,x,xt,y) = WeakForms.fsi_uvp_residual_Ωf(strategy,x,xt,y,μ_f,ρ_f,E_m,ν_m)
    jac_FSI_Ωf(t,x,xt,dx,y) = WeakForms.fsi_uvp_jacobian_Ωf(strategy,x,xt,dx,y,μ_f,ρ_f,E_m,ν_m)
    jac_t_FSI_Ωf(t,x,xt,dxt,y) = WeakForms.fsi_uvp_jacobian_t_Ωf(strategy,x,xt,dxt,y,ρ_f)
    res_FSI_Ωs(t,x,xt,y) = WeakForms.fsi_uvp_residual_Ωs(strategy,x,xt,y,ρ_s,E_s,ν_s)
    jac_FSI_Ωs(t,x,xt,dx,y) = WeakForms.fsi_uvp_jacobian_Ωs(strategy,x,xt,dx,y,ρ_s,E_s,ν_s)
    jac_t_FSI_Ωs(t,x,xt,dxt,y) = WeakForms.fsi_uvp_jacobian_t_Ωs(strategy,x,xt,dxt,y,ρ_s)
    res_FSI_Γi(x,y) = WeakForms.fsi_uvp_residual_Γi(strategy,x,y,n_Γi,E_m,ν_m)
    jac_FSI_Γi(x,dx,y) = WeakForms.fsi_uvp_jacobian_Γi(strategy,x,dx,y,n_Γi,E_m,ν_m)
    t_FSI_Ωf = FETerm(res_FSI_Ωf, jac_FSI_Ωf, jac_t_FSI_Ωf, trian_fluid, quad_fluid)
    t_FSI_Ωs = FETerm(res_FSI_Ωs, jac_FSI_Ωs, jac_t_FSI_Ωs, trian_solid, quad_solid)
    t_FSI_Γi = FETerm(res_FSI_Γi,jac_FSI_Γi,trian_Γi,quad_Γi)
    op_FSI = TransientFEOperator(X_FSI,Y_FSI,t_FSI_Ωf,t_FSI_Ωs)#,t_FSI_Γi)

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
        writePVD(filePath, trian_fluid, [(xh, 0.0)])
    end

    # Solve FSI problem
    @timeit "FSI problem" begin
        println("Solving FSI problem")
		    xh0  = interpolate(X_FSI(0.0),xh)
		    nls = NLSolver(
				    #GmresSolver(preconditioner=ilu,τ=1.0e-6),
				    #GmresSolver(preconditioner=AMGPreconditioner{SmoothedAggregation}),
				    show_trace = true,
				    method = :newton,
				    #linesearch = HagerZhang(),
				    linesearch = BackTracking(),
				    ftol = 1.0e-6,
            iterations = 50
		    )
		    odes =  ThetaMethod(nls, dt, 0.5)
		    solver = TransientFESolver(odes)
		    sol_FSI = solve(solver, op_FSI, xh0, t0, tf)
		    writePVD(filePath, trian_fluid, sol_FSI, append=true)
    end

end

function get_boundary_conditions(strategy::WeakForms.MeshStrategy{:linearElasticity},u_in,u_noSlip)
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
function get_boundary_conditions(strategy::WeakForms.MeshStrategy{:neoHookean},u_in,u_noSlip)
    _strategy = WeakForms.MeshStrategy{:linearElasticity}()
    get_boundary_conditions(_strategy,u_in,u_noSlip)
end
function get_boundary_conditions(strategy::WeakForms.MeshStrategy{:biharmonic},u_in,u_noSlip)
    _strategy = WeakForms.MeshStrategy{:linearElasticity}()
    get_boundary_conditions(_strategy,u_in,u_noSlip)
end

function get_FE_spaces(strategy::WeakForms.MeshStrategy{:linearElasticity},model,model_fluid,order,bconds)
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
        X_FSI = MultiFieldFESpace([Uu_FSI,Uv_FSI,P])
    )
end

function get_FE_spaces(strategy::WeakForms.MeshStrategy{:neoHookean},model,model_fluid,order,bconds)
    _strategy = WeakForms.MeshStrategy{:linearElasticity}()
    get_FE_spaces(_strategy,model,model_fluid,order,bconds)
end

function get_FE_spaces(strategy::WeakForms.MeshStrategy{:biharmonic},model,model_fluid,order,bconds)
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
        X_FSI = MultiFieldFESpace([Uw_FSI,Uu_FSI,Uv_FSI,P])
    )
end
