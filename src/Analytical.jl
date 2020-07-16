function execute(problem::Problem{:analytical};kwargs...)

    # Problem setting (Default FSI-2)
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
    # Time stepping
    t0 = _get_kwarg(:t0,kwargs,0.0)
    tf = _get_kwarg(:tf,kwargs,0.1)
    dt = _get_kwarg(:dt,kwargs,0.1)

    # Mesh strategy
    strategyName = _get_kwarg(:strategy,kwargs,"linearElasticity")
    strategy = WeakForms.MeshStrategy{Symbol(strategyName)}()

    # Define BC functions
    println("Defining Boundary conditions")
    # global ut(x, t) = VectorValue( x[1]^2*x[2] , -x[1]*x[2]^2) * t
    # global vt(x, t) = VectorValue( x[1]^2*x[2] , -x[1]*x[2]^2)
    # global ut(t::Real) = x -> ut(x,t)
    # global vt(t::Real) = x -> vt(x,t)
    # ∂tut(t) = x -> VectorValue(ForwardDiff.derivative(t -> get_array(ut(x,t)),t))
    # ∂tvt(t) = x -> VectorValue(ForwardDiff.derivative(t -> get_array(vt(x,t)),t))
    # ∂tut(x,t) = ∂tut(t)(x)
    # ∂tvt(x,t) = ∂tvt(t)(x)
    # ∂t(::typeof(ut)) = ∂tut
    # ∂t(::typeof(vt)) = ∂tvt

    #u_in(x, t) = VectorValue(1.5 * Um * x[2] * (H - x[2]) / ((H / 2)^2), 0.0)
    #u_noSlip(x, t) = VectorValue(0.0, 0.0)
    #ut_in(t::Real) = x -> u_in(x, t)
    #ut_noSlip(t::Real) = x -> u_noSlip(x, t)
    #∂tu_in(t) = x -> VectorValue(0.0, 0.0)
    #∂tu_in(x, t) = ∂tu_in(t)(x)
    #∂t(::typeof(ut_in)) = ∂tu_in
    #∂t(::typeof(ut_noSlip)) = ∂tu_in
    bconds = get_boundary_conditions(problem,strategy,ut,vt)

    # Discrete model
    println("Defining discrete model")
    domain = (-1,1,-1,1)
    partition = (n_m,n_m)
    model = CartesianDiscreteModel(domain,partition)
    trian = Triangulation(model)
    R = 0.5
    function is_in(coords)
		    n = length(coords)
		    x = (1/n)*sum(coords)
		    d = x[1]^2 + x[2]^2 - R^2
		    d < 0
    end
    oldcell_to_coods = get_cell_coordinates(trian)
    oldcell_to_is_in = collect1d(apply(is_in,oldcell_to_coods))
    incell_to_cell = findall(oldcell_to_is_in)
    outcell_to_cell = findall(collect(Bool, .! oldcell_to_is_in))
    model_solid = DiscreteModel(model,incell_to_cell)
    model_fluid = DiscreteModel(model,outcell_to_cell)

    # Triangulations
    println("Defining triangulations")
    trian_solid = Triangulation(model_solid)
    trian_fluid = Triangulation(model_fluid)

    # Construct boundary triangulation at the fluid-solid interface
    topo = get_grid_topology(model_fluid)
    labeling = get_face_labeling(model_fluid)
    D = num_cell_dims(model_fluid)
    fluid_boundary_mask = collect(Bool,get_isboundary_face(topo,D-1))
    fluid_outer_boundary_mask = get_face_mask(labeling,"boundary",D-1)
    fluid_interface_mask = collect(Bool,fluid_boundary_mask .!= fluid_outer_boundary_mask)
    trian_Γi = BoundaryTriangulation(model_fluid,fluid_interface_mask)
    n_Γi = get_normal_vector(trian_Γi)

    println(get_face_entity(labeling,D-1))

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

    # Test FE Spaces
    println("Defining FE spaces")
    Y_ST, X_ST, Y_NS, X_NS = get_FE_spaces(problem,strategy,model,model_fluid,order,bconds)

    # # Stokes problem for initial solution
    # println("Defining Stokes operator")
    # res_ST(x,y) = WeakForms.stokes_residual(strategy,x,y,μ_f)
    # jac_ST(x,dx,y) = WeakForms.stokes_residual(strategy,dx,y,μ_f)
    # t_ST_Ωf = FETerm(res_ST, jac_ST, trian_fluid, quad_fluid)
    # op_ST = FEOperator(X_ST,Y_ST,t_ST_Ωf)

    # # Setup equation parameters
    # fsi_f_params = Dict("μ"=>μ_f,
    #                     "ρ"=>ρ_f,
    #                     "E"=>E_m,
    #                     "ν"=>ν_m,
    #                     "vol"=>vol_fluid)
    # fsi_s_params = Dict("ρ"=>ρ_s,
    #                     "E"=>E_s,
    #                     "ν"=>ν_s,
    #                     "vol"=>vol_solid)
    # fsi_Γi_params = Dict("n"=>n_Γi,
    #                      "E"=>E_m,
    #                      "ν"=>ν_m,
    #                      "vol"=>vol_Γi)

    # # FSI problem
    # println("Defining FSI operator")
    # res_FSI_Ωf(t,x,xt,y) = WeakForms.fsi_residual_Ωf(strategy,x,xt,y,fsi_f_params)
    # jac_FSI_Ωf(t,x,xt,dx,y) = WeakForms.fsi_jacobian_Ωf(strategy,x,xt,dx,y,fsi_f_params)
    # jac_t_FSI_Ωf(t,x,xt,dxt,y) = WeakForms.fsi_jacobian_t_Ωf(strategy,x,xt,dxt,y,fsi_f_params)
    # res_FSI_Ωs(t,x,xt,y) = WeakForms.fsi_residual_Ωs(strategy,x,xt,y,fsi_s_params)
    # jac_FSI_Ωs(t,x,xt,dx,y) = WeakForms.fsi_jacobian_Ωs(strategy,x,xt,dx,y,fsi_s_params)
    # jac_t_FSI_Ωs(t,x,xt,dxt,y) = WeakForms.fsi_jacobian_t_Ωs(strategy,x,xt,dxt,y,fsi_s_params)
    # res_FSI_Γi(x,y) = WeakForms.fsi_residual_Γi(strategy,x,y,fsi_Γi_params)
    # jac_FSI_Γi(x,dx,y) = WeakForms.fsi_jacobian_Γi(strategy,x,dx,y,fsi_Γi_params)
    # t_FSI_Ωf = FETerm(res_FSI_Ωf, jac_FSI_Ωf, jac_t_FSI_Ωf, trian_fluid, quad_fluid)
    # t_FSI_Ωs = FETerm(res_FSI_Ωs, jac_FSI_Ωs, jac_t_FSI_Ωs, trian_solid, quad_solid)
    # t_FSI_Γi = FETerm(res_FSI_Γi,jac_FSI_Γi,trian_Γi,quad_Γi)
    # op_FSI = TransientFEOperator(X_FSI,Y_FSI,t_FSI_Ωf,t_FSI_Ωs,t_FSI_Γi)

    # # Setup output files
    # folderName = "fsi-results"
    # fileName = "fields"
    # if !isdir(folderName)
    #     mkdir(folderName)
    # end
    # filePath = join([folderName, fileName], "/")

    # # Solve Stokes problem
    # @timeit "ST problem" begin
    #     println("Solving Stokes problem")
    #     xh = solve(op_ST)
    #     writePVD(filePath, trian_fluid, [(xh, 0.0)])
    # end

    # # Solve FSI problem
    # @timeit "FSI problem" begin
    #     println("Solving FSI problem")
		#     xh0  = interpolate(X_FSI(0.0),xh)
		#     nls = NLSolver(
		# 		    #GmresSolver(preconditioner=ilu,τ=1.0e-6),
		# 		    #GmresSolver(preconditioner=AMGPreconditioner{SmoothedAggregation}),
		# 		    show_trace = true,
		# 		    method = :newton,
		# 		    #linesearch = HagerZhang(),
		# 		    linesearch = BackTracking(),
		# 		    ftol = 1.0e-6,
    #         iterations = 50
		#     )
		#     odes =  ThetaMethod(nls, dt, 0.5)
		#     solver = TransientFESolver(odes)
		#     sol_FSI = solve(solver, op_FSI, xh0, t0, tf)
		#     writePVD(filePath, trian_fluid, sol_FSI, append=true)
    # end

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
        ST_Vu_values = [u,u],
        ST_Vv_values = [v,v],
    )
end
