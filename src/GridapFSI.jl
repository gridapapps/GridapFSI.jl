module GridapFSI

using TimerOutputs

export main

include("WeakForms.jl")
include("FSIDrivers.jl")

using GridapFSI.FSIDrivers

# Generic struct for multiple dispatch
#struct Problem{Kind} end
# struct Geometry{Kind} end

# # Include problem-specific implementations
# include("ElasticFlag.jl")

# # ### Problem parameters
# const Um = 1.0
# const H = 0.41
# const ⌀ = 0.1
# const ρ = 1.0
# const t0 = 0.0
# const Re = 100.0

# # ### Boundary conditions
# u_in(x, t) = VectorValue(1.5 * Um * x[2] * (H - x[2]) / ((H / 2)^2), 0.0)
# u_noSlip(x, t) = VectorValue(0.0, 0.0)
# u_in(t::Real) = x -> u_in(x, t)
# u_noSlip(t::Real) = x -> u_noSlip(x, t)
# ∂tu_in(t) = x -> VectorValue(0.0, 0.0)
# ∂tu_in(x, t) = ∂tu_in(t)(x)
# ∂t(::typeof(u_in)) = ∂tu_in
# ∂t(::typeof(u_noSlip)) = ∂tu_in

# # ### Material properties
# function lame_parameters(E,ν)
# 		λ = (E*ν)/((1+ν)*(1-2*ν))
# 		μ = E/(2*(1+ν))
# 		(λ, μ)
# end
# const (λ_s,μ_s) = lame_parameters(1.4e6, 0.4)
# const ρ_s = 1.0e4
# const ρ_f = 1.0e3
# const μ_f = ρ_f * Um * ⌀ / Re
# const E_m = 1.0
# const ν_m = -0.1

# # Output function
# function writePVD(filePath::String, trian::Triangulation, sol; append=false)
#     outfiles = paraview_collection(filePath, append=append) do pvd
#         for (i, (xh, t)) in enumerate(sol)
#             uh = restrict(xh.blocks[1],trian)
#             vh = restrict(xh.blocks[2],trian)
# 						ph = restrict(xh.blocks[3],trian)
#             pvd[t] = createvtk(
#                 trian,
#                 filePath * "_$t.vtu",
#                 cellfields = ["uh" => uh, "vh" => vh, "ph" => ph]
#             )
#         end
#     end
# end

# include("Laws.jl")
# include("WeakForms.jl")

function main(;problemName::String="elasticFlag",kwargs...)

    reset_timer!()

    # Setup problem
    problem = Problem{Symbol(problemName)}()

    # Execute driver
    execute(problem; kwargs...)


    # # Discrete model

    # modelName = geometryName * meshName
    # model = DiscreteModelFromFile("../models/"*modelName*".json")
    # model_solid = DiscreteModel(model,"solid")
    # model_fluid = DiscreteModel(model,"fluid")  # # Define functions
    # # params = get_parameters(problem)

    # # Output boundary conditions
    # bconds = Dict(
    #     # Tags
    #     "FSI_Vu_tags" => ["inlet", "noSlip", "cylinder","fixed","outlet"],
    #     "FSI_Vv_tags" => ["inlet", "noSlip", "cylinder","fixed"],
    #     "ST_Vu_tags" => ["inlet", "noSlip", "cylinder","interface","outlet"],
    #     "ST_Vv_tags" => ["inlet", "noSlip", "cylinder","interface"],
    #     # Values,
    #     "FSI_Vu_values" => [u_noSlip, u_noSlip, u_noSlip, u_noSlip, u_noSlip],
    #     "FSI_Vv_values" => [u_in, u_noSlip, u_noSlip, u_noSlip],
    #     "ST_Vu_values" => [u_noSlip(0.0), u_noSlip(0.0), u_noSlip(0.0), u_noSlip(0.0), u_noSlip(0.0)],
    #     "ST_Vv_values" => [u_in(0.0), u_noSlip(0.0), u_noSlip(0.0), u_noSlip(0.0)],
    # )

    # # Triangulations and quadratures
    # trian = Triangulation(model)
    # trian_solid = Triangulation(model_solid)
    # trian_fluid = Triangulation(model_fluid)
    # trian_Γi = BoundaryTriangulation(model_fluid,"interface")
    # k = 2
    # degree = 2*k
    # bdegree = 2*k
    # quad_solid = CellQuadrature(trian_solid,degree)
    # quad_fluid = CellQuadrature(trian_fluid,degree)
    # quad_Γi = CellQuadrature(trian_Γi,bdegree)
    # n_Γi = get_normal_vector(trian_Γi)

    # # Compute cell area (auxiliar quantity for mesh motion eq.)
    # vol = cell_measure(trian_fluid,trian)
    # vol_fluid = reindex(vol,trian_fluid)
    # vol_Γi = reindex(vol,trian_Γi)

    # # Test FE Spaces
    # Vu_FSI = TestFESpace(
    #     model=model,
    #     valuetype=VectorValue{2,Float64},
    #     reffe=:Lagrangian,
    #     order=k,
    #     conformity =:H1,
    #     dirichlet_tags=bconds["FSI_Vu_tags"])
    # Vv_FSI = TestFESpace(
    #     model=model,
    #     valuetype=VectorValue{2,Float64},
    #     reffe=:Lagrangian,
    #     order=k,
    #     conformity =:H1,
    #     dirichlet_tags=bconds["FSI_Vv_tags"])
    # Vu_ST = TestFESpace(
    #     model=model_fluid,
		#     #model=model,
    #     valuetype=VectorValue{2,Float64},
    #     reffe=:Lagrangian,
    #     order=k,
    #     conformity =:H1,
    #     dirichlet_tags=bconds["ST_Vu_tags"])
    # Vv_ST = TestFESpace(
    #     model=model_fluid,
		#     #model=model,
    #     valuetype=VectorValue{2,Float64},
    #     reffe=:Lagrangian,
    #     order=k,
    #     conformity =:H1,
    #     dirichlet_tags=bconds["ST_Vv_tags"])
    # Q = TestFESpace(
		#     model=model_fluid,
		#     #model=model,
		#     valuetype=Float64,
		#     order=k-1,
		#     reffe=:Lagrangian,
		#     conformity=:C0)

    # # Trial FE Spaces
    # Uu_ST = TrialFESpace(Vu_ST,bconds["ST_Vu_values"])
    # Uv_ST = TrialFESpace(Vv_ST,bconds["ST_Vv_values"])
    # Uu_FSI = TransientTrialFESpace(Vu_FSI,bconds["FSI_Vu_values"])
    # Uv_FSI = TransientTrialFESpace(Vv_FSI,bconds["FSI_Vv_values"])
    # P = TrialFESpace(Q)

    # # Multifield FE Spaces
    # Y_ST = MultiFieldFESpace([Vu_ST,Vv_ST,Q])
    # X_ST = MultiFieldFESpace([Uu_ST,Uv_ST,P])
    # Y_FSI = MultiFieldFESpace([Vu_FSI,Vv_FSI,Q])
    # X_FSI = MultiFieldFESpace([Uu_FSI,Uv_FSI,P])

    # # Stokes Bilinear forms
    # function aST_ϕ_f(x, y)
    #     u, v, p = x
    #     ϕ, φ, q = y
		#     (∇(ϕ) ⊙ ∇(u))
    # end
    # function aST_φ_f(x, y)
    #     u, v, p = x
    #     ϕ, φ, q = y
		#     visc_term = ε(φ) ⊙ ( 2*μ_f*ε(v) )
		#     pres_term = (∇⋅φ) * p

		#     visc_term + pres_term
    # end
    # function aST_q_f(x, y)
    #     u, v, p = x
    #     ϕ, φ, q = y
		#     q * (∇⋅v)
    # end

    # # FSI Bilinear forms
    # function aFSI_ϕ_f(x, xt, y)
    #     u, v, p = x
    #     ut, vt, pt = xt
    #     ϕ, φ, q = y
		#     (∇(ϕ) ⊙ σ_m(α(J(∇(u))),ε(ut)))
    # end
    # function aFSI_ϕ_Γi(x,xt,y)
		#     u,v,p = x
    #     ut, vt, pt = xt
    #     ϕ,φ,q = y
    #     - (ϕ ⋅  (n_Γi⋅σ_m(α(J(∇(u))),ε(ut))) )
    # end
    # function aFSI_φ_f(x, xt, y)
    #     u, v, p = x
    #     ut, vt, pt = xt
    #     ϕ, φ, q = y
		#     temp_term = φ ⋅ ( J(∇(u)) * ρ_f * vt )
		#     conv_term = φ ⋅ ( J(∇(u)) * ρ_f * conv(	Finv(∇(u))⋅(v-ut), ∇(v)) )
		#     visc_term = ( ∇(φ) ⊙ ( J(∇(u)) * σ_dev(∇(v),Finv(∇(u))) ⋅ FinvT(∇(u))) )
		#     pres_term = (∇⋅φ) * J(∇(u)) * p * tr(FinvT(∇(u)))
		#     temp_term + conv_term + visc_term + pres_term
    # end
    # function aFSI_q_f(x, y)
    #     u, v, p = x
    #     ϕ, φ, q = y
		#     #q * (∇⋅v)
		#     q * (J(∇(u))*(∇(v)⊙FinvT(∇(u))))

    # end
    # function aFSI_ϕ_s(x, xt, y)
    #     u,v,p = x
		#     ut, vt,pt = xt
    #     ϕ,φ,q = y
		#     (ϕ⋅ut) + 0.0*(u⋅ϕ) - (ϕ⋅v) 
    # end
    # function aFSI_φ_s(x, xt, y)
    #     u,v,p = x
		#     ut, vt,pt = xt
    #     ϕ,φ,q = y
		#     (φ⋅(ρ_s*vt)) + 0.0*(φ⋅(ρ_s*v)) + (∇(φ) ⊙ (F(∇(u))⋅S(∇(u))))
    # end

    # # FSI Jacobians
    # function daFSI_du_ϕ_f(x, xt, dx, y)
		#     u, v, p = x
		#     ut, vt, pt = xt
    #     du, dv, dp = dx
    #     ϕ, φ, q = y
		#     (∇(ϕ) ⊙ dσ_m(dα(J(∇(u)),dJ(∇(u),∇(du))),ε(ut)))
    # end
    # function daFSI_dut_ϕ_f(x, dxt, y)
		#     u, v, p = x
    #     dut, dvt, dpt = dxt
    #     ϕ, φ, q = y
		#     (∇(ϕ) ⊙ σ_m(α(J(∇(u))),ε(dut))) 
    # end
    # function daFSI_du_ϕ_Γi(x,xt,dx,y)
		#     u,v,p = x
    #     ut, vt, pt = xt
		#     du,dv,dp = dx
    #     ϕ,φ,q = y
    #     - (ϕ ⋅  (n_Γi⋅dσ_m(dα(J(∇(u)),dJ(∇(u),∇(du))),ε(ut))) )
    # end
    # function daFSI_dut_ϕ_Γi(x,dxt,y)
		#     u,v,p = x
		#     dut,dvt,dpt = dxt
    #     ϕ,φ,q = y
		#     - (ϕ ⋅  (n_Γi⋅σ_m(α(J(∇(u))),ε(dut))) ) 
    # end
    # function daFSI_du_φ_f(x, xt, dx, y)
    #     u, v, p = x
    #     ut, vt, pt = xt
    #     du, dv, dp = dx
    #     ϕ, φ, q = y
		#     temp_term = φ ⋅ ( dJ(∇(u),∇(du)) * ρ_f * vt )
		#     conv_term = φ ⋅ ( ( dJ(∇(u),∇(du)) * ρ_f * conv( Finv(∇(u))⋅(v-ut), ∇(v)) ) +
		# 									    ( J(∇(u)) * ρ_f * conv(	dFinv(∇(u),∇(du))⋅(v-ut), ∇(v)) ) )
		#     visc_term = ∇(φ) ⊙ ( dJ(∇(u),∇(du)) * σ_dev(∇(v),Finv(∇(u))) ⋅ FinvT(∇(u)) +
		# 										     J(∇(u)) * σ_dev(∇(v),dFinv(∇(u),∇(du))) ⋅ FinvT(∇(u)) +
		# 										     J(∇(u)) * σ_dev(∇(v),Finv(∇(u))) ⋅ dFinvT(∇(u),∇(du)) )
		#     pres_term = (∇⋅φ) * p * ( dJ(∇(u),∇(du)) * tr(FinvT(∇(u))) +
		# 													    J(∇(u)) * tr(dFinvT(∇(u),∇(du))) )
		#     temp_term + conv_term + visc_term + pres_term
    # end
    # function daFSI_dv_φ_f(x, xt, dx, y)
    #     u, v, p = x
    #     ut, vt, pt = xt
    #     du, dv, dp = dx
    #     ϕ, φ, q = y
		#     conv_term = φ ⋅ ( J(∇(u)) * ρ_f * dconv( Finv(∇(u))⋅dv, ∇(dv), Finv(∇(u))⋅(v-ut) , ∇(v)) )
		#     visc_term = ( ∇(φ) ⊙ ( J(∇(u)) * σ_dev(∇(dv),Finv(∇(u))) ⋅ FinvT(∇(u))) )
		#     conv_term + visc_term
    # end
    # function daFSI_dp_φ_f(x, dx, y)
    #     u, v, p = x
    #     du, dv, dp = dx
    #     ϕ, φ, q = y
		#     pres_term = (∇⋅φ) * J(∇(u)) * dp * tr(FinvT(∇(u)))
    # end
    # function daFSI_dut_φ_f(x, dxt, y)
    #     u, v, p = x
    #     dut, dvt, dpt = dxt
    #     ϕ, φ, q = y
		#     conv_term = - φ ⋅ ( J(∇(u)) * ρ_f * conv(	Finv(∇(u))⋅dut, ∇(v)) )
    # end
    # function daFSI_dvt_φ_f(x, dxt, y)
    #     u, v, p = x
    #     dut, dvt, dpt = dxt
    #     ϕ, φ, q = y
		#     temp_term = φ ⋅ ( J(∇(u)) * ρ_f * dvt )
    # end
    # function daFSI_du_q_f(x, dx, y)
    #     u, v, p = x
    #     du, dv, dp = dx
    #     ϕ, φ, q = y
		#     q * ( dJ(∇(u),∇(du))*(∇(v)⊙FinvT(∇(u))) + J(∇(u))*(∇(v)⊙dFinvT(∇(u),∇(du))) )
    # end
    # function daFSI_dv_q_f(x, dx, y)
    #     u, v, p = x
    #     du, dv, dp = dx
    #     ϕ, φ, q = y
		#     q * ( J(∇(u))*(∇(dv)⊙FinvT(∇(u))) )
    # end
    # function daFSI_ϕ_s(x, dx, y)
    #     u,v,p = x
    #     du,dv,dp = dx
    #     ϕ,φ,q = y
		#     0.0*(du⋅ϕ) - (ϕ⋅dv) 
    # end
    # function daFSI_φ_s(x, dx, y)
    #     u,v,p = x
    #     du,dv,dp = dx
    #     ϕ,φ,q = y
		#     0.0*(φ⋅(ρ_s*dv)) + (∇(φ) ⊙ ( dF(∇(du))⋅S(∇(u)) + (F(∇(u))⋅dS(∇(u),∇(du))) ) )
    # end
    # function daFSI_dt_s(x, dxt, y)
    #     u, v, p = x
    #     dut, dvt, dpt = dxt
    #     ϕ, φ, q = y
		#     ϕ⋅dut + (φ⋅(ρ_s*dvt))
    # end


    # # Stokes FE Operator
    # res_ST_f(x,y) = aST_ϕ_f(x,y) + aST_φ_f(x,y) + aST_q_f(x,y)
    # jac_ST_f(x,dx,y) = aST_ϕ_f(dx,y) + aST_φ_f(dx,y) + aST_q_f(dx,y)
    # t_ST_Ωf = FETerm(res_ST_f, jac_ST_f, trian_fluid, quad_fluid)
    # op_ST = FEOperator(X_ST,Y_ST,t_ST_Ωf)

    # # FSI FE Operator
    # res_FSI_f(t,x,xt,y) =	aFSI_ϕ_f(x,xt,y) + aFSI_φ_f(x,xt,y) + aFSI_q_f(x,y)
    # jac_FSI_f(t,x,xt,dx,y) = 	daFSI_du_ϕ_f(x,xt,dx,y) + daFSI_du_φ_f(x,xt,dx,y) + daFSI_dv_φ_f(x,xt,dx,y) + daFSI_dp_φ_f(x,dx,y) + daFSI_du_q_f(x,dx,y) + daFSI_dv_q_f(x,dx,y)
    # #jac_FSI_f(t,x,xt,dx,y) = 	daFSI_du_φ_f(x,xt,dx,y) + daFSI_dv_φ_f(x,xt,dx,y) + daFSI_dp_φ_f(x,dx,y) + daFSI_du_q_f(x,dx,y) + daFSI_dv_q_f(x,dx,y)
    # jac_t_FSI_f(t,x,xt,dxt,y) =	 daFSI_dut_ϕ_f(x,dxt,y) + daFSI_dut_φ_f(x,dxt,y) + daFSI_dvt_φ_f(x,dxt,y)
    # #jac_t_FSI_f(t,x,xt,dxt,y) =	 daFSI_dut_φ_f(x,dxt,y) + daFSI_dvt_φ_f(x,dxt,y)
    # res_FSI_s(t,x,xt,y) = aFSI_ϕ_s(x,xt,y) + aFSI_φ_s(x,xt,y)
    # jac_FSI_s(t,x,xt,dx,y) = daFSI_ϕ_s(x,dx,y) + daFSI_φ_s(x,dx,y)
    # jac_t_FSI_s(t,x,xt,dxt,y) = daFSI_dt_s(x,dxt,y)
    # res_FSI_fΓi(t,x,xt,y) = aFSI_ϕ_Γi(x,xt,y)
    # jac_FSI_fΓi(t,x,xt,dx,y) = 	daFSI_du_ϕ_Γi(x,xt,dx,y)
    # jac_t_FSI_fΓi(t,x,xt,dxt,y) = 	daFSI_dut_ϕ_Γi(x,dxt,y)
    # res_FSI_fΓi(x,y) = aFSI_ϕ_Γi(x,y)
    # jac_FSI_fΓi(x,dx,y) = 	daFSI_du_ϕ_Γi(x,dx,y)

    # t_FSI_Ωf = FETerm(res_FSI_f, jac_FSI_f, jac_t_FSI_f, trian_fluid, quad_fluid)
    # t_FSI_Ωs = FETerm(res_FSI_s, jac_FSI_s, jac_t_FSI_s, trian_solid, quad_solid)
    # t_FSI_Γi = FETerm(res_FSI_fΓi,jac_FSI_fΓi,jac_t_FSI_fΓi,trian_Γi,quad_Γi)
    # #tFSI_Γi = FETerm(res_FSI_fΓi,jac_FSI_fΓi,trian_Γi,quad_Γi)
    # op_FSI = TransientFEOperator(X_FSI,Y_FSI,t_FSI_Ωf,t_FSI_Ωs,t_FSI_Γi)

    # folderName = "fsi-results"
    # fileName = "fields"
    # if !isdir(folderName)
    #     mkdir(folderName)
    # end
    # filePath = join([folderName, fileName], "/")

    # # Solve Stokes problem
    # @timeit "ST problem" begin
    #     xh = solve(op_ST)
    #     writePVD(filePath, trian_fluid, [(xh, 0.0)])
    # end

    # # Solve FSI problem
    # @timeit "FSI problem" begin
		#     xh0  = interpolate(X_FSI(0.0),xh)
		#     nls = NLSolver(
		# 		    #GmresSolver(preconditioner=ilu,τ=1.0e-6),
		# 		    #GmresSolver(preconditioner=AMGPreconditioner{SmoothedAggregation}),
		# 		    show_trace = true,
		# 		    method = :newton,
		# 		    #linesearch = HagerZhang(),
		# 		    linesearch = BackTracking(),
		# 		    ftol = 1.0e-6
		#     )
		#     odes =  ThetaMethod(nls, 0.1, 0.5)
		#     solver = TransientFESolver(odes)
		#     sol_FSI = solve(solver, op_FSI, xh0, 0.0, 30.0)
		#     writePVD(filePath, trian_fluid, sol_FSI, append=true)
    # end

    print_timer()
    println()

end

end # module
