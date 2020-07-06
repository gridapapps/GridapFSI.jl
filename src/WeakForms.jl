module WeakForms

using Gridap

# Laws
include("Laws.jl")

# Material properties
# ===================
function lame_parameters(E,ν)
		λ = (E*ν)/((1+ν)*(1-2*ν))
		μ = E/(2*(1+ν))
		(λ, μ)
end

# Stokes
# ======
function stokes_vp_residual(x,y,μ)
    v, p = x
    φ, q = y
    (ε(φ) ⊙ σ_dev(μ,ε(v))) - ((∇⋅φ) * p) + (q * (∇⋅v))
end

function stokes_uvp_residual(x,y,μ)
    u, v, p = x
    ϕ, φ, q = y
    (ε(ϕ) ⊙ σ_dev(μ,ε(u))) + stokes_vp_residual([v,p],[φ,q],μ)
end

# FSI
# ===
# Residuals
function a_uvp_ϕ_Ωf(x,y,E,ν)
    u, v, p = x
    ϕ, φ, q = y
    (λ_m,μ_m) = lame_parameters(E,ν)
    α(u) = α_m(J(∇(u)))
    λ(u) = α(u)*λ_m
    μ(u) = α(u)*μ_m
    (ε(ϕ) ⊙ σ_m(λ(u),μ(u),ε(u)))
end
function a_uvp_φ_Ωf(x, xt, y, μ, ρ)
    u, v, p = x
    ut, vt, pt = xt
    ϕ, φ, q = y
		temp_term = φ ⋅ ( J(∇(u)) * ρ * vt )
		conv_term = φ ⋅ ( J(∇(u)) * ρ * conv(	Finv(∇(u))⋅(v-ut), ∇(v)) )
		visc_term = ( ∇(φ) ⊙ ( J(∇(u)) * σ_dev(μ,∇(v),Finv(∇(u))) ⋅ FinvT(∇(u))) )
	  pres_term = - (∇⋅φ) * J(∇(u)) * p * tr(FinvT(∇(u)))
		temp_term + conv_term + visc_term + pres_term
end
function a_uvp_q_Ωf(x, y)
    u, v, p = x
    ϕ, φ, q = y
		q * (J(∇(u))*(∇(v)⊙FinvT(∇(u))))
end
function a_uvp_ϕ_Ωs(x, xt, y)
    u,v,p = x
		ut, vt,pt = xt
    ϕ,φ,q = y
		(ϕ⋅ut) + 0.0*(u⋅ϕ) - (ϕ⋅v) 
end
function a_uvp_φ_Ωs(x, xt, y, ρ, E, ν)
    u,v,p = x
		ut, vt,pt = xt
    ϕ,φ,q = y
    (λ,μ) = lame_parameters(E,ν)
		(φ⋅(ρ*vt)) + 0.0*(φ⋅(ρ*v)) + (∇(φ) ⊙ (F(∇(u))⋅S(∇(u))))
end
function a_uvp_ϕ_Γi(x,y,n,E,ν)
    u, v, p = x
    ϕ, φ, q = y
    (λ_m,μ_m) = lame_parameters(E,ν)
    α(u) = α_m(J(∇(u)))
    λ(u) = α(u)*λ_m
    μ(u) = α(u)*μ_m
    - (ϕ ⋅  (n⋅σ_m(λ(u),μ(u),ε(u))) )
end

# Jacobians
function da_uvp_du_ϕ_Ωf(x, dx, y, E, ν)
		u, v, p = x
    du, dv, dp = dx
    ϕ, φ, q = y
    (λ_m,μ_m) = lame_parameters(E,ν)
    α(u) = α_m(J(∇(u)))
    λ(u) = α(u)*λ_m
    μ(u) = α(u)*μ_m
    dα(u,du) = dα_m(J(∇(u)),dJ(∇(u),∇(du)))
    dλ(u,du) = dα(u,du)*λ_m
    dμ(u,du) = dα(u,du)*μ_m
    (ε(ϕ) ⊙ dσ_m(λ(u),dλ(u,du),μ(u),dμ(u,du),ε(u),ε(du)))
end
function da_uvp_du_φ_Ωf(x, xt, dx, y, μ, ρ)
    u, v, p = x
    ut, vt, pt = xt
    du, dv, dp = dx
    ϕ, φ, q = y
		temp_term = φ ⋅ ( dJ(∇(u),∇(du)) * ρ * vt )
		conv_term = φ ⋅ ( ( dJ(∇(u),∇(du)) * ρ * conv( Finv(∇(u))⋅(v-ut), ∇(v)) ) +
											( J(∇(u)) * ρ * conv(	dFinv(∇(u),∇(du))⋅(v-ut), ∇(v)) ) )
		visc_term = ∇(φ) ⊙ ( dJ(∇(u),∇(du)) * σ_dev(μ,∇(v),Finv(∇(u))) ⋅ FinvT(∇(u)) +
												 J(∇(u)) * σ_dev(μ,∇(v),dFinv(∇(u),∇(du))) ⋅ FinvT(∇(u)) +
												 J(∇(u)) * σ_dev(μ,∇(v),Finv(∇(u))) ⋅ dFinvT(∇(u),∇(du)) )
		pres_term = - (∇⋅φ) * p * ( dJ(∇(u),∇(du)) * tr(FinvT(∇(u))) +
															J(∇(u)) * tr(dFinvT(∇(u),∇(du))) )
		temp_term + conv_term + visc_term + pres_term
end
function da_uvp_dv_φ_Ωf(x, xt, dx, y, μ, ρ)
    u, v, p = x
    ut, vt, pt = xt
    du, dv, dp = dx
    ϕ, φ, q = y
		conv_term = φ ⋅ ( J(∇(u)) * ρ * dconv( Finv(∇(u))⋅dv, ∇(dv), Finv(∇(u))⋅(v-ut) , ∇(v)) )
		visc_term = ( ∇(φ) ⊙ ( J(∇(u)) * σ_dev(μ,∇(dv),Finv(∇(u))) ⋅ FinvT(∇(u))) )
		conv_term + visc_term
end
function da_uvp_dp_φ_Ωf(x, dx, y)
    u, v, p = x
    du, dv, dp = dx
    ϕ, φ, q = y
		pres_term = - (∇⋅φ) * J(∇(u)) * dp * tr(FinvT(∇(u)))
end
function da_uvp_dut_φ_Ωf(x, dxt, y, ρ)
    u, v, p = x
    dut, dvt, dpt = dxt
    ϕ, φ, q = y
		conv_term = - φ ⋅ ( J(∇(u)) * ρ * conv(	Finv(∇(u))⋅dut, ∇(v)) )
end
function da_uvp_dvt_φ_Ωf(x, dxt, y, ρ)
    u, v, p = x
    dut, dvt, dpt = dxt
    ϕ, φ, q = y
		temp_term = φ ⋅ ( J(∇(u)) * ρ * dvt )
end
function da_uvp_du_q_Ωf(x, dx, y)
    u, v, p = x
    du, dv, dp = dx
    ϕ, φ, q = y
		q * ( dJ(∇(u),∇(du))*(∇(v)⊙FinvT(∇(u))) + J(∇(u))*(∇(v)⊙dFinvT(∇(u),∇(du))) )
end
function da_uvp_dv_q_Ωf(x, dx, y)
    u, v, p = x
    du, dv, dp = dx
    ϕ, φ, q = y
		q * ( J(∇(u))*(∇(dv)⊙FinvT(∇(u))) )
end
function da_uvp_dx_ϕ_Ωs(x, dx, y)
    u,v,p = x
    du,dv,dp = dx
    ϕ,φ,q = y
		0.0*(du⋅ϕ) - (ϕ⋅dv) 
end
function da_uvp_dx_φ_Ωs(x, dx, y, ρ, E, ν)
    u,v,p = x
    du,dv,dp = dx
    ϕ,φ,q = y
    (λ,μ) = lame_parameters(E,ν)
		0.0*(φ⋅(ρ*dv)) + (∇(φ) ⊙ ( dF(∇(du))⋅S(∇(u)) + (F(∇(u))⋅dS(∇(u),∇(du))) ) )
end
function da_uvp_dxt_Ωs(x, dxt, y, ρ)
    u, v, p = x
    dut, dvt, dpt = dxt
    ϕ, φ, q = y
		ϕ⋅dut + (φ⋅(ρ*dvt))
end
function da_uvp_du_ϕ_Γi(x,dx,y, n, E, ν)
		u, v, p = x
    du, dv, dp = dx
    ϕ, φ, q = y
    (λ_m,μ_m) = lame_parameters(E,ν)
    α(u) = α_m(J(∇(u)))
    λ(u) = α(u)*λ_m
    μ(u) = α(u)*μ_m
    dα(u,du) = dα_m(J(∇(u)),dJ(∇(u),∇(du)))
    dλ(u,du) = dα(u,du)*λ_m
    dμ(u,du) = dα(u,du)*μ_m
    - (ϕ ⋅  (n⋅dσ_m(λ(u),dλ(u,du),μ(u),dμ(u,du),ε(u),ε(du))) )
end

fsi_uvp_residual_Ωf(x,xt,y,μ,ρ,E,ν) = 
    a_uvp_ϕ_Ωf(x,y,E,ν) +
    a_uvp_φ_Ωf(x,xt,y,μ,ρ) +
    a_uvp_q_Ωf(x, y)

fsi_uvp_jacobian_Ωf(x,xt,dx,y,μ,ρ,E,ν) = 
    da_uvp_du_ϕ_Ωf(x,dx,y,E,ν) +
    da_uvp_du_φ_Ωf(x,xt,dx,y,μ,ρ) +
    da_uvp_dv_φ_Ωf(x,xt,dx,y,μ,ρ) +
    da_uvp_dp_φ_Ωf(x,dx,y) +
    da_uvp_du_q_Ωf(x,dx,y) +
    da_uvp_dv_q_Ωf(x,dx,y)

fsi_uvp_jacobian_t_Ωf(x,xt,dxt,y,ρ) =
    da_uvp_dut_φ_Ωf(x,dxt,y,ρ) +
    da_uvp_dvt_φ_Ωf(x,dxt,y,ρ)

fsi_uvp_residual_Ωs(x,xt,y,ρ,E,ν) =
    a_uvp_ϕ_Ωs(x,xt,y) +
    a_uvp_φ_Ωs(x,xt,y,ρ,E,ν)

fsi_uvp_jacobian_Ωs(x,xt,dx,y,ρ,E,ν) =
    da_uvp_dx_ϕ_Ωs(x,dx,y) +
    da_uvp_dx_φ_Ωs(x,dx,y,ρ,E,ν) 

fsi_uvp_jacobian_t_Ωs(x,xt,dxt,y,ρ) =
    da_uvp_dxt_Ωs(x,dxt,y,ρ)

fsi_uvp_residual_Γi(x,y,n,E,ν) = a_uvp_ϕ_Γi(x,y,n,E,ν)
fsi_uvp_jacobian_Γi(x,dx,y,n,E,ν) = da_uvp_du_ϕ_Γi(x,dx,y,n,E,ν)

end

# # FSI Bilinear forms


# # FSI Jacobians
# function daFSI_dut_ϕ_f(x, dxt, y)
# 		u, v, p = x
#     dut, dvt, dpt = dxt
#     ϕ, φ, q = y
# 		(∇(ϕ) ⊙ σ_m(α(J(∇(u))),ε(dut))) 
# end
# function daFSI_du_ϕ_Γi(x,xt,dx,y)
# 		u,v,p = x
#     ut, vt, pt = xt
# 		du,dv,dp = dx
#     ϕ,φ,q = y
#     - (ϕ ⋅  (n_Γi⋅dσ_m(dα(J(∇(u)),dJ(∇(u),∇(du))),ε(ut))) )
# end
# function daFSI_dut_ϕ_Γi(x,dxt,y)
# 		u,v,p = x
# 		dut,dvt,dpt = dxt
#     ϕ,φ,q = y
# 		- (ϕ ⋅  (n_Γi⋅σ_m(α(J(∇(u))),ε(dut))) ) 
# end


