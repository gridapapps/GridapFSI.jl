module WeakForms

using Gridap
using LinearAlgebra: tr, inv, det

export MeshStrategy

struct MeshStrategy{Kind} end

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

function stokes_uvp_residual(strategy::MeshStrategy,x,y,μ)
    u, v, p = x
    ϕ, φ, q = y
    (ε(ϕ) ⊙ σ_dev(μ,ε(u))) + stokes_vp_residual([v,p],[φ,q],μ)
end

function stokes_uvp_residual(strategy::MeshStrategy{:biharmonic},x,y,μ)
    w, u, v, p = x
    ψ, ϕ, φ, q = y
    _strategy = MeshStrategy{:linearElasticity}()
    (ε(ψ) ⊙ σ_dev(μ,ε(w))) + stokes_uvp_residual(_strategy,[u,v,p],[ϕ,φ,q],μ)
end

# FSI
# ===
# Residuals
function a_uvp_ϕ_Ωf(strategy::MeshStrategy{:linearElasticity},x,y,E,ν)
    u, v, p = x
    ϕ, φ, q = y
    (λ_m,μ_m) = lame_parameters(E,ν)
    α(u) = α_m(J(∇(u)))
    λ(u) = α(u)*λ_m
    μ(u) = α(u)*μ_m
    (ε(ϕ) ⊙ σ_m(λ(u),μ(u),ε(u)))
end
function a_uvp_ϕ_Ωf(strategy::MeshStrategy{:neoHookean},x,y,E,ν)
    u, v, p = x
    ϕ, φ, q = y
    (λ_m,μ_m) = lame_parameters(E,ν)
    (dE(∇(ϕ),∇(u)) ⊙ S_NH(∇(u)))
end
function a_wuvp_ψ_Ωf(strategy::MeshStrategy{:biharmonic},x,y,vol)
    w, u, v, p = x
    ψ, ϕ, φ, q = y
    α = 1.0e-8
    vol*α*(ψ⋅w) + α*(∇(ψ) ⊙ ∇(u))
end
function a_wuvp_ϕ_Ωf(strategy::MeshStrategy{:biharmonic},x,y,vol)
    w, u, v, p = x
    ψ, ϕ, φ, q = y
    α = 1.0e-8
    vol*α*(∇(ϕ) ⊙ ∇(w))
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
		(φ⋅(ρ*vt)) + 0.0*(φ⋅(ρ*v)) + (∇(φ) ⊙ (F(∇(u))⋅S_SV(∇(u))))
end
function a_uvp_ϕ_Γi(strategy::MeshStrategy{:linearElasticity},x,y,n,E,ν)
    u, v, p = x
    ϕ, φ, q = y
    (λ_m,μ_m) = lame_parameters(E,ν)
    α(u) = α_m(J(∇(u)))
    λ(u) = α(u)*λ_m
    μ(u) = α(u)*μ_m
    - (ϕ ⋅  (n⋅σ_m(λ(u),μ(u),ε(u))) )
end
function a_uvp_ϕ_Γi(strategy::MeshStrategy{:neoHookean},x,y,n,E,ν)
    u, v, p = x
    ϕ, φ, q = y
    (λ_m,μ_m) = lame_parameters(E,ν)
    - (ϕ ⋅  (n⋅S_NH(∇(u))) )
end
function a_wuvp_ψ_Γi(strategy::MeshStrategy{:biharmonic},x,y,n,vol)
    w, u, v, p = x
    ψ, ϕ, φ, q = y
    α = 1.0e-8
    - vol * α * (ψ ⋅  (n⋅∇(u)))
end
function a_wuvp_ϕ_Γi(strategy::MeshStrategy{:biharmonic},x,y,n,vol)
    w, u, v, p = x
    ψ, ϕ, φ, q = y
    α = 1.0e-8
    - vol * α * (ϕ ⋅  (n⋅∇(w)))
end

# Jacobians
function da_uvp_du_ϕ_Ωf(strategy::MeshStrategy{:linearElasticity},x, dx, y, E, ν)
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
function da_uvp_du_ϕ_Ωf(strategy::MeshStrategy{:neoHookean},x, dx, y, E, ν)
		u, v, p = x
    du, dv, dp = dx
    ϕ, φ, q = y
    (λ_m,μ_m) = lame_parameters(E,ν)
    ( dE(∇(ϕ),∇(u)) ⊙ dS_NH(∇(du),∇(u)) ) + ( ∇(ϕ) ⊙ ( S_NH(∇(u))⋅∇(du) ) )
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
		0.0*(φ⋅(ρ*dv)) + (∇(φ) ⊙ ( dF(∇(du))⋅S_SV(∇(u)) + (F(∇(u))⋅dS_SV(∇(u),∇(du))) ) )
end
function da_uvp_dxt_Ωs(x, dxt, y, ρ)
    u, v, p = x
    dut, dvt, dpt = dxt
    ϕ, φ, q = y
		ϕ⋅dut + (φ⋅(ρ*dvt))
end
function da_uvp_du_ϕ_Γi(strategy::MeshStrategy{:linearElasticity},x,dx,y, n, E, ν)
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
function da_uvp_du_ϕ_Γi(strategy::MeshStrategy{:neoHookean},x,dx,y, n, E, ν)
		u, v, p = x
    du, dv, dp = dx
    ϕ, φ, q = y
    (λ_m,μ_m) = lame_parameters(E,ν)
    - (ϕ ⋅  (n⋅dS_NH(∇(du),∇(u))) )
end

fsi_uvp_residual_Ωf(strategy::MeshStrategy,x,xt,y,μ,ρ,E,ν,vol) = 
    a_uvp_ϕ_Ωf(strategy,x,y,E,ν) +
    a_uvp_φ_Ωf(x,xt,y,μ,ρ) +
    a_uvp_q_Ωf(x, y)

fsi_uvp_jacobian_Ωf(strategy::MeshStrategy,x,xt,dx,y,μ,ρ,E,ν,vol) = 
    da_uvp_du_ϕ_Ωf(strategy,x,dx,y,E,ν) +
    da_uvp_du_φ_Ωf(x,xt,dx,y,μ,ρ) +
    da_uvp_dv_φ_Ωf(x,xt,dx,y,μ,ρ) +
    da_uvp_dp_φ_Ωf(x,dx,y) +
    da_uvp_du_q_Ωf(x,dx,y) +
    da_uvp_dv_q_Ωf(x,dx,y)

fsi_uvp_jacobian_t_Ωf(strategy::MeshStrategy,x,xt,dxt,y,ρ) =
    da_uvp_dut_φ_Ωf(x,dxt,y,ρ) +
    da_uvp_dvt_φ_Ωf(x,dxt,y,ρ)

function fsi_uvp_residual_Ωf(strategy::MeshStrategy{:biharmonic},x,xt,y,μ,ρ,E,ν,vol)
    w, u, v, p = x
    wt, ut, vt, pt = xt
    ψ, ϕ, φ, q = y
    a_wuvp_ψ_Ωf(strategy,x,y,vol) +
    a_wuvp_ϕ_Ωf(strategy,x,y,vol) +
    a_uvp_φ_Ωf([u,v,p],[ut,vt,pt],[ϕ,φ,q],μ,ρ) +
    a_uvp_q_Ωf([u,v,p],[ϕ,φ,q])
end

function fsi_uvp_jacobian_Ωf(strategy::MeshStrategy{:biharmonic},x,xt,dx,y,μ,ρ,E,ν,vol)
    w, u, v, p = x
    wt, ut, vt, pt = xt
    dw, du, dv, dp = dx
    ψ, ϕ, φ, q = y
    a_wuvp_ψ_Ωf(strategy,dx,y,vol) +
    a_wuvp_ϕ_Ωf(strategy,dx,y,vol) +
    da_uvp_du_φ_Ωf([u,v,p],[ut,vt,pt],[du,dv,dp],[ϕ,φ,q],μ,ρ) +
    da_uvp_dv_φ_Ωf([u,v,p],[ut,vt,pt],[du,dv,dp],[ϕ,φ,q],μ,ρ) +
    da_uvp_dp_φ_Ωf([u,v,p],[du,dv,dp],[ϕ,φ,q]) +
    da_uvp_du_q_Ωf([u,v,p],[du,dv,dp],[ϕ,φ,q]) +
    da_uvp_dv_q_Ωf([u,v,p],[du,dv,dp],[ϕ,φ,q])
end

function fsi_uvp_jacobian_t_Ωf(strategy::MeshStrategy{:biharmonic},x,xt,dxt,y,ρ)
    w, u, v, p = x
    dwt, dut, dvt, dpt = dxt
    ψ, ϕ, φ, q = y
    da_uvp_dut_φ_Ωf([u,v,p],[dut,dvt,dpt],[ϕ,φ,q],ρ) +
    da_uvp_dvt_φ_Ωf([u,v,p],[dut,dvt,dpt],[ϕ,φ,q],ρ)
end

fsi_uvp_residual_Ωs(strategy::MeshStrategy,x,xt,y,ρ,E,ν,vol) =
    a_uvp_ϕ_Ωs(x,xt,y) +
    a_uvp_φ_Ωs(x,xt,y,ρ,E,ν)

fsi_uvp_jacobian_Ωs(strategy::MeshStrategy,x,xt,dx,y,ρ,E,ν,vol) =
    da_uvp_dx_ϕ_Ωs(x,dx,y) +
    da_uvp_dx_φ_Ωs(x,dx,y,ρ,E,ν) 

fsi_uvp_jacobian_t_Ωs(strategy::MeshStrategy,x,xt,dxt,y,ρ) =
    da_uvp_dxt_Ωs(x,dxt,y,ρ)

function fsi_uvp_residual_Ωs(strategy::MeshStrategy{:biharmonic},x,xt,y,ρ,E,ν,vol)
    w, u, v, p = x
    wt, ut, vt, pt = xt
    ψ, ϕ, φ, q = y
    a_wuvp_ψ_Ωf(strategy,x,y,vol) +
    a_uvp_ϕ_Ωs([u,v,p],[ut,vt,pt],[ϕ,φ,q]) +
    a_uvp_φ_Ωs([u,v,p],[ut,vt,pt],[ϕ,φ,q],ρ,E,ν)
end

function fsi_uvp_jacobian_Ωs(strategy::MeshStrategy{:biharmonic},x,xt,dx,y,ρ,E,ν,vol)
    w, u, v, p = x
    dw, du, dv, dp = dx
    ψ, ϕ, φ, q = y
    a_wuvp_ψ_Ωf(strategy,dx,y,vol) +
    da_uvp_dx_ϕ_Ωs([u,v,p],[du,dv,dp],[ϕ,φ,q]) +
    da_uvp_dx_φ_Ωs([u,v,p],[du,dv,dp],[ϕ,φ,q],ρ,E,ν)
end

function fsi_uvp_jacobian_t_Ωs(strategy::MeshStrategy{:biharmonic},x,xt,dxt,y,ρ)
    w, u, v, p = x
    dwt, dut, dvt, dpt = dxt
    ψ, ϕ, φ, q = y
    da_uvp_dxt_Ωs([u,v,p],[dut,dvt,dpt],[ϕ,φ,q],ρ)
end

fsi_uvp_residual_Γi(strategy::MeshStrategy,x,y,n,E,ν,vol) = a_uvp_ϕ_Γi(strategy,x,y,n,E,ν)
fsi_uvp_jacobian_Γi(strategy::MeshStrategy,x,dx,y,n,E,ν,vol) = da_uvp_du_ϕ_Γi(strategy,x,dx,y,n,E,ν)
fsi_uvp_residual_Γi(strategy::MeshStrategy{:biharmonic},x,y,n,E,ν,vol) = a_wuvp_ϕ_Γi(strategy,x,y,n,vol)#a_wuvp_ψ_Γi(strategy,x,y,n,vol) + a_wuvp_ϕ_Γi(strategy,x,y,n,vol)
fsi_uvp_jacobian_Γi(strategy::MeshStrategy{:biharmonic},x,dx,y,n,E,ν,vol) = a_wuvp_ϕ_Γi(strategy,dx,y,n,vol)#a_wuvp_ψ_Γi(strategy,dx,y,n,vol) + a_wuvp_ϕ_Γi(strategy,dx,y,n,vol)

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


