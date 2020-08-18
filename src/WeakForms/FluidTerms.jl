# Stokes
# ======
function a_ST_vp(x,y,μ)
  v, p = x
  φ, q = y
  (ε(φ) ⊙ σ_dev(μ,ε(v))) - ((∇⋅φ) * p) + (q * (∇⋅v))
end
function l_ST_vp(y,μ,f)
  φ, q = y
  (φ⋅f)
end

# FSI
# ===
# Residuals
function a_FSI_ϕ_Ωf(strategy::MeshStrategy{:laplacian},x,y,α)
  u, v, p = x
  ϕ, φ, q = y
  α * (∇(ϕ) ⊙ ∇(u))
end
# Residuals
function a_FSI_ϕ_Ωf(strategy::MeshStrategy{:linearElasticity},x,y,E,ν)
  u, v, p = x
  ϕ, φ, q = y
  (λ_m,μ_m) = lame_parameters(E,ν)
  α(u) = 1.0#α_m(J(∇(u)))
  λ(u) = α(u)*λ_m
  μ(u) = α(u)*μ_m
  (ε(ϕ) ⊙ σ_m(λ(u),μ(u),ε(u)))
end
function a_FSI_ϕ_Ωf(strategy::MeshStrategy{:neoHookean},x,y,E,ν)
  u, v, p = x
  ϕ, φ, q = y
  (λ_m,μ_m) = lame_parameters(E,ν)
  (dE(∇(ϕ),∇(u)) ⊙ S_NH(∇(u),λ_m,μ_m))
end
function a_FSI_ψ_Ωf(strategy::MeshStrategy{:biharmonic},x,y,α)
  w, u, v, p = x
  ψ, ϕ, φ, q = y
  α * ( (ψ ⋅ w) + (∇(ψ) ⊙ ∇(u)) )
end
function a_FSI_ϕ_Ωf(strategy::MeshStrategy{:biharmonic},x,y,α)
  w, u, v, p = x
  ψ, ϕ, φ, q = y
  α * (∇(ϕ) ⊙ ∇(w))
end
function l_FSI_ϕ_Ωf(strategy::MeshStrategy,y,f,t)
  ϕ, φ, q = y
  (ϕ⋅f(t))
end
function l_FSI_ψ_Ωf(strategy::MeshStrategy{:biharmonic},y,f,t)
  ψ, ϕ, φ, q = y
  (ψ⋅f(t))
end
function a_FSI_φ_Ωf(x, xt, y, μ, ρ)
  u, v, p = x
  ut, vt, pt = xt
  ϕ, φ, q = y
  temp_term = φ ⋅ ( J(∇(u)) * ρ * vt )
  conv_term = φ ⋅ ( J(∇(u)) * ρ * conv(	Finv(∇(u))⋅(v-ut), ∇(v)) )
  visc_term = ( ∇(φ) ⊙ Pf_dev(μ,u,v) )
  pres_term = (∇⋅φ) * Pf_vol(u,p)
  temp_term + conv_term + visc_term + pres_term
end
function l_FSI_φ_Ωf(y,f,t)
  ϕ, φ, q = y
  (φ⋅f(t))
end
function a_FSI_q_Ωf(x, y)
  u, v, p = x
  ϕ, φ, q = y
  q * (J(∇(u))*(∇(v)⊙FinvT(∇(u))))
end
# Jacobians
function da_FSI_du_ϕ_Ωf(strategy::MeshStrategy{:linearElasticity},x, dx, y, E, ν)
  u, v, p = x
  du, dv, dp = dx
  ϕ, φ, q = y
  (λ_m,μ_m) = lame_parameters(E,ν)
  α(u) = 1.0#α_m(J(∇(u)))
  λ(u) = α(u)*λ_m
  μ(u) = α(u)*μ_m
  dα(u,du) = 0.0#dα_m(J(∇(u)),dJ(∇(u),∇(du)))
  dλ(u,du) = dα(u,du)*λ_m
  dμ(u,du) = dα(u,du)*μ_m
  (ε(ϕ) ⊙ dσ_m(λ(u),dλ(u,du),μ(u),dμ(u,du),ε(u),ε(du)))
end
function da_FSI_du_ϕ_Ωf(strategy::MeshStrategy{:neoHookean},x, dx, y, E, ν)
  u, v, p = x
  du, dv, dp = dx
  ϕ, φ, q = y
  (λ_m,μ_m) = lame_parameters(E,ν)
  ( dE(∇(ϕ),∇(u)) ⊙ dS_NH(∇(u),∇(du),λ_m,μ_m) ) + ( ∇(ϕ) ⊙ ( S_NH(∇(u),λ_m,μ_m)⋅∇(du) ) )
end
function da_FSI_dx_ψ_Ωf(strategy::MeshStrategy{:biharmonic},x,dx,y,α)
  w, u, v, p = x
  dw, du, dv, dp = dx
  ψ, ϕ, φ, q = y
  α * (  (ψ⋅dw) + (∇(ψ) ⊙ ∇(du)) )
end
function da_FSI_dx_ϕ_Ωf(strategy::MeshStrategy{:biharmonic},x,dx,y,α)
  w, u, v, p = x
  dw, du, dv, dp = dx
  ψ, ϕ, φ, q = y
  α * ( (∇(ϕ) ⊙ ∇(dw)) )
end
function da_FSI_du_φ_Ωf(x, xt, dx, y, μ, ρ)
  u, v, p = x
  ut, vt, pt = xt
  du, dv, dp = dx
  ϕ, φ, q = y
  temp_term = φ ⋅ ( dJ(∇(u),∇(du)) * ρ * vt )
  conv_term = φ ⋅ ( ( dJ(∇(u),∇(du)) * ρ * conv( Finv(∇(u))⋅(v-ut), ∇(v)) ) +
                    ( J(∇(u)) * ρ * conv(	dFinv(∇(u),∇(du))⋅(v-ut), ∇(v)) ) )
  visc_term = ∇(φ) ⊙ dPf_dev_du(μ,u,du,v)
  pres_term = (∇⋅φ) * dPf_vol_du(u,du,p)
  temp_term + conv_term + visc_term + pres_term
end
function da_FSI_dv_φ_Ωf(x, xt, dx, y, μ, ρ)
  u, v, p = x
  ut, vt, pt = xt
  du, dv, dp = dx
  ϕ, φ, q = y
  conv_term = φ ⋅ ( J(∇(u)) * ρ * dconv( Finv(∇(u))⋅dv, ∇(dv), Finv(∇(u))⋅(v-ut) , ∇(v)) )
  visc_term = ( ∇(φ) ⊙ Pf_dev(μ,u,dv) )
  conv_term + visc_term
end
function da_FSI_dp_φ_Ωf(x, dx, y)
  u, v, p = x
  du, dv, dp = dx
  ϕ, φ, q = y
  pres_term = (∇⋅φ) * Pf_vol(u,dp)
end
function da_FSI_dut_φ_Ωf(x, dxt, y, ρ)
  u, v, p = x
  dut, dvt, dpt = dxt
  ϕ, φ, q = y
  conv_term = - φ ⋅ ( J(∇(u)) * ρ * conv(	Finv(∇(u))⋅dut, ∇(v)) )
end

function da_FSI_dvt_φ_Ωf(x, dxt, y, ρ)
  u, v, p = x
  dut, dvt, dpt = dxt
  ϕ, φ, q = y
  temp_term = φ ⋅ ( J(∇(u)) * ρ * dvt )
end
function da_FSI_du_q_Ωf(x, dx, y)
  u, v, p = x
  du, dv, dp = dx
  ϕ, φ, q = y
  q * ( dJ(∇(u),∇(du))*(∇(v)⊙FinvT(∇(u))) + J(∇(u))*(∇(v)⊙dFinvT(∇(u),∇(du))) )
end
function da_FSI_dv_q_Ωf(x, dx, y)
  u, v, p = x
  du, dv, dp = dx
  ϕ, φ, q = y
  q * ( J(∇(u))*(∇(dv)⊙FinvT(∇(u))) )
end
