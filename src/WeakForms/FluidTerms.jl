# Stokes
# ======
a_ST((v,p),(φ,q),μ,dΩ) = ∫( ε(φ) ⊙ (σᵥ_Ωf∘(μ,ε(v))) - (∇⋅φ) * p + q * (∇⋅v) )dΩ
l_ST((φ,q),f,dΩ) = ∫( φ⋅f )dΩ


# Navier-Stokes
# =============
function a_NS((v,p),(vt,pt),(φ,q),μ,ρ,dΩ)
  ∫( ρ*(ϕ ⋅ vt) + ρ*( conv∘(v,∇(v))) + ε(φ) ⊙ (σᵥ_Ωf∘(μ,ε(v))) - (∇⋅φ) * p + q * (∇⋅v) )dΩ
end
function da_NS_dx((v,p),(dv,dp),(φ,q),μ,ρ,dΩ)
  ∫( 0.0*ρ*(ϕ ⋅ dv) + ρ*( dconv∘(dv,∇(dv),v,∇(v))) + ε(φ) ⊙ (σᵥ_Ωf∘(μ,ε(dv))) - (∇⋅φ) * dp + q * (∇⋅dv) )dΩ
end
function da_NS_dxt((dvt,dpt),(φ,q),ρ,dΩ)
  ∫( ρ*(ϕ ⋅ dvt) )dΩ
end
l_NS((φ,q),f,dΩ) = ∫( φ⋅f )dΩ


# FSI
# ===
# LHS terms
function a_FSI_ϕ_Ωf(strategy::MeshStrategy{:laplacian},(u,v,p),(ϕ,φ,q),α,dΩ)
  ∫( α * (∇(ϕ) ⊙ ∇(u)) )dΩ
end
function a_FSI_ϕ_Ωf(strategy::MeshStrategy{:linearElasticity},(u,v,p),(ϕ,φ,q),λ,μ,dΩ)
  ∫( ε(ϕ) ⊙ (σₘ∘(λ,μ,ε(u))) )dΩ
end
function a_FSI_ψϕ_Ωf(strategy::MeshStrategy{:biharmonic},(w,u,v,p),(ψ,ϕ,φ,q),α,dΩ)
  ∫( α * ( (ψ ⋅ w) + (∇(ψ) ⊙ ∇(u)) + (∇(ϕ) ⊙ ∇(w)) ) )dΩ
end
function a_FSI_φ_Ωf((u,v,p),(ut,vt,pt),(ϕ,φ,q),μ,ρ,dΩ)
  ∫( φ ⋅ ( J∘∇(u) * ρ * vt ) +
     φ ⋅ ( J∘∇(u) * ρ * conv∘( (Finv∘∇(u))⋅(v-ut), ∇(v) ) ) +
     ∇(φ) ⊙ Pᵥ_Ωf(μ,u,v) +
     (∇⋅φ) * Pₚ_Ωf(u,p) )dΩ
end
function a_FSI_q_Ωf((u,v,p),(ϕ,φ,q),dΩ)
  ∫( q * (J∘∇(u) * (∇(v) ⊙ FinvT∘∇(u))) )dΩ
end
function a_FSI_ΓfD((v,p),(vt,pt),(φ,q),t,vD,n,μ,γ,h,dΓ)
  ∫( 0.0*(ϕ⋅u) + γ*μ/h*(ϕ⋅(ut-vD(t))) + γ*μ/h*(φ⋅(v-vD(t))) +
    -(φ ⋅ (n⋅Pᵥ_Ωf(μ,u,v))) - (φ ⋅ n)*Pₚ_Ωf(u,p) +
    (n⋅Pᵥ_Ωf(μ,u,φ)) ⋅ (v-vD(t)) + (n*Pₚ_Ωf(u,q)) ⋅ (v-vD(t)) )dΓ
end

# LHS linearized terms
function da_FSI_φ_Ωf_du((u,v,p),(ut,vt,pt),(du,dv,dp),(ϕ,φ,q),μ,ρ,dΩ)
  ∫( φ ⋅ ( dJ∘(∇(u),∇(du)) * ρ * vt ) +
     φ ⋅ ( ( dJ∘(∇(u),∇(du)) * ρ * conv∘( Finv∘∇(u)⋅(v-ut), ∇(v)) ) +
           ( J∘∇(u) * ρ * conv∘( dFinv∘(∇(u),∇(du))⋅(v-ut), ∇(v)) ) ) +
     ∇(φ) ⊙ dPᵥ_Ωf_du(μ,u,du,v) +
     (∇⋅φ) * dPₚ_Ωf_du(u,du,p) )dΩ
end
function da_FSI_φ_Ωf_dv((u,v,p),(ut,vt,pt),(du,dv,dp),(ϕ,φ,q),μ,ρ,dΩ)
  ∫( φ ⋅ ( J∘∇(u) * ρ * dconv( (Finv∘∇(u))⋅dv, ∇(dv), (Finv∘∇(u))⋅(v-ut) , ∇(v)) ) +
     ∇(φ) ⊙ Pᵥ_Ωf(μ,u,dv) )dΩ
end
function da_FSI_φ_Ωf_dp((u,v,p),(du,dv,dp),(ϕ,φ,q),dΩ)
  ∫( (∇⋅φ) * Pₚ_Ωf_dp(u,dp) )dΩ
end
function da_FSI_φ_Ωf_dut((u,v,p),(dut,dvt,dpt),(ϕ,φ,q),ρ,dΩ)
  ∫( - φ ⋅ ( (J∘∇(u)) * ρ * conv∘(	(Finv∘∇(u))⋅dut, ∇(v)) ) )dΩ
end
function da_FSI_φ_Ωf_dvt((u,v,p),(dut,dvt,dpt),(ϕ,φ,q),ρ,dΩ)
  ∫( φ ⋅ ( (J∘∇(u)) * ρ * dvt ) )dΩ
end
function da_FSI_q_Ωf_du((u,v,p),(du,dv,dp),(ϕ,φ,q),dΩ)
  ∫( q * ( (dJ∘(∇(u),∇(du))) * (∇(v) ⊙ (FinvT∘∇(u))) + (J∘∇(u)) * (∇(v) ⊙ (dFinvT∘(∇(u),∇(du)))) ) )dΩ
end
function da_FSI_q_Ωf_dv((u,v,p),(du,dv,dp),(ϕ,φ,q),dΩ)
  ∫( q * ( (J∘∇(u)) * (∇(dv) ⊙ (FinvT∘∇(u))) ) )dΩ
end
function da_FSI_ΓfD_dx((u,v,p),(ut,vt,pt),(du,dv,dp),(ϕ,φ,q),n,μ,γ,h,dΓ)
  dP_tensor(u,du,v,dv,p) = dPᵥ_Ωf_du(μ,u,du,v) + Pᵥ_Ωf_dv(μ,u,dv) + dPₚ_Ωf_du(u,du,p)
  dP_scalar(u,dp) = Pₚ_Ωf_dp(u,dp)
  ∫( 0.0*(ϕ⋅du) + γ*μ/h*(φ⋅dv) +
     - φ ⋅ ( n⋅dP_tensor(u,du,v,dv,p) + n*dP_scalar(u,dp) ) +
     (n⋅Pᵥ_Ωf(μ,u,φ)) ⋅ dv + (n⋅dPᵥ_Ωf_du(μ,u,du,φ)) ⋅ v +
     (n*Pₚ_Ωf(u,q)) ⋅ dv + (n⋅dPₚ_Ωf_du(u,du,q)) ⋅ v )dΓ
end
function da_FSI_ΓfD_dxt((dut,dvt,dpt),(ϕ,φ,q),μ,γ,h,dΓ)
  ∫( γ*μ/h*(ϕ⋅dut) )dΓ
end

# RHS terms
function l_FSI_ϕ_Ωf(strategy::MeshStrategy,(ϕ,φ,q),f,t,dΩ)
  ∫( ϕ⋅f(t) )dΩ
end
function l_FSI_ψ_Ωf(strategy::MeshStrategy{:biharmonic},(ψ,ϕ,φ,q),f,t,dΩ)
  ∫( ψ⋅f(t) )dΩ
end
function l_FSI_φ_Ωf((ϕ,φ,q),f,t,dΩ)
  ∫( φ⋅f(t) )dΩ
end
