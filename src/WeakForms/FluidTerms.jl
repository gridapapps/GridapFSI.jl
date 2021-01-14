# Stokes
# ======
#a_ST((v,p),(φ,q),μ,dΩ) = ∫( ε(φ) ⊙ σᵥ_Ωf(μ,v) - (∇⋅φ) * p + q * (∇⋅v) )dΩ
a_ST((v,p),(φ,q),μ,dΩ) = ∫( ∇(φ) ⊙ ∇(v) - (∇⋅φ) * p + q * (∇⋅v) )dΩ
l_ST((φ,q),f,dΩ) = ∫( φ⋅f )dΩ


# Navier-Stokes
# =============
function a_NS((v,p),(vt,pt),(φ,q),μ,ρ,dΩ)
  ∫( ρ*(ϕ ⋅ vt) + ρ*( conv∘(v,∇(v))) + ε(φ) ⊙ σᵥ_Ωf(μ,v) - (∇⋅φ) * p + q * (∇⋅v) )dΩ
end
function da_NS_dx((v,p),(dv,dp),(φ,q),μ,ρ,dΩ)
  ∫( 0.0*ρ*(ϕ ⋅ dv) + ρ*( dconv∘(dv,∇(dv),v,∇(v))) + ε(φ) ⊙ (σᵥ_Ωf(μ)∘(ε(dv))) - (∇⋅φ) * dp + q * (∇⋅dv) )dΩ
end
function da_NS_dxt((dvt,dpt),(φ,q),ρ,dΩ)
  ∫( ρ*(ϕ ⋅ dvt) )dΩ
end
l_NS((φ,q),f,dΩ) = ∫( φ⋅f )dΩ


# Navier-Stokes ALE
# =================
# LHS terms
function a_NS_ALE((u,v,p),(ut,vt,pt),(ϕ,φ,q),μ,ρ,dΩ)
  ∫( φ ⋅ ( J∘∇(u) * ρ * vt ) +
     φ ⋅ ( J∘∇(u) * ρ * conv∘( (Finv∘∇(u))⋅(v-ut), ∇(v) ) ) +
     ∇(φ) ⊙ Pᵥ_Ωf(μ,u,v) +
     (∇⋅φ) * Pₚ_Ωf(u,p) +
     q * (J∘∇(u) * (∇(v) ⊙ FinvT∘∇(u))) )dΩ
end
function a_NS_ALE_ΓD((u,v,p),(ut,vt,pt),(ϕ,φ,q),t,vD,n,μ,γ,h,dΓ)
  ∫( 0.0*(ϕ⋅u) + γ*μ/h*(ϕ⋅(ut-vD(t))) + γ*μ/h*(φ⋅(v-vD(t))) +
    -(φ ⋅ (n⋅Pᵥ_Ωf(μ,u,v))) - (φ ⋅ n)*Pₚ_Ωf(u,p) +
    (n⋅Pᵥ_Ωf(μ,u,φ)) ⋅ (v-vD(t)) + (n*Pₚ_Ωf(u,q)) ⋅ (v-vD(t)) )dΓ
end

# LHS linearized terms
function da_NS_ALE_dx((u,v,p),(ut,vt,pt),(du,dv,dp),(ϕ,φ,q),μ,ρ,dΩ)
  ∫( φ ⋅ ( dJ∘(∇(u),∇(du)) * ρ * vt ) +
     φ ⋅ ( ( dJ∘(∇(u),∇(du)) * ρ * conv∘( Finv∘∇(u)⋅(v-ut), ∇(v)) ) +
           ( J∘∇(u) * ρ * conv∘( dFinv∘(∇(u),∇(du))⋅(v-ut), ∇(v)) ) ) +
     ∇(φ) ⊙ dPᵥ_Ωf_du(μ,u,du,v) +
     (∇⋅φ) * dPₚ_Ωf_du(u,du,p) +
     φ ⋅ ( J∘∇(u) * ρ * dconv( (Finv∘∇(u))⋅dv, ∇(dv), (Finv∘∇(u))⋅(v-ut) , ∇(v)) ) +
     ∇(φ) ⊙ Pᵥ_Ωf(μ,u,dv) +
     (∇⋅φ) * Pₚ_Ωf_dp(u,dp) +
     q * ( (dJ∘(∇(u),∇(du))) * (∇(v) ⊙ (FinvT∘∇(u))) +
           (J∘∇(u)) * (∇(v) ⊙ (dFinvT∘(∇(u),∇(du)))) ) +
     q * ( (J∘∇(u)) * (∇(dv) ⊙ (FinvT∘∇(u))) ) )dΩ
end
function da_NS_ALE_dxt((u,v,p),(dut,dvt,dpt),(ϕ,φ,q),ρ,dΩ)
  ∫( - φ ⋅ ( (J∘∇(u)) * ρ * conv∘(	(Finv∘∇(u))⋅dut, ∇(v)) ) +
     φ ⋅ ( (J∘∇(u)) * ρ * dvt ) )dΩ
end
function da_NS_ALE_ΓD_dx((u,v,p),(ut,vt,pt),(du,dv,dp),(ϕ,φ,q),n,μ,γ,h,dΓ)
  dP_tensor(u,du,v,dv,p) = dPᵥ_Ωf_du(μ,u,du,v) + Pᵥ_Ωf_dv(μ,u,dv) + dPₚ_Ωf_du(u,du,p)
  dP_scalar(u,dp) = Pₚ_Ωf_dp(u,dp)
  ∫( 0.0*(ϕ⋅du) + γ*μ/h*(φ⋅dv) +
     - φ ⋅ ( n⋅dP_tensor(u,du,v,dv,p) + n*dP_scalar(u,dp) ) +
     (n⋅Pᵥ_Ωf(μ,u,φ)) ⋅ dv + (n⋅dPᵥ_Ωf_du(μ,u,du,φ)) ⋅ v +
     (n*Pₚ_Ωf(u,q)) ⋅ dv + (n⋅dPₚ_Ωf_du(u,du,q)) ⋅ v )dΓ
end
function da_NS_ALE_ΓD_dxt((dut,dvt,dpt),(ϕ,φ,q),μ,γ,h,dΓ)
  ∫( γ*μ/h*(ϕ⋅dut) )dΓ
end

# RHS terms
function l_NS_ALE((ϕ,φ,q),f,t,dΩ)
  ∫( φ⋅f(t) )dΩ
end

# Monolithic mesh motion
# ======================
function a_mesh(strategy::MeshStrategy{:laplacian},(u,v,p),(ϕ,φ,q),α,dΩ)
  ∫( α * (∇(ϕ) ⊙ ∇(u)) )dΩ
end
function a_mesh(strategy::MeshStrategy{:linearElasticity},(u,v,p),(ϕ,φ,q),λ,μ,dΩ)
  ∫( ε(ϕ) ⊙ (σₘ(λ,μ)∘ε(u)) )dΩ
end
function a_mesh(strategy::MeshStrategy{:biharmonic},(w,u,v,p),(ψ,ϕ,φ,q),α₁,α₂,dΩ)
  ∫( α₁ * ( - (ψ ⋅ w) + (∇(ψ) ⊙ ∇(u)) ) + α₂ * ( (∇(ϕ) ⊙ ∇(w)) ) )dΩ
end
function a_mesh_Γi(strategy::MeshStrategy{:laplacian},(u,v,p),(ϕ,φ,q),n,α,dΓ)
  ∫( - α * (ϕ ⋅ (n.⁺⋅∇(u))) )dΓ
end
function a_mesh_Γi(strategy::MeshStrategy{:linearElasticity},(u,v,p),(ϕ,φ,q),n,λ,μ,dΓ)
  ∫( - (ϕ ⋅  (n.⁺⋅σₘ(λ,μ)∘ε(u)) ) )dΓ
end
function a_mesh_Γi(strategy::MeshStrategy{:biharmonic},(w,u,v,p),(ψ,ϕ,φ,q),n,α₁,α₂,dΓ)
  ∫( - α₁ * (ψ ⋅  (n.⁺⋅∇(u))) - α₂ * (ϕ ⋅  (n.⁺⋅∇(w))) )dΓ
end
function l_mesh(strategy::MeshStrategy,(ϕ,φ,q),f,t,dΩ)
  ∫( ϕ⋅f(t) )dΩ
end
function l_mesh(strategy::MeshStrategy{:biharmonic},(ψ,ϕ,φ,q),f,t,dΩ)
  ∫( ψ⋅f(t) )dΩ
end
