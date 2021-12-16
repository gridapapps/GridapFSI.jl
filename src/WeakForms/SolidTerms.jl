# Primal Finite Elasticity (Saint-Venant)
# ========================
function a_PFE((u,v),(ϕ,φ),ρ,λ,μ,dΩ)
  ∫( ϕ⋅∂t(u) - ϕ⋅v +
     φ⋅(ρ*∂t(v)) + (∇(φ) ⊙ Pₛᵥ_Ωs(λ,μ,u)) )dΩ
end
function da_PFE_dx((u,v),(du,dv),(ϕ,φ),ρ,λ,μ,dΩ)
  ∫( 0.0*(du⋅ϕ) - (ϕ⋅dv) +
     0.0*(φ⋅(ρ*dv)) + (∇(φ) ⊙ dPₛᵥ_Ωs_du(λ,μ,u,du)) )dΩ
end
function da_PFE_dxt((dut,dvt),(ϕ,φ),ρ,dΩ)
  ∫( dut⋅ϕ + φ⋅(ρ*dvt) )dΩ
end
function l_PFE((ϕ,φ),f₁,f₂,t,dΩ)
  ∫( ϕ⋅f₁(t) + φ⋅f₂(t) )dΩ
end
