# Primal Finite Elasticity
# ========================
function a_PFE((u,v),(ut,vt),(ϕ,φ),λ,μ,dΩ)
  ∫( ϕ⋅ut - ϕ⋅v +
     φ⋅(ρ*vt) + (∇(φ) ⊙ Pₛᵥ_Ωs(λ,μ,u)) )dΩ
end
function da_PFE_dx((u,v),(du,dv),(ϕ,φ),ρ,λ,μ,dΩ)
  ∫( 0.0*(du⋅ϕ) - (ϕ⋅dv) +
     0.0*(φ⋅(ρ*dv)) + (∇(φ) ⊙ dPₛᵥ_Ωs_du(λ,μ,u)) )dΩ
end
function da_PFE_dxt((dut,dvt),(ϕ,φ),ρ,dΩ)
  ∫( dut⋅ϕ + φ⋅(ρ*dvt) )dΩ
end
function l_PFE((ϕ,φ),f₁,f₂,t,dΩ)
  ∫( ϕ⋅f₁(t) + φ⋅f₂(t) )dΩ
end
