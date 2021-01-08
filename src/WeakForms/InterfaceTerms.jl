# Fluid-Solid Interface
# =====================
function a_FSI_weak_Γi((uf_Γ,vf_Γ,pf_Γ,us_Γ,vs_Γ),(ϕf_Γ,φf_Γ,qf_Γ,ϕs_Γ,φs_Γ),n,μ,γ,h,dt,dΓ)
  uf, vf, pf, us, vs = uf_Γ.⁺, vf_Γ.⁺, pf_Γ.⁺, us_Γ.⁻, vs_Γ.⁻
  ϕf, φf, qf, ϕs, φs = ϕf_Γ.⁺, φf_Γ.⁺, qf_Γ.⁺, ϕs_Γ.⁻, φs_Γ.⁻
  ∫( γ*μ/h/dt*(ϕf-ϕs)⋅(uf-us) + γ*μ/h*(φf-φs)⋅(vf-vs) +
    - (φf-φs) ⋅ ( n.⁺⋅Pᵥ_Ωf(μ,uf,vf) + n.⁺*Pₚ_Ωf(uf,pf) ) +
    (n.⁺⋅Pᵥ_Ωf(μ,uf,φf)) ⋅ (vf-vs) + (n.⁺*Pₚ_Ωf(uf,qf)) ⋅ (vf-vs) )dΓ
end
function da_FSI_weak_Γi_dx(x,dx,y,n,μ,γ,h,dt,dΓ)
  uf_Γ, vf_Γ, pf_Γ, us_Γ, vs_Γ = x
  duf_Γ, dvf_Γ, dpf_Γ, dus_Γ, dvs_Γ = dx
  ϕf_Γ, φf_Γ, qf_Γ, ϕs_Γ, φs_Γ = y
  uf, vf, pf, us, vs = uf_Γ.⁺, vf_Γ.⁺, pf_Γ.⁺, us_Γ.⁻, vs_Γ.⁻
  duf, dvf, dpf, dus, dvs = duf_Γ.⁺, dvf_Γ.⁺, dpf_Γ.⁺, dus_Γ.⁻, dvs_Γ.⁻
  ϕf, φf, qf, ϕs, φs = ϕf_Γ.⁺, φf_Γ.⁺, qf_Γ.⁺, ϕs_Γ.⁻, φs_Γ.⁻
  dP_tensor(uf,duf,vf,dvf,pf) = dPᵥ_Ωf_du(μ,uf,duf,vf) + Pᵥ_Ωf_dv(μ,uf,dvf) + dPₚ_Ωf_du(uf,duf,pf)
  dP_scalar(uf,dpf) = dPₚ_Ωf_vol(uf,dpf)
  ∫( γ*μ/h/dt*(ϕf-ϕs)⋅(duf-dus) + γ*μ/h*(φf-φs)⋅(dvf-dvs) +
     - (φf-φs) ⋅ ( n.⁺⋅dP_tensor(uf,duf,vf,dvf,pf) + n.⁺*dP_scalar(uf,dpf) ) +
     (n⋅Pᵥ_Ωf(μ,uf,φf)) ⋅ (dvf-dvs) + (n⋅dPᵥ_Ω_du(μ,uf,duf,φf)) ⋅ (vf-vs) +
     (n*Pₚ_Ωf(uf,qf)) ⋅ (dvf-dvs) + (n⋅dPₚ_Ωf_du(uf,duf,qf)) ⋅ (vf-vs) )dΓ
end
