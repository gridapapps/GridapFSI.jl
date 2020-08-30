function a_FSI_ϕ_Γi(strategy::MeshStrategy{:laplacian},x,y,n,α)
  u, v, p = x
  ϕ, φ, q = y
  - α * (ϕ ⋅ (n⋅∇(u)))
end
function a_FSI_ϕ_Γi(strategy::MeshStrategy{:linearElasticity},x,y,n,E,ν)
  u, v, p = x
  ϕ, φ, q = y
  (λ_m,μ_m) = lame_parameters(E,ν)
  α(u) = 1.0#α_m(J(∇(u)))
  λ(u) = α(u)*λ_m
  μ(u) = α(u)*μ_m
  - (ϕ ⋅  (n⋅σ_m(λ(u),μ(u),ε(u))) )
end
function a_FSI_ϕ_Γi(strategy::MeshStrategy{:neoHookean},x,y,n,E,ν)
  u, v, p = x
  ϕ, φ, q = y
  (λ_m,μ_m) = lame_parameters(E,ν)
  - (ϕ ⋅  (n ⋅ S_NH(∇(u),λ_m,μ_m) ) )
end
function a_FSI_ψ_Γi(strategy::MeshStrategy{:biharmonic},x,y,n,α)
  w, u, v, p = x
  ψ, ϕ, φ, q = y
  - α * (ψ ⋅  (n⋅∇(u)))
end
function a_FSI_ϕ_Γi(strategy::MeshStrategy{:biharmonic},x,y,n,α)
  w, u, v, p = x
  ψ, ϕ, φ, q = y
  - α * (ϕ ⋅  (n⋅∇(w)))
end
@law jump_law(af,as) = af.inward - as.outward
function a_FSI_Nitsche_ϕ_Γi(x,y,n,μ,γ,h)
uf_Γ, vf_Γ, pf_Γ, us_Γ, vs_Γ = x
ϕf_Γ, φf_Γ, qf_Γ, ϕs_Γ, φs_Γ = y
uf, vf, pf, us, vs = uf_Γ.⁺, vf_Γ.⁺, pf_Γ.⁺, us_Γ.⁻, vs_Γ.⁻
ϕf, φf, qf, ϕs, φs = ϕf_Γ.⁺, φf_Γ.⁺, qf_Γ.⁺, ϕs_Γ.⁻, φs_Γ.⁻
#γ*μ/h*(ϕf-ϕs)⋅(uf-us)
γ*μ/h*(φf-φs)⋅(vf-vs)
- ((φf-φs) ⋅ (n⋅Pf_dev(μ,uf_Γ,vf_Γ).⁺))
+ ((φf-φs) ⋅ (n⋅Pf_vol(uf_Γ,pf_Γ).⁺))
+ ((n⋅Pf_dev(μ,uf_Γ,φf_Γ).⁺) ⋅ (vf-vs))
- ((n⋅Pf_vol(uf_Γ,qf_Γ).⁺) ⋅ (vf-vs))
# γ*μ/h*(ϕf.inward⋅uf.inward) - γ*μ/h*(ϕf.inward⋅us.outward) - γ*μ/h*(ϕs.outward⋅uf.inward) + γ*μ/h*(ϕs.outward⋅uf.inward)
# + γ*μ/h*(φf.inward⋅vf.inward) - γ*μ/h*(φf.inward⋅vs.outward) - γ*μ/h*(φs.outward⋅vf.inward) + γ*μ/h*(φs.outward⋅vf.inward)
# - φf.inward ⋅ (n⋅Pf_dev(μ,uf,vf).inward) + φs.outward ⋅ (n⋅Pf_dev(μ,uf,vf).inward)
# + φf.inward ⋅ (n*Pf_vol(uf,pf).inward) - φs.outward ⋅ (n*Pf_vol(uf,pf).inward)
# + ((n⋅Pf_dev(μ,uf,φf).inward) ⋅ vf.inward) - (n⋅Pf_dev(μ,uf,φf).inward) ⋅ vs.outward
# - (n*Pf_vol(uf,qf).inward) ⋅ vf.inward + (n*Pf_vol(uf,qf).inward) ⋅ vs.outward
#(γ*μ/h*(jump_law(ϕf,ϕs)⋅jump_law(uf,us)))  + (γ*μ/h*(jump_law(φf,φs)⋅jump_law(vf,vs))
#- (jump_law(φf,φs) ⋅ (n⋅Pf_dev(μ,uf.inward,vf.inward))) + (jump_law(φf,φs)⋅(n*p))
#- ((n⋅Pf_dev(μ,ϕf.inward,φf.inward)) ⋅ jump_law(vf,vs)) - ((n*q) ⋅ jump_law(vf,vs))
end


function da_FSI_du_ϕ_Γi(strategy::MeshStrategy{:linearElasticity},x,dx,y, n, E, ν)
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
  - (ϕ ⋅  (n⋅dσ_m(λ(u),dλ(u,du),μ(u),dμ(u,du),ε(u),ε(du))) )
end
function da_FSI_du_ϕ_Γi(strategy::MeshStrategy{:neoHookean},x,dx,y, n, E, ν)
  u, v, p = x
  du, dv, dp = dx
  ϕ, φ, q = y
  (λ_m,μ_m) = lame_parameters(E,ν)
  - (ϕ ⋅  (n⋅dS_NH(∇(u),∇(du),λ_m,μ_m)) )
end
function da_FSI_dx_ψ_Γi(strategy::MeshStrategy{:biharmonic},x,dx,y,n,α)
  w, u, v, p = x
  dw, du, dv, dp = dx
  ψ, ϕ, φ, q = y
  - α * ( ψ ⋅  (n⋅∇(du)) )
end
function da_FSI_dx_ϕ_Γi(strategy::MeshStrategy{:biharmonic},x,dx,y,n,α)
  w, u, v, p = x
  dw, du, dv, dp = dx
  ψ, ϕ, φ, q = y
  - α * ( ϕ ⋅  (n⋅∇(dw)) )
end
function da_FSI_Nitsche_ϕ_Γi(x,dx,y,n,μ,γ,h)
  uf_Γ, vf_Γ, pf_Γ, us_Γ, vs_Γ = x
  duf_Γ, dvf_Γ, dpf_Γ, dus_Γ, dvs_Γ = dx
  ϕf_Γ, φf_Γ, qf_Γ, ϕs_Γ, φs_Γ = y
  uf, vf, pf, us, vs = uf_Γ.⁺, vf_Γ.⁺, pf_Γ.⁺, us_Γ.⁻, vs_Γ.⁻
  duf, dvf, dpf, dus, dvs = duf_Γ.⁺, dvf_Γ.⁺, dpf_Γ.⁺, dus_Γ.⁻, dvs_Γ.⁻
  ϕf, φf, qf, ϕs, φs = ϕf_Γ.⁺, φf_Γ.⁺, qf_Γ.⁺, ϕs_Γ.⁻, φs_Γ.⁻
  #γ*μ/h*(ϕf-ϕs)⋅(duf-dus)
  γ*μ/h*(φf-φs)⋅(dvf-dvs)
  - ((φf-φs) ⋅ (n⋅dPf_dev_du(μ,uf_Γ,duf_Γ,vf_Γ).⁺))
  - ((φf-φs) ⋅ (n⋅Pf_dev(μ,uf_Γ,dvf_Γ).⁺))
  + ((φf-φs) ⋅ (n⋅dPf_vol_du(uf_Γ,duf_Γ,pf_Γ).⁺))
  + ((φf-φs) ⋅ (n⋅Pf_vol(uf_Γ,dpf_Γ).⁺))
  + ((n⋅dPf_dev_du(μ,uf_Γ,duf_Γ,φf_Γ).⁺) ⋅ (vf-vs))
  + ((n⋅Pf_dev(μ,uf_Γ,φf_Γ).⁺) ⋅ (dvf-dvs))
  - ((n⋅dPf_vol_du(uf_Γ,duf_Γ,qf_Γ).⁺) ⋅ (vf-vs))
  - ((n⋅Pf_vol(uf_Γ,qf_Γ).⁺) ⋅ (dvf-dvs))
  end
