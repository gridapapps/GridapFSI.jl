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
function a_FSI_Nitsche_ϕ_Γi(x,y,n,μ,γ,h,dt)
uf_Γ, vf_Γ, pf_Γ, us_Γ, vs_Γ = x
ϕf_Γ, φf_Γ, qf_Γ, ϕs_Γ, φs_Γ = y
uf, vf, pf, us, vs = uf_Γ.⁺, vf_Γ.⁺, pf_Γ.⁺, us_Γ.⁻, vs_Γ.⁻
ϕf, φf, qf, ϕs, φs = ϕf_Γ.⁺, φf_Γ.⁺, qf_Γ.⁺, ϕs_Γ.⁻, φs_Γ.⁻
penalty = γ*μ/h/dt*(ϕf-ϕs)⋅(uf-us) + γ*μ/h*(φf-φs)⋅(vf-vs)
integration_by_parts = - (φf-φs) ⋅ ( n⋅Pf_dev(μ,uf,vf) + n*Pf_vol(uf,pf) )
nitsche = (n⋅Pf_dev(μ,uf,φf)) ⋅ (vf-vs) + (n*Pf_vol(uf,qf)) ⋅ (vf-vs)
penalty + integration_by_parts + nitsche
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
function da_FSI_Nitsche_ϕ_Γi(x,dx,y,n,μ,γ,h,dt)
  uf_Γ, vf_Γ, pf_Γ, us_Γ, vs_Γ = x
  duf_Γ, dvf_Γ, dpf_Γ, dus_Γ, dvs_Γ = dx
  ϕf_Γ, φf_Γ, qf_Γ, ϕs_Γ, φs_Γ = y
  uf, vf, pf, us, vs = uf_Γ.⁺, vf_Γ.⁺, pf_Γ.⁺, us_Γ.⁻, vs_Γ.⁻
  duf, dvf, dpf, dus, dvs = duf_Γ.⁺, dvf_Γ.⁺, dpf_Γ.⁺, dus_Γ.⁻, dvs_Γ.⁻
  ϕf, φf, qf, ϕs, φs = ϕf_Γ.⁺, φf_Γ.⁺, qf_Γ.⁺, ϕs_Γ.⁻, φs_Γ.⁻
  dP_tensor(uf,duf,vf,dvf,pf) = dPf_dev_du(μ,uf,duf,vf) + Pf_dev(μ,uf,dvf) + dPf_vol_du(uf,duf,pf)
  dP_scalar(uf,dpf) = Pf_vol(uf,dpf)
  penalty = γ*μ/h/dt*(ϕf-ϕs)⋅(duf-dus) + γ*μ/h*(φf-φs)⋅(dvf-dvs)
  integration_by_parts = - (φf-φs) ⋅ ( n⋅dP_tensor(uf,duf,vf,dvf,pf) + n*dP_scalar(uf,dpf) )
  nitsche = (n⋅Pf_dev(μ,uf,φf)) ⋅ (dvf-dvs) + (n⋅dPf_dev_du(μ,uf,duf,φf)) ⋅ (vf-vs) + (n*Pf_vol(uf,qf)) ⋅ (dvf-dvs) + (n⋅dPf_vol_du(uf,duf,qf)) ⋅ (vf-vs)
  penalty + integration_by_parts + nitsche
  end
