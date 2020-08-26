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
uf, vf, pf, us, vs = x
ϕf, φf, qf, ϕs, φs = y
p = pf.inward
q = qf.inward
γ*μ/h*(ϕf.inward⋅uf.inward) - γ*μ/h*(ϕf.inward⋅us.outward) - γ*μ/h*(ϕs.outward⋅uf.inward) + γ*μ/h*(ϕs.outward⋅uf.inward)
+ γ*μ/h*(φf.inward⋅vf.inward) - γ*μ/h*(φf.inward⋅vs.outward) - γ*μ/h*(φs.outward⋅vf.inward) + γ*μ/h*(φs.outward⋅vf.inward)
- φf.inward ⋅ (n⋅Pf_dev(μ,uf,vf).inward) + φs.outward ⋅ (n⋅Pf_dev(μ,uf,vf).inward)
+ φf.inward ⋅ (n*Pf_vol(uf,pf).inward) - φs.outward ⋅ (n*Pf_vol(uf,pf).inward)
+ ((n⋅Pf_dev(μ,uf,φf).inward) ⋅ vf.inward) - (n⋅Pf_dev(μ,uf,φf).inward) ⋅ vs.outward
- (n*Pf_vol(uf,qf).inward) ⋅ vf.inward + (n*Pf_vol(uf,qf).inward) ⋅ vs.outward
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
  uf, vf, pf, us, vs = x
  duf, dvf, dpf, dus, dvs = dx
  ϕf, φf, qf, ϕs, φs = y
  du_jump = duf.inward - dus.outward
  dv_jump = dvf.inward - dvs.outward
  dp = dpf.inward
  ϕ_jump = ϕf.inward - ϕs.outward
  φ_jump = φf.inward - φs.outward
  q = qf.inward
  (γ*μ/h*(ϕ_jump⋅du_jump))   + (γ*μ/h*(φ_jump⋅dv_jump))
  - (φ_jump ⋅ (n⋅dPf_dev(μ,uf.inward,duf.inward,vf.inward)))
  - (φ_jump ⋅ (n⋅Pf_dev(μ,uf.inward,dvf.inward))) + (φ_jump⋅(n*dp))
  - ((n⋅Pf_dev(μ,ϕf.inward,φf.inward)) ⋅ dv_jump) - ((n*q) ⋅ dv_jump)
  end
