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
function a_FSI_Nitsche_ϕ_Γi(x,y,n,γ,h)
uf, vf, pf, us, vs = x
ϕf, φf, qf, ϕs, φs = y
u_jump = uf.inward - us.outward
v_jump = vf.inward - vs.outward
p = pf.inward
q = qf.inward
ϕ_jump = ϕf.inward - ϕs.outward
φ_jump = φf.inward - φs.outward
(γ*μ/h*(ϕ_jump⋅u_jump))   + (γ*μ/h*(φ_jump⋅v_jump))
- (φ_jump ⋅ (n⋅Pf_dev(μ,uf.inward,vf.inward))) + (φ_jump⋅(n*p))
- ((n⋅Pf_dev(μ,ϕf.inward,φf.inward)) ⋅ v_jump) - ((n*q) ⋅ v_jump)
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
