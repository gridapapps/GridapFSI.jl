function a_FSI_ϕ_Ωs(x, xt, y)
  u,v,p = x
  ut, vt,pt = xt
  ϕ,φ,q = y
  (ϕ⋅ut) - (ϕ⋅v)
end
function l_FSI_ϕ_Ωs(y,f,t)
  ϕ,φ,q = y
  (ϕ⋅f(t))
end
function a_FSI_φ_Ωs(x, xt, y, ρ, Es, νs)
  u,v,p = x
  ut, vt,pt = xt
  ϕ,φ,q = y
  (λ,μ) = lame_parameters(Es,νs)
  (φ⋅(ρ*vt)) + (∇(φ) ⊙ Ps(λ,μ,u))
end
function l_FSI_φ_Ωs(y,f,t)
  ϕ,φ,q = y
  (φ⋅f(t))
end


function da_FSI_dx_ϕ_Ωs(x, dx, y)
  u,v,p = x
  du,dv,dp = dx
  ϕ,φ,q = y
  0.0*(du⋅ϕ) - (ϕ⋅dv)
end
function da_FSI_dx_φ_Ωs(x, dx, y, ρ, E, ν)
  u,v,p = x
  du,dv,dp = dx
  ϕ,φ,q = y
  (λ,μ) = lame_parameters(E,ν)
  0.0*(φ⋅(ρ*dv)) + (∇(φ) ⊙ ( dF(∇(du))⋅S_SV(∇(u),λ,μ) ) ) + ( ∇(φ) ⊙ (F(∇(u))⋅dS_SV(∇(u),∇(du),λ,μ)) )
end
function da_FSI_dxt_Ωs(x, dxt, y, ρ)
  u, v, p = x
  dut, dvt, dpt = dxt
  ϕ, φ, q = y
  ϕ⋅dut + (φ⋅(ρ*dvt))
end
