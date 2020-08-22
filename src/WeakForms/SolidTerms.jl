function a_FSI_ϕ_Ωs(x, xt, y)
  u,v = x
  ut, vt = xt
  ϕ,φ = y
  (ϕ⋅ut) - (ϕ⋅v)
end
function l_FSI_ϕ_Ωs(y,f,t)
  ϕ,φ = y
  (ϕ⋅f(t))
end
function a_FSI_φ_Ωs(x, xt, y, ρ, Es, νs)
  u,v = x
  ut, vt = xt
  ϕ,φ = y
  (λ,μ) = lame_parameters(Es,νs)
  (φ⋅(ρ*vt)) + (∇(φ) ⊙ Ps(λ,μ,u))
end
function l_FSI_φ_Ωs(y,f,t)
  ϕ,φ = y
  (φ⋅f(t))
end


function da_FSI_dx_ϕ_Ωs(x, dx, y)
  u,v = x
  du,dv = dx
  ϕ,φ = y
  0.0*(du⋅ϕ) - (ϕ⋅dv)
end
function da_FSI_dx_φ_Ωs(x, dx, y, ρ, E, ν)
  u,v = x
  du,dv = dx
  ϕ,φ = y
  (λ,μ) = lame_parameters(E,ν)
  0.0*(φ⋅(ρ*dv)) + (∇(φ) ⊙ ( dF(∇(du))⋅S_SV(∇(u),λ,μ) ) ) + ( ∇(φ) ⊙ (F(∇(u))⋅dS_SV(∇(u),∇(du),λ,μ)) )
end
function da_FSI_dxt_Ωs(x, dxt, y, ρ)
  u, v = x
  dut, dvt = dxt
  ϕ, φ = y
  ϕ⋅dut + (φ⋅(ρ*dvt))
end
