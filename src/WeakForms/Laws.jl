# Map operations
I(A) = one(A)
F(∇u) = ∇u + I(∇u)
J(∇u) = det(F(∇u))
Finv(∇u) = inv(F(∇u))
FinvT(∇u) = (Finv(∇u)')
E(∇u) = 0.5 * ((F(∇u)')⋅F(∇u) - I(∇u))

# Map operations derivatives
dF(∇du) = ∇du
dJ(∇u,∇du) = J(∇u)*tr(Finv(∇u)⋅dF(∇du))
dFinv(∇u,∇du) = -Finv(∇u) ⋅ dF(∇du) ⋅ Finv(∇u)
dFinvT(∇u,∇du) = (dFinv(∇u,∇du)')
dE(∇u,∇du) = 0.5 * ((dF(∇du)')⋅F(∇u) + (F(∇u)')⋅dF(∇du))

# Fluid constitutive laws
# Cauchy stress
σᵥ_Ωf(μ,u) = 2.0*μ*ε(u)
σᵥ_Ωf(μ) = (∇v,Finv) -> μ*(∇v⋅Finv + (Finv')⋅(∇v'))
# First Piola-Kirchhoff stress
Pᵥ_Ωf(μ,u,v) = (J∘∇(u)) * (σᵥ_Ωf(μ)∘(∇(v),Finv∘∇(u)))' ⋅ (FinvT∘∇(u))
function Pₚ_Ωf(u,p)
  typeof(- (J∘∇(u)) * p * tr((FinvT∘∇(u))))
  - (J∘∇(u)) * p * tr((FinvT∘∇(u)))
end
# First Piola-Kirchhoff stress Jacobian
dPᵥ_Ωf_du(μ,u,du,v) = (dJ∘(∇(u),∇(du))) * (σᵥ_Ωf(μ)∘(∇(v),(Finv∘∇(u))))' ⋅ (FinvT∘∇(u)) +
                      (J∘∇(u)) * (σᵥ_Ωf(μ)∘(∇(v),(dFinv∘(∇(u),∇(du)))))' ⋅ (FinvT∘∇(u)) +
                      (J∘∇(u)) * (σᵥ_Ωf(μ)∘(∇(v),(Finv∘∇(u))))' ⋅ (dFinvT∘(∇(u),∇(du)))
dPₚ_Ωf_du(u,du,p) = - p * ( (dJ∘(∇(u),∇(du))) * tr((FinvT∘∇(u))) +
                           (J∘∇(u)) * tr((dFinvT∘(∇(u),∇(du)))) )
dPᵥ_Ωf_dv(μ,u,dv) = Pᵥ_Ωf(μ,u,dv)
dPₚ_Ωf_dp(u,dp) = Pₚ_Ωf(u,dp)
# Convective term
conv(c,∇v) = (∇v') ⋅ c
dconv(dc,∇dv,c,∇v) = conv(c,∇dv) + conv(dc,∇v)

# Solid constitutive laws
# Second Piola-Kirchhoff stress (Saint-Venant)
Sₛᵥ_Ωs(λ,μ,u) = 2*μ*(E∘∇(u)) + λ*tr((E∘∇(u)))*(I∘∇(u))
dSₛᵥ_Ωs_du(λ,μ,u,du) = 2*μ*(dE∘(∇(u),∇(du))) + λ*tr((dE∘(∇(u),∇(du))))*(I∘∇(u))
# First Piola-Kirchhoff stress (Saint-Venant)
Pₛᵥ_Ωs(λ,μ,u) = (F∘∇(u)) ⋅ Sₛᵥ_Ωs(λ,μ,u)
dPₛᵥ_Ωs_du(λ,μ,u,du) = (dF∘∇(du)) ⋅ Sₛᵥ_Ωs(λ,μ,u) + (F∘∇(u)) ⋅ dSₛᵥ_Ωs_du(λ,μ,u,du)
# C(F) = (F')⋅F
# function S_NH(∇u,λ,μ)
#   Cinv = invLaw(C(F(∇u)))
#   μ*(I(∇u)-Cinv) + λ*logLaw(J(∇u))*Cinv
# end
# function dS_NH(∇u,∇du,λ,μ)
#   Cinv = invLaw(C(F(∇u)))
#   _dE = dE(∇du,∇u)
#   λ*(Cinv⊙_dE)*Cinv + 2*(μ-λ*logLaw(J(∇u)))*Cinv⋅_dE⋅(Cinv')
# end

# Mesh constitutive laws
αₘ(J) = 1.0e-5 / J
σₘ(λ,μ) = ε -> λ*tr(ε)*I(ε) + 2.0*μ*ε
dαₘ(J,dJ) = - 1.0e-5 * 1.0 / (J*J) * dJ
dσₘ(λ,dλ,μ,dμ,ε,dε) = σ_m(λ,μ,dε) + σ_m(dλ,dμ,ε)
