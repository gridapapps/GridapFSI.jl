# Fluid laws
@law σ_dev(μ,∇v,Finv) = μ*(∇v⋅Finv + (Finv')⋅(∇v'))
@law σ_dev(μ,ε) = 2.0*μ*ε
@law conv(c,∇v) = (∇v') ⋅ c
@law dconv(dc,∇dv,c,∇v) = conv(c,∇dv) + conv(dc,∇v)

# Mesh laws
@law α_m(J) = 1.0e-6 # 1.0e-6 / J
@law σ_m(λ,μ,ε) = λ*tr(ε)*one(ε) + 2.0*μ*ε
@law dα_m(J,dJ) = 0.0 #- 1.0e-6 * 2.0 / (J*J*J) * dJ
@law dσ_m(λ,dλ,μ,dμ,ε,dε) = σ_m(λ,μ,dε) + σ_m(dλ,dμ,ε)

# Solid laws
@law E(∇u) = 0.5 * ((F(∇u)')⋅F(∇u) - one(F(∇u)))
@law dE(∇u,∇du) = 0.5 * ((dF(∇du)')⋅F(∇u) + (F(∇u)')⋅dF(∇du))
#@law S(λ,μ,∇u) = 2*μ*E(∇u) + λ*tr(E(∇u))*one(E(∇u))
@law S_SV(∇u) = 2*0.5e6*E(∇u) + 2.0e6*tr(E(∇u))*one(E(∇u))
#@law dS(λ,μ,∇u,∇du) = 2*μ*dE(∇u,∇du) + λ*tr(dE(∇u,∇du))*one(E(∇u))
@law dS_SV(∇u,∇du) = 2*0.5e6*dE(∇u,∇du) + 2.0e6*tr(dE(∇u,∇du))*one(E(∇u))
C(F) = (F')⋅F
@law function S_NH(∇u)
    (λ,μ) = lame_parameters(1.0e-10,0.49)
    Cinv = inv(C(F(∇u)))
    μ*(one(∇u)-Cinv) + λ*log(J(F(∇u)))*Cinv
end
@law function dS_NH(∇du,∇u)
    (λ,μ) = lame_parameters(1.0e-10,0.49)
    Cinv = inv(C(F(∇u)))
    _dE = dE(∇du,∇u)
	  λ*(Cinv⊙_dE)*Cinv + 2*(μ-λ*log(J(F(∇u))))*Cinv⋅_dE⋅(Cinv')
end

# Maps
@law F(∇u) = ∇u + one(∇u)
@law J(∇u) = det(F(∇u))
@law Finv(∇u) = inv(F(∇u))
@law FinvT(∇u) = (Finv(∇u)')

# Map derivatives
@law dF(∇du) = ∇du
@law dJ(∇u,∇du) = J(F(∇u))*tr(inv(F(∇u))⋅dF(∇du))
@law dFinv(∇u,∇du) = -Finv(∇u) ⋅ dF(∇du) ⋅ Finv(∇u)
@law dFinvT(∇u,∇du) = (dFinv(∇u,∇du)')

# # ### Constitutive laws
# #@law σ_m(∇ut,Finv) = λ_m*tr(∇ut⋅Finv)*one(Finv) + μ_m*(∇ut⋅Finv + (Finv')⋅(∇ut')) 
# #@law σ_m(∇ut,Finv) = β_m*tr(∇ut⋅Finv)*one(Finv) + (∇ut⋅Finv + (Finv')⋅(∇ut')) 
# #@law σ_m2(α,∇ut,Finv) = 16.0*α*μ_s*tr(∇ut⋅Finv)*one(Finv) + α*μ_s*(∇ut⋅Finv + (Finv')⋅(∇ut')) 
# @law Sm(∇u) = 2*μ_m*E(∇u) + λ_m*tr(E(∇u))*one(E(∇u))
# @law α(J) = 1.0e-6/(J*J)
# #@law λ_m(J) = (E_m*ν_m) / ((1.0+ν_m)*(1.0-2.0*ν_m))
# #@law μ_m(J) = E_m / (2.0*(1.0+ν_m))

# # Derivatives:

# @law dSm(∇u,∇du) = 2*μ_m*dE(∇u,∇du) + λ_m*tr(dE(∇u,∇du))*one(E(∇u))
# @law dλ_m(dα) = dα * (E_m*ν_m) / ((1.0+ν_m)*(1.0-2.0*ν_m))
# @law dμ_m(dα) = dα * E_m / (2.0*(1.0+ν_m))
