# Fluid laws
@law σ_dev(μ,∇v,Finv) = μ*(∇v⋅Finv + (Finv')⋅(∇v'))
@law σ_dev(μ,ε) = 2.0*μ*ε
@law conv(c,∇v) = (∇v') ⋅ c
@law dconv(dc,∇dv,c,∇v) = conv(c,∇dv) + conv(dc,∇v)

# Mesh laws
@law α_m(J) = 1.0e-5 / J
@law σ_m(λ,μ,ε) = λ*tr(ε)*one(ε) + 2.0*μ*ε
@law dα_m(J,dJ) = - 1.0e-5 * 1.0 / (J*J) * dJ
@law dσ_m(λ,dλ,μ,dμ,ε,dε) = σ_m(λ,μ,dε) + σ_m(dλ,dμ,ε)

# Solid laws
@law I(A) = one(A)
@law E(∇u) = 0.5 * ((F(∇u)')⋅F(∇u) - I(∇u))
@law dE(∇u,∇du) = 0.5 * ((dF(∇du)')⋅F(∇u) + (F(∇u)')⋅dF(∇du))
S_SV(∇u,λ,μ) =  2*μ*E(∇u) + λ*tr(E(∇u))*I(∇u)
dS_SV(∇u,∇du,λ,μ) = 2*μ*dE(∇u,∇du) + λ*tr(dE(∇u,∇du))*I(∇u)
C(F) = (F')⋅F
function S_NH(∇u,λ,μ)
    Cinv = inv(C(F(∇u)))
    μ*(one(∇u)-Cinv) + λ*log(J(F(∇u)))*Cinv
end
function dS_NH(∇u,∇du,λ,μ)
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
