function execute(problem::Problem{:oscillator}; kwargs...)

  # Parameters
  A₀ = 1.5
  A₁ = 8.5
  A₂ = -11.1
  A₃ = 2.6
  B₀ = 4.2
  B₁ = 11.3
  B₂ = -68.7
  B₃ = 50.5

  # Physics Properties
  k = 1.0
  m = 1.0
  D = 0.1
  L = 1.0
  ρ = 1.0e3
  b = 1.0
  V = 1.0

  # Parameters
  q₁ = 2
  St = 0.2
  Cx0 = 2.0
  Cy0 = 0.3
  A = 12
  ε = 0.3

  # Auxiliar variables
  mₐ = π*ρ*D^2*L/4
  mʳ = m/mₐ
  ωₙ = sqrt(k/(m+mₐ))
  ωₛ = St*V/(2*π*D)
  Ωₙ = ωₙ/ωₛ
  ζ = b/(2*sqrt((m+mₐ)*k))

function oscillator(x,p,t)
  y, q, dy, dq = x
  y_dot = dy
  q_dot = dq
  Cvy= (-2^π*St*dy*Cx0 + Cy0/q₁*q) * sqrt(1+4*π^2*St^2*dy^2)
  dy_dot = ρ*D^2*L/(m+mₐ)*1.0/(8*π^2*St^2)*Cvy - 2*ζ*Ωₙ*dy - Ωₙ^2*y
  dq_dot = A*dy_dot - ε*(q^2-1)*dq - q
  dx = [y_dot; q_dot; dy_dot; dq_dot]
end

prob = ODEProblem(oscillator,[0.0;2.0;0.0;0.0],(0.0,100.0))
sol=DifferentialEquations.solve(prob)

end
