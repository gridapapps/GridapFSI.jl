function execute(problem::PotentialFlowProblem{:analytical};kwargs...)

  # Parameters
  L = 2*π
  H = 1.0
  n = 4
  order = 2
  g = 9.81
  ξ = 0.1
  λ = L/2
  k = 2*π/L
  h = L/n
  ω = √(g*k*tanh(k*H))
  t₀ = 0.0
  tf = 8*π
  Δt = h/(10*λ*ω)
  γ = 0.5
  β = 0.25

  # Exact solution
  ϕₑ(x,t) = ω/k * ξ * (cosh(k*(x[2]))) / sinh(k*H) * sin(k*x[1] - ω*t)
  ηₑ(x,t) = ξ * cos(k*x[1] - ω*t)
  ϕₑ(t::Real) = x -> ϕₑ(x,t)
  ηₑ(t::Real) = x -> ηₑ(x,t)

  # Domain
  domain = (0,L,0,H)
  partition = (n,n)
  model = CartesianDiscreteModel(domain,partition;isperiodic=(true,false))

  # Boundaries
  labels = get_face_labeling(model)
  add_tag_from_tags!(labels,"bottom",[1,2,5])
  add_tag_from_tags!(labels,"free_surface",[3,4,6])

  # Triangulation
  Ω = Triangulation(model)
  Γ = BoundaryTriangulation(model,tags="free_surface")
  dΩ = Measure(Ω,2*order)
  dΓ = Measure(Γ,2*order)

  # FE spaces
  reffe = ReferenceFE(lagrangian,Float64,order)
  V = TestFESpace(model,reffe,conformity=:H1)
  U = TransientTrialFESpace(V)

  # Weak form
  m(ϕₜₜ,w) = ∫( 1/g*ϕₜₜ*w )dΓ
  c(ϕₜ,w) = ∫( 0.0*ϕₜ*w )dΓ
  a(ϕ,w) = ∫( ∇(ϕ)⋅∇(w) )dΩ
  b(w) = ∫( 0.0*w )dΓ
  op = TransientConstantFEOperator(m,c,a,b,U,V)

  # Solver
  ls = LUSolver()
  odes = Newmark(ls,Δt,γ,β)
  solver = TransientFESolver(odes)

  # Initial solution
  ϕ₀ = interpolate_everywhere(ϕₑ(0.0),U(0.0))
  ∂ϕ₀_∂t = interpolate_everywhere(∂t(ϕₑ)(0.0),U(0.0))
  ∂ϕ₀_∂tt = interpolate_everywhere(∂tt(ϕₑ)(0.0),U(0.0))

  # Solution
  sol_t = solve(solver,op,ϕ₀,∂ϕ₀_∂t,∂ϕ₀_∂tt,t₀,tf)

  # Post-process
  #η(∂tϕₙ) = -β*Δt/γ*∂tϕₙ/g
  l2(w) = √(∑(∫(w*w)dΩ))
  E_kin(w) = 0.5*∑( ∫(∇(w)⋅∇(w))dΩ )
  #E_pot(w) = g*0.5*∑( ∫(η(w)⋅η(w))dΓ )
  Eₑ = 0.5*g*ξ^2*L

  next = Base.iterate(sol_t)
  while next != nothing
    (ϕₙ,tₙ), state = next
    U_aux, ode_state = state
    uf,u0,v0,a0,tf,cache = ode_state
    Eₜ = E_kin(ϕₙ) #+ E_pot(v0)
    println(Eₜ/Eₑ)
    next = Base.iterate(sol_t,state)
  end


  println("xxxx")

end