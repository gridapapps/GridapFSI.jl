function execute(problem::PotentialFlowProblem{:beam};kwargs...)

  # Parameters
  L = 2*π
  H = 1.0
  xb₀ = 0.8*π
  xb₁ = 1.2*π
  n = 200
  order = 2
  g = 9.81
  ξ = 0.01
  λ = L/2
  k = 2*π/L
  h = L/n
  ω = √(g*k*tanh(k*H))
  t₀ = 0.0
  tf = 2*π
  Δt = 1/(50*ω)#h/(10*λ*ω)
  θ = 0.5

  # Exact solution
  ϕₑ(x,t) = ω/k * ξ * (cosh(k*(x[2]))) / sinh(k*H) * sin(k*x[1] - ω*t)
  ηₑ(x,t) = ξ * cos(k*x[1] - ω*t)
  ϕₑ(t::Real) = x -> ϕₑ(x,t)
  ηₑ(t::Real) = x -> ηₑ(x,t)

  # Interior Domain
  domain = (0,L,0,H)
  partition = (n,n/10)
  model_Ω = simplexify(CartesianDiscreteModel(domain,partition))#;isperiodic=(true,false)))
  Ω = Triangulation(model_Ω)

  # # Surface Domain
  # labels_Ω = get_face_labeling(model_Ω)
  # add_tag_from_tags!(labels_Ω,"bottom",[1,2,5])
  # add_tag_from_tags!(labels_Ω,"top",[3,4,6])
  # bgface_to_mask = get_face_mask(labels_Ω,[3,4,6],1)
  # model_Γ =BoundaryDiscreteModel(Polytope{1},model_Ω,bgface_to_mask)

  # Surface domains
  function is_beam(coords)
    n = length(coords)
		x = (1/n)*sum(coords)
    (xb₀ <= x[1] <= xb₁ ) * ( x[2] ≈ H )
  end

  # Boundary domain
  labels = get_face_labeling(model_Ω)
  bgface_to_mask = get_face_mask(labels,[3,4,6],1)
  Γface_to_bgface = findall(bgface_to_mask)
  model_Γ = BoundaryDiscreteModel(Polytope{1},model_Ω,Γface_to_bgface)
  Γ = Triangulation(model_Γ)

  # Get beam - free surface triangulations
  Γface_coords = get_cell_coordinates(Γ)
  Γface_mask = lazy_map(is_beam,Γface_coords)
  Γbface_Γface = findall(Γface_mask)
  Γfface_Γface = findall(!,Γface_mask)
  Γb = BoundaryTriangulation(model_Ω,view(Γface_to_bgface,Γbface_Γface))
  Γf = BoundaryTriangulation(model_Ω,view(Γface_to_bgface,Γfface_Γface))

  # writevtk(Γb,"Gb")
  # writevtk(Γf,"Gf")

  # Quadratures
  dΩ = Measure(Ω,2*order)
  dΓb = Measure(Γb,2*order)
  dΓf = Measure(Γf,2*order)

  # FE spaces
  reffe = ReferenceFE(lagrangian,Float64,order)
  V_Ω = TestFESpace(model_Ω,reffe,conformity=:H1)
  V_Γ = TestFESpace(model_Γ,reffe,conformity=:H1)
  U_Ω = TransientTrialFESpace(V_Ω)
  U_Γ = TransientTrialFESpace(V_Γ)
  X = TransientMultiFieldFESpace([U_Ω,U_Γ])
  Y = MultiFieldFESpace([V_Ω,V_Γ])

  # Weak form
  α = 2/Δt
  #m((ϕtt,ηtt),(w,v)) = ∫( v*ηtt )dΓb
  #c((ϕt,ηt),(w,v)) = ∫( 0.5*(α/g*w + v)*ϕt - w*ηt )dΓb + ∫( 0.5*(α/g*w + v)*ϕt - w*ηt )dΓf
  m((ϕt,ηt),(w,v)) = ∫( 0.5*(α/g*w + v)*ϕt - w*ηt )dΓb + ∫( 0.5*(α/g*w + v)*ϕt - w*ηt )dΓf
  a((ϕ,η),(w,v)) = ∫( ∇(ϕ)⋅∇(w) )dΩ + ∫( Δ(v)*Δ(η) + 0.5*(α/g*w + v)*g*η )dΓb + ∫( 0.5*(α/g*w + v)*g*η )dΓf
  b((w,v)) = ∫( 0.5*(α/g*w + v)*(-0.5) )dΓb
  #op = TransientConstantFEOperator(m,c,a,b,X,Y)
  op = TransientConstantFEOperator(m,a,b,X,Y)

  # Solver
  ls = LUSolver()
  #odes = Newmark(ls,Δt,0.5,0.25)
  odes = ThetaMethod(ls,Δt,0.5)
  solver = TransientFESolver(odes)

  # Initial solution
  x₀ = interpolate_everywhere([ϕₑ(0.0),ηₑ(0.0),],X(0.0))
  v₀ = interpolate_everywhere([∂t(ϕₑ)(0.0),∂t(ηₑ)(0.0)],X(0.0))
  a₀ = interpolate_everywhere([∂tt(ϕₑ)(0.0),∂tt(ηₑ)(0.0)],X(0.0))

  # Solution
  #sol_t = solve(solver,op,(x₀,v₀,a₀),t₀,tf)
  sol_t = solve(solver,op,x₀,t₀,tf)

  # Post-process
  l2_Ω(w) = √(∑(∫(w*w)dΩ))
  l2_Γb(v) = √(∑(∫(v*v)dΓb))
  l2_Γf(v) = √(∑(∫(v*v)dΓf))
  E_kin(w) = 0.5*∑( ∫(∇(w)⋅∇(w))dΩ )
  E_pot_b(v) = g*0.5*∑( ∫(v*v)dΓb )
  E_pot_f(v) = g*0.5*∑( ∫(v*v)dΓf )
  Eₑ = 0.5*g*ξ^2*L

  folderName = "PFlow-results"
  if !isdir(folderName)
    mkdir(folderName)
  end
  filePath_Ω = folderName * "/fields_O"
  filePath_Γb = folderName * "/fields_Gb"
  filePath_Γf = folderName * "/fields_Gf"
  pvd_Ω = paraview_collection(filePath_Ω, append=false)
  pvd_Γb = paraview_collection(filePath_Γb, append=false)
  pvd_Γf = paraview_collection(filePath_Γf, append=false)
  for ((ϕn,ηn),tn) in sol_t
    # E = E_kin(ϕn) + E_pot_b(ηn) + E_pot_f(ηn)
    # error_ϕ = l2_Ω(ϕn-ϕₑ(tn))
    # error_η = l2_Γb(ηn-ηₑ(tn)) + l2_Γf(ηn-ηₑ(tn))
    # println(E/Eₑ," ", error_ϕ," ",error_η)

    pvd_Ω[tn] = createvtk(Ω,filePath_Ω * "_$tn.vtu",cellfields = ["phi" => ϕn],order=2)
    pvd_Γb[tn] = createvtk(Γb,filePath_Γb * "_$tn.vtu",cellfields = ["eta" => ηn],order=2)
    pvd_Γf[tn] = createvtk(Γf,filePath_Γf * "_$tn.vtu",cellfields = ["eta" => ηn],order=2)
  end
  vtk_save(pvd_Ω)
  vtk_save(pvd_Γb)
  vtk_save(pvd_Γf)
end
