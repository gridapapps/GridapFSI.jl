include("FSI_FESpaces.jl")
include("ElasticFlag.jl")
include("Analytical.jl")
include("Oscillator.jl")
include("ElasticFlagAggFEM.jl")

# Output function
function writePVD(filePath::String, trian::Triangulation, sol; append=false)
    outfiles = paraview_collection(filePath, append=append) do pvd
        for (i, (xh, t)) in enumerate(sol)
            println("STEP: $i, TIME: $t")
            println("============================")
            uh = xh[1]
            vh = xh[2]
						ph = xh[3]
            pvd[t] = createvtk(
                trian,
                filePath * "_$t.vtu",
                cellfields = ["uh" => uh, "vh" => vh, "ph" => ph]
            )
        end
    end
end

function get_FSI_triangulations(models,coupling)
  trian = Triangulation(models[:Ω])
  trian_s = Triangulation(models[:Ωs])
  trian_f = Triangulation(models[:Ωf])
  function Γi_triangulation(coupling)
    if typeof(coupling) == WeakForms.Coupling{:weak}
      InterfaceTriangulation(models[:Ωf],models[:Ωs])
    else
      BoundaryTriangulation(models[:Ωf],tags="interface")
    end
  end
  trian_Γi = Γi_triangulation(coupling)
  Dict(:Ω=>trian, :Ωs=>trian_s, :Ωf=>trian_f, :Γi=>trian_Γi)
end

function get_FSI_measures(triangulations,order)
  degree = 2*order
  bdegree = 2*order
  dΩ  = Measure(triangulations[:Ω],degree)
  dΩs = Measure(triangulations[:Ωs],degree)
  dΩf = Measure(triangulations[:Ωf],degree)
  dΓi = Measure(triangulations[:Γi],bdegree)
  Dict(:Ω=>dΩ, :Ωs=>dΩs, :Ωf=>dΩf, :Γi=>dΓi)
end

function get_Stokes_operator(X,Y,strategy,dΩ,μ,f)
  res(x,y) = WeakForms.stokes_residual(strategy,x,y,μ,f,dΩ)
  jac(x,dx,y) = WeakForms.stokes_jacobian(strategy,dx,y,μ,dΩ)
  op = FEOperator(res,jac,X,Y)
end

function get_FSI_operator(X,Y,coupling,strategy,Tₕ,dTₕ,params)
  m_params, f_params, s_params, Γi_params = params

  # Compute cell area (auxiliar quantity for mesh motion eq.)
  α_m = m_params[:α]
  if( m_params[:w_strategy] == "volume")
    # volf = get_cell_measure(Tₕ[:Ωf])
    # vols = get_cell_measure(Tₕ[:Ωs])
    # α_Ωf = α_m * reindex(volf,Tₕ[:Ωf])
    # α_Ωs = α_m * reindex(vols,Tₕ[:Ωs])
    α_Ωf = α_m * get_cell_measure(Tₕ[:Ωf])
    α_Ωs = α_m * get_cell_measure(Tₕ[:Ωs])
    # α_Ωs = α_Ωf
    if ( typeof(coupling) == WeakForms.Coupling{:weak} )
       α_Γi = lazy_map(Reindex(α_Ωf),get_cell_to_bgcell(Tₕ[:Γi].⁺))
    else
      α_Γi = lazy_map(Reindex(α_Ωf),get_cell_to_bgcell(Tₕ[:Γi]))
    end
  else
    α_Ωf = α_m; α_Ωs = α_m; α_Γi = α_m
  end

  # Compute interface element size (if weak coupling)
  if ( typeof(coupling) == WeakForms.Coupling{:weak} )
    dim = num_cell_dims(Tₕ[:Γi])
    h_Γfs = get_array(∫(1)dTₕ[:Γi])
    hΓᵢ = CellField( lazy_map(h->(h.^(-dim)),h_Γfs), Tₕ[:Γi])
    # trian_boundary_Γi = get_left_boundary(Tₕ[:Γi])
    # hΓᵢ = reindex(cell_measure(trian_boundary_Γi,Tₕ[:Ω]),trian_boundary_Γi)
  else
    hΓᵢ = 0.0
  end

  # Interface normal vector
  n_Γi = get_normal_vector(Tₕ[:Γi])

  # Complete params
  push!(f_params, :α=>α_Ωf)
  push!(f_params, :E=>m_params[:E])
  push!(f_params, :ν=>m_params[:ν])
  push!(s_params, :α=>α_Ωs)
  push!(Γi_params, :α=>α_Γi)
  push!(Γi_params, :E=>m_params[:E])
  push!(Γi_params, :ν=>m_params[:ν])
  push!(Γi_params, :n=>n_Γi)
  push!(Γi_params, :h=>hΓᵢ)

  # Define operator
  function res(t,(x,xt),y)
    WeakForms.fluid_residual_Ω(strategy,coupling,t,x,xt,y,f_params,dTₕ[:Ωf]) +
    WeakForms.solid_residual_Ω(strategy,coupling,t,x,xt,y,s_params,dTₕ[:Ωs]) +
    WeakForms.fsi_residual_Γi(strategy,coupling,x,y,Γi_params,dTₕ[:Γi])
  end
  function jac(t,(x,xt),dx,y)
    WeakForms.fluid_jacobian_Ω(strategy,coupling,t,x,xt,dx,y,f_params,dTₕ[:Ωf]) +
    WeakForms.solid_jacobian_Ω(strategy,coupling,t,x,xt,dx,y,s_params,dTₕ[:Ωs]) +
    WeakForms.fsi_jacobian_Γi(strategy,coupling,x,dx,y,Γi_params,dTₕ[:Γi])
  end
  function jac_t(t,(x,xt),dxt,y)
    WeakForms.fluid_jacobian_t_Ω(strategy,coupling,t,x,xt,dxt,y,f_params,dTₕ[:Ωf]) +
    WeakForms.solid_jacobian_t_Ω(strategy,coupling,t,x,xt,dxt,y,s_params,dTₕ[:Ωs])
  end
  op = TransientFEOperator(res,jac,jac_t,X,Y)

end
