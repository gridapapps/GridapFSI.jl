include("FSI_FESpaces.jl")
include("ElasticFlag.jl")
include("Analytical.jl")
include("Oscillator.jl")

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

function get_FSI_measures(Tₕ,order)
  degree = 2*order
  bdegree = 2*order
  dΩ  = Measure(Tₕ[:Ω],degree)
  dΩs = Measure(Tₕ[:Ωs],degree)
  dΩf = Measure(Tₕ[:Ωf],degree)
  dΓi = Measure(Tₕ[:Γi],bdegree)
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
  if m_params[:w_strategy] == "volume"
    α_Ωf = α_m * get_cell_measure(Tₕ[:Ωf])
    α_Ωs = α_m * get_cell_measure(Tₕ[:Ωs],Tₕ[:Ω])  # not sure why I need to reindex here...
    if ( typeof(coupling) == WeakForms.Coupling{:weak} )
      α_Γi = α_m * get_cell_measure(Tₕ[:Γi].⁺)
    else
      α_Γi = α_m * get_cell_measure(Tₕ[:Γi])
    end
  else
    α_Ωf = α_m; α_Ωs = α_m; α_Γi = α_m
  end

  # Compute interface element size (if weak coupling)
  if ( typeof(coupling) == WeakForms.Coupling{:weak} )
    dim = num_cell_dims(Tₕ[:Γi])
    h_Γfs = get_array(∫(1)dTₕ[:Γi])
    hΓᵢ = CellField( lazy_map(h->(h.^(-dim)),h_Γfs), Tₕ[:Γi])
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
  function res(t,x,y)
    residual_Ωf(strategy,coupling,t,x,y,f_params,dTₕ[:Ωf]) +
    residual_Ωs(strategy,coupling,t,x,y,s_params,dTₕ[:Ωs]) +
    residual_Γi(strategy,coupling,x,y,Γi_params,dTₕ[:Γi])
  end
  function jac(t,x,dx,y)
    jacobian_Ωf(strategy,coupling,t,x,dx,y,f_params,dTₕ[:Ωf]) +
    jacobian_Ωs(strategy,coupling,t,x,dx,y,s_params,dTₕ[:Ωs]) +
    jacobian_Γi(strategy,coupling,x,dx,y,Γi_params,dTₕ[:Γi])
  end
  function jac_t(t,x,dxt,y)
    jacobian_t_Ωf(strategy,coupling,t,x,dxt,y,f_params,dTₕ[:Ωf]) +
    jacobian_t_Ωs(strategy,coupling,t,x,dxt,y,s_params,dTₕ[:Ωs])
  end
  op = TransientFEOperator(res,jac,jac_t,X,Y)

end
