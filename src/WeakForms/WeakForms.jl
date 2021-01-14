module WeakForms

using Gridap
using LinearAlgebra: tr, inv, det

export MeshStrategy
export Coupling

struct MeshStrategy{Kind} end
struct Coupling{Kind} end

# Material properties
# ===================
function lame_parameters(E,ν)
  λ = (E*ν)/((1+ν)*(1-2*ν))
  μ = E/(2*(1+ν))
  (λ, μ)
end

# Laws
include("Laws.jl")
# Terms
include("FluidTerms.jl")
include("SolidTerms.jl")
include("InterfaceTerms.jl")

# Stokes
"""
Stokes residual:
R([u,p],[v,q]=a([u,p],[v,q])-l([v,q]))
"""
function stokes_residual(x,y,μ::Real,f,dΩ)
  a_ST(x,y,μ,dΩ) #- l_ST(y,f,dΩ)
end
function stokes_residual(strategy::MeshStrategy,(u,v,p),(ϕ,φ,q),μ::Real,f,dΩ)
  ∫(ϕ⋅u)dΩ + stokes_residual((v,p),(φ,q),μ,f,dΩ)
end
function stokes_residual(strategy::MeshStrategy{:biharmonic},(w,u,v,p),(ψ,ϕ,φ,q),μ::Real,f,dΩ)
  ∫(ϕ⋅u + ψ⋅w)dΩ + stokes_residual((v,p),(φ,q),μ,f,dΩ)
end
function stokes_jacobian(dx,y,μ::Real,dΩ)
  a_ST(dx,y,μ,dΩ)
end
function stokes_jacobian(strategy::MeshStrategy,(du,dv,dp),(ϕ,φ,q),μ::Real,dΩ)
  ∫(ϕ⋅du)dΩ + stokes_jacobian((dv,dp),(φ,q),μ,dΩ)
end
function stokes_jacobian(strategy::MeshStrategy{:biharmonic},(dw,du,dv,dp),(ψ,ϕ,φ,q),μ::Real,dΩ)
  ∫(ϕ⋅du + ψ⋅dw)dΩ + stokes_jacobian((dv,dp),(φ,q),μ,dΩ)
end

# FSI (fluid)
# ===========
# Residual
function fluid_residual_Ω(st::MeshStrategy{:linearElasticity},t,x,xt,y,params,dΩ)
  λₘ,μₘ = lame_parameters(params[:E],params[:ν])
  μ = params[:μ]; ρ = params[:ρ]; fᵤ = params[:fu]; fᵥ = params[:fv]
  a_mesh(st,x,y,λₘ,μₘ,dΩ) + a_NS_ALE(x,xt,y,μ,ρ,dΩ) - l_mesh(st,y,fᵤ,t,dΩ) - l_NS_ALE(y,fᵥ,t,dΩ)
end
function fluid_residual_Ω(st::MeshStrategy{:laplacian},t,x,xt,y,params,dΩ)
  αₘ = params[:α]; μ = params[:μ]; ρ = params[:ρ]; fᵤ = params[:fu]; fᵥ = params[:fv]
  a_mesh(st,x,y,αₘ,dΩ) + a_NS_ALE(x,xt,y,μ,ρ,dΩ) - l_mesh(st,y,fᵤ,t,dΩ) - l_NS_ALE(y,fᵥ,t,dΩ)
end
function fluid_residual_Ω(st::MeshStrategy{:biharmonic},t,x,xt,y,params,dΩ)
  _x = x[2:4]
  _xt = xt[2:4]
  _y = y[2:4]
  αₘ = params[:α]; μ = params[:μ]; ρ = params[:ρ]; fᵤ = params[:fu]; fᵥ = params[:fv]
  a_mesh(st,x,y,αₘ,dΩ) + a_NS_ALE(_x,_xt,_y,μ,ρ,dΩ) - l_mesh(st,y,fᵤ,t,dΩ) - l_NS_ALE(_y,fᵥ,t,dΩ)
end
# Spatial jacobian
function fluid_jacobian_Ω(st::MeshStrategy{:linearElasticity},x,xt,dx,y,params,dΩ)
  λₘ,μₘ = lame_parameters(params[:E],params[:ν])
  μ = params[:μ]; ρ = params[:ρ]
  a_mesh(st,dx,y,λₘ,μₘ,dΩ) +  da_NS_ALE_dx(x,xt,dx,y,μ,ρ,dΩ)
end
function fluid_jacobian_Ω(st::MeshStrategy{:laplacian},x,xt,dx,y,params,dΩ)
  αₘ = params[:α]; μ = params[:μ]; ρ = params[:ρ]
  a_mesh(st,dx,y,αₘ,dΩ) +  da_NS_ALE_dx(x,xt,dx,y,μ,ρ,dΩ)
end
function fluid_jacobian_Ω(st::MeshStrategy{:biharmonic},x,xt,dx,y,params,dΩ)
  _x = x[2:4]
  _xt = xt[2:4]
  _dx = dx[2:4]
  _y = y[2:4]
  αₘ = params[:α]; μ = params[:μ]; ρ = params[:ρ]
  a_mesh(st,dx,y,αₘ,dΩ) +  da_NS_ALE_dx(_x,_xt,_dx,_y,μ,ρ,dΩ)
end
# Temporal Jacobian
function fluid_jacobian_t_Ω(strategy::MeshStrategy,x,xt,dxt,y,params,dΩ)
  ρ = params[:ρ]
  da_NS_ALE_dxt(x,dxt,y,ρ,dΩ)
end
function fluid_jacobian_t_Ω(strategy::MeshStrategy{:biharmonic},x,xt,dxt,y,params,dΩ)
  _x = x[2:4]
  _dxt = dxt[2:4]
  _y = y[2:4]
  ρ = params[:ρ]
  da_NS_ALE_dxt(_x,_dxt,_y,ρ,dΩ)
end


# FSI (solid)
# ===========
# Residual
function solid_residual_Ω(st::MeshStrategy,t,x,xt,y,params,dΩ)
  _x = x[1:2]
  _xt = xt[1:2]
  _y = y[1:2]
  λ,μ = lame_parameters(params[:E],params[:ν])
  ρ = params[:ρ]; fᵤ = params[:fu]; fᵥ = params[:fv]
  a_PFE(_x,_xt,_y,ρ,λ,μ,dΩ) - l_PFE(_y,fᵤ,fᵥ,t,dΩ)
end
function solid_residual_Ω(st::MeshStrategy{:biharmonic},t,x,xt,y,params,dΩ)
  _x = x[1:2]
  _xt = xt[1:2]
  _y = y[1:2]
  λ,μ = lame_parameters(params[:E],params[:ν])
  α = params[:α]; ρ = params[:ρ]; fₘ = params[:fu]; fᵤ = params[:fu]; fᵥ = params[:fv]
  a_mesh(st,x,y,α,dΩ) + a_PFE(_x,_xt,_y,ρ,λ,μ,dΩ) - l_mesh(st,y,fₘ,t,dΩ) - l_PFE(_y,fᵤ,fᵥ,t,dΩ)
end
# Spatial Jacobian
function solid_jacobian_Ω(st::MeshStrategy,t,x,xt,dx,y,params,dΩ)
  _x = x[1:2]
  _xt = xt[1:2]
  _dx = dx[1:2]
  _y = y[1:2]
  λ,μ = lame_parameters(params[:E],params[:ν])
  ρ = params[:ρ]
  da_PFE_dx(_x,_dx,_y,ρ,λ,μ,dΩ)
end
function solid_jacobian_Ω(st::MeshStrategy{:biharmonic},t,x,xt,dx,y,params,dΩ)
  _x = x[2:3]
  _xt = xt[2:3]
  _dx = dx[2:3]
  _y = y[2:3]
  λ,μ = lame_parameters(params[:E],params[:ν])
  α = params[:α]; ρ = params[:ρ]
  a_mesh(st,dx,y,α,dΩ) + da_PFE_dx(_x,_dx,_y,ρ,λ,μ,dΩ)
end
# Temporal Jacobian
function solid_jacobian_t_Ω(st::MeshStrategy,t,x,xt,dxt,y,params,dΩ)
  _dxt = dxt[1:2]
  _y = y[1:2]
  ρ = params[:ρ]
  da_PFE_dxt(_dxt,_y,ρ,dΩ)
end
function solid_jacobian_t_Ω(st::MeshStrategy{:biharmonic},t,x,xt,dxt,y,params,dΩ)
  _dxt = dxt[2:3]
  _y = y[2:3]
  ρ = params[:ρ]
  da_PFE_dxt(_dxt,_y,ρ,dΩ)
end


# FSI (interface)
# ===============
# Residual
function fsi_residual_Γi(st::MeshStrategy{:linearElasticity},x,y,params,dΓ)
  λₘ,μₘ = lame_parameters(params[:E],params[:ν])
  n = params[:n]
  a_mesh_Γi(st,x,y,n,λₘ,μₘ,dΓ)
end
function fsi_residual_Γi(st::MeshStrategy{:laplacian},x,y,params,dΓ)
  n = params[:n]; α = params[:α]
  a_mesh_Γi(st,x,y,n,α,dΓ)
end
function fsi_residual_Γi(st::MeshStrategy{:biharmonic},x,y,params,dΓ)
  n = params[:n]; α = params[:α]
  a_mesh_Γi(st,x,y,n,α,α,dΓ)
end
# Jacobian
function fsi_jacobian_Γi(st::MeshStrategy{:linearElasticity},x,dx,y,params,dΓ)
  λₘ,μₘ = lame_parameters(params[:E],params[:ν])
  n = params[:n]
  a_mesh_Γi(st,dx,y,n,λₘ,μₘ,dΓ)
end
function fsi_jacobian_Γi(st::MeshStrategy{:laplacian},x,dx,y,params,dΓ)
  n = params[:n]; α = params[:α]
  a_mesh_Γi(st,dx,y,n,α,dΓ)
end
function fsi_jacobian_Γi(st::MeshStrategy{:biharmonic},x,dx,y,params,dΓ)
  n = params[:n]; α = params[:α]
  a_mesh_Γi(st,dx,y,n,α,α,dΓ)
end


# Coupling management
# ===================
get_fluid_vars_Ω(::MeshStrategy,::Coupling{:strong},x) = x
get_fluid_vars_Ω(::MeshStrategy,::Coupling{:weak},x) = x[1:3]
get_fluid_vars_Ω(::MeshStrategy{:biharmonic},::Coupling{:weak},x) = x[1:4]
get_fluid_vars_Γi(::MeshStrategy,::Coupling{:weak},x) = (x[1].⁺,x[2].⁺,x[3].⁺)
get_fluid_vars_Γi(::MeshStrategy{:biharmonic},::Coupling{:weak},x) = (x[1].⁺,x[2].⁺,x[3].⁺,x[4].⁺)
get_solid_vars_Ω(::MeshStrategy,::Coupling{:strong},x) = x
get_solid_vars_Ω(::MeshStrategy,::Coupling{:weak},x) = x[4:5]
get_solid_vars_Ω(::MeshStrategy{:biharmonic},::Coupling{:weak},x) = x[1,5:6]
# Fluid
function fluid_residual_Ω(st::MeshStrategy,c::Coupling,t,x,xt,y,params,dΩ)
  xf = get_fluid_vars_Ω(st,c,x)
  xtf = get_fluid_vars_Ω(st,c,xt)
  yf = get_fluid_vars_Ω(st,c,y)
  fluid_residual_Ω(st,t,xf,xtf,yf,params,dΩ)
end
function fluid_jacobian_Ω(st::MeshStrategy,::Coupling,t,x,xt,dx,y,params,dΩ)
  xf = get_fluid_vars_Ω(st,c,x)
  xtf = get_fluid_vars_Ω(st,c,xt)
  dxf = get_fluid_vars_Ω(st,c,dx)
  yf = get_fluid_vars_Ω(st,c,y)
  fluid_residual_Ω(st,t,xf,xtf,dxf,yf,params,dΩ)
end
function fluid_jacobian_t_Ω(st::MeshStrategy,::Coupling,t,x,xt,dxt,y,params,dΩ)
  xf = get_fluid_vars_Ω(st,c,x)
  xtf = get_fluid_vars_Ω(st,c,xt)
  dxtf = get_fluid_vars_Ω(st,c,dxt)
  yf = get_fluid_vars_Ω(st,c,y)
  fluid_residual_Ω(st,t,xf,xtf,dxtf,yf,params,dΩ)
end
# Solid
function solid_residual_Ω(st::MeshStrategy,c::Coupling,t,x,xt,y,params,dΩ)
  xf = get_solid_vars_Ω(st,c,x)
  xtf = get_solid_vars_Ω(st,c,xt)
  yf = get_solid_vars_Ω(st,c,y)
  solid_residual_Ω(st,t,xf,xtf,yf,params,dΩ)
end
function solid_jacobian_Ω(st::MeshStrategy,c::Coupling,t,x,xt,dx,y,params,dΩ)
  xf = get_solid_vars_Ω(st,c,x)
  xtf = get_solid_vars_Ω(st,c,xt)
  dxf = get_solid_vars_Ω(st,c,dx)
  yf = get_solid_vars_Ω(st,c,y)
  solid_residual_Ω(st,t,xf,xtf,dxf,yf,params,dΩ)
end
function solid_jacobian_t_Ω(st::MeshStrategy,c::Coupling,t,x,xt,dxt,y,params,dΩ)
  xf = get_solid_vars_Ω(st,c,x)
  xtf = get_solid_vars_Ω(st,c,xt)
  dxtf = get_solid_vars_Ω(st,c,dxt)
  yf = get_solid_vars_Ω(st,c,y)
  solid_residual_Ω(st,t,xf,xtf,dxtf,yf,params,dΩ)
end
# Interface
function fsi_residual_Γi(st::MeshStrategy,c::Coupling{:weak},x,y,params,dΓ)
  x_Γf = get_fluid_vars_Γi(st,c,x)
  y_Γf = get_fluid_vars_Γi(st,c,y)
  fsi_residual_Γi(st,x_Γf,y_Γf,params,dΓ) +
  a_FSI_weak_Γi(x,y,params[:n],params[:μ],params[:γ],params[:h],params[:dt],dΓ)
end
function fsi_jacobian_Γi(st::MeshStrategy,c::Coupling{:weak},x,dx,y,params,dΓ)
  x_Γf = get_fluid_vars_Γi(st,c,x)
  dx_Γf = get_fluid_vars_Γi(st,c,dx)
  y_Γf = get_fluid_vars_Γi(st,c,y)
  fsi_jacobian_Γi(strategy,x_Γf,dx_Γf,y_Γf,params,dΓ) +
  da_FSI_weak_Γi_dx(x,dx,y,params[:n],params[:μ],params[:γ],params[:h],params[:dt],dΓ)
end
# Fuild boundary
function fluid_residual_Γ(st::MeshStrategy,t,x,xt,y,params,dΓ)
  xf = get_fluid_vars_Ω(st,c,x)
  xtf = get_fluid_vars_Ω(st,c,xt)
  yf = get_fluid_vars_Ω(st,c,y)
  fsi_residual_Γi(st,xf,yf,params,dΓ) +
  a_NS_ALE_ΓD(xf,xtf,yf,t,params[:vD],params[:n],params[:μ],params[:γ],params[:h],dΓ)
end
function fluid_jacobian_Γ(st::MeshStrategy,t,x,xt,y,params,dΓ)
  xf = get_fluid_vars_Ω(st,c,x)
  xtf = get_fluid_vars_Ω(st,c,xt)
  dxf = get_fluid_vars_Ω(st,c,dx)
  yf = get_fluid_vars_Ω(st,c,y)
  fsi_jacobian_Γi(st,xf,dxf,yf,params,dΓ) +
  da_NS_ALE_ΓD_dx(xf,xtf,dxf,yf,params[:n],params[:μ],params[:γ],params[:h],dΓ)
end
function fluid_jacobian_t_Γ(st::MeshStrategy,t,x,xt,y,params,dΓ)
  xf = get_fluid_vars_Ω(st,c,x)
  xtf = get_fluid_vars_Ω(st,c,xt)
  dxtf = get_fluid_vars_Ω(st,c,dxt)
  yf = get_fluid_vars_Ω(st,c,y)
  da_NS_ALE_ΓD_dxt(dxtf,yf,params[:μ],params[:γ],params[:h],dΓ)
end

end
