module WeakForms

using Gridap
using GridapODEs.TransientFETools: ∂t
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
  a_ST(x,y,μ,dΩ) - l_ST(y,f,dΩ)
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
function stokes_residual_Γd(x,y,n,μ::Real,γ::Real,h,vD,dΓ)
  a_ST_Γd(x,y,n,μ,γ,h,dΓ) - l_ST_Γd(y,n,μ,γ,h,vD,dΓ)
end


# FSI (fluid)
# ===========
# Residual
function residual_Ωf(st::MeshStrategy{:linearElasticity},t,x,y,params,dΩ)
  λₘ,μₘ = lame_parameters(params[:E],params[:ν])
  μ = params[:μ]; ρ = params[:ρ]; fᵤ = params[:fu]; fᵥ = params[:fv]
  a_mesh(st,x,y,λₘ,μₘ,dΩ) + a_NS_ALE(x,y,μ,ρ,dΩ) - l_mesh(st,y,fᵤ,t,dΩ) - l_NS_ALE(y,fᵥ,t,dΩ)
end
function residual_Ωf(st::MeshStrategy{:laplacian},t,x,y,params,dΩ)
  αₘ = params[:α]; μ = params[:μ]; ρ = params[:ρ]; fᵤ = params[:fu]; fᵥ = params[:fv]
  a_mesh(st,x,y,αₘ,dΩ) + a_NS_ALE(x,y,μ,ρ,dΩ) - l_mesh(st,y,fᵤ,t,dΩ) - l_NS_ALE(y,fᵥ,t,dΩ)
end
function residual_Ωf(st::MeshStrategy{:biharmonic},t,x,y,params,dΩ)
  _x = (x[2],x[3],x[4])
  _y = (y[2],y[3],y[4])
  αₘ = params[:α]; μ = params[:μ]; ρ = params[:ρ]; fᵤ = params[:fu]; fᵥ = params[:fv]
  a_mesh(st,x,y,αₘ,αₘ,dΩ) + a_NS_ALE(_x,_y,μ,ρ,dΩ) - l_mesh(st,y,fᵤ,t,dΩ) - l_NS_ALE(_y,fᵥ,t,dΩ)
end
# Spatial jacobian
function jacobian_Ωf(st::MeshStrategy{:linearElasticity},x,dx,y,params,dΩ)
  λₘ,μₘ = lame_parameters(params[:E],params[:ν])
  μ = params[:μ]; ρ = params[:ρ]
  a_mesh(st,dx,y,λₘ,μₘ,dΩ) +  da_NS_ALE_dx(x,dx,y,μ,ρ,dΩ)
end
function jacobian_Ωf(st::MeshStrategy{:laplacian},x,dx,y,params,dΩ)
  αₘ = params[:α]; μ = params[:μ]; ρ = params[:ρ]
  a_mesh(st,dx,y,αₘ,dΩ) +  da_NS_ALE_dx(x,dx,y,μ,ρ,dΩ)
end
function jacobian_Ωf(st::MeshStrategy{:biharmonic},x,dx,y,params,dΩ)
  _x = (x[2],x[3],x[4])
  _dx = (dx[2],dx[3],dx[4])
  _y = (y[2],y[3],y[4])
  αₘ = params[:α]; μ = params[:μ]; ρ = params[:ρ]
  a_mesh(st,dx,y,αₘ,αₘ,dΩ) +  da_NS_ALE_dx(_x,_dx,_y,μ,ρ,dΩ)
end
# Temporal Jacobian
function jacobian_t_Ωf(strategy::MeshStrategy,x,dxt,y,params,dΩ)
  ρ = params[:ρ]
  da_NS_ALE_dxt(x,dxt,y,ρ,dΩ)
end
function jacobian_t_Ωf(strategy::MeshStrategy{:biharmonic},x,dxt,y,params,dΩ)
  _x = (x[2],x[3],x[4])
  _dxt = (dxt[2],dxt[3],dxt[4])
  _y = (y[2],y[3],y[4])
  ρ = params[:ρ]
  da_NS_ALE_dxt(_x,_dxt,_y,ρ,dΩ)
end


# FSI (solid)
# ===========
# Residual
residual_Ωs(t,x,y,ρ,λ,μ,fᵤ,fᵥ,dΩ) = a_PFE(x,y,ρ,λ,μ,dΩ) - l_PFE(y,fᵤ,fᵥ,t,dΩ)
function residual_Ωs(st::MeshStrategy,t,x,y,params,dΩ)
  _x = (x[1],x[2])
  _y = (y[1],y[2])
  λ,μ = lame_parameters(params[:E],params[:ν])
  ρ = params[:ρ]; fᵤ = params[:fu]; fᵥ = params[:fv]
  residual_Ωs(t,_x,_y,ρ,λ,μ,fᵤ,fᵥ,dΩ)
end
function residual_Ωs(st::MeshStrategy{:biharmonic},t,x,y,params,dΩ)
  _x = (x[2],x[3])
  _y = (y[2],y[3])
  λ,μ = lame_parameters(params[:E],params[:ν])
  α = params[:α]; ρ = params[:ρ]; fₘ = params[:fu]; fᵤ = params[:fu]; fᵥ = params[:fv]
  a_mesh(st,x,y,α,α,dΩ) - l_mesh(st,y,fₘ,t,dΩ) + residual_Ωs(t,_x,_y,ρ,λ,μ,fᵤ,fᵥ,dΩ)
end
residual_Ωs_weak(st::MeshStrategy,t,x,y,params,dΩ) = residual_Ωs(st,t,x,y,params,dΩ)
function residual_Ωs_weak(st::MeshStrategy{:biharmonic},t,x,y,params,dΩ)
  _x = (x[2],x[3])
  _y = (y[2],y[3])
  λ,μ = lame_parameters(params[:E],params[:ν])
  ρ = params[:ρ]; fᵤ = params[:fu]; fᵥ = params[:fv]
  residual_Ωs(t,_x,_y,ρ,λ,μ,fᵤ,fᵥ,dΩ)
end
# Spatial Jacobian
jacobian_Ωs(x,dx,y,ρ,λ,μ,dΩ) = da_PFE_dx(x,dx,y,ρ,λ,μ,dΩ)
function jacobian_Ωs(st::MeshStrategy,x,dx,y,params,dΩ)
  _x = (x[1],x[2])
  _dx = (dx[1],dx[2])
  _y = (y[1],y[2])
  λ,μ = lame_parameters(params[:E],params[:ν])
  ρ = params[:ρ]
  jacobian_Ωs(_x,_dx,_y,ρ,λ,μ,dΩ)
end
function jacobian_Ωs(st::MeshStrategy{:biharmonic},x,dx,y,params,dΩ)
  _x = (x[2],x[3])
  _dx = (dx[2],dx[3])
  _y = (y[2],y[3])
  λ,μ = lame_parameters(params[:E],params[:ν])
  α = params[:α]; ρ = params[:ρ]
  a_mesh(st,dx,y,α,α,dΩ) + jacobian_Ωs(_x,_dx,_y,ρ,λ,μ,dΩ)
end
jacobian_Ωs_weak(st::MeshStrategy,x,dx,y,params,dΩ) = jacobian_Ωs(st,x,dx,y,params,dΩ)
function jacobian_Ωs_weak(st::MeshStrategy{:biharmonic},x,dx,y,params,dΩ)
  _x = (x[2],x[3])
  _dx = (dx[2],dx[3])
  _y = (y[2],y[3])
  λ,μ = lame_parameters(params[:E],params[:ν])
  α = params[:α]; ρ = params[:ρ]
  jacobian_Ωs(_x,_dx,_y,ρ,λ,μ,dΩ)
end
# Temporal Jacobian
function jacobian_t_Ωs(st::MeshStrategy,x,dxt,y,params,dΩ)
  _dxt = (dxt[1],dxt[2])
  _y = (y[1],y[2])
  ρ = params[:ρ]
  da_PFE_dxt(_dxt,_y,ρ,dΩ)
end
function jacobian_t_Ωs(st::MeshStrategy{:biharmonic},x,dxt,y,params,dΩ)
  _dxt = (dxt[2],dxt[3])
  _y = (y[2],y[3])
  ρ = params[:ρ]
  da_PFE_dxt(_dxt,_y,ρ,dΩ)
end


# FSI (interface)
# ===============
# Residual
function residual_Γi(st::MeshStrategy{:linearElasticity},x,y,params,dΓ)
  λₘ,μₘ = lame_parameters(params[:E],params[:ν])
  n = params[:n]
  a_mesh_Γi(st,x,y,n,λₘ,μₘ,dΓ)
end
function residual_Γi(st::MeshStrategy{:laplacian},x,y,params,dΓ)
  n = params[:n]; α = params[:α]
  a_mesh_Γi(st,x,y,n,α,dΓ)
end
function residual_Γi(st::MeshStrategy{:biharmonic},x,y,params,dΓ)
  n = params[:n]; α = params[:α]
  a_mesh_Γi(st,x,y,n,α,α,dΓ)
end
# Jacobian
function jacobian_Γi(st::MeshStrategy{:linearElasticity},x,dx,y,params,dΓ)
  λₘ,μₘ = lame_parameters(params[:E],params[:ν])
  n = params[:n]
  a_mesh_Γi(st,dx,y,n,λₘ,μₘ,dΓ)
end
function jacobian_Γi(st::MeshStrategy{:laplacian},x,dx,y,params,dΓ)
  n = params[:n]; α = params[:α]
  a_mesh_Γi(st,dx,y,n,α,dΓ)
end
function jacobian_Γi(st::MeshStrategy{:biharmonic},x,dx,y,params,dΓ)
  n = params[:n]; α = params[:α]
  a_mesh_Γi(st,dx,y,n,α,α,dΓ)
end


# Coupling management
# ===================
get_fluid_vars_Ω(::MeshStrategy,::Coupling{:strong},x) = x
get_fluid_vars_Ω(::MeshStrategy,::Coupling{:weak},x) = (x[1],x[2],x[3])
get_fluid_vars_Ω(::MeshStrategy{:biharmonic},::Coupling{:weak},x) = (x[1],x[2],x[3],x[4])
get_fluid_vars_Γi(::MeshStrategy,::Coupling{:strong},x) = x
get_fluid_vars_Γi(::MeshStrategy,::Coupling{:weak},x) = (x[1].⁺,x[2].⁺,x[3].⁺)
get_fluid_vars_Γi(::MeshStrategy{:biharmonic},::Coupling{:weak},x) = (x[1].⁺,x[2].⁺,x[3].⁺,x[4].⁺)
get_solid_vars_Ω(::MeshStrategy,::Coupling{:strong},x) = x
get_solid_vars_Ω(::MeshStrategy,::Coupling{:weak},x) = (x[4],x[5])
get_solid_vars_Ω(::MeshStrategy{:biharmonic},::Coupling{:weak},x) = (x[1],x[5],x[6])
get_interface_vars_Γi(::MeshStrategy,::Coupling{:weak},x) = x
get_interface_vars_Γi(::MeshStrategy{:biharmonic},::Coupling{:weak},x) = (x[2],x[3],x[4],x[5],x[6])
# Fluid
function residual_Ωf(st::MeshStrategy,c::Coupling,t,x,y,params,dΩ)
  xf = get_fluid_vars_Ω(st,c,x)
  yf = get_fluid_vars_Ω(st,c,y)
  residual_Ωf(st,t,xf,yf,params,dΩ)
end
function jacobian_Ωf(st::MeshStrategy,c::Coupling,t,x,dx,y,params,dΩ)
  xf = get_fluid_vars_Ω(st,c,x)
  dxf = get_fluid_vars_Ω(st,c,dx)
  yf = get_fluid_vars_Ω(st,c,y)
  jacobian_Ωf(st,xf,dxf,yf,params,dΩ)
end
function jacobian_t_Ωf(st::MeshStrategy,c::Coupling,t,x,dxt,y,params,dΩ)
  xf = get_fluid_vars_Ω(st,c,x)
  dxtf = get_fluid_vars_Ω(st,c,dxt)
  yf = get_fluid_vars_Ω(st,c,y)
  jacobian_t_Ωf(st,xf,dxtf,yf,params,dΩ)
end
# Solid
function residual_Ωs(st::MeshStrategy,c::Coupling,t,x,y,params,dΩ)
  xf = get_solid_vars_Ω(st,c,x)
  yf = get_solid_vars_Ω(st,c,y)
  residual_Ωs(st,t,xf,yf,params,dΩ)
end
function residual_Ωs(st::MeshStrategy,c::Coupling{:weak},t,x,y,params,dΩ)
  xf = get_solid_vars_Ω(st,c,x)
  yf = get_solid_vars_Ω(st,c,y)
  residual_Ωs_weak(st,t,xf,yf,params,dΩ)
end
function jacobian_Ωs(st::MeshStrategy,c::Coupling,t,x,dx,y,params,dΩ)
  xf = get_solid_vars_Ω(st,c,x)
  dxf = get_solid_vars_Ω(st,c,dx)
  yf = get_solid_vars_Ω(st,c,y)
  jacobian_Ωs(st,xf,dxf,yf,params,dΩ)
end
function jacobian_Ωs(st::MeshStrategy,c::Coupling{:weak},t,x,dx,y,params,dΩ)
  xf = get_solid_vars_Ω(st,c,x)
  dxf = get_solid_vars_Ω(st,c,dx)
  yf = get_solid_vars_Ω(st,c,y)
  jacobian_Ωs_weak(st,xf,dxf,yf,params,dΩ)
end
function jacobian_t_Ωs(st::MeshStrategy,c::Coupling,t,x,dxt,y,params,dΩ)
  xf = get_solid_vars_Ω(st,c,x)
  dxtf = get_solid_vars_Ω(st,c,dxt)
  yf = get_solid_vars_Ω(st,c,y)
  jacobian_t_Ωs(st,xf,dxtf,yf,params,dΩ)
end
# Interface
function residual_Γi(st::MeshStrategy,c::Coupling,x,y,params,dΓ)
  x_Γf = get_fluid_vars_Γi(st,c,x)
  y_Γf = get_fluid_vars_Γi(st,c,y)
  residual_Γi(st,x_Γf,y_Γf,params,dΓ)
end
function jacobian_Γi(st::MeshStrategy,c::Coupling,x,dx,y,params,dΓ)
  x_Γf = get_fluid_vars_Γi(st,c,x)
  dx_Γf = get_fluid_vars_Γi(st,c,dx)
  y_Γf = get_fluid_vars_Γi(st,c,y)
  jacobian_Γi(st,x_Γf,dx_Γf,y_Γf,params,dΓ)
end
function residual_Γi(st::MeshStrategy,c::Coupling{:weak},x,y,params,dΓ)
  x_Γf = get_fluid_vars_Γi(st,c,x)
  y_Γf = get_fluid_vars_Γi(st,c,y)
  x_Γi = get_interface_vars_Γi(st,c,x)
  y_Γi = get_interface_vars_Γi(st,c,y)
  residual_Γi(st,x_Γf,y_Γf,params,dΓ) +
  a_FSI_weak_Γi(x_Γi,y_Γi,params[:n],params[:μ],params[:γ],params[:h],params[:dt],dΓ)
end
function jacobian_Γi(st::MeshStrategy,c::Coupling{:weak},x,dx,y,params,dΓ)
  x_Γf = get_fluid_vars_Γi(st,c,x)
  dx_Γf = get_fluid_vars_Γi(st,c,dx)
  y_Γf = get_fluid_vars_Γi(st,c,y)
  x_Γi = get_interface_vars_Γi(st,c,x)
  y_Γi = get_interface_vars_Γi(st,c,y)
  dx_Γi = get_interface_vars_Γi(st,c,dx)
  jacobian_Γi(st,x_Γf,dx_Γf,y_Γf,params,dΓ) +
  da_FSI_weak_Γi_dx(x_Γi,dx_Γi,y_Γi,params[:n],params[:μ],params[:γ],params[:h],params[:dt],dΓ)
end
# Fuild boundary
function residual_Γf(st::MeshStrategy,t,x,y,params,dΓ)
  xf = get_fluid_vars_Ω(st,c,x)
  yf = get_fluid_vars_Ω(st,c,y)
  residual_Γi(st,xf,yf,params,dΓ) +
  a_NS_ALE_ΓD(xf,yf,t,params[:n],params[:μ],params[:γ],params[:h],dΓ) -
  l_NS_ALE_ΓD(yf,t,params[:vD],params[:n],params[:μ],params[:γ],params[:h],dΓ)
end
function jacobian_Γf(st::MeshStrategy,t,x,y,params,dΓ)
  xf = get_fluid_vars_Ω(st,c,x)
  dxf = get_fluid_vars_Ω(st,c,dx)
  yf = get_fluid_vars_Ω(st,c,y)
  jacobian_Γi(st,xf,dxf,yf,params,dΓ) +
  da_NS_ALE_ΓD_dx(xf,dxf,yf,params[:n],params[:μ],params[:γ],params[:h],dΓ)
end
function jacobian_t_Γf(st::MeshStrategy,t,x,y,params,dΓ)
  xf = get_fluid_vars_Ω(st,c,x)
  dxtf = get_fluid_vars_Ω(st,c,dxt)
  yf = get_fluid_vars_Ω(st,c,y)
  da_NS_ALE_ΓD_dxt(dxtf,yf,params[:μ],params[:γ],params[:h],dΓ)
end

end
