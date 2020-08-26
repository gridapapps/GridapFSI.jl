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
function stokes_residual(x,y,μ::Real,f)
    a_ST_vp(x,y,μ) - l_ST_vp(y,μ,f)
end
function stokes_residual(strategy::MeshStrategy,x,y,μ::Real,f)
    u, v, p = x
    ϕ, φ, q = y
    (ε(ϕ) ⊙ σ_dev(μ,ε(u))) + stokes_residual([v,p],[φ,q],μ::Real,f)
end
function stokes_residual(strategy::MeshStrategy{:biharmonic},x,y,μ::Real,f)
    w, u, v, p = x
    ψ, ϕ, φ, q = y
    _strategy = MeshStrategy{:linearElasticity}()
    (ε(ψ) ⊙ σ_dev(μ,ε(w))) + stokes_residual(_strategy,[u,v,p],[ϕ,φ,q],μ,f)
end
function stokes_jacobian(dx,y,μ::Real)
    a_ST_vp(dx,y,μ)
end
function stokes_jacobian(strategy::MeshStrategy,dx,y,μ::Real)
    du, dv, dp = dx
    ϕ, φ, q = y
    (ε(ϕ) ⊙ σ_dev(μ,ε(du))) + stokes_jacobian([dv,dp],[φ,q],μ)
end
function stokes_jacobian(strategy::MeshStrategy{:biharmonic},dx,y,μ::Real)
    dw, du, dv, dp = dx
    ψ, ϕ, φ, q = y
    _strategy = MeshStrategy{:linearElasticity}()
    (ε(ψ) ⊙ σ_dev(μ,ε(dw))) + stokes_jacobian(_strategy,[du,dv,dp],[ϕ,φ,q],μ)
end


# FSI (fluid)
fsi_residual_Ωf(strategy::MeshStrategy,t,x,xt,y,params) =
    a_FSI_ϕ_Ωf(strategy,x,y,params["E"],params["ν"]) +
    a_FSI_φ_Ωf(x,xt,y,params["μ"],params["ρ"]) +
    a_FSI_q_Ωf(x, y) -
    l_FSI_ϕ_Ωf(strategy,y,params["fu"],t) -
    l_FSI_φ_Ωf(y,params["fv"],t)
fsi_jacobian_Ωf(strategy::MeshStrategy,x,xt,dx,y,params) =
    da_FSI_du_ϕ_Ωf(strategy,x,dx,y,params["E"],params["ν"]) +
    da_FSI_du_φ_Ωf(x,xt,dx,y,params["μ"],params["ρ"]) +
    da_FSI_dv_φ_Ωf(x,xt,dx,y,params["μ"],params["ρ"]) +
    da_FSI_dp_φ_Ωf(x,dx,y) +
    da_FSI_du_q_Ωf(x,dx,y) +
    da_FSI_dv_q_Ωf(x,dx,y)
fsi_jacobian_t_Ωf(strategy::MeshStrategy,x,xt,dxt,y,params) =
    da_FSI_dut_φ_Ωf(x,dxt,y,params["ρ"]) +
    da_FSI_dvt_φ_Ωf(x,dxt,y,params["ρ"])

fsi_residual_Ωf(strategy::MeshStrategy{:laplacian},t,x,xt,y,params) =
    a_FSI_ϕ_Ωf(strategy,x,y,params["α"]) +
    a_FSI_φ_Ωf(x,xt,y,params["μ"],params["ρ"]) +
    a_FSI_q_Ωf(x, y) -
    l_FSI_ϕ_Ωf(strategy,y,params["fu"],t) -
    l_FSI_φ_Ωf(y,params["fv"],t)
fsi_jacobian_Ωf(strategy::MeshStrategy{:laplacian},x,xt,dx,y,params) =
    a_FSI_ϕ_Ωf(strategy,dx,y,params["α"]) +
    da_FSI_du_φ_Ωf(x,xt,dx,y,params["μ"],params["ρ"]) +
    da_FSI_dv_φ_Ωf(x,xt,dx,y,params["μ"],params["ρ"]) +
    da_FSI_dp_φ_Ωf(x,dx,y) +
    da_FSI_du_q_Ωf(x,dx,y) +
    da_FSI_dv_q_Ωf(x,dx,y)

function fsi_residual_Ωf(strategy::MeshStrategy{:biharmonic},t,x,xt,y,params)
    w, u, v, p = x
    wt, ut, vt, pt = xt
    ψ, ϕ, φ, q = y
    a_FSI_ψ_Ωf(strategy,x,y,params["α"]) +
    a_FSI_ϕ_Ωf(strategy,x,y,params["α"]) +
    a_FSI_φ_Ωf([u,v,p],[ut,vt,pt],[ϕ,φ,q],params["μ"],params["ρ"]) +
    a_FSI_q_Ωf([u,v,p],[ϕ,φ,q]) -
    l_FSI_ϕ_Ωf(strategy,[ϕ,φ,q],params["fu"],t) -
    l_FSI_φ_Ωf([ϕ,φ,q],params["fv"],t)
end
function fsi_jacobian_Ωf(strategy::MeshStrategy{:biharmonic},x,xt,dx,y,params)
    w, u, v, p = x
    wt, ut, vt, pt = xt
    dw, du, dv, dp = dx
    ψ, ϕ, φ, q = y
    da_FSI_dx_ψ_Ωf(strategy,x,dx,y,params["α"]) +
    da_FSI_dx_ϕ_Ωf(strategy,x,dx,y,params["α"]) +
    da_FSI_du_φ_Ωf([u,v,p],[ut,vt,pt],[du,dv,dp],[ϕ,φ,q],params["μ"],params["ρ"]) +
    da_FSI_dv_φ_Ωf([u,v,p],[ut,vt,pt],[du,dv,dp],[ϕ,φ,q],params["μ"],params["ρ"]) +
    da_FSI_dp_φ_Ωf([u,v,p],[du,dv,dp],[ϕ,φ,q]) +
    da_FSI_du_q_Ωf([u,v,p],[du,dv,dp],[ϕ,φ,q]) +
    da_FSI_dv_q_Ωf([u,v,p],[du,dv,dp],[ϕ,φ,q])
end
function fsi_jacobian_t_Ωf(strategy::MeshStrategy{:biharmonic},x,xt,dxt,y,params)
    w, u, v, p = x
    dwt, dut, dvt, dpt = dxt
    ψ, ϕ, φ, q = y
    da_FSI_dut_φ_Ωf([u,v,p],[dut,dvt,dpt],[ϕ,φ,q],params["ρ"]) +
        da_FSI_dvt_φ_Ωf([u,v,p],[dut,dvt,dpt],[ϕ,φ,q],params["ρ"])
end


# FSI (solid)
function fsi_residual_Ωs(strategy::MeshStrategy,t,x,xt,y,params)
  u, v, p = x
  ut, vt, pt = xt
  ϕ, φ, q = y
  a_FSI_ϕ_Ωs([u,v],[ut,vt],[ϕ,φ]) +
  a_FSI_φ_Ωs([u,v],[ut,vt],[ϕ,φ],params["ρ"],params["E"],params["ν"]) -
  l_FSI_ϕ_Ωs([ϕ,φ],params["fu"],t) -
  l_FSI_φ_Ωs([ϕ,φ],params["fv"],t)
end
function fsi_jacobian_Ωs(strategy::MeshStrategy,x,xt,dx,y,params)
  u, v, p = x
  ut, vt, pt = xt
  du, dv, dp = dx
  ϕ, φ, q = y
  da_FSI_dx_ϕ_Ωs([u,v],[du,dv],[ϕ,φ]) +
  da_FSI_dx_φ_Ωs([u,v],[du,dv],[ϕ,φ],params["ρ"],params["E"],params["ν"])
end
function fsi_jacobian_t_Ωs(strategy::MeshStrategy,x,xt,dxt,y,params)
  u, v, p = x
  dut, dvt, dpt = dxt
  ϕ, φ, q = y
  da_FSI_dxt_Ωs([u,v],[dut,dvt],[ϕ,φ],params["ρ"])
end

function fsi_residual_Ωs(strategy::MeshStrategy{:biharmonic},t,x,xt,y,params)
    w, u, v, p = x
    wt, ut, vt, pt = xt
    ψ, ϕ, φ, q = y
    a_FSI_ψ_Ωf(strategy,x,y,params["α"]) +
    a_FSI_ϕ_Ωs([u,v],[ut,vt],[ϕ,φ]) +
    a_FSI_φ_Ωs([u,v],[ut,vt],[ϕ,φ],params["ρ"],params["E"],params["ν"]) -
    l_FSI_ϕ_Ωs([ϕ,φ],params["fu"],t) -
    l_FSI_φ_Ωs([ϕ,φ],params["fv"],t)
end
function fsi_jacobian_Ωs(strategy::MeshStrategy{:biharmonic},x,xt,dx,y,params)
    w, u, v, p = x
    dw, du, dv, dp = dx
    ψ, ϕ, φ, q = y
    da_FSI_dx_ψ_Ωf(strategy,x,dx,y,params["α"]) +
    da_FSI_dx_ϕ_Ωs([u,v],[du,dv],[ϕ,φ]) +
    da_FSI_dx_φ_Ωs([u,v],[du,dv],[ϕ,φ],params["ρ"],params["E"],params["ν"])
end
function fsi_jacobian_t_Ωs(strategy::MeshStrategy{:biharmonic},x,xt,dxt,y,params)
    w, u, v, p = x
    dwt, dut, dvt, dpt = dxt
    ψ, ϕ, φ, q = y
    da_FSI_dxt_Ωs([u,v],[dut,dvt],[ϕ,φ],params["ρ"])
end

# FSI (interface)
fsi_residual_Γi(strategy::MeshStrategy,x,y,params) = a_FSI_ϕ_Γi(strategy,x,y,params["n"],params["E"],params["ν"])
fsi_jacobian_Γi(strategy::MeshStrategy,x,dx,y,params) = da_FSI_du_ϕ_Γi(strategy,x,dx,y,params["n"],params["E"],params["ν"])
fsi_residual_Γi(strategy::MeshStrategy{:laplacian},x,y,params) = a_FSI_ϕ_Γi(strategy,x,y,params["n"],params["α"])
fsi_jacobian_Γi(strategy::MeshStrategy{:laplacian},x,dx,y,params) = a_FSI_ϕ_Γi(strategy,dx,y,params["n"],params["α"])
fsi_residual_Γi(strategy::MeshStrategy{:biharmonic},x,y,params) = a_FSI_ϕ_Γi(strategy,x,y,params["n"],params["α"])
fsi_jacobian_Γi(strategy::MeshStrategy{:biharmonic},x,dx,y,params) = da_FSI_dx_ϕ_Γi(strategy,x,dx,y,params["n"],params["α"])

# Coupling management
fsi_residual_Ωf(strategy,coupling::Coupling,t,x,xt,y,params) = fsi_residual_Ωf(strategy,t,x,xt,y,params)
fsi_jacobian_Ωf(strategy,coupling::Coupling,x,xt,dx,y,params) = fsi_jacobian_Ωf(strategy,x,xt,dx,y,params)
fsi_jacobian_t_Ωf(strategy,coupling::Coupling,x,xt,dxt,y,params) = fsi_jacobian_t_Ωf(strategy,x,xt,dxt,y,params)
fsi_residual_Ωs(strategy,coupling::Coupling,t,x,xt,y,params) = fsi_residual_Ωs(strategy,t,x,xt,y,params)
fsi_jacobian_Ωs(strategy,coupling::Coupling,x,xt,dx,y,params) = fsi_jacobian_Ωs(strategy,x,xt,dx,y,params)
fsi_jacobian_t_Ωs(strategy,coupling::Coupling,x,xt,dxt,y,params) = fsi_jacobian_t_Ωs(strategy,x,xt,dxt,y,params)
fsi_residual_Γi(strategy,coupling::Coupling,x,y,params) = fsi_residual_Γi(strategy,x,y,params)
fsi_jacobian_Γi(strategy,coupling::Coupling,x,dx,y,params) = fsi_jacobian_Γi(strategy,x,dx,y,params)

function fsi_residual_Ωf(strategy::MeshStrategy{:biharmonic},coupling::Coupling{:weak},t,x,xt,y,params)
  w_f, u_f, v_f, p, u_s, v_s = x
  wt_f, ut_f, vt_f, pt, ut_s, vt_s = xt
  ψ_f, ϕ_f, φ_f, q, ϕ_s, φ_s = y
  x_f = w_f, u_f, v_f, p
  xt_f = wt_f, ut_f, vt_f, pt
  y_f = ψ_f, ϕ_f, φ_f, q
  fsi_residual_Ωf(strategy,t,x_f,xt_f,y_f,params)
end
function fsi_jacobian_Ωf(strategy::MeshStrategy{:biharmonic},coupling::Coupling{:weak},x,xt,dx,y,params)
  w_f, u_f, v_f, p, u_s, v_s = x
  wt_f, ut_f, vt_f, pt, ut_s, vt_s = xt
  dw_f, du_f, dv_f, dp, du_s, dv_s = dx
  ψ_f, ϕ_f, φ_f, q, ϕ_s, φ_s = y
  x_f = w_f, u_f, v_f, p
  xt_f = wt_f, ut_f, vt_f, pt
  dx_f = dw_f, du_f, dv_f, dp
  y_f = ψ_f, ϕ_f, φ_f, q
  fsi_jacobian_Ωf(strategy,x_f,xt_f,dx_f,y_f,params)
end
function fsi_jacobian_t_Ωf(strategy::MeshStrategy{:biharmonic},coupling::Coupling{:weak},x,xt,dxt,y,params)
  w_f, u_f, v_f, p, u_s, v_s = x
  wt_f, ut_f, vt_f, pt, ut_s, vt_s = xt
  dwt_f, dut_f, dvt_f, dpt, dut_s, dvt_s = dxt
  ψ_f, ϕ_f, φ_f, q, ϕ_s, φ_s = y
  x_f = w_f, u_f, v_f, p
  xt_f = wt_f, ut_f, vt_f, pt
  dxt_f = dwt_f, dut_f, dvt_f, dpt
  y_f = ψ_f, ϕ_f, φ_f, q
  fsi_jacobian_t_Ωf(strategy,x_f,xt_f,dxt_f,y_f,params)
end
function fsi_residual_Ωs(strategy::MeshStrategy{:biharmonic},coupling::Coupling{:weak},t,x,xt,y,params)
  w_f, u_f, v_f, p, u_s, v_s = x
  wt_f, ut_f, vt_f, pt, ut_s, vt_s = xt
  ψ_f, ϕ_f, φ_f, q, ϕ_s, φ_s = y
  x_s = w_f, u_s, v_s, p
  xt_s = wt_f, ut_s, vt_s, pt
  y_s = ψ_f, ϕ_s, φ_s, q
  fsi_residual_Ωs(strategy,t,x_s,xt_s,y_s,params)
end
function fsi_jacobian_Ωs(strategy::MeshStrategy{:biharmonic},coupling::Coupling{:weak},x,xt,dx,y,params)
  w_f, u_f, v_f, p, u_s, v_s = x
  wt_f, ut_f, vt_f, pt, ut_s, vt_s = xt
  dw_f, du_f, dv_f, dp, du_s, dv_s = dx
  ψ_f, ϕ_f, φ_f, q, ϕ_s, φ_s = y
  x_s = w_f, u_s, v_s, p
  dx_s = dw_f, du_s, dv_s, dp
  xt_s = wt_f, ut_s, vt_s, pt
  y_s = ψ_f, ϕ_s, φ_s, q
  fsi_jacobian_Ωs(strategy,x_s,xt_s,dx_s,y_s,params)
end
function fsi_jacobian_t_Ωs(strategy::MeshStrategy{:biharmonic},coupling::Coupling{:weak},x,xt,dxt,y,params)
  w_f, u_f, v_f, p, u_s, v_s = x
  wt_f, ut_f, vt_f, pt, ut_s, vt_s = xt
  dwt_f, dut_f, dvt_f, dpt, dut_s, dvt_s = dxt
  ψ_f, ϕ_f, φ_f, q, ϕ_s, φ_s = y
  x_s = w_f, u_s, v_s, p
  dxt_s = dwt_f, dut_s, dvt_s, dpt
  xt_s = wt_f, ut_s, vt_s, pt
  y_s = ψ_f, ϕ_s, φ_s, q
  fsi_jacobian_t_Ωs(strategy,x_s,xt_s,dxt_s,y_s,params)
end
function fsi_residual_Γi(strategy::MeshStrategy{:biharmonic},coupling::Coupling{:weak},x,y,params)
  w_f, u_f, v_f, p, u_s, v_s = x
  ψ_f, ϕ_f, φ_f, q, ϕ_s, φ_s = y
  x_Γf = w_f.inward, u_f.inward, v_f.inward, p.inward
  y_Γf = ψ_f.inward, ϕ_f.inward, φ_f.inward, q.inward
  fsi_residual_Γi(strategy,x_Γf,y_Γf,params) +
  a_FSI_Nitsche_ϕ_Γi([u_f,v_f,p,u_s,v_s],[ϕ_f, φ_f, q, ϕ_s, φ_s],params["n"],params["μ"],params["γ"],params["h"])
end
function fsi_jacobian_Γi(strategy::MeshStrategy{:biharmonic},coupling::Coupling{:weak},x,dx,y,params)
  w_f, u_f, v_f, p, u_s, v_s = x
  dw_f, du_f, dv_f, dp, du_s, dv_s = dx
  ψ_f, ϕ_f, φ_f, q, ϕ_s, φ_s = y
  x_Γf = w_f.inward, u_f.inward, v_f.inward, p.inward
  dx_Γf = dw_f.inward, du_f.inward, dv_f.inward, dp.inward
  y_Γf = ψ_f.inward, ϕ_f.inward, φ_f.inward, q.inward
  fsi_jacobian_Γi(strategy,x_Γf,dx_Γf,y_Γf,params) +
  da_FSI_Nitsche_ϕ_Γi([u_f,v_f,p,u_s,v_s],[du_f,dv_f,dp,du_s,dv_s],[ϕ_f, φ_f, q, ϕ_s, φ_s],params["n"],params["μ"],params["γ"],params["h"])
end
end
