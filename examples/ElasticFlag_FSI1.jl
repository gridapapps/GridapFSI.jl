using GridapFSI

output = main(
  Um=0.2,
  Re=20,
  rho_s=1.0e-3,
  strategy="biharmonic",
  alpha_m=1.0e-6,
  alpha_m_weight="volume",
  model="../models/elasticFlagFine.json",
  dt=1.0,
  tf=25.0,
  theta=0.6
)
