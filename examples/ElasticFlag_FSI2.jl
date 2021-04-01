using GridapFSI

output = main(
  strategy="biharmonic",
  alpha_m=1.0e-6,
  alpha_m_weight="volume",
  model="models/elasticFlagFine.json",
  dt=0.01,
  tf=10.0,
  theta=0.6
)
