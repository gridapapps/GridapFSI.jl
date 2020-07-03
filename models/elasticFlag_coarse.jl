using Gridap
using Gridap.Io
using GridapGmsh

model = GmshDiscreteModel("elasticFlag_coarse.msh")

writevtk(model,"elasticFlag_coarse")

fn = "elasticFlag_coarse.json"
to_json_file(model,fn)
