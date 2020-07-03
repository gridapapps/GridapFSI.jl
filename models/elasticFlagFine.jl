using Gridap
using Gridap.Io
using GridapGmsh

model = GmshDiscreteModel("elasticFlagFine.msh")

writevtk(model,"elasticFlagFine")

fn = "elasticFlagFine.json"
to_json_file(model,fn)
