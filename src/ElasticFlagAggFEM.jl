function execute(problem::Problem{:elasticFlagAggFem}; kwargs...)

  ## Build geometry
  # Cilinder
  R = 0.05
  C = Point(0.2,0.2)
  cylinder = disk(R,x0=C)
  # Flag
  C1 = Point(0.2,0.19)
  d1 = VectorValue(0.4,0.0)
  d2 = VectorValue(0.0,0.2)
  flag = quadrilateral(x0=C1,d1=d1,d2=d2)
  # Fluid and solid domains
  solid_geo = setdiff(flag,cylinder,name="solid")
  not_fluid = union(cylinder,flag)
  fluid_geo = !(not_fluid,name="fluid")

  # Background model
  pmin = Point(0.0,0.0)
  pmax = Point(2.2,0.41)
  partition = (100,20)
  bgmodel = CartesianDiscreteModel(pmin,pmax,partition)

  # Fluid and solid models
  cutgeo = cut(bgmodel,union(fluid_geo,solid_geo))
  model_fluid = DiscreteModel(cutgeo,"fluid")
  model_solid = DiscreteModel(cutgeo,"solid")
  writevtk(model_fluid,"model_fluid")
  writevtk(model_solid,"model_solid")

end
