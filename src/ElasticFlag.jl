function get_parameters(problem::Problem{:elasticFlag})
    parameters = [
        # Problem Setting
        "Um" => 1.0,
        "H" => 0.41,
        "⌀" => 0.1,
        "t0" => 0.0,
        "Re" => 100.0,

        # Solid properties
        "E_s" => 1.4e6,
        "ν_s" => 0.4,
        "ρ_s" => 1.0e4,

        # Fluid properties
        "μ_f" => 1.0,
        "ρ_f" => 1.0e3,

        # Mesh properties
        "E_m" => 1.0,
        "ν_m" => 0.1,
    ]
end


