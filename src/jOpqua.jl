module jOpqua

include("simulation/constants.jl")

include("parameters/parameters.jl")

include("simulation/pathogen.jl")
include("simulation/immunity.jl")
include("simulation/host.jl")
include("simulation/class.jl")
include("simulation/population.jl")
include("simulation/simulation.jl")

include("simulation/model.jl")

include("analysis/data.jl")
include("analysis/plots.jl")

greet() = print("Hello World!")

end # module jOpqua
