module jOpqua

# Groundwork
include("utils.jl")
include("simulation/flexle/src/sampler.jl")
include("simulation/flexle/src/interface.jl")

# Structure
include("simulation/constants.jl")
include("simulation/structs.jl")

# Extensions
include("extensions/immunity/immunity.jl")
include("extensions/intrahost/intrahost.jl")

# Input
include("parameters/default.jl")
include("parameters/setup.jl")

# Simulation
include("simulation/initializers.jl")
include("simulation/choice.jl")
include("simulation/weights.jl")
include("simulation/events.jl")
include("simulation/simulation.jl")

# Output
include("analysis/data.jl")
include("analysis/plots.jl")

end # module jOpqua
