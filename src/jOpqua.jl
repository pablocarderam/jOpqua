module jOpqua

# Groundwork
include("utils.jl")

# Structure
include("simulation/constants.jl")
include("simulation/structs.jl")

# Weight/rate computation
include("simulation/weights.jl")

# Extensions
include("extensions/immunity/immunity.jl")
include("extensions/intrahost/intrahost.jl")

# Input
include("parameters/default.jl")
include("parameters/setup.jl")

# Simulation
include("simulation/initializers.jl")
include("simulation/choice.jl")
include("simulation/events.jl")
include("simulation/simulation.jl")

# Output
include("analysis/data.jl")
include("analysis/plots.jl")

end # module jOpqua
