module jOpqua

# Groundwork
include("utils.jl")

# Structure
include("simulation/constants.jl")
include("simulation/structs.jl")

# Input
include("parameters/parameters.jl")

# Core
include("simulation/initializers.jl")
include("simulation/flexlev.jl")
include("simulation/choice.jl")
include("simulation/weights.jl")
include("simulation/events.jl")
include("simulation/simulation.jl")

# Output
include("analysis/data.jl")
include("analysis/plots.jl")

# Extensions
include("extensions/immunity/immunity.jl")
include("extensions/intrahost/intrahost.jl")

end # module jOpqua
