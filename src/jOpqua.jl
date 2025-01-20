module jOpqua

include("utils.jl")

include("simulation/constants.jl")
include("simulation/structs.jl")

include("parameters/parameters.jl")

include("simulation/initializers.jl")
include("simulation/choice.jl")
include("simulation/weights.jl")
include("simulation/events.jl")
include("simulation/simulation.jl")

include("analysis/data.jl")
include("analysis/plots.jl")

end # module jOpqua
