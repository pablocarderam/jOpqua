module jOpqua

include("simulation/constants.jl")
include("simulation/structs.jl")

include("simulation/ImmunityType.jl")
include("parameters/Parameters.jl")

include("simulation/Pathogen.jl")
include("simulation/Immunity.jl")
include("simulation/Host.jl")
include("simulation/Class.jl")
include("simulation/Population.jl")
include("simulation/Simulation.jl")

include("simulation/Model.jl")

include("analysis/data.jl")
include("analysis/plots.jl")

end # module jOpqua
