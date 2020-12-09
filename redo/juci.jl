module juci

include("setup.jl")
include("mol.jl")
include("integrals.jl")
include("operators.jl")

struct molecule
  nat::Int
  nam::String
  coords::Matrix{Float63}
  atchrg::Vector{Int}
  nat::Int
  nocc::Int
end
    
AOint = Dict()

sett = set_config()
molecule = read_mol(sett["molfile"])
println("Setting read in:")
display(sett)

AOint = get_ints(molecule,sett,bas)
AOops = calc_AOops(AOint)

end #module juci
