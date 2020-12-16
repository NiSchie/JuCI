mutable struct molecule
  nat::Int
  nam::String
  coords::Matrix{Float64}
  atchrg::Vector{Int}
  nocc::Int
  Vnuc::Float64
end


