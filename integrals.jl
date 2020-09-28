module Integrals

using Lints

mutable struct IntegralData{T}

  cache::Dict{String,Array}
  basnam::Dict{String,String}
  basis::Dict{String,Lints.BasisSetAllocated}

end #struct IntegralData





end #module Integrals
