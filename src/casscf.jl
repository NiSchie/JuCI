module casscf

include("juci.jl")
include("rhf.jl")

function run_casscf()
  #setup calculation, get AO integrals, start orbitals (AOOps["C"])
  sett = juci.set_config()
  mol = juci.read_mol(sett["molfile"])
  AOint,dims = juci.get_ints(mol,sett)
  AOops = rhf.core(mol,sett,AOint,dims)

  strl = gen_strl(sett)

end #function casscf

function gen_strl(sett)

  nact = parse(Int,sett["nact"])
  nel  = Int(parse(Int,sett["nel"])/2)

  str =[1:1:nel;]

  #build graph field (to get dimenstion and adress weights)

  #loop over all determinants, calculate <det|a^+_p a_q |det> to get couplings

end #function gen_strl

function permute_string(str)
  

end #function permute_string

end #module casscf
