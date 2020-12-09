using Printf

function read_mol(molfile)

  @printf("Reading in the molecule\n")
  if isfile(sett["molfile"])
    lines = readlines(sett["molfile"])
  else
    throw("Molecule file not found")
  end

      nat = parse(Int,lines[1])
  pos   = zeros(Float64,3,nat)
  charg = Array{Int}(undef,nat)
  
  for i = 1:nat
    coord = [0.0, 0.0, 0.0]
    splitline = split(lines[i+2]," ")
    charg[i] = get(PhysConsts.atlist,splitline[1],0)
    filter!(e->e!="",splitline)
    popat!(splitline,1)
    coord .= parse.(Float64,splitline)
    pos[:,i] = collect(coord) .* PhysConsts.angstrom_to_bohr
  end
  @printf("Coordinates in bohr:\n")
  display(pos)
  
  nocc = Int(sum(charg)/2)
  mol = molecule(nat,lines[2],pos,charg,nocc)
  
  Vnuc = 0.0
  for i in 1:nat, j in 1:(i-1)
    global Vnuc
    x = mol.coords[1,i] - mol.coords[1,j]
    y = mol.coords[2,i] - mol.coords[2,j]
    z = mol.coords[3,i] - mol.coords[3,j]
    Vnuc += (mol.atchrg[i]*mol.atchrg[j])/(sqrt(x*x + y*y + z*z))
  end
  @printf("\n\nVnuc: %.5f\n\n",Vnuc)

  return mol