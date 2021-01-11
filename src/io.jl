module io

using DelimitedFiles

function write_arr(fn,arr,ldlm)
  if isfile(fn) 
    rm(fn)
  end
  f = open(fn,"w")
  if(ldlm) 
    writedlm(f,arr)
  else
    write(f,arr)
  end #if
  close(f)
end #write_arr_fmt

end #module io
