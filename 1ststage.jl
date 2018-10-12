using PyPlot
using MUMPS


# Read dumped matrix
# Format is m, n integers setting the sizeof
# and m*n entries
function openMat(filename)
  if (isfile(filename))
    fp=open(filename)
    m=read(fp, Int32)
    n=read(fp, Int32)
    A=read(fp, Float64, m, n)
    close(fp)
    println("Reading matrix from ", filename, " ", m, "x", n)
  else
    println("Reading empty matrix from ", filename)
    A=Array{Float64}(0,2)
  end
  return A
end

# Read first stage matrix and compute solution
function firststage()
  println("Reading matrix M")
  A=openMat("1ststageM.dmp")
  # show structure
  # spy(A,precision=0)
  # show()
  # Reading RHS  
  println("Reading RHS")
  rhs=openMat("1ststageRHS.dmp")
  # Reading solution computed in PIPS  
  println("Reading solution")
  sol=openMat("1ststageSol.dmp")
  # Solve system in Julia  
  println("Solving system")
  nz = 0
  for i in A
    if i != 0 
      nz = nz + 1
    end
  end
  println("Non zeros: ", countnz(A))
  println("nz: ", nz)
  x=solveMUMPS(sparse(A),rhs,1)
  #x=\(A,rhs)
  
  # Compare both solutions  
  println("Solution from file:")
  println(norm(sol))
  println("Computed solution:")
  println(norm(x))
  
  println("Condition number: ", cond(A,Inf))
end

firststage()
