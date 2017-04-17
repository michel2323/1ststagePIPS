using PyPlot

function openMat(filename)
  fp=open(filename)
  m=read(fp, Int32)
  n=read(fp, Int32)
  A=read(fp, Float64, m, n)
  close(fp)
  println(filename," ", m," ",n)
  return A
end

function firststage()
  println("Reading matrix M")
  A=openMat("1ststageM.dmp")
  spy(A)
  show()
  
  println("Reading RHS")
  rhs=openMat("1ststageRHS.dmp")
  
  println("Reading solution")
  sol=openMat("1ststageSol.dmp")
  
  println("Solving system")
  x=\(A,rhs)
  
  println("Solution from file:")
  println(sol)
  println("Computed solution:")
  println(x)
  
  println("Condition number: ", cond(A))
end

function buildMaster()
  Q=openMat("globalQ00_1.dmp")
  spy(Q)
  show()
  B=openMat("globalB00_1.dmp")
  spy(B)
  show()
  C=openMat("globalC00_1.dmp")
  spy(B)
  show()
end
function buildW(child)
  filename="globalQs" * "$child" * "_1.dmp"
  println(filename)
  Q=openMat(filename)
  spy(Q)
  show()
  filename="globalBs" * "$child" * "_1.dmp"
  B=openMat(filename)
  spy(B)
  show()
  filename="globalDs" * "$child" * "_1.dmp"
  D=openMat(filename)
  spy(D)
  show()
end

ioff()
#firststage()
#M=buildMaster()
W=buildW(0)

