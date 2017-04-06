using PyPlot

function openMat(filename)
  fp=open(filename)
  m=read(fp, Int32)
  n=read(fp, Int32)
  A=read(fp, Float64, m, n)
  close(fp)
  println(m,n)
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

function globalmat()
  A=openMat("globalA0.dmp")
  spy(A)
  show()
  A=openMat("globalA1.dmp")
  spy(A)
  show()
  C=openMat("globalC0.dmp")
  spy(C)
  show()
  C=openMat("globalC1.dmp")
  spy(C)
  show()
  R=openMat("globalR0.dmp")
  spy(R)
  show()
  R=openMat("globalR1.dmp")
  spy(R)
  show()
end

ioff()
#firststage()
globalmat()

