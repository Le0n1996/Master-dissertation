#data = read.csv(amzn, sep="")
#data

start.time <- Sys.time()

t = 0.0001*0:10000
D = 0
for (g in 1:10){
  b = 0
  w = 0.01*rnorm(10000, mean = 0, sd = 1)
  for (i in 1:10000){
    s = 0
    for (j in 1:i){
      s = s + w[j]
    }
    b = c(b, s)
  }
  b = b - t*b[10001]
  D = c(D, max(abs(b)))
}

l = length(D[-1])
print(0.5*(sort(D[-1])[floor(0.05*l)]+sort(D[-1])[ceiling(0.05*l)]))
print(0.5*(sort(D[-1])[floor(0.1*l)]+sort(D[-1])[ceiling(0.1*l)]))
print(0.5*(sort(D[-1])[floor(0.25*l)]+sort(D[-1])[ceiling(0.25*l)]))
print(0.5*(sort(D[-1])[floor(0.5*l)]+sort(D[-1])[ceiling(0.5*l)]))
print(0.5*(sort(D[-1])[floor(0.75*l)]+sort(D[-1])[ceiling(0.75*l)]))
print(0.5*(sort(D[-1])[floor(0.90*l)]+sort(D[-1])[ceiling(0.90*l)]))
print(0.5*(sort(D[-1])[floor(0.95*l)]+sort(D[-1])[ceiling(0.95*l)]))
print(0.5*(sort(D[-1])[floor(0.99*l)]+sort(D[-1])[ceiling(0.99*l)]))

end.time <- Sys.time()
time.taken <- end.time - start.time
print(time.taken)


l = length(Night[-1])
print(0.5*(sort(Night[-1])[floor(0.05*l)]+sort(Night[-1])[ceiling(0.05*l)]))
print(0.5*(sort(Night[-1])[floor(0.1*l)]+sort(Night[-1])[ceiling(0.1*l)]))
print(0.5*(sort(Night[-1])[floor(0.25*l)]+sort(Night[-1])[ceiling(0.25*l)]))
print(0.5*(sort(Night[-1])[floor(0.5*l)]+sort(Night[-1])[ceiling(0.5*l)]))
print(0.5*(sort(Night[-1])[floor(0.75*l)]+sort(Night[-1])[ceiling(0.75*l)]))
print(0.5*(sort(Night[-1])[floor(0.90*l)]+sort(Night[-1])[ceiling(0.90*l)]))
print(0.5*(sort(Night[-1])[floor(0.95*l)]+sort(Night[-1])[ceiling(0.95*l)]))
print(0.5*(sort(Night[-1])[floor(0.99*l)]+sort(Night[-1])[ceiling(0.99*l)]))


start.time <- Sys.time()

t = 0.0001*0:10000
D = 0
for (g in 1:10000){
  b = 0
  s = 0
  w = 1/sqrt(10000)*rnorm(10000, mean = 0, sd = 1)
  s = cumsum(w)
  b = c(b, s)
  b = b - t*b[10001]
  D = c(D, max(abs(b)))
}
l = length(D[-1])
D = sort(D[-1])

print(0.5*(D[floor(0.05*l)]+D[ceiling(0.05*l)]))
print(0.5*(D[floor(0.1*l)]+D[ceiling(0.1*l)]))
print(0.5*(D[floor(0.25*l)]+D[ceiling(0.25*l)]))
print(0.5*(D[floor(0.5*l)]+D[ceiling(0.5*l)]))
print(0.5*(D[floor(0.75*l)]+D[ceiling(0.75*l)]))
print(0.5*(D[floor(0.90*l)]+D[ceiling(0.90*l)]))
print(0.5*(D[floor(0.95*l)]+D[ceiling(0.95*l)]))
print(0.5*(D[floor(0.99*l)]+D[ceiling(0.99*l)]))

end.time <- Sys.time()
time.taken <- end.time - start.time
print(time.taken)