using Dates
using HDF5
using CSV

try
  cd("C:/Users/avtei/OneDrive - University of Strathclyde/ML/scikit-learn/mordred/Judred")
catch
  try
    cd("C:/Users/rudolf/OneDrive - University of Strathclyde/ML/scikit-learn/mordred/Judred")
  catch
    cd("/users/rkb19187/Desktop/judred")
  end
end

include("peptideutils.jl")

ticks() = round(Int64, time() * 1000)

if length(ARGS) == 1
  L = parse(Int8, ARGS[1])
else
  L = 3
end

pep = fill("A", (L,))


letters_1 = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"]
#sp2 carbons in side-chain only
SP2 =       [0,    0,   1,   1,   6,   0,   3,   0,   0,   0,   0,   1,   0,   1,   1,   0,   0,   0,   8,   6]
SP3 =       [1,    1,   1,   2,   1,   0,   1,   4,   4,   4,   3,   1,   3,   2,   3,   1,   2,   3,   1,   1]
NH2 =       [0,    0,   0,   0,   0,   0,   0,   0,   1,   0,   0,   1,   0,   1,   2,   0,   0,   0,   0,   0]
MW =        [89.10, 121.16, 133.11, 147.13, 165.19, 75.07, 155.16, 131.18, 146.19, 131.18, 149.21, 132.12, 115.13, 146.15, 174.20, 105.09, 119.12, 117.15, 204.23, 181.19]
S =         [0,    1,   0,   0,   0,   0,   0,   0,   0,   0,   1,   0,   0,   0,   0,   0,   0,   0,   0,   0]
charge =    [0,    0,  -1,  -1,  0,    0,   0,   0,   1,   0,   0,   0,   0,   0,   1,   0,   0,   0,   0,   0]
# ASP ARG and LYS are as charged side chains
Gwif = [0.17, -0.24, 1.23, 2.02, -1.13, 0.01, 0.17, -0.31, 0.99, -0.56, -0.23, 0.42, 0.45, 0.58, 0.81, 0.13, 0.14, 0.07, -1.85, -0.94, ] #kcal / mol
Gwoct = [0.5, -0.02, 3.64, 3.63, -1.71, 1.15, 0.11, -1.12, 2.8, -1.25, -0.67, 0.85, 0.14, 0.77, 1.81, 0.46, 0.25, -0.46, -2.09, -0.71, ] #kcal / mol
#Tien et al. 2013 (theory)
MaxASA =    [129, 167, 193, 223, 240, 104, 224, 197, 236, 201, 224, 195, 159, 225, 274, 155, 172, 174, 285, 263]
# Zimmerman J.M., Eliezer N., Simha R. J. Theor. Biol. 21:170-201(1968).
#RotRatio =  [0.2, 0.333, 0.375, 0.444, 0.25, 0.25, 0.273, 0.375, 0.556, 0.375, 0.5, 0.375, 0.125, 0.444, 0.455, 0.333, 0.286, 0.286, 0.188, 0.231]
bulky =     [11.50, 13.46, 11.68, 13.57, 19.80, 3.4, 13.69, 21.40, 15.71, 21.4, 16.25, 12.82, 17.43, 14.45, 14.28, 9.47, 15.77, 21.57, 21.67, 18.03]
OH =        [0,  0,   0,   0,    0,     0,    0,   0,   0,  0,   0,  0,  0,    0,   0,  1,  1,    0,   0,  1]
#nBase =     [0,   0,  0,   0,    0,     0,    0,   0,   1,  0,   0,  0,  0,    0,   3,    0,  0,  0,   0,  0]

peptides = Array{String}(undef, 20^L)
#data = zeros(Float32, L, 6)
APs_data = zeros(Float32, L)
features = ["SPRatio", "NH2", "MW", "S", "LogP WW", "Z", "MaxASA", "RotRatio", "Bulkiness", "OH"]
data = zeros(Float32, length(features), 20^L)

function descriptors(pep, data_index)
  #SP2, NH2, MW, S, LogP WW Z
  MW_H2O = 18.01528 #g/mol
  local_SP3 = 0
  for letter in pep
    letter = string(letter)
    index = findall(x->x==letter, letters_1)
    index = index[1]
    data[1, data_index] += SP2[index]
    local_SP3 += SP3[index]
    data[2, data_index] += NH2[index]
    data[3, data_index] += MW[index]
    data[4, data_index] += S[index]
    data[5, data_index] += (Gwif[index] - Gwoct[index])
    data[6, data_index] += charge[index]
    data[7, data_index] += MaxASA[index]
    #data[8, data_index] = RotRatio[index]
    data[9, data_index] += bulky[index]
    data[10, data_index] += OH[index]
  end
  data[8, data_index] = data[1, data_index] / local_SP3
  #data[8, data_index] = data[8, data_index] / length(pep)
  #data[N,6] = data[N,6] / length(pep)
  #data[N,9] = data[N,9] / length(pep)
  data[3, data_index] -= MW_H2O * (length(pep)-1)
end

st = ticks() # ms
println("Running")
for n in range(1, step=1, stop=20^L)

  for i in range(1, step=1, stop=L)
    i = L - i +1
    if pep[i] == "Y"
      pep[i] = "A"
    else
      index = findall(x->x==pep[i], letters_1)[1]
      pep[i] = letters_1[index+1]
      break
    end
  end
  peptides[n] = peptideutils.translate1to3(join(pep))
  descriptors(pep, n)
  #break
end

took = (ticks() - st) / 1000

perpep = took / 20^L
penta = (perpep * (20^5)) / 60
hexa = (perpep * (20^6)) / 60
hepta = (perpep * (20^7)) / 60 / 60
octa = (perpep * (20^8)) / 60 / 60 / 24

if L == 2
  fid=h5open("dipeptides.hdf5","w")
elseif L == 3
    fid=h5open("tripeptides.hdf5","w")
elseif L == 4
  fid=h5open("tetrapeptides.hdf5","w")
elseif L == 5
  fid=h5open("pentapeptides.hdf5","w")
elseif L == 6
  fid=h5open("hexapeptides.hdf5","w")
elseif L == 7
  fid=h5open("heptapeptides.hdf5","w")
elseif L == 8
  fid=h5open("octapeptides.hdf5","w")
end

fid["peptides"] = peptides
fid["data"] = data
fid["features"] = features
close(fid)


println("ofsize(data): ", (sizeof(data)/10^6)+(sizeof(features)/10^6)+(sizeof(peptides)/10^6), " MB")
println("on ONE cpu core")
println("Took ", took, " seconds")
println("For pentapeptides this would take: ", penta, " mins ")
println("For hexapeptides this would take: ", hexa, " mins ")
println("For heptapeptides this would take: ", hepta, " hours")
println("For octapeptides this would take: ", octa, " days")
