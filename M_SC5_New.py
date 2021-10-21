# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# import items for the Mij matrix Below
# import the acid table
table_file0 = open("cas9table0.txt")
lines0 = table_file0.readlines()
table_file0.close()

# import the cas9 amino acid sequence
seq_file0 = open("aminosequence0.txt")
seq_string0 = seq_file0.read()
seq_file0.close()

# import items for the Hij matrix Below
# import the hamiltonian chart
table_file = open("HamiltonianChart.txt")
lines = table_file.readlines()
table_file.close()

# import the amino acid list
seq_file = open("list.txt")
seq_string = seq_file.read()
seq_file.close()

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# print items for the Mij matrix below
# print out the list of amino acids
print(" ")
print("acids")
acids0 = lines0[0].split()
print(acids0)

# print out the cas9 amino acid sequence
print(" ")
print("cas9 sequence")
seq_list0 = seq_string0.split()
print(seq_list0)

# print items for the Hij matrix Below
# print out the list of amino acids
print(" ")
print("acids")
seq_list = seq_string.split()
print(seq_list)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# parameters for the Hij matrix below
n = len(seq_list)
seq_index = []

# parameters for the Mij matrix below
m = len(seq_list0)
seq_index0 = []

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# matrices corresponding to the Mij matrix below
# hydrogen bond probability matrix
acid_table0 = []

# Mij matrix below
Mij = []

# matrices corresponding to the Hij matrix below
# matrix for the types of hydrogen bonds
hamiltonian_table = []

# matrix for the energy levels
Eij = []

# hamiltonian matrix
Hij = []

# tensor product for acid_table0 and Eij
Aij = []

# matrix for the Kronecker delta
Rij = []

# tensor product for Mij and Kij
Dij = []

# Energy corresponding to the sidechains
Hsc = []

Es = []

S = 0

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# construct the hydrogen bond probability matrix
for line0 in lines0[1:]:
    
    row0 = line0.split()[1:]
    numbers0 = list(map(float,row0))
    acid_table0.append(numbers0)

for item0 in seq_list0:
    seq_index0.append(acids0.index(item0))

# construct the Mij matrtix
for i in range(m):
    Mij.append([])
    for j in range(m):
        prob = acid_table0[seq_index0[i]][seq_index0[j]]
        Mij[i].append(prob)

# construct the matrix for the types of hydrogen bonds
for line in lines[1:]:
    
    row = line.split()[1:]
    letters = list(map(str,row))
    hamiltonian_table.append(letters)

# construct the matrix for the energy levels
for i in range(n):
    Eij.append([])
    for j in range(n):
        prob = hamiltonian_table[i][j]
        if prob == "N":
            prob2 = -1.64013e-22
            Eij[i].append(prob2)
        elif prob == "O":
            prob2 = -2.09727e-22
            Eij[i].append(prob2)
        elif prob == "P":
            prob2 = 0.0
            Eij[i].append(prob2)
        else:
            prob2 = 0.0
            Eij[i].append(prob2)

# tensor product of acid_table0 and Eij
for i in range(n):
    Aij.append([])
    for j in range(n):
            prob = acid_table0[i][j]*Eij[i][j]
            Aij[i].append(prob)

# kronecker delta matrix
for i in range(m):
        Rij.append([])
        for j in range(m):
            ith_index = i
            jth_index = j
            Rij[i].append("d(r" + str(ith_index) + " " +"- r" + str(jth_index) + " " + "- r)" )

# hamiltonian matrix
for i in range(m):
    Hij.append([])
    for j in range(m):
        prob = Eij[seq_index0[i]][seq_index0[j]]
        Hij[i].append(prob)

# tensor product of Mij and Kij
for i in range(m):
    Dij.append([])
    for j in range(m):
            prob = Mij[i][j]
            Dij[i].append(prob)

# steric hinderance correction
for i in range(m):
    Es.append([])
    for j in range(m):
            prob = str(2.26006e-22) + str(Mij[i][j]) + Rij[i][j]  
            Es[i].append(prob)

# Energy corresponding to the sidechains
for i in range(m):
    Hsc.append([])
    for j in range(m):
        prob = Hij[i][j] * Dij[i][j]
        if i == j + 1 or j == i + 1:
            Hsc[i].append(prob)
        else:
            Hsc[i].append(0)

for i in range(m):
    for j in range(m):
        S = S + Hsc[i][j]

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
print(Hsc)
print(S)
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 