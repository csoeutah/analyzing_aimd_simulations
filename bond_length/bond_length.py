import pandas as pd
import math
import matplotlib.pyplot as plt
import sys

#*XDATCAR.xlsx file as input to the script
input_file = sys.argv[1]
dataset = pd.read_excel(input_file)
dataset = dataset.fillna(0)

#lattice parameters in angstrom
a = (dataset.iloc[[0], :]).values.tolist()[0][0]
b = (dataset.iloc[[1], :]).values.tolist()[0][1]
c = (dataset.iloc[[2], :]).values.tolist()[0][2]

#atoms and number of atoms
atom = (dataset.iloc[[3], :]).values.tolist()[0]
number_of_atom = (dataset.iloc[[4], :]).values.tolist()[0]

#coordinates
X = (dataset.iloc[5:,[0]]).values.tolist()
Y = (dataset.iloc[5:,[1]]).values.tolist()
Z = (dataset.iloc[5:,[2]]).values.tolist()

#raising errors for invalid input file
if (len(X) != len(Y)):
    print("Check Input File. Coordinates missing/extra!")
    sys.exit()

if (len(X) != len(Z)):
    print("Check Input File. Coordinates missing/extra!")
    sys.exit()

#frames in dynamic simualtion
number_of_frames = math.floor(len(X)/sum(number_of_atom))

#raising errors for missing frame data
if ((len(X)/sum(number_of_atom)) > number_of_frames):
    print("Missing data in the last frame!")
    sys.exit()

#atom numbers as input to the script for tracking the bond distance
atom_number_one = int(sys.argv[2]) - 1
atom_number_two = int(sys.argv[3]) - 1

#raising errors for invalid atom indexes
if (atom_number_one > sum(number_of_atom) - 1):
    print("Atom One Index Out Of Range")
    sys.exit()
elif (atom_number_one < 0):
    print("Invalid Atom One Index")
    sys.exit()
if (atom_number_two > sum(number_of_atom) - 1):
    print("Atom Two Index Out Of Range")
    sys.exit()
elif (atom_number_two < 0):
    print("Invalid Atom Two Index")
    sys.exit()

#initializing distances
periodic_distance = [[[0 for k in range(3)] for j in range(3)] for i in range(3)]
distance = []
avg_distance = []

#calculating distances
for frames in range(number_of_frames):
    x_atom1, y_atom1, z_atom1 = X[frames*sum(number_of_atom) + atom_number_one][0], Y[frames*sum(number_of_atom) + atom_number_one][0], Z[frames*sum(number_of_atom) + atom_number_one][0]
    x_atom2, y_atom2, z_atom2 = X[frames*sum(number_of_atom) + atom_number_two][0], Y[frames*sum(number_of_atom) + atom_number_two][0], Z[frames*sum(number_of_atom) + atom_number_two][0]
    for m in range(-1, 2):
        for n in range(-1, 2):
            for o in range(-1, 2):
                periodic_distance[m+1][n+1][o+1] = math.dist((x_atom1*a, y_atom1*b, z_atom1*c), ((x_atom2 + m)*a, (y_atom2 + n)*b, (z_atom2 + o)*c))
    distance += [min(min(min(periodic_distance)))]
    avg_distance += [sum(distance)/len(distance)]

#plotting and printing output
fig = plt.figure(figsize=(10, 4))
plt.plot(distance, color = "darkgrey")
plt.plot(avg_distance, "lightcoral", label = "Average", linewidth = 3)
plt.xlabel('Number of Frames')
plt.xlim([0, number_of_frames])
plt.ylabel(f"Bond Distance between Atom {atom_number_one + 1} and Atom {atom_number_two + 1}" +'in $\mathrm{\AA}$')
plt.legend(loc="upper right")
output_file = sys.argv[4]
plt.savefig(output_file)
