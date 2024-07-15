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

#atoms, number of atoms, and list of atoms
atom = (dataset.iloc[[3], :]).values.tolist()[0]
number_of_atom = (dataset.iloc[[4], :]).values.tolist()[0]
atom_list = sum([[s] * n for s, n in zip(atom, number_of_atom)], [])

#coordinates
X = (dataset.iloc[5:,[0]]).values.tolist()
Y = (dataset.iloc[5:,[1]]).values.tolist()
Z = (dataset.iloc[5:,[2]]).values.tolist()

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
    
#reference atom as input for the script
index_reference_atom = int(sys.argv[2]) - 1

#raising errors for invalid atom indexes
if (index_reference_atom > sum(number_of_atom) - 1):
    print("Reference Atom Index Out Of Range")
    sys.exit()
elif (index_reference_atom < 0):
    print("Invalid Reference Atom Index")
    sys.exit()

#atom dictionary having vander Waal radius of different atoms
element_dict = {
    "H" : 1.20,
    "C" : 1.70,
    "N" : 1.55,
    "O" : 1.52,
    "F" : 1.47,
    "Si" : 2.10, 
    "P" : 1.80,
    "S" : 1.80,
    "Cl" : 1.75,
}

#initializing distances
periodic_distance = [[[0 for k in range(3)] for j in range(3)] for i in range(3)]
distance = [[0 for j in range(sum(number_of_atom))] for i in range(number_of_frames)]

#calculating distances
for frames in range(number_of_frames):
    x_reference_atom, y_reference_atom, z_reference_atom = X[frames*sum(number_of_atom) + index_reference_atom][0], Y[frames*sum(number_of_atom) + index_reference_atom][0], Z[frames*sum(number_of_atom) + index_reference_atom][0]
    for index_target_atom in range(sum(number_of_atom)):
        x_target_atom, y_target_atom, z_target_atom = X[frames*sum(number_of_atom) + index_target_atom][0], Y[frames*sum(number_of_atom) + index_target_atom][0], Z[frames*sum(number_of_atom) + index_target_atom][0]
        for m in range(-1, 2):
            for n in range(-1, 2):  
                for o in range(-1, 2):
                    periodic_distance[m+1][n+1][o+1] = math.dist((x_reference_atom*a, y_reference_atom*b, z_reference_atom*c), ((x_target_atom + m)*a, (y_target_atom + n)*b, (z_target_atom + o)*c))
        distance[frames][index_target_atom] = min(min(min(periodic_distance)))

total_distance = [0 for i in range(sum(number_of_atom))]
for frames in range(number_of_frames):
    total_distance = [sum(x) for x in zip(total_distance, distance[frames])]
avg_distance = [x/number_of_frames for x in total_distance]

#defining the grid size
start = element_dict.get(atom_list[index_reference_atom])
end = max(avg_distance)
num_points = int(sys.argv[3])
def grid(start, end, num_points):
    step = (end - start) / (num_points - 1)
    points = [start + i * step for i in range(num_points)]
    return points
dx = grid(start, end, num_points)

#initializing intersection volumes
volume = [0 for i in range(len(dx))]
volume_fraction = [0 for i in range(len(dx))]

#calculating intersection volumes
pi = math.pi
for i in range(0, len(dx)):
    V = 0
    for j in range(0, len(atom_list)):
        atom1 = atom_list[index_reference_atom]
        r1 = dx[i]
        atom2 = atom_list[j]
        r2 = element_dict.get(atom2)
        dist = avg_distance[j]
        if dist == 0:
            V += 0
        elif r1 + r2 <= dist:
            V += 0
        elif abs(r1-r2) >= dist:
            if r1 < r2:
                V += 4*pi*r1*r1*r1/3
            else:
                V += 4*pi*r2*r2*r2/3
        else:
            V += pi*((r1+r2-dist)**2)*(dist**2+2*dist*(r1+r2)+6*r1*r2-3*(r1**2+r2**2))/(12*dist)
    volume[i] = V
    volume_fraction[i] = V/((4*pi*dx[i]*dx[i]*dx[i])/3)

#volume fraction distribution
fig = plt.figure(figsize=(14, 6))
fig = plt.subplot(1, 2, 1)
plt.plot(dx, volume, '.')
plt.xlabel('Distance from the reference atom')
plt.ylabel('Average Volume Distribution in $\mathrm{\AA}^3$')
fig = plt.subplot(1, 2, 2)
plt.plot(dx, volume_fraction, '.', c = 'salmon')
plt.xlabel('Distance from the reference atom')
plt.ylabel('Average Volume Fraction Distribution')
output_file = sys.argv[4]
plt.savefig(output_file)

