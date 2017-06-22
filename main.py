# Finite Element analysis of bar supported at one end, subject to an axial force
# 2017, June 22 @ 7:28PM
# c0dew0rth.blogspot.com

import gauss_seidel
import matplotlib.pyplot as plt
import numpy as np
import copy

# LINEAR VARIATION OF DIAMETER
E = 200000  # in N/mm2
force = float(input("Force(N): "))  # in Newtons
length = float(input("Length(m): ")) * 1000  # in mm
num_elements = int(input("Number of elements: "))
diaI = float(input("Diameter of first end(mm): "))  # in mm
diaII = float(input("Diameter of second end(mm): "))  # in mm
elementaldias = []
for i in range(0, num_elements):  # create list of elemental diameters
    x_i = i * length / num_elements  # x of first node of this element
    x_iplus1 = (i + 1) * length / num_elements  # x of second node of this element
    del_i = x_i / length * (diaII - diaI)  # dia corrxn of first node of element
    del_iplus1 = x_iplus1 / length * (diaII - diaI)  # dia corrxn of second node of element
    del_avg = (del_i + del_iplus1) / 2  # average dia corrxn of this element
    dia = diaI + del_avg  # dia of this element
    elementaldias.append(dia)
print("Elemental diameters created:")
print(elementaldias)
k = [[0 for i in range(0, num_elements + 1)] for i in
     range(0, num_elements + 1)]  # stiffness matrix, order=number of nodes=num_elements+1
for i in range(0, num_elements):  # create stiffness matrix
    k_i = E * (3.14 * pow(elementaldias[i], 2) / 4) / (length / num_elements)
    if i == 0:
        k[i][i] = k_i
    else:
        k[i][i] = k_old + k_i
    k[i][i + 1] = -k_i
    k[i + 1][i] = -k_i
    k[i + 1][i + 1] = k_i
    k_old = k_i
F = [0 for i in range(0, num_elements + 1)]  # force vector
F[num_elements] = force  # assign input force to the last node
# Application of Boundary condition (u1=0)
# i.e. eliminating first row of F and first row and first column of k matrix
F_cpy = copy.deepcopy(F)
k_cpy = copy.deepcopy(k)
del F_cpy[0]
del k_cpy[0]  # delete first row
for i in range(0, num_elements):
    del k_cpy[i][0]  # delete first column member of the row
# By now, k_cpy matrix is an order smaller than k due to the application of boundary condition
u = gauss_seidel.gauss_seidel_solve(num_elements, k_cpy, F_cpy)  # list(np.linalg.solve(k_cpy,F_cpy))
u.insert(0, 0)  # prepend the boundary condition displacement i.e. 0 to the solution
print("Nodal displacements: ")
print(u)
X = np.linspace(length / num_elements, length, num_elements + 1)
Y = np.array(u)
forces = []
for i in range(0, num_elements + 1):
    sum = 0
    for j in range(0, num_elements + 1):
        sum += k[i][j] * u[j]
    forces.append(sum)
print("Nodal forces: ")
print(forces)
plt.plot(X, Y)
plt.show()

# #LINEAR VARIATION OF AREA
# E = 200000  # in N/mm2
# force = float(input("Force(N): "))  # in Newtons
# length = float(input("Length(m): ")) * 1000  # in mm
# num_elements = int(input("Number of elements: "))
# diaI = float(input("Diameter of first end(mm): "))  # in mm
# diaII = float(input("Diameter of second end(mm): "))  # in mm
# elementalareas = []
# for i in range(0, num_elements):  # create list of elemental diameters
#     x_i = i * length / num_elements  # x of first node of this element
#     x_iplus1 = (i + 1) * length / num_elements  # x of second node of this element
#     del_i = x_i / length * ((diaII ** 2 - diaI ** 2) * 3.14 / 4)  # area corrxn of first node of element
#     del_iplus1 = x_iplus1 / length * ((diaII ** 2 - diaI ** 2) * 3.14 / 4)  # area corrxn of second node of element
#     del_avg = (del_i + del_iplus1) / 2  # average area corrxn of this element
#     area = (diaI ** 2) * 3.14 / 4 + del_avg  # dia of this element
#     elementalareas.append(area)
# print("Elemental areas created:")
# print(elementalareas)
# k = [[0 for i in range(0, num_elements + 1)] for i in
#      range(0, num_elements + 1)]  # stiffness matrix, order=number of nodes=num_elements+1
# for i in range(0, num_elements):  # create stiffness matrix
#     k_i = E * (elementalareas[i]) / (length / num_elements)
#     if i == 0:
#         k[i][i] = k_i
#     else:
#         k[i][i] = k_i_old + k_i
#     k[i][i + 1] = -k_i
#     k[i + 1][i] = -k_i
#     k[i + 1][i + 1] = k_i
#     k_i_old = k_i
# F = [0 for i in range(0, num_elements + 1)]  # force vector
# F[num_elements] = force  # assign input force to the last node
# # Application of Boundary condition (u1=0)
# # i.e. eliminating first row of F and first row and first column of k matrix
# F_cpy = copy.deepcopy(F)
# k_cpy = copy.deepcopy(k)
# del F_cpy[0]
# del k_cpy[0]  # delete first row
# for i in range(0, num_elements):
#     del k_cpy[i][0]  # delete first column member of the row
# # By now, k_cpy matrix is an order smaller than k due to the application of boundary condition
# u = list(np.linalg.solve(k_cpy, F_cpy))  # gauss_seidel.gauss_seidel_solve(num_elements, k_cpy, F_cpy)
# u.insert(0, 0)  # prepend the boundary condition displacement i.e. 0 to the solution
# print("Nodal displacements: ")
# print(u)
# X = np.linspace(length / num_elements, length, num_elements + 1)
# Y = np.array(u)
# forces = []
# for i in range(0, num_elements + 1):
#     sum = 0
#     for j in range(0, num_elements + 1):
#         sum += k[i][j] * u[j]
#     forces.append(sum)
# print("Nodal forces: ")
# print(forces)
# plt.plot(X, Y)
# plt.show()
