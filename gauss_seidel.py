
# Example:
# a = [[8, -3, 2], [4, 11, -1], [6, 3, 12]]
# b = [20, 33, 35]
# numberofunknowns = 3
#
def gauss_seidel_solve(numberofunknowns,a,b):
    x = [0 for i in range(numberofunknowns)]  # initial guess starting from x1
    x_new = [0 for i in range(numberofunknowns)]  # used in loop
    error = 1  # any number above tolerance will do
    tolerance = 0.000001

    while error > tolerance:
        for i in range(numberofunknowns):
            remaindersum = 0
            for j in range(numberofunknowns):
                if i != j:
                    remaindersum += a[i][j] * x_new[j]
            x_new[i] = (b[i] - remaindersum) / float(a[i][i])
        error = max([abs(x_new[i] - x[i]) for i in range(numberofunknowns)])
        x = list(x_new)
    return x
