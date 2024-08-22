import numpy as np
import matplotlib.pyplot as plt

def solve_beam_finite_element(n_elements, bc_start, bc_end, force_coefficients):
    def shape_functions(co,cod,codd,coddd,h):
        # Coefficients of shape functions
        co = [[0.5, -0.75, 0, 0.25],[-0.125*h, 0.125*h, 0.125*h, -0.125*h],
              [0.5, 0.75, 0, -0.25],[0.125*h, 0.125*h, -0.125*h, -0.125*h]]
        cod = [[-0.75, 0, 0.75], [0.125*h, 0.25*h, -0.375*h],
               [0.75, 0, -0.75],[0.125*h, -0.25*h, -0.375*h]]
        codd = [[0, 1.5], [0.25*h, -0.75*h], [0, -1.5], [-0.25*h, -0.75*h]]
        coddd = [[1.5], [-0.75*h], [-1.5], [-0.75*h]]
        return co, cod, codd, coddd

    def integrate(re):
        # Integration approximating coefficients
        coi = np.array([[0.23862, 0.46791],[-0.23862, 0.46791],[0.66121, 0.36076],
                        [-0.66121, 0.36076],[0.93247, 0.17132],[-0.93247, 0.17132]])
        s = 0
        for i in range(0, 6):
            su = 0
            for j in range(0, int(re.shape[0])):
                su = su + (re[j] * coi[i][0]**j)
            s = s + su * coi[i][1]
        return s

    def element_matrix(f):
      # Polynomial multiplication of coefficients in an elemental matrix
      for i in range(0, 4):
          element_forces[i][0] = (h / 2) * integrate(np.polynomial.polynomial.polymul(co[i], f[0]))
      return element_forces

    def get_value(co, x, p, d):
        # Mapping of element position to global matrix position
        if d == 0:
            s = 0
            for i in range(0, p+1):
                s = s + co[i] * x**i
        elif d == 1:
            s = 0
            for i in range(0, p):
                s = s + co[i] * x**i
        elif d == 2:
            s = 0
            for i in range(0, p-1):
                s = s + co[i] * x**i
        else:
            s = 0
            for i in range(0, p-2):
                s = s + co[i] * x**i
        return s

    h = 1 / n_elements
    co = np.zeros((5, 5))
    cod = np.zeros((5, 5))
    codd = np.zeros((5, 5))
    coddd = np.zeros((5, 5))
    co, cod, codd, coddd = shape_functions(co, cod, codd, coddd, h)

    if bc_start == 2:
        force1 = float(input("Enter force= "))
        moment1 = float(input("Enter moment= "))
    else:
        displacement1 = 0
        slope1 = 0
    if bc_end == 2:
        force2 = 5
        moment2 = 5
    else:
        displacement2 = float(input("Enter displacement= "))
        slope2 = float(input("Enter slope= "))

    node_locations = np.zeros((n_elements, 2))
    for i in range(0, n_elements):
        for j in range(0, 2):
            if (i == j == 0):
                continue
            elif (j % 2 == 0):
                node_locations[i][j] = node_locations[i - 1][j + 1]
            else:
                node_locations[i][j] = node_locations[i][j - 1] + h

    stiffness_matrix = np.zeros((2 * n_elements + 2, 2 * n_elements + 2))
    forces = np.zeros((2 * n_elements + 2, 1))
    loads = np.zeros((2 * n_elements + 2, 1))
    element_stiffness = np.zeros((4, 4))
    element_stiffness += np.array([[6, -3*h, -6, -3*h], [-3*h, 2*h*h, 3*h, h*h],
                    [-6, 3*h, 6, 3*h], [-3*h, h*h, 3*h, 2*h*h]])
    element_stiffness = (1000 / (3 * h ** 3)) * element_stiffness

    for i in range(0, n_elements):
        element_forces = np.zeros((4, 1))
        sunod = (node_locations[i][0] + node_locations[i][1]) / 2
        forces_of_element = np.array([[force_coefficients[0][0] + force_coefficients[0][1] * sunod +
                                       force_coefficients[0][2] * sunod * 2,
                                       force_coefficients[0][1] * (h / 2) + force_coefficients[0][2] * h * sunod,
                                       force_coefficients[0][2] * (h / 4)]])
        element_forces = element_matrix(forces_of_element)
        if (i == 0):
            stiffness_matrix[:4, :4] += element_stiffness
            forces[:4, 0:1] += element_forces
        else:
            stiffness_matrix[2 * i:2 * (i + 2), 2 * i:2 * (i + 2)] += element_stiffness
            forces[2 * i:2 * (i + 2), 0:1] += element_forces

    stiffness_matrix_final = stiffness_matrix

    if (bc_start == 1):
        for i in range(1, 2 * n_elements + 2):
            forces[i] = forces[i] - displacement1 * stiffness_matrix_final[i][0]
        for i in range(0, 2 * n_elements + 2):
            for j in range(0, n_elements * 2 + 2):
                if (i == 0 or j == 0):
                    stiffness_matrix_final[i][j] = 0
        stiffness_matrix_final[0][0] = 1
        forces[0][0] = displacement1
        loads[0][0] = 0
        for i in range(2, 2 * n_elements + 2):
            forces[i] = forces[i] - slope1 * stiffness_matrix_final[i][1]
        for i in range(0, 2 * n_elements + 2):
            for j in range(0, 2 * n_elements + 2):
                if (i == 1 or j == 1):
                    stiffness_matrix_final[i][j] = 0
        stiffness_matrix_final[1][1] = 1
        forces[1][0] = slope1
        loads[1][0] = 0

    if (bc_end == 1):
        for i in range(0, 2 * n_elements + 1):
            forces[i] = forces[i] - slope2 * stiffness_matrix_final[i][2 * n_elements + 1]
        for i in range(0, n_elements * 2 + 2):
            for j in range(0, 2 * n_elements + 2):
                if (i == 2 * n_elements + 1 or j == n_elements * 2 + 1):
                    stiffness_matrix_final[i][j] = 0
        stiffness_matrix_final[n_elements * 2 + 1][n_elements * 2 + 1] = 1
        forces[n_elements * 2 + 1][0] = slope2
        loads[n_elements * 2 + 1][0] = 0
        for i in range(0, 2 * n_elements):
            forces[i] = forces[i] - displacement2 * stiffness_matrix_final[i][2 * n_elements]
        for i in range(0, n_elements * 2 + 2):
            for j in range(0, 2 * n_elements + 2):
                if (i == 2 * n_elements or j == n_elements * 2):
                    stiffness_matrix_final[i][j] = 0
        stiffness_matrix_final[n_elements * 2][n_elements * 2] = 1
        forces[n_elements * 2][0] = displacement2
        loads[n_elements * 2][0] = 0

    if (bc_start == 2):
        loads[0][0] = force1
        loads[1][0] = -moment1
    if (bc_end == 2):
        loads[2 * n_elements][0] = force2
        loads[2 * n_elements + 1][0] = -moment2

    U = np.linalg.inv(stiffness_matrix_final) @ (forces + loads)

    x = np.linspace(0, 1, 1000)
    y_deflection = []
    y_slope = []
    moment = []
    shear_force = []
    stress = []

    for i in range(0, n_elements):
        for j in range(0, 1000):
            if x[j] >= node_locations[i][0] and x[j] <= node_locations[i][1]:
                aux = (2 * x[j] - (node_locations[i][0] + node_locations[i][1])) / h
                w = 0
                wslope = 0
                mom = 0
                sh = 0
                st = 0
                b = 0
                for m in range(i * 2, i * 2 + 4):
                    w = w + (U[m][0] * get_value(co[b], aux, 3, 0))
                    wslope = wslope + (U[m][0] * get_value(cod[b], aux, 3, 1)) * (2 / h)
                    mom = mom + (500 / 3) * (U[m][0] * get_value(codd[b], aux, 3, 2)) * (2 / h) ** 2
                    sh = sh + (-500 / 3) * (U[m][0] * get_value(coddd[b], aux, 3, 3)) * (2 / h) ** 3
                    st = st + (-1000) * (U[m][0] * get_value(codd[b], aux, 3, 2)) * (2 / h) ** 2
                    b = b + 1
                y_deflection.append(w)
                y_slope.append(wslope)
                moment.append(mom)
                shear_force.append(sh)
                stress.append(st)

    y_exact_deflection = []
    y_exact_slope = []
    i = 0
    while i <= 1:
        y_exact_deflection.append((3 / 500) * ((i ** 4) / 24 - (i ** 3) / 3 + (5 * i * i)/ 4))
        y_exact_slope.append(((i ** 3) / 6 - i * i + 5 * i / 2) * (3 / 500))
        i = i + 0.001

    error = 0
    for i in range(0, 999):
        error = error + ((y_exact_deflection[i] - y_deflection[i]) ** 2)

    rmse_error = np.sqrt(error / 1000) * 100
    print("RMSE error% = ", rmse_error, "%")

    if (np.size(x) < np.size(y_deflection)):
        for i in range(0, (np.size(y_deflection) - np.size(x))):
            y_deflection.pop()
            y_slope.pop()
    elif (np.size(y_deflection) < np.size(x)):
        for i in range(0, (np.size(x) - np.size(y_deflection))):
            x = x[:-1]
            y_exact_deflection = y_exact_deflection[:-1]
            y_exact_slope = y_exact_slope[:-1]

    return y_deflection, y_exact_deflection, y_slope, y_exact_slope, moment, shear_force, stress, x

n_elements_list = [1, 4, 10, 50, 100]
force_coefficients = np.zeros((1, 3))
y_deflection_final = []
y_exact_deflection_final = []
y_slope_final = []
y_exact_slope_final = []
moment_final = []
shear_force_final = []
stress_final = []
x_final = []
force_coefficients[0][0] = 1

for i in range(0, 5):
    y_deflection, y_exact_deflection, y_slope, y_exact_slope, moment, shear_force, stress, x = \
        solve_beam_finite_element(n_elements_list[i], 1, 2, force_coefficients)
    y_deflection_final.append(y_deflection)
    y_exact_deflection_final.append(y_exact_deflection)
    y_slope_final.append(y_slope)
    y_exact_slope_final.append(y_exact_slope)
    moment_final.append(moment)
    shear_force_final.append(shear_force)
    stress_final.append(stress)
    x_final.append(x)

plt.figure()
plt.plot(x_final[0], y_deflection_final[0], label="1 element", color='blue')
plt.plot(x_final[1], y_deflection_final[1], label="4 element", color='green')
plt.plot(x_final[2], y_deflection_final[2], label="10 element", color='red')
plt.plot(x_final[3], y_deflection_final[3], label="50 element", color='orange')
plt.plot(x_final[4], y_deflection_final[4], label="100 element", color='purple')
plt.plot(x_final[0], y_exact_deflection_final[0], label="Exact", color='black', linestyle='dashed')
plt.legend()
plt.title("Deflection with x")
plt.ylabel("Deflection")
plt.xlabel("x")
plt.show()

plt.figure()
plt.plot(x_final[0], y_slope_final[0], label="1 element", color='red')
plt.plot(x_final[1], y_slope_final[1], label="4 element", color='purple')
plt.plot(x_final[2], y_slope_final[2], label="10 element", color='orange')
plt.plot(x_final[3], y_slope_final[3], label="50 element", color='blue')
plt.plot(x_final[4], y_slope_final[4], label="100 element", color='green')
plt.plot(x_final[0], y_exact_slope_final[0], label="Exact", color='black', linestyle='dashed')
plt.legend()
plt.title("Slope with x")
plt.ylabel("Slope")
plt.xlabel("x")
plt.show()

plt.figure()
plt.plot(x_final[0], moment_final[0], label="1 element", color='red')
plt.plot(x_final[1], moment_final[1], label="4 element", color='purple')
plt.plot(x_final[2], moment_final[2], label="10 element", color='orange')
plt.plot(x_final[3], moment_final[3], label="50 element", color='blue')
plt.plot(x_final[4], moment_final[4], label="100 element", color='green')
plt.plot(x_final[4], moment_final[4], label="Exact", color='black', linestyle='dashed')
plt.legend()
plt.title("Moment with x")
plt.ylabel("Moment")
plt.xlabel("x")
plt.show()

plt.figure()
plt.plot(x_final[0], shear_force_final[0], label="1 element", color='red')
plt.plot(x_final[1], shear_force_final[1], label="4 element", color='purple')
plt.plot(x_final[2], shear_force_final[2], label="10 element", color='orange')
plt.plot(x_final[3], shear_force_final[3], label="50 element", color='blue')
plt.plot(x_final[4], shear_force_final[4], label="100 element", color='green')
plt.plot(x_final[4], shear_force_final[4], label="Exact", color='black', linestyle='dashed')
plt.legend()
plt.title("Shear force with x")
plt.ylabel("Shear force")
plt.xlabel("x")
plt.show()

plt.figure()
plt.plot(x_final[0], stress_final[0], label="1 element", color='red')
plt.plot(x_final[1], stress_final[1], label="4 element", color='purple')
plt.plot(x_final[2], stress_final[2], label="10 element", color='orange')
plt.plot(x_final[3], stress_final[3], label="50 element", color='blue')
plt.plot(x_final[4], stress_final[4], label="100 element", color='green')
plt.plot(x_final[4], stress_final[4], label="Exact", color='black', linestyle='dashed')
plt.legend()
plt.title("Stress with x")
plt.ylabel("Stress (MPa)")
plt.xlabel("x")
plt.show()
