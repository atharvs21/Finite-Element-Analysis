import numpy as np
import matplotlib.pyplot as plt

def finite_element_solver(n_elements, degree, selection, bc_left, bc_right, a_coeff, c_coeff, f_coeff):
    def shape_function(degree, coordinates, coordinates_derivative, selection):
        if selection == 1:
            if degree == 1:
                coordinates = [[0.5, -0.5],
                               [0.5, 0.5]]
                coordinates_derivative = [[-0.5],
                                          [0.5]]
            elif degree == 2:
                coordinates = [[0, -0.5, 0.5],
                               [1, 0, -1],
                               [0, 0.5, 0.5]]
                coordinates_derivative = [[-0.5, 1],
                                          [0, -2],
                                          [0.5, 1]]
            elif degree == 3:
                coordinates = [[-0.0625, 0.0625, 0.5625, -0.5625],
                               [0.5625, -1.6875, -0.5625, 1.6875],
                               [0.5625, 1.6875, -0.5625, -1.6875],
                               [-0.0625, -0.0625, 0.5625, 0.5625]]
                coordinates_derivative = [[0.0625, 1.125, -1.6875],
                                          [-1.6875, -1.125, 5.0625],
                                          [1.6875, -1.125, -5.0625],
                                          [-0.0625, 1.125, 1.6875]]
            elif degree == 4:
                coordinates = [[0, 0.1667, -0.1667, -0.6667, 0.6667],
                               [0, -1.3333, 2.6667, 1.3333, -2.6667],
                               [1, 0, -5, 0, 4],
                               [0, 1.3333, 2.6667, -1.3333, -2.6667],
                               [0, -0.1667, -0.1667, 0.6667, 0.6667]]
                coordinates_derivative = [[0.1667, -0.3333, -2, 2.6667],
                                          [-1.3333, 5.3333, 4, -10.6667],
                                          [0, -10, 0, 16],
                                          [1.3333, 5.3333, -4, -10.6667],
                                          [-0.1667, -0.3333, 2, 2.6667]]
        else:
            if degree == 1:
                coordinates = [[0.5, -0.5],
                               [0.5, 0.5]]
                coordinates_derivative = [[-0.5],
                                          [0.5]]
            elif degree == 2:
                coordinates = [[0.5, -0.5, 0, 0, 0],
                               [-0.6124, 0, 0.6124, 0, 0],
                               [0.5, 0.5, 0, 0, 0]]
                coordinates_derivative = [[-0.5, 0],
                                          [0, 1.2247],
                                          [0.5, 0]]
            elif degree == 3:
                coordinates = [[0.5, -0.5, 0, 0, 0],
                               [0, -0.7906, 0, 0.7906, 0],
                               [-0.6124, 0, 0.6124, 0, 0],
                               [0.5, 0.5, 0, 0, 0]]
                coordinates_derivative = [[-0.5, 0, 0, 0],
                                          [-0.7906, 0, 2.3717, 0],
                                          [0, 1.2247, 0, 0],
                                          [0.5, 0, 0, 0]]
            elif degree == 4:
                coordinates = [[0.5, -0.5, 0, 0, 0],
                               [0.2338, 0, -1.4031, 0, 1.1693],
                               [0, -0.7906, 0, 0.7906, 0],
                               [-0.6124, 0, 0.6124, 0, 0],
                               [0.5, 0.5, 0, 0, 0]]
                coordinates_derivative = [[-0.5, 0, 0, 0],
                                          [0, -2.8062, 0, 4.6771],
                                          [-0.7906, 0, 2.3717, 0],
                                          [0, 1.2247, 0, 0],
                                          [0.5, 0, 0, 0]]
        return coordinates, coordinates_derivative

    def integration_rule(re):
        coordinates_integration = np.array([[0.23862, 0.46791], [-0.23862, 0.46791], [0.66121, 0.36076],
                                            [-0.66121, 0.36076], [0.93247, 0.17132], [-0.93247, 0.17132]])
        s = 0
        for i in range(0, 6):
            su = 0
            for j in range(0, int(re.shape[0])):
                su = su + (re[j] * coordinates_integration[i][0] ** j)
            s = s + su * coordinates_integration[i][1]
        return s

    def element_matrix(ae, c, f, p):
        for i in range(0, p + 1):
            for j in range(0, p + 1):
                ke[i][j] = (2 / element_length) * integration_rule(
                    np.polynomial.polynomial.polymul(np.polynomial.polynomial.polymul(cod[i], cod[j]), ae[0]))
                ge[i][j] = (element_length / 2) * integration_rule(
                    np.polynomial.polynomial.polymul(np.polynomial.polynomial.polymul(co[i], co[j]), c[0]))
            fe[i][0] = (element_length / 2) * integration_rule(
                np.polynomial.polynomial.polymul(co[i], f[0]))
        return ke, ge, fe

    def get_value(co, x, p, d):
        if d == 0:
            s = 0
            for i in range(0, p + 1):
                s = s + co[i] * x ** i
        else:
            s = 0
            for i in range(0, p):
                s = s + co[i] * x ** i
        return s

    co = np.zeros((5, 5))
    cod = np.zeros((5, 5))
    co, cod = shape_function(degree, co, cod, selection)

    if bc_left == 2:
        force_left = float(input("Enter force= "))
    elif bc_left == 1:
        displacement_left = 0
    else:
        spring_constant_left = float(input("Enter spring constant= "))
        deviation_left = float(input("Enter spring deviation= "))

    if bc_right == 2:
        force_right = 0
    elif bc_right == 1:
        displacement_right = float(input("Enter displacement= "))
    else:
        spring_constant_right = float(input("Enter spring constant= "))
        deviation_right = float(input("Enter spring deviation= "))

    element_length = 1 / n_elements
    node_location = np.zeros((n_elements, 2))

    for i in range(0, n_elements):
        for j in range(0, 2):
            if i == j == 0:
                continue
            elif j % 2 == 0:
                node_location[i][j] = node_location[i - 1][j + 1]
            else:
                node_location[i][j] = node_location[i][j - 1] + element_length

    stiffness_matrix = np.zeros((n_elements * degree + 1, n_elements * degree + 1))
    gradient_matrix = np.zeros((n_elements * degree + 1, n_elements * degree + 1))
    force_matrix = np.zeros((n_elements * degree + 1, 1))
    load_matrix = np.zeros((n_elements * degree + 1, 1))

    for i in range(0, n_elements):
        ke = np.zeros((degree + 1, degree + 1))
        ge = np.zeros((degree + 1, degree + 1))
        fe = np.zeros((degree + 1, 1))

        midpoint = (node_location[i][0] + node_location[i][1]) / 2

        ae_coefficient = np.array([[a_coeff[0][0] + a_coeff[0][1] * midpoint + a_coeff[0][2] * midpoint ** 2,
                                    a_coeff[0][1] * (element_length / 2) + a_coeff[0][2] * element_length * midpoint,
                                    a_coeff[0][2] * (element_length / 4)]])
        c_coefficient = np.array([[c_coeff[0][0] + c_coeff[0][1] * midpoint + c_coeff[0][2] * midpoint ** 2,
                                    c_coeff[0][1] * (element_length / 2) + c_coeff[0][2] * element_length * midpoint,
                                    c_coeff[0][2] * (element_length / 4)]])
        f_coefficient = np.array([[f_coeff[0][0] + f_coeff[0][1] * midpoint + f_coeff[0][2] * midpoint ** 2,
                                    f_coeff[0][1] * (element_length / 2) + f_coeff[0][2] * element_length * midpoint,
                                    f_coeff[0][2] * (element_length / 4)]])

        ke, ge, fe = element_matrix(ae_coefficient, c_coefficient, f_coefficient, degree)

        if i == 0:
            stiffness_matrix[:degree + 1, :degree + 1] += ke
            gradient_matrix[:degree + 1, :degree + 1] += ge
            force_matrix[:degree + 1, 0:1] += fe
        else:
            stiffness_matrix[i * degree:i * degree + degree + 1, i * degree:i * degree + degree + 1] += ke
            gradient_matrix[i * degree:i * degree + degree + 1, i * degree:i * degree + degree + 1] += ge
            force_matrix[i * degree:i * degree + degree + 1, 0:1] += fe

    stiffness_plus_gradient = stiffness_matrix + gradient_matrix

    if bc_left == 1:
        for i in range(1, n_elements * degree + 1):
            force_matrix[i] = force_matrix[i] - displacement_left * stiffness_plus_gradient[i][0]
        for i in range(0, n_elements * degree + 1):
            for j in range(0, n_elements * degree + 1):
                if i == 0 or j == 0:
                    stiffness_plus_gradient[i][j] = 0
        stiffness_plus_gradient[0][0] = 1
        force_matrix[0][0] = displacement_left
        load_matrix[0][0] = 0

    if bc_right == 1:
        for i in range(0, n_elements * degree):
            force_matrix[i] = force_matrix[i] - displacement_right * stiffness_plus_gradient[i][n_elements * degree]
        for i in range(0, n_elements * degree + 1):
            for j in range(0, n_elements * degree + 1):
                if i == n_elements * degree or j == n_elements * degree:
                    stiffness_plus_gradient[i][j] = 0
        stiffness_plus_gradient[n_elements * degree][n_elements * degree] = 1
        force_matrix[n_elements * degree][0] = displacement_right
        load_matrix[n_elements * degree][0] = 0

    if bc_left == 2:
        load_matrix[0][0] = force_left

    if bc_right == 2:
        load_matrix[n_elements * degree][0] = force_right

    if bc_left == 3:
        stiffness_plus_gradient[0][0] += spring_constant_left
        load_matrix[0][0] = spring_constant_left * deviation_left

    if bc_right == 3:
        stiffness_plus_gradient[n_elements * degree][n_elements * degree] += spring_constant_right
        load_matrix[n_elements * degree][0] = spring_constant_right * deviation_right

    displacement_solution = np.linalg.inv(stiffness_plus_gradient) @ (force_matrix + load_matrix)

    x_values = np.linspace(0, 1, 1000)
    y_values = []
    y_derivative_values = []

    for i in range(0, n_elements):
        for j in range(0, 1000):
            if x_values[j] >= node_location[i][0] and x_values[j] <= node_location[i][1]:
                auxiliary = (2 * x_values[j] - (node_location[i][0] + node_location[i][1])) / element_length
                fr = 0
                frd = 0
                b = 0
                for m in range(i * degree, i * degree + degree + 1):
                    fr = fr + (displacement_solution[m][0] * get_value(co[b], auxiliary, degree, 0))
                    frd = frd + (displacement_solution[m][0] * get_value(cod[b], auxiliary, degree, 1)) * (2 / element_length)
                    b = b + 1
                y_values.append(fr)
                y_derivative_values.append(frd)

    exact_values = []
    exact_derivative_values = []
    i = 0

    while i < 1:
        exact_values.append(-i ** 2 / 2 + (force_right + 1) * i + displacement_left)
        exact_derivative_values.append(-i + (force_right + 1))
        i = i + 0.001

    if np.size(y_values) < np.size(exact_values):
        for i in range(0, (np.size(exact_values) - np.size(y_values))):
            exact_values.pop()
            x_values = x_values[:-1]
    elif np.size(exact_values) < np.size(y_values):
        for i in range(0, (np.size(y_values) - np.size(exact_values))):
            y_values.pop()

    return y_values, y_derivative_values, exact_values, exact_derivative_values, x_values


number_of_elements = [1, 2, 5, 10, 100]
orders = [1, 2]
selection = 1
ae_coefficient = np.zeros((1, 3))
c_coefficient = np.zeros((1, 3))
f_coefficient = np.zeros((1, 3))
ae_order = 0

for i in range(0, ae_order + 1):
    ae_coefficient[0][i] = 1

c_order = 0

for i in range(0, c_order + 1):
    c_coefficient[0][i] = 0

f_order = 0

for i in range(0, f_order + 1):
    f_coefficient[0][i] = 1

bc_left = 1
bc_right = 2

for order in orders:
    y_values, y_derivative_values, exact_values, exact_derivative_values, x_values = finite_element_solver(1, order, selection, bc_left, bc_right, ae_coefficient, c_coefficient, f_coefficient)
    err=[exact_values[i]-exact_values[i] for i in range(len(exact_values))]
    line1, = plt.plot(x_values, err, label="Exact solution")
    err=[exact_values[i]-y_values[i] for i in range(len(exact_values))]
    line2, = plt.plot(x_values, err, label="1 element")
    y_values, y_derivative_values, exact_values, exact_derivative_values, x_values = finite_element_solver(2, order, selection, bc_left, bc_right, ae_coefficient, c_coefficient, f_coefficient)
    err=[exact_values[i]-y_values[i] for i in range(len(exact_values))]
    line3, = plt.plot(x_values, err, label="2 element")
    y_values, y_derivative_values, exact_values, exact_derivative_values, x_values = finite_element_solver(5, order, selection, bc_left, bc_right, ae_coefficient, c_coefficient, f_coefficient)
    err=[exact_values[i]-y_values[i] for i in range(len(exact_values))]
    line4, = plt.plot(x_values, err, label="5 element")
    y_values, y_derivative_values, exact_values, exact_derivative_values, x_values = finite_element_solver(10, order, selection, bc_left, bc_right, ae_coefficient, c_coefficient, f_coefficient)
    err=[exact_values[i]-y_values[i] for i in range(len(exact_values))]
    line5, = plt.plot(x_values, err, label="10 element")
    y_values, y_derivative_values, exact_values, exact_derivative_values, x_values = finite_element_solver(100, order, selection, bc_left, bc_right, ae_coefficient, c_coefficient, f_coefficient)
    err=[exact_values[i]-y_values[i] for i in range(len(exact_values))]
    line6, = plt.plot(x_values, err, label="100 element")
    plt.legend(handles=[line1, line2, line3, line4, line5, line6])
    plt.title("Plot of error between exact solution and FEM solution")
    plt.xlabel("x")
    plt.ylabel("Difference in values(Error)")
    plt.show()

for order in orders:
    line1, = plt.plot(x_values, exact_derivative_values, label="Exact solution")
    y_values, y_derivative_values, exact_values, exact_derivative_values, x_values = finite_element_solver(1, order, selection, bc_left, bc_right, ae_coefficient, c_coefficient, f_coefficient)
    line2, = plt.plot(x_values, y_derivative_values, label="1 element")
    y_values, y_derivative_values, exact_values, exact_derivative_values, x_values = finite_element_solver(2, order, selection, bc_left, bc_right, ae_coefficient, c_coefficient, f_coefficient)
    line3, = plt.plot(x_values, y_derivative_values, label="2 element")
    y_values, y_derivative_values, exact_values, exact_derivative_values, x_values = finite_element_solver(5, order, selection, bc_left, bc_right, ae_coefficient, c_coefficient, f_coefficient)
    line4, = plt.plot(x_values, y_derivative_values, label="5 element")
    y_values, y_derivative_values, exact_values, exact_derivative_values, x_values = finite_element_solver(10, order, selection, bc_left, bc_right, ae_coefficient, c_coefficient, f_coefficient)
    line5, = plt.plot(x_values, y_derivative_values, label="10 element")
    y_values, y_derivative_values, exact_values, exact_derivative_values, x_values = finite_element_solver(100, order, selection, bc_left, bc_right, ae_coefficient, c_coefficient, f_coefficient)
    line6, = plt.plot(x_values, y_derivative_values, label="100 element")
    plt.legend(handles=[line1, line2, line3, line4, line5, line6])
    plt.title("Plot of derivatives")
    plt.xlabel("x")
    plt.ylabel("Force")
    plt.show()

errors = []
log_n = []
energies = []

for n in number_of_elements:
    y_values, y_derivative_values, exact_values, exact_derivative_values, x_values = finite_element_solver(n, 1, selection, bc_left, bc_right, ae_coefficient, c_coefficient, f_coefficient)
    error = np.log(np.linalg.norm(np.array(exact_values) - np.array(y_values)))
    errors.append(error)
    log_n.append(np.log(n))
    total_energy = 0

    for i in range(0, 998):
        fx1 = y_derivative_values[i] ** 2
        fx2 = y_derivative_values[i + 1] ** 2
        total_energy = total_energy + 0.5 * (fx1 + fx2) * (x_values[i + 1] - x_values[i])

    energies.append(total_energy)

print(log_n, errors)

line1, = plt.plot(number_of_elements, energies, label="linear", marker=".")
errors = []
log_n = []
energies = []
exact_energies = []

for n in number_of_elements:
    y_values, y_derivative_values, exact_values, exact_derivative_values, x_values = finite_element_solver(n, 2, selection, bc_left, bc_right, ae_coefficient, c_coefficient, f_coefficient)
    error = np.log(np.linalg.norm(np.array(exact_values) - np.array(y_values)))
    errors.append(error)
    log_n.append(np.log(n))
    total_energy = 0

    for i in range(0, 998):
        fx1 = y_derivative_values[i] ** 2
        fx2 = y_derivative_values[i + 1] ** 2
        total_energy = total_energy + 0.5 * (fx1 + fx2) * (x_values[i + 1] - x_values[i])

    energies.append(total_energy)

    total_exact_energy = 0

    for i in range(0, 998):
        fx1 = exact_derivative_values[i] ** 2
        fx2 = exact_derivative_values[i + 1] ** 2
        total_exact_energy = total_exact_energy + 0.5 * (fx1 + fx2) * (x_values[i + 1] - x_values[i])

    exact_energies.append(total_exact_energy)

line2, = plt.plot(number_of_elements, energies, label="quadratic", marker=".")
plt.legend(handles=[line1, line2])
plt.xlabel("Number of elements")
plt.ylabel("Energy")
plt.title("Strain energy by approximating degree")
plt.show()

line1, = plt.plot(number_of_elements, energies, label="FEM solution", marker=".")
line2, = plt.plot(number_of_elements, exact_energies, label="Exact solution", marker=".")
plt.legend(handles=[line1, line2])
plt.xlabel("Number of elements")
plt.ylabel("Energy")
plt.title("Strain Energy of FEM solution and Exact solution")
plt.show()