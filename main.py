import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def surface(x, y):
    return x**2 + y**2

def normal_vector(x, y):
    df_dx = 2 * x
    df_dy = 2 * y
    n = np.array([-df_dx, -df_dy, 1])
    n /= np.linalg.norm(n)
    return n

def reflect_velocity(v_before, n, k=1):
    v_after = v_before - 2 * np.dot(v_before, n) * n * np.sqrt(k)
    return v_after

def simulate_bounces(x0, y0, z0, v0x, v0y, v0z, g, num_bounces, k=1):
    t = 0
    bounces = 0

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    impact_times = []
    energy_data = []

    x_values = [x0]
    y_values = [y0]
    z_values = [z0]

    while bounces < num_bounces:
        x = x0 + v0x * t
        y = y0 + v0y * t
        z = z0 + v0z * t - 0.5 * g * t**2

        if z <= surface(x, y):
            a = 0.5 * g
            b = v0z
            c = z0 - surface(x0, y0)
            discriminant = b**2 - 4 * a * c
            t_impact = (-b - np.sqrt(discriminant)) / (2 * a)
            impact_times.append(t + t_impact)

            v0z = -v0z
            t = 0
            bounces += 1
            print(f"Bounce {bounces} - Impact Point: ({x}, {y}, {z}), Time: {impact_times[-1]} seconds")

            m = 1
            v_before = np.array([v0x, v0y, v0z - g * t])
            v_before_length = np.linalg.norm(v_before)

            v_before_normalized = v_before / v_before_length
            n = normal_vector(x, y)
            v_after = reflect_velocity(v_before_normalized, n, k)
            v0x, v0y, v0z = v_after * v_before_length

            # Energia ale coś źle liczy
            Ekin = 0.5 * m * (v_before_length ** 2)
            Epot = m * g * z
            Etotal = Ekin + Epot
            energy_data.append((Ekin, Epot, Etotal))

            # Update cord
            x0, y0, z0 = x, y, z

        x_values.append(x)
        y_values.append(y)
        z_values.append(z)
        t += 0.01

    ax.plot(x_values, y_values, z_values, c='r', linewidth=0.5)

    ax.set_xlim([-20, 20])
    ax.set_ylim([-20, 20])
    ax.set_zlim([0, 500])

    x_surface = np.linspace(-20, 20, 100)
    y_surface = np.linspace(-20, 20, 100)
    X, Y = np.meshgrid(x_surface, y_surface)
    Z = surface(X, Y)
    ax.plot_surface(X, Y, Z, alpha=0.5, cmap='viridis')

    plt.show()

    return energy_data

if __name__ == "__main__":
    x0, y0, z0 = -5, 0, 200
    v0x, v0y, v0z = 2, 4, 0
    g = 10
    num_bounces = 5 
    k = 0.8  

    energy_data = simulate_bounces(x0, y0, z0, v0x, v0y, v0z, g, num_bounces, k)

    print("Energy Data:")
    print("Bounce | Kinetic Energy | Potential Energy | Total Energy")
    for i, (Ekin, Epot, Etotal) in enumerate(energy_data):
        print(f"{i + 1} | {Ekin} | {Epot} | {Etotal}")
