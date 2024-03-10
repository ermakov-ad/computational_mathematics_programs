import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

path = '/home/sq-wm/Документы/Comp_math_8_sem/'

def animate_file(filename):
    
    def init():
        line.set_data([], [])
        return line,

    def animate(i):
        y = np.float_(data[i].split(', '))
        x = np.linspace(0, 2.0, len(y))
        line.set_data(x, y)
        return line,

    with open(path + str(filename) + '.txt', 'r') as fp:
        data = fp.read().split('\n')

    fig = plt.figure()
    ax = plt.axes(xlim=(0, 2.0), ylim=(-0.1, 1.1))
    ax.set_title(str(filename))
    line, = ax.plot([], [], lw=3)
    plt.grid(True)

    anim = FuncAnimation(fig, animate, init_func=init, frames=len(data), interval=100, blit=True)
    anim.save(path + str(filename) + '.gif', writer='imagemagick', fps=60)

with open(path + 'transport_eq_hybrid_scheme.txt', 'r') as fp:
    data_hybrid = fp.read().split('\n')

with open(path + 'transport_equation_analitical_solution.txt', 'r') as fp:
    data_analitic = fp.read().split('\n')

with open(path + 'transport_eq_point_O.txt', 'r') as fp:
    data_O = fp.read().split('\n')

with open(path + 'transport_eq_point_X.txt', 'r') as fp:
    data_X = fp.read().split('\n')

with open(path + 'transport_eq_point_Y.txt', 'r') as fp:
    data_Y = fp.read().split('\n')

with open(path + 'transport_eq_point_Z.txt', 'r') as fp:
    data_Z = fp.read().split('\n')

with open(path + 'transport_eq_point_H_2_order.txt', 'r') as fp:
    data_H = fp.read().split('\n')

with open(path + 'transport_eq_point_T_3_order.txt', 'r') as fp:
    data_T = fp.read().split('\n')

animate_file('transport_eq_hybrid_scheme')
animate_file('transport_eq_point_O')
animate_file('transport_eq_point_X')
animate_file('transport_eq_point_Y')
animate_file('transport_eq_point_Z')
animate_file('transport_eq_point_H_2_order')
animate_file('transport_eq_point_T_3_order')

plt.axes(xlim=(0.0, 2.0), ylim=(-0.1, 1.1))
x = np.linspace(0, 2.0, len(data_hybrid[0].split(', ')))
plt.plot(x, np.float_(data_hybrid[0].split(', ')))
plt.plot(x, np.float_(data_analitic[0].split(', ')))
plt.plot(x, np.float_(data_O[len(data_O) - 1].split(', ')))
plt.plot(x, np.float_(data_X[len(data_X) - 1].split(', ')))
plt.plot(x, np.float_(data_Y[len(data_Y) - 1].split(', ')))
plt.plot(x, np.float_(data_Z[len(data_Z) - 1].split(', ')))
plt.plot(x, np.float_(data_H[len(data_H) - 1].split(', ')))
plt.plot(x, np.float_(data_T[len(data_T) - 1].split(', ')))
plt.plot(x, np.float_(data_hybrid[len(data_hybrid) - 1].split(', ')))

names_list = ['start figure', 'analitical solution', 'scheme point O (1 order)', 'scheme point X (1 order)', 
              'scheme point Y (1 order)', 'scheme point Z (1 order)', 'scheme point H (2 order)', 'scheme point T (3 order)', 'hybrid scheme']

plt.grid(True)
plt.title('schemes, CFL = 0.2')
plt.xlabel('L')
plt.ylabel('U')
plt.legend(names_list)
plt.show()
