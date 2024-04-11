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

with open(path + 'breakup_of_the_gap_hybrid_p.txt', 'r') as fp:
    data_hybrid_p = fp.read().split('\n')

with open(path + 'breakup_of_the_gap_hybrid_u.txt', 'r') as fp:
    data_hybrid_u = fp.read().split('\n')

with open(path + 'breakup_of_the_gap_analitical_solution_p.txt', 'r') as fp:
    data_analitical_p = fp.read().split('\n')

with open(path + 'breakup_of_the_gap_analitical_solution_u.txt', 'r') as fp:
    data_analitical_u = fp.read().split('\n')

with open(path + 'breakup_of_the_gap_point_O_p.txt', 'r') as fp:
    data_O_p = fp.read().split('\n')

with open(path + 'breakup_of_the_gap_point_O_u.txt', 'r') as fp:
    data_O_u = fp.read().split('\n')

with open(path + 'breakup_of_the_gap_point_H_p.txt', 'r') as fp:
    data_H_p = fp.read().split('\n')

with open(path + 'breakup_of_the_gap_point_H_u.txt', 'r') as fp:
    data_H_u = fp.read().split('\n')

#animate_file('transport_eq_hybrid_scheme')
#animate_file('transport_eq_point_O')
#animate_file('transport_eq_point_X')
#animate_file('transport_eq_point_Y')
#animate_file('transport_eq_point_Z')
#animate_file('transport_eq_point_H_2_order')
#animate_file('transport_eq_point_T_3_order')

x = np.linspace(-1.0, 1.0, len(data_analitical_p[0].split(', ')))

names_list = ['start figure', 'analitical solution', 'hybrid scheme']
fig, axs = plt.subplots(2, 1)
axs[0].plot(x, np.float_(data_hybrid_p[0].split(', ')))
axs[0].plot(x, np.float_(data_analitical_p[0].split(', ')))
axs[0].plot(x, np.float_(data_hybrid_p[len(data_hybrid_p) - 1].split(', ')))
axs[0].grid(True, linestyle='--')
axs[0].legend(names_list)
axs[0].set_title('hybrid scheme, p')
axs[0].set_xlim((-0.5, 0.5))
axs[0].set_xlabel('L')
axs[0].set_ylabel('p')

axs[1].plot(x, np.float_(data_hybrid_u[0].split(', ')))
axs[1].plot(x, np.float_(data_analitical_u[0].split(', ')))
axs[1].plot(x, np.float_(data_hybrid_u[len(data_hybrid_u) - 1].split(', ')))
axs[1].grid(True, linestyle='--')
axs[1].legend(names_list)
axs[1].set_title('hybrid scheme, u')
axs[1].set_xlim((-0.5, 0.5))
axs[1].set_xlabel('L')
axs[1].set_ylabel('u')
plt.show()

names_list = ['start figure', 'analitical solution', 'scheme 2 order (H point)']
fig, axs = plt.subplots(2, 1)
axs[0].plot(x, np.float_(data_H_p[0].split(', ')))
axs[0].plot(x, np.float_(data_analitical_p[0].split(', ')))
axs[0].plot(x, np.float_(data_H_p[len(data_H_p) - 1].split(', ')))
axs[0].grid(True, linestyle='--')
axs[0].legend(names_list)
axs[0].set_title('scheme 2 order (H point), p')
axs[0].set_xlim((-0.5, 0.5))
axs[0].set_xlabel('L')
axs[0].set_ylabel('p')

axs[1].plot(x, np.float_(data_H_u[0].split(', ')))
axs[1].plot(x, np.float_(data_analitical_u[0].split(', ')))
axs[1].plot(x, np.float_(data_H_u[len(data_H_u) - 1].split(', ')))
axs[1].grid(True, linestyle='--')
axs[1].legend(names_list)
axs[1].set_title('scheme 2 order (H point), u')
axs[1].set_xlim((-0.5, 0.5))
axs[1].set_xlabel('L')
axs[1].set_ylabel('u')
plt.show()

names_list = ['start figure', 'analitical solution', 'scheme 1 order (O point)']
fig, axs = plt.subplots(2, 1)
axs[0].plot(x, np.float_(data_O_p[0].split(', ')))
axs[0].plot(x, np.float_(data_analitical_p[0].split(', ')))
axs[0].plot(x, np.float_(data_O_p[len(data_O_p) - 1].split(', ')))
axs[0].grid(True, linestyle='--')
axs[0].legend(names_list)
axs[0].set_title('scheme 1 order (O point), p')
axs[0].set_xlim((-0.5, 0.5))
axs[0].set_xlabel('L')
axs[0].set_ylabel('p')

axs[1].plot(x, np.float_(data_O_u[0].split(', ')))
axs[1].plot(x, np.float_(data_analitical_u[0].split(', ')))
axs[1].plot(x, np.float_(data_O_u[len(data_O_u) - 1].split(', ')))
axs[1].grid(True, linestyle='--')
axs[1].legend(names_list)
axs[1].set_title('scheme 2 order (O point), u')
axs[1].set_xlim((-0.5, 0.5))
axs[1].set_xlabel('L')
axs[1].set_ylabel('u')
plt.show()

#x = np.linspace(0, 2.0, len(data_hybrid[0].split(', ')))
#plt.plot(x, np.float_(data_hybrid[0].split(', ')))
#plt.plot(x, np.float_(data_analitic[0].split(', ')))
#plt.plot(x, np.float_(data_O[len(data_O) - 1].split(', ')))
#plt.plot(x, np.float_(data_X[len(data_X) - 1].split(', ')))
#plt.plot(x, np.float_(data_Y[len(data_Y) - 1].split(', ')))
#plt.plot(x, np.float_(data_Z[len(data_Z) - 1].split(', ')))
#plt.plot(x, np.float_(data_H[len(data_H) - 1].split(', ')))
#plt.plot(x, np.float_(data_T[len(data_T) - 1].split(', ')))
#plt.plot(x, np.float_(data_hybrid[len(data_hybrid) - 1].split(', ')))
#plt.plot(x, np.float_(data_H1[len(data_H1) - 1].split(', ')))
#plt.plot(x, np.float_(data_H2[len(data_H2) - 1].split(', ')))
#plt.plot(x, np.float_(data_hybrid_4[len(data_hybrid_4) - 1].split(', ')))
#plt.plot(x, np.float_(data_hybrid_7[len(data_hybrid_7) - 1].split(', ')))
#plt.plot(x, np.float_(data_hybrid_3_and_1[len(data_hybrid_3_and_1) - 1].split(', ')))

#names_list = ['start figure', 'analitical solution', 'scheme point O (1 order)', 'scheme point X (1 order)', 
#              'scheme point Y (1 order)', 'scheme point Z (1 order)', 'scheme point H (2 order)', 'scheme point T (3 order)', 'hybrid scheme']
#names_list = ['start figure', 'analitical solution', 'hybrid scheme from 3 points (3 order, 2 order, 1 order)']

#plt.grid(True)

#plt.legend(names_list)
#plt.show()
