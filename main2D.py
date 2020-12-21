from utils import *
import matplotlib.animation as animation
import matplotlib.pyplot as plt
import random

def step_simulation(q, q_dot, old_p, density, dt, num_axis, grid_res, pix_size, g):
    avection(q, q_dot, dt)
    np.clip(q, 0, pix_size * grid_res, q)
    external_force(q_dot, g, dt)
    staggered_grids = project_to_staggered_grid(q, q_dot, num_axis, pix_size, grid_res)
    for i, grid in enumerate(staggered_grids):
        #filter out the unwanted speed component
        filter_speed(grid, i, num_axis)
    flattened_staggered_grids = flatten(staggered_grids)
    RHS = B @ flattened_staggered_grids
    new_p = conj_grad(A, old_p, RHS)
    staggered_updater = get_speed_updater(D, new_p, dt, density, num_axis, grid_res, pix_size)
    speed_updater = unproject_from_staggered_grid(q, staggered_updater, num_axis, pix_size)
    q_dot += speed_updater
    return q, q_dot

def step_simulation2(q, q_dot, old_p, density, dt, num_axis, grid_res, pix_size, g):
    avection(q, q_dot, dt)
    np.clip(q, 0, pix_size * grid_res, q)
    external_force(q_dot, g, dt)
    staggered_grids = project_to_staggered_grid(q, q_dot, num_axis, pix_size, grid_res)
    for i, grid in enumerate(staggered_grids):
        #filter out the unwanted speed component
        filter_speed(grid, i, num_axis)
    flattened_staggered_grids = flatten(staggered_grids)
    RHS = B @ flattened_staggered_grids
    new_p = conj_grad(A, old_p, RHS)
    staggered_updater = get_speed_updater(D, new_p, dt, density, num_axis, grid_res, pix_size)
    for g, u in zip(staggered_grids, staggered_updater):
        g += u
    q_dot = unproject_from_staggered_grid(q, staggered_grids, num_axis, pix_size)
    return q, q_dot

num_particles = 5000
box_size = 1
grid_res = 100
pix_size = box_size / grid_res
num_axis = 2
dt = 0.01
num_axis = 2
g = np.asarray([0, -0.98])
density = 1

q = np.zeros((num_particles, num_axis))
q[:,1] = (np.random.rand(num_particles) - 0.5) * 0.3 + 0.5
q *= box_size
q_dot = np.zeros((num_particles, num_axis))
q_dot[:,0] = np.random.rand(num_particles)

B = construct_B(pix_size, grid_res, num_axis)
D = construct_D(pix_size, grid_res, num_axis)
A = B @ D

p = np.random.rand(grid_res ** num_axis)

q, qdot = step_simulation(q, q_dot, p, density, dt, num_axis, grid_res, pix_size, g)

fig = plt.figure()
ax = plt.axes(xlim=(0, box_size), ylim=(0, box_size))
particles, = ax.plot(q[:5,0], q[:5,1], 'o', ms=5)

def animate(i):
    q[:i*5], qdot[:i*5] = step_simulation2(q[:i*5], qdot[:i*5], p, density, dt, num_axis, grid_res, pix_size, g)
    particles.set_data(q[:i*5,0], q[:i*5,1])
    return particles

plt.xlim(0, box_size)
plt.ylim(0, box_size)
anim = animation.FuncAnimation(fig, animate, interval=10)
plt.show()
