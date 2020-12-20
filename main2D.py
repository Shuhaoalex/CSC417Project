from utils import *
import matplotlib.animation as animation
import matplotlib.pyplot as plt
import random

num_particles = 5000
box_size = 1
grid_res = 100
pix_size = box_size / grid_res
num_axis = 2

q = np.random.rand((num_particles, 2))
q_dot = np.zeros((num_particles, 2))
