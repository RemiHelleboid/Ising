import numpy as np
import matplotlib.pyplot as plt

# try:
#     plt.style.use(['science', 'high-vis'])
# except:
#     plt.style.use(['seaborn-paper'])


class Ising_2d:
    def __init__(self, grid_size) -> None:
        self.configuration = np.zeros((grid_size, grid_size))
        self.neighbors = np.array([[-1, 0], [1, 0], [0, -1], [0, 1]])
        self.grid_size = grid_size
        self.temperature = 1.0
        self.iteration = 0

    def randomize_configuration(self, probability):
        for i in range(self.grid_size):
            for j in range(self.grid_size):
                self.configuration[i, j] = 1 if (
                    np.random.rand() < probability) else -1

    def get_total_energy(self):
        energy = 0
        for i in range(self.grid_size):
            for j in range(self.grid_size):
                energy += -self.configuration[i, j] * (
                    self.configuration[(i + 1) % self.grid_size, j] +
                    self.configuration[(i - 1) % self.grid_size, j] +
                    self.configuration[i, (j + 1) % self.grid_size] +
                    self.configuration[i, (j - 1) % self.grid_size])
        return energy / 2.0

    def metropolis_iter(self):
        for i in range(self.grid_size):
            for j in range(self.grid_size):
                rand_index_x = np.random.randint(0, self.grid_size)
                rand_index_y = np.random.randint(0, self.grid_size)
                spin_rd = self.configuration[rand_index_x, rand_index_y]
                delta_e = 2.0 * spin_rd * (self.configuration[(rand_index_x + 1) % self.grid_size, rand_index_y] +
                                           self.configuration[(rand_index_x - 1) % self.grid_size, rand_index_y] +
                                           self.configuration[rand_index_x, (rand_index_y + 1) % self.grid_size] +
                                           self.configuration[rand_index_x, (rand_index_y - 1) % self.grid_size])
                if delta_e <= 0 or np.random.rand() < np.exp(-delta_e / self.temperature):
                    self.configuration[rand_index_x, rand_index_y] *= -1
        self.iteration += 1

    def get_magnetization(self):
        return np.sum(self.configuration)

    def plot_configuration(self, ax):
        ax.imshow(self.configuration, cmap='rainbow', interpolation='gaussian')


def main():
    grid_size = 150
    temperature = 1.0e-3
    max_iteration = 1000
    list_magnetization = []
    list_total_energy = []

    ising = Ising_2d(grid_size)
    ising.randomize_configuration(0.5)
    list_magnetization.append(ising.get_magnetization())
    list_total_energy.append(ising.get_total_energy())

    fig, axs = plt.subplots(1, 3, figsize=(15, 5))
    axs[0].set_title('Configuration')
    axs[1].set_title('Magnetization')
    axs[2].set_title('Total Energy')
    fig.tight_layout()
    for iter in range(max_iteration):
        ising.metropolis_iter()
        # list_magnetization.append(ising.get_magnetization())
        # list_total_energy.append(ising.get_total_energy())
        ising.plot_configuration(axs[0])
        # axs[1].plot(list_magnetization, lw=2.0, c='b')
        # axs[2].plot(list_total_energy, lw=2.0, c='r')
        # axs[1].set_title('Magnetization')
        # axs[2].set_title('Total Energy')
        plt.pause(0.00001)
        # axs[0].clear()
        # axs[1].clear()
        # axs[2].clear()
    plt.show()


if __name__ == '__main__':
    main()
