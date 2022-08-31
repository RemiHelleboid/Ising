import numpy as np
import matplotlib.pyplot as plt
from argparse import ArgumentParser
import matplotlib.animation as animation


def get_configuration(filename):
    X, Y, Spins = np.loadtxt(filename, unpack=True, delimiter=",", skiprows=1)
    nx = int(np.sqrt(len(X)))
    ny = int(np.sqrt(len(Y)))

    Spins = Spins.reshape((nx, ny))

    return Spins


def get_files_spins_from_dir(dirname):
    import os
    import glob
    files = glob.glob(os.path.join(dirname, "*_0*.csv"))
    files.sort()
    return files


def main(dirname, nb_configurations):
    files = get_files_spins_from_dir(dirname)
    fig, axs = plt.subplots(1, figsize=(10, 10))
    axs.set_axis_off()
    axs.set_xmargin(0.0)
    axs.set_ymargin(0.0)

    InitialConfig = get_configuration(files[0])
    print("InitialConfig.shape", InitialConfig.shape)
    im = axs.imshow(InitialConfig, interpolation='bicubic', cmap='plasma')
    axs.set_title(f"Iteration 0", fontsize=20)
    fig.tight_layout()

    def init():
        im.set_data(InitialConfig)
        return [im]

    def animate(i):
        data = get_configuration(files[i])
        im.set_array(data)
        axs.set_title(f"Iteration {i}", fontsize=20)
        return [im]

    ani = animation.FuncAnimation(fig, animate, frames=nb_configurations,
                                  init_func=init, blit=False, interval=1.0/30, repeat=False)
    plt.show()
    print("Saving animation to Ising_2d.mp4")
    ani.save("Ising_2d.mp4", writer="ffmpeg", fps=30)
    print("Done !")


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument(
        "-d", "--directory", help="directory to parse", default=".", dest="directory")
    parser.add_argument(
        "-N", "--NbConfig", help="Number of config to parse", dest="N")

    args = parser.parse_args()
    DIRECTORY = args.directory
    FILES = get_files_spins_from_dir(DIRECTORY)
    NB_CONFIG = int(args.N) if args.N else len(FILES)

    main(DIRECTORY, NB_CONFIG)
