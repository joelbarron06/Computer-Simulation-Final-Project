from src.system import SolarSystem
import matplotlib.pyplot as plt
import os

def main():

    system = SolarSystem()
    system.read_parameters("base_parameters.json")
    system.run_simulation()
    fig, ax, ani = system.animation()

    plt.show()

if __name__ == "__main__":
    main()
