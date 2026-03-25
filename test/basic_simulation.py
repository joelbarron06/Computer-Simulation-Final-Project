from src.system import SolarSystem
import matplotlib.pyplot as plt
import os

def main():

    system = SolarSystem()
    system.read_parameters("parameters_solar.json")
    system.run_simulation()
    fig, ax = system.animation()

    BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    plt.savefig(os.path.join(BASE_DIR, "figures", "basic_animation"))
    plt.show()

if __name__ == "__main__":
    main()
