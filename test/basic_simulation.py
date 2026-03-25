from src.system import SolarSystem
import matplotlib.pyplot as plt

def main():

    system = SolarSystem()
    system.read_parameters("parameters_solar.json")
    system.run_simulation()
    fig, ax = system.animation()

    plt.show()

if __name__ == "__main__":
    main()
