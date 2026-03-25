from src.system import SolarSystem
import matplotlib.pyplot as plt

def main():
    """
    Runs simulation and animation using paraeterms of basic configuration
    """
    system = SolarSystem()
    system.read_parameters("base_parameters.json")
    system.run_simulation()
    fig, ax, ani = system.animation()

    plt.show()

if __name__ == "__main__":
    main()
    
