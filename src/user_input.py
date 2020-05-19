"""
Class to store user input
"""
class Inputs():
    def __init__(self):
        # Use init method to set default values
        # EoS
        self.gamma = 1.4

        # Artificial Viscosity
        self.cl = 0.5
        self.cq = 0.75

        # Mesh Extents - Used only for hard-coded spherical Sod
        self.x0 = 0.0
        self.x1 = 1.0
        self.y0 = 0.0
        self.y1 = 1.0

        # Timestep, etc.
        self.t0 = 0.0
        self.tend = 0.205
        self.dtinit = 0.0001
        self.dtmax = 0.0001
        self.growth = 1.02

        # Mesh resolution
        self.nregions = 1
        self.meshx = 50
        self.meshy = 50

        # Debug settings
        self.debug_step_count = 0   # Set to zero to disable debugging


        self.material_list = {
            1: 1,
            2: 2
        }

    def read_input(self, inputfile):
        import configparser
        import sys

        core_sections = ["Q", "SPHSOD MESH", "CONTROL", "DEBUG"]

        parser = configparser.ConfigParser()
        parser.read(inputfile)

        # Read core data
        for section in parser.sections():
            if section.upper() in core_sections:
                for option in parser.options(section):
                    if hasattr(self, option):
                        value = parser.get(section, option)
                        current = getattr(self, option)
                        value = type(current)(value)

                        setattr(self, option, value)
                    else:
                        print(
                            "Error: Unexpected input option {}"
                            " in section {}".format(option, section)
                        )
                        sys.exit()

            elif section.upper() == "MATERIAL":
                for option in parser.options(section):
                    if option == "material_numbers":
                        values = parser.get(section, option).split(",")

                        self.material_list = {}
                        index = 1
                        for matno in [int(v) for v in values]:
                            self.material_list[index] = matno
                            index += 1

                        print(self.material_list)
                    else:
                        print(
                            "Error: Unexpected input option {}"
                            " in section {}".format(option, section)
                        )
                        sys.exit()



if __name__ == "__main__":
    inputs = Inputs()
    inputs.read_input("input.in")
