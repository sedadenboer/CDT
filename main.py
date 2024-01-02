# main.py
#
# Author: Seda den Boer
# Date: 02-01-2024
#
# Description:
# Main file for the project. This file will be used to run the project.

from code_cdt.classes.universe import Universe

def main():
    universe = Universe(10, 20)
    print(universe)

if __name__ == "__main__":
    main()
