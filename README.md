# Causal Dynamical Triangulations in 1+1D and 2+1D

This repository contains code and documentation for simulating and studying causal dynamical triangulations (CDT) in both 1+1D and 2+1D spacetime dimensions.

## Introduction

Causal Dynamical Triangulations is an approach to quantum gravity, where spacetime is approximated by a simplicial complex. The dynamics of spacetime is described through the evolution of this complex, which is constructed by gluing together simple building blocks known as simplices. In (1+1)D and (2+1)D, these simplices are triangles and tetrahedra respectively.

## Repository Structure

```
CDT/
├── src/
│   │
│   ├── 1+1/
│   │   │
│   │   ├── classes/
│   │   │   ├── __init__.py
│   │   │   ├── bag.py
│   │   │   ├── pool.py
│   │   │   ├── simulation.py
│   │   │   ├── triangle.py
│   │   │   ├── universe.py
│   │   │   └── vertex.py
│   │   │
│   │   └── main.py
│   │
│   ├── 2+1/
│   │   │
│   │   ├── classes/
│   │   │   │
│   │   │   ├── helper_functions/
│   │   │   │   ├── __init__.py
│   │   │   │   ├── generate_initial.py
│   │   │   │   └── helpers.py
│   │   │   │
│   │   │   ├── initial_universes/
│   │   │   │   └── initial_t3.CDT
│   │   │   │
│   │   │   ├── __init__.py
│   │   │   ├── bag.py
│   │   │   ├── check_mt.py
│   │   │   ├── observable.py
│   │   │   ├── pool.py
│   │   │   ├── simulation.py
│   │   │   ├── simulation_mt.py
│   │   │   ├── tetra.py
│   │   │   ├── universe.py
│   │   │   ├── universe_mt.py
│   │   │   └── vertex.py
│   │   │
│   │   ├── __init__.py
│   │   └── main.py
│   │
│   └── __init__.py
│
├── .gitignore
├── README.md
├── requirements.txt
└── __init__.py
```


## Usage

To use the code in this repository, follow the steps below:

1. **Install the dependencies**:
    Ensure you have Python installed. Then install the required Python packages using pip:
    ```bash
    pip install -r requirements.txt
    ```

3. **Running simulations**:
    - For (1+1)D simulations, navigate to the `src/1+1/` directory and run the main script:
        ```bash
        cd src/1+1/
        python main.py
        ```
    - For (2+1)D simulations, navigate to the `src/2+1/` directory and run the main script:
        ```bash
        cd src/2+1/
        python main.py
        ```

4. **Configuration**:
    - You can configure the simulation parameters by editing the configuration files or passing parameters directly to the simulation scripts.

## Acknowledgements

Reference paper and codebase:
- [Simulating CDT quantum gravity (Brunekreef, Loll, and Görlich, 2024)](https://www.sciencedirect.com/science/article/pii/S0010465524000936)

## References

For more information about Causal Dynamical Triangulations, refer to the following resources:
- [Quantum Gravity from Causal Dynamical Triangulations: A Review (Loll, 2019)](https://arxiv.org/abs/1905.08669)
- [Causal Dynamical Triangulations Wikipedia](https://en.wikipedia.org/wiki/Causal_dynamical_triangulation)
- [CDT Research Papers](https://search.arxiv.org/?query=causal+dynamical+triangulations&in=)

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for more details.
