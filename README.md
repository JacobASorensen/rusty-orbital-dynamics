# rusty-orbital-dynamics
my goal is to make an accurate orbital dynamics simulator in rust, callable in python
![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)

### Building the library

```
maturin develop  # for local development
maturin build    # to build a wheel
```

## Milestones

### Implement First Lagrange Point Test Case
Write test case simulating James Web Space Telescope for 21 days to verify correct-ness of code
#### in Rust:
- Implement N-body ODE system solver
- Implement RKF4(5) integrator 
#### In Python: 
- Implement class to pull data from JPL Horizons (this can be scavenged from my Orbital Mechanics repo)
- Make Rust integrator callable in python