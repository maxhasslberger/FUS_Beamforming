# Fast Optimization for Multi-Target Focused Ultrasound in Linear Regimes
A fast, constraint-based optimization framework for multi-target focused ultrasound using phased array transducers and k-wave-based acoustic propagation models.

This repository contains the code used to generate the results reported in the associated publication Rapid optimization of focused ultrasound for complex targeting with phased array transducers and precomputed propagation operators (DOI: 10.1088/1361-6560/ae3afe).

## Dependencies
- [k-wave](http://www.k-wave.org)
- [kArray tools](http://www.k-wave.org/downloads/kWaveArray_alpha_0.3.zip)
- Matlab Optimization Toolbox
- Matlab Deep Learning Toolbox
- Matlab Deep Learning Toolbox Converter for TensorFlow Models
- ([Matlab2tikz](https://www.mathworks.com/matlabcentral/fileexchange/22022-matlab2tikz-matlab2tikz))

## Quick Start
After installing all dependencies, run `simulationApp.mlapp` to open the GUI. Follow the tabs at the top to define medium properties, phased array transducer geometry, position, and orientation, target area(s), and optimizer parameters. An (outdated, but still helpful) [Video](https://drive.google.com/file/d/1dEiv4D78owbD0Kw8CSpThkDKQ44hQXB1/view?usp=sharing) provides a demo, sonicating the right amygdala, while applying off-target constraints. 

Please refer to the `README` files in the respective subfolders for details on the required input file formats.

## Final Notes
This codebase is intended for research and prototyping purposes and has not been optimized for production deployment.

Please don't hesitate to [contact me](mailto:hasslberger.max@gmail.com) with questions, recommendations, or feedback on this work!
