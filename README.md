## Phase Separation Coupling between Membrane and Cytoplasm

This project is aiming to study some fundamental phase behaviors of cell membrane and cytoplasm under the coupling effect between each other. The whole system should include three parts-outer solvent, cell membrane and the cytoplasm, and there exists both thermodynamic and hydrodynamic coupling effects between out solvent and cell membrane, and between cell membrane and cytoplasm. At the early stage, we ignore the out solvent’s effect for simplicity, and only consider the thermodynamic coupling effect between cell membrane and the cytoplasm. 

Given the length scale we study, we use phase field method to simulation kinetic procedure happening both in cell membrane and cytoplasm.

### Cell Membrane
The true cell membrane has a lipid bilayer structure, with proteins embedded as well. Given the problem we study, we treat it as a simple 2-dimensional spherical sheet with no internal structure, and use the simple 2 components model. In numerical scheme, we use Yin-Yang grids to divide the whole sheet into two parts. The two parts communicate with each other at the internal boundaries.

### Cytoplasm/Inner Solvent
The true cytoplasm has a very complex structure and various different components. Given the problem we study, we only consider it as a 3-dimension sphere with 2 components. In numerical scheme, the whole cytoplasm can be divided into three parts-Yin, Yang and core. These three parts talk to each other. We also use the no-flux boundary condition at the surface of the cytoplasm.