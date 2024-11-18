# Structure-Preserving Numerical Methods for Fokker-Planck Equations

[![License: MIT](https://img.shields.io/badge/License-MIT-success.svg)](https://opensource.org/licenses/MIT)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10955446.svg)](https://doi.org/10.5281/zenodo.10955446)


This repository contains information and code to reproduce the results presented in the
article
```bibtex
@article{bartel2024structure,
  title={Structure-Preserving Numerical Methods for {F}okker-{P}lanck
         Equations},
  author={Bartel, Hanna and Lampert, Joshua and Ranocha, Hendrik},
  journal={Proceedings in Applied Mathematics and Mechanics},
  pages={e202400007},
  year={2024},
  month={11},
  doi={10.1002/pamm.202400007},
  eprint={2404.07641},
  eprinttype={arxiv},
  eprintclass={math.NA}
}
```

If you find these results useful, please cite the article mentioned above. If you
use the implementations provided here, please **also** cite this repository as
```bibtex
@misc{bartel2024structureRepro,
  title={Reproducibility repository for
         "{S}tructure-Preserving Numerical Methods
           for {F}okker-{P}lanck Equations"},
  author={Bartel, Hanna and Lampert, Joshua and Ranocha, Hendrik},
  year={2024},
  howpublished={\url{https://github.com/JoshuaLampert/2024\_fokker\_planck}},
  doi={10.5281/zenodo.10955446}
}
```

## Abstract

A common way to numerically solve Fokker-Planck equations is the Chang-Cooper method
in space combined with one of the Euler methods in time. However, the explicit Euler method is not
unconditionally positive, leading to severe restrictions on the time step to ensure positivity.
On the other hand, the implicit Euler method is robust but nonlinearly implicit. Instead, we
propose to combine the Chang-Cooper method with unconditionally positive Patankar-type time
integration methods, since they are unconditionally positive, robust for stiff problems,
only linearly implicit, and also higher-order accurate. We describe the combined approach,
analyse it, and present a relevant numerical example demonstrating advantages compared to schemes
proposed in the literature.


## Numerical experiments

To reproduce the numerical experiments presented in this article, you need
to install [Julia](https://julialang.org/). The numerical experiments presented
in this article were performed using Julia v1.10.1.

First, you need to download this repository, e.g., by cloning it with `git`
or by downloading an archive via the GitHub interface. Then, you need to start
Julia in the `code` directory of this repository and follow the instructions
described in the `README.md` file therein.


## Authors

- Hanna Bartel (University of Hamburg, Germany)
- Joshua Lampert (University of Hamburg, Germany)
- [Hendrik Ranocha](https://ranocha.de) (Johannes Gutenberg University Mainz, Germany)


## License

The code in this repository is published under the MIT license, see the
`LICENSE` file.


## Disclaimer

Everything is provided as is and without warranty. Use at your own risk!
