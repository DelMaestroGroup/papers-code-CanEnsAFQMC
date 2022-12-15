[![Paper](https://img.shields.io/badge/paper-arXiv%3AXXXX.YYYYY-B31B1B.svg)](https://arxiv.org/abs/XXXX.YYYYY)
[![DOI](https://zenodo.org/badge/214220909.svg)](https://zenodo.org/badge/latestdoi/214220909)

# A Stable, Recursive Auxiliary Field Quantum Monte Carlo Algorithm in the Canonical Ensemble: Applications to Thermometry and the Hubbard Model

Tong Shen, Hatem Barghathi, Jiangyong Yu, Adrian Del Maestro, Brenda M. Rubenstein

[arXiv:XXXX.YYYYY](https://arxiv.org/abs/XXXX.YYYYY)

### Abstract
Many experimentally-accessible, finite-sized interacting quantum systems are most appropriately described by the canonical ensemble of statistical mechanics. Conventional numerical simulation methods either approximate them as being coupled to a particle bath, or use projective algorithms which may suffer from non-optimal scaling with system size or large algorithmic prefactors.  In this paper, we introduce a highly stable, recursive Auxiliary Field Quantum Monte Carlo approach that can directly simulate systems in the canonical ensemble.  We apply the method to the fermionic Hubbard model in one and two spatial dimensions in a regime known to exhibit a significant "sign" problem and find improved performance over existing approaches including rapid convergence to ground state expectation values.  The effects of excitations above the ground state are quantified using an estimator-agnostic approach including studying the temperature dependence of the purity and overlap fidelity of the canonical and grand canonical density matrices.  As an important application, we show that thermometry approaches often exploited in ultra-cold atoms that employ an analysis of the velocity distribution in the grand canonical ensemble may be subject to errors leading to an under-estimation of extracted temperatures with respect to the Fermi temperature.

### Description
This repository includes links, code, scripts, and data to generate the figures in a paper.

### Requirements
The data in this project was generated via averaging Monte Carlo sample data.  Everything included in the [data](https://github.com/DelMaestroGroup/papers-code-CanEnsAFQMC/tree/main/data) directory was generated from the raw data via:

* [Dependency Name](https://dependencelink)

### Support
The creation of these materials was supported in part by the NSF CTMC CAREER Award 2046744 and the NSF under Grant No. DMR-2041995.

[<img width="100px" src="https://www.nsf.gov/images/logos/NSF_4-Color_bitmap_Logo.png">](http://www.nsf.gov/awardsearch/showAward?AWD_ID=1553991)

### Figures

#### Figure 04: Purity Comparison between Ensembles
<img src="https://raw.githubusercontent.com/DelMaestroGroup/papers-code-CanEnsAFQMC/main/figures/Purity_Lx6Ly6.svg" width="400px">

#### Figure 05: Fidelity Demonstrating Convergence to the Ground State
<img src="https://raw.githubusercontent.com/DelMaestroGroup/papers-code-CanEnsAFQMC/main/figures/Fidelity_Lx6Ly6.svg" width="400px">

#### Figure 07: Comparison of Charge Structure Factors
<img src="https://raw.githubusercontent.com/DelMaestroGroup/papers-code-CanEnsAFQMC/main/figures/Cq_Lx6Ly6_U2_N46.svg" width="400px">

