---
title: "GasDispersion.jl"
excerpt: "A julia package for dispersion modeling of chemical releases"
header:
  overlay_image: /images/gaussian_dispersion_example_files/veeterzy-unsplash-header.jpg
  teaser: /images/projects/gasdispersion.jpg
  caption: "Photo by [veeterzy](https://unsplash.com/@veeterzy?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText) on [Unsplash](https://unsplash.com/)"
gallery:
  - url: https://nbviewer.org/github/aefarrell/GasDispersion.jl/blob/main/examples/01_gaussian_plume_models.ipynb
    image_path: /images/projects/gasdispersion_gaussian.svg
    alt: "gaussian plume dispersion"
  - url: https://nbviewer.org/github/aefarrell/GasDispersion.jl/blob/main/examples/02_gaussian_puff_models.ipynb
    image_path: /images/projects/gasdispersion_puff.gif
    alt: "gaussian puff dispersion"
  - url: https://nbviewer.org/github/aefarrell/GasDispersion.jl/blob/main/examples/03_britter_mcquaid_plume_models.ipynb
    image_path: /images/projects/gasdispersion_brittermcquaid.jpg
    alt: "Britter-McQuaid dense gas dispersion"
  - url: https://nbviewer.org/github/aefarrell/GasDispersion.jl/blob/main/examples/05_turbulent_jet_models.ipynb
    image_path: /images/projects/gasdispersion_simplejet.svg
    alt: "Simple turbulent jet dispersion"
comments: false
---

[![LICENSE](https://img.shields.io/badge/license-MIT-lightgrey.svg)](https://github.com/aefarrell/GasDispersion.jl/blob/main/LICENSE) [![Documentation](https://img.shields.io/badge/docs-dev-blue)](https://aefarrell.github.io/GasDispersion.jl/dev/) [![Build Status](https://github.com/aefarrell/GasDispersion.jl/workflows/CI/badge.svg)](https://github.com/aefarrell/GasDispersion.jl/actions) [![Coverage](https://codecov.io/gh/aefarrell/GasDispersion.jl/branch/main/graph/badge.svg?token=PB3LOR80K2)](https://codecov.io/gh/aefarrell/GasDispersion.jl)

[GasDispersion.jl](https://github.com/aefarrell/GasDispersion.jl) aims to bring together several models for dispersion modeling of chemical releases with a consistent interface. Currently it implements gaussian dispersion modeling for simple plumes and puffs.


{% include gallery caption="Examples of GasDispersion.jl in action." %}


---
