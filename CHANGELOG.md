# Changelog

All notable changes to this project will be documented in this file. The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/), and this project **does not** adhere to [Semantic Versioning](https://semver.org/spec/v2.0.0.html). We refer to [GitHub issues](https://github.com/CsatiZoltan/CristalX/issues) by their numbers. If there is no issue associated to the change, the commit hashes implementing the change are linked.



## [Unreleased]






## [1.1.0] - 2021-03-04

### Added

- Grains embedded into other grains are now identified in the splinegon creation algorithm. [4cef7ec](https://github.com/CsatiZoltan/CristalX/commit/4cef7ecbbe4c3fc5cddb9947b29b3d0b454eecb8), [6c4ad0a](https://github.com/CsatiZoltan/CristalX/commit/6c4ad0aff2c4dae49fb648656a94b14436925bed)
- Computation of the infinitesimal and Green-Lagrange strain tensors for DIC displacement field. [7eaa0b5](https://github.com/CsatiZoltan/CristalX/commit/7eaa0b523c0bc7e7a585f955aa8c31206cff751a), [9e5466e](https://github.com/CsatiZoltan/CristalX/commit/9e5466ebe1f4ab19defb9e91ff2ed8a961b7044b)
- Computation of the equivalent von Mises strain. [2b9b0be](https://github.com/CsatiZoltan/CristalX/commit/2b9b0be36cf38f5d1683821573646a65ca541160)
- Characterization of the strain localization (intergranular/intragranular) in the microstructure. [a74f954](https://github.com/CsatiZoltan/CristalX/commit/a74f95424e08d4915ff833e335a2863103da5ca4), [13576f7](https://github.com/CsatiZoltan/CristalX/commit/13576f7ed85f23b9680eecad030229fc4636e50d)
- The README file shows how to cite our paper. [#30](https://github.com/CsatiZoltan/CristalX/issues/30)
- Added a new segmented microstructure ([6b539ba](https://github.com/CsatiZoltan/CristalX/commit/6b539ba356bb2e27f6dca4fbc299ca90635d3ab2)) and the corresponding DIC measurements ([3c845e7](https://github.com/CsatiZoltan/CristalX/commit/3c845e7e180a75acb2ea27404d75c308d35ca786)).

### Deprecated

- The content of `dic.DIC.plot_strain` will be replaced with that of `dic.plot_strain` in version 1.3.0. The latter function will be removed. [b8633c2](https://github.com/CsatiZoltan/CristalX/commit/b8633c202cff0b8bc7b698cc771473b8d40055a9)




## [1.0.1] - 2020-11-19

### Added

- Document on describing the versioning scheme we will follow. [e74474c](https://github.com/CsatiZoltan/CristalX/commit/e74474cb3d1448a1784281f0aa934d3789197662)
- Added a changelog to the project. [b48821b](https://github.com/CsatiZoltan/CristalX/commit/b48821b5e2df5923ba19e5895c20a3ec248815e9), [f6ca9c1](https://github.com/CsatiZoltan/CristalX/commit/f6ca9c1da56a77af5696c1262d78f1184c6b1deb)

### Removed

- The documentation no longer shows the recent git commits. [1fd3a50](https://github.com/CsatiZoltan/CristalX/commit/1fd3a50930074b0ff79744bf6069af4a211ec0ca)




## [1.0.0] - 2020-11-16

This is the initial release of *CristalX*.



[unreleased]: https://github.com/CsatiZoltan/CristalX/compare/v1.1.0...HEAD
[1.1.0]: https://github.com/CsatiZoltan/CristalX/compare/v1.0.1...v1.1.0
[1.0.1]: https://github.com/CsatiZoltan/CristalX/compare/v1.0.0...v1.0.1
[1.0.0]: https://github.com/CsatiZoltan/CristalX/compare/981dbcd...v1.0.0