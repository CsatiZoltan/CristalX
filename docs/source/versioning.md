# Versioning



The versioning of *CristalX* starts with 1.0.0 and it somewhat follows the rules of [Semantic Versioning 2.0.0](https://semver.org/), with the syntax MAJOR.MINOR.PATCH, where MAJOR introduces significant changes, MINOR comes with smaller changes, and PATCH provides a fixture or a tiny improvement either in the code or in the documentation.



## Why not SemVer?

Semantic Versioning ([SemVer](https://semver.org/)) is a widely used versioning scheme, applicable for **public APIs**. Its purpose is to be rigorous on how to indicate when a bug fix, new features or incompatible changes in the public API are introduced. There is an excellent discussion about it [here](https://gist.github.com/jashkenas/cbd2b088e20279ae2c8e) and a detailed guide [here](https://www.jering.tech/articles/semantic-versioning-in-practice). Many criticize it for not being indicative about the rate of important changes. E.g. 1.8.5 --> 1.9.0 may include dozens of relevant improvements, while 1.9.5 --> 2.0.0 may merely be a simple clean-up that changes the public API. However, [SemVer was never meant to be used for software version numbers](https://gist.github.com/jashkenas/cbd2b088e20279ae2c8e#gistcomment-3448638).

[*CristalX* is not a library](design), but a collection of tools that operate at a high level. Most of its functions are exposed to the user, except the ones marked with a leading underscore or two leading underscores. However, the public API would be a subset of these functions: e.g. the functions of the *utils* module are exposed but they are mostly intended to be used by other modules, not by the user. As *CristalX* is not a library but rather a tool for rapid prototyping, semantic versioning would not make much sense.



## Why not CalVer?

Calendar Versioning ([CalVer](https://calver.org/)) is an alternative to SemVer. *CristalX* may not come with regular changes in the future or it will get updates at irregular intervals. We do not want to give the impression that relevant changes are introduced linearly in time. Moreover, the changelog and the time stamp in the git commits clearly show when new releases are published.



## Our versioning

As there is no unconditionally best versioning system, we came up with our own, which seems to fit well for research code like *CristalX*. The starting point is SemVer with some differences. Our versioning is intended for humans. The MAJOR version in increased only for substantially new features. This is, of course, subjective but we want to avoid large MAJOR version numbers as e.g. in Firefox. The guideline to follow is that the novelty of a feature makes the MAJOR version increase, not the number of additions. Here are some examples. Assume that we are at 1.2.3. We fixed a set of related bugs in the code: 1.2.3 --> 1.2.4. Then we implemented several functionalities to speed up the code: 1.2.4 --> 1.3.0. A new module was created and several others were modified that allowed us to carry out groundbreaking research: 1.3.0 --> 2.0.0. An existing algorithm was improved to handle the corner cases: 2.0.0 --> 2.1.0.

Similarly to SemVer, when one of the digits increases, the ones right to it are set to zero (e.g. 1.0.1 --> 1.0.2,  1.0.3 --> 1.1.0,  1.8.6 --> 2.0.0).

We try to keep backward compatibility. Insignificant changes are postponed, and are included as part of a MINOR release. If you often feel the need to introduce changes to the function signatures, rather **add** new functions and give **deprecation notices** than remove or modify existing ones. This does not lead to code bloat in the long run because deprecated syntax is removed from time to time.



### Connection with the Git workflow

The version numbers are reflected in the tag names. Your normal Git workflow stays the same: commit modifications and push them to the remote repository. When you want to mark a commit yet to be included in a certain version, type

```bash
git tag -m "Concise message" v<version_number>
```

where `<version_number>` has the form MAJOR.MINOR.PATCH. Note that it is preceded by the `v` letter, conventionally used for tags. [As an example](https://github.com/CsatiZoltan/CristalX/tags):

```bash
git tag -m "Initial release." v1.0.0
```

Keep the tag message short: the detailed changes since the previous version are collected in the changelog.

The tag is pushed by

```bash
git push --tags
```

Whenever you publish a tag, and hence update the changelog, also [create a release](https://github.com/CsatiZoltan/CristalX/releases) for that tag on GitHub. Copy the changes the new version brings from the changelog to the description of the release. In this regard, we follow the way of [JabRef](https://github.com/JabRef/jabref/releases). The zipped size of *CristalX* is quite small, so [size constraint will not be a problem on GitHub](https://docs.github.com/en/free-pro-team@latest/github/managing-large-files/distributing-large-binaries).

