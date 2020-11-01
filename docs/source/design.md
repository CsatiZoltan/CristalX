# Program design

*CristalX* is not a black-box library such as BLAS, neither is a GUI-based application intended for end-users. It is rather an easy-to-use and extensible set of Python codes that provide the basic functionalities that scientists can extend based on their needs. The following ideas were kept in mind while writing and maintaining *CristalX*.

- Driven by actual needs

Only implement features that are currently used. Adding extra features requires more testing, possibly more dependencies and therefore code bloat, and increases the cognitive load of the user. Instead, the emphasis is on creating a stable minimum core library that can be easily extended according to users’ demands. Consequently, application code is separated from the core modules.

- Build on well-established packages

We rely on the scientific Python stack: *NumPy* for array manipulations, *SciPy* for interpolation and some other computations, *Matplotlib* for visualization and *scikit-image* for image processing. This ensures interoperability with other scientific codes and that our software is hopefully bug-free.

- Minimize the dependencies
  

Rapid prototyping is essential in scientific code development and Python is an excellent choice to satisfy this requirement. At the same time, relying on fast libraries ensures that the computations are reasonably fast. The libraries mentioned in the previous point are easy to install, often already pre-installed in certain Python distributions.

- High-quality documentation

Future contributors will benefit from the rich documentation. Python doctests are extensively used, serving both as test cases and as examples of usage. The docstrings conform to the *numpydoc* style guide.

- We strive for decoupling the modules

Although part of the *grains* package, if the modules are independent, they can be reused in other projects too just by copy-pasting the required functions.

- Do not overuse classes

In the prototyping phase, prefer using free functions to methods. As an idea evolves, you will naturally find data and algorithms that belong together, and can refactor free functions into member functions of a class.

- Do not use deep hierarchies

Initially, stay away from excessive nesting to avoid fragmenting the code base. If the project grows big, you can still refactor the code by introducing deeper hierarchies. Deep nesting causes unnecessary cognitive load and it also makes the code more verbose at the caller's site. Compare

   ```python
   import package1.package2.module
   ```

   with

   ```python
   import package.module
   ```

- Gradually refactor code

As more and more features are added to the project, we will often find that similar tasks emerge in different contexts. It is a good time to think about how they can be generalized and to reconsider your model. This way, you will come up with utility functions best put into *utils.py*.

- Start writing code only after careful thinking

It is no point in writing code before you completely understand your problem domain you want to model. It is more efficient to build abstractions in your head or on paper, then to split it into modular chunks, and only after that start coding.


- Write the documentation *before* the code

If you document the function parameters and the return values in advance, as well as construct a doctest, you are enforced to think about the problem deeply and to create a good interface. Moreover, it guarantees that the documentation is not missing (what you would anyway have to write at some point, so why not at the beginning?).

- Give doctest-compatible examples

You hit two birds with the same stone: provide an example for the user and get some confidence that your code works as intended (at least for the particular example). As mentioned in the previous point, write them *before* the actual code implementation. Doctests do not replace careful testing.

- Keep the documentation as part of your code

The problem with wiki pages is that they are version controlled in a different Git repository. It makes it longer to change a hosting service (e.g. moving from GitHub to GitLab), you need to maintain two repositories, cannot change the documentation and the code in the same commit, and you have to rely on the rendering capabilities of the hosting service (e.g. GitHub cannot render math). It is therefore better to keep the documentation as part of your project in a dedicated directory (`docs/` in our case), use a documentation generator (*Sphinx* in our case) and host it online (on *Read the Docs* in our case).
