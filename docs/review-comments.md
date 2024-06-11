# Review Comments 💁

## Project Structure 🚀

The project is in good structure already. I like it!
You did a good job. File-naming is consistent, nice small kebab-case folder names. Already in good shape with BPA Repository Guidelines I think, (also there is much more I like to see, later =)). Overall continue what I see and beeing clean and precise etc. ❤️

- `has_sample` variables in `.zarr` files implies `bool` values. Maybe consider
  `sample(s)`. Check also for other use cases of this pattern.

- Maybe add a `.modos` extension to `ZARR` root folder (representing the file
  format). This gives the whole thing a much nicer semantic meaning. Yes its a file-format =).

- Maybe put all top-level files `.bcf` etc. in `ZARR` files into `root/assets` (e.g. `root=<repo-root>/data/ex`). To make the toplevel of the file `.modos` clean. =)

- I would introduce the following folder structure on top-level
  to have it even more clean.

  - move `modos` -> `src/modos` : The `modos` source code.
  - move `data` -> `examples/`
  - create `tools` : Tools and scripts needed in the future for this repository. Tooling around certain development aspects.
  - move `deploy` -> `tools/deployment`: Move the folder maybe there.
  - move `nix` -> `tools/nix`.

  **or** if you want a more mono-repository structure:

  - `components/modos`: the python component providing the python source code, its a folder which acts as a full python project.

  - `components/modos-backend` : the stuff which is needed to deploy the `backend`, e.g. `htsget` etc...

  This mono-repo structure makes it easier to add more stuff later, if we need another micro-service or a database etc.
  But maybe overkill for now!

- Provide a `.devcontainer` such that `code` users can simply run everything inside this container. I can make a short PR on this providing the minimal setup (ubuntu with nix installed, nice shell, etc). To make it more proper:
  _We need a `python` dev-container maintained centrally to pull this off properly. I can provide all this but needs discussion first in BPA and SDSC in general how you did things etc._

## Code Review 🐊

TODO: ...