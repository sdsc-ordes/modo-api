
:::{image} img/modo1.png
   :width: 200
   :alt: modo logo
:::


# Welcome to modo-api's documentation!

modo-api (__Multi-Omics Digital Object__) is a python library and command line tool to process, store and serve multi-omics data together with their metadata.
It allows remote storage and access using [zarr](https://github.com/zarr-developers/zarr-python)'s S3 compatibility and integrates [htsget](https://github.com/ga4gh/htsget) to stream and access genomics data.

::::{grid} 3

:::{grid-item-card} {octicon}`rocket;2em`  Quick start
A user guide to get started and try modo-api

:::{button-ref} intro/quickstart
:ref-type: myst
:expand:
:color: dark
:click-parent:

Try it out!
:::

:::


:::{grid-item-card} {octicon}`project;2em` Tutorials
A list of tutorials and how-to user guides.

:::{button-ref} tutorials/tutorials
:ref-type: myst
:expand:
:color: dark
:click-parent:

Guided tutorials!
:::

:::

:::{grid-item-card} {octicon}`report;2em` API reference
A reference guide to all python api classes and functions.

:::{button-ref} autoapi/index
:ref-type: myst
:expand:
:color: dark
:click-parent:

To the python API!
:::

:::

:::{grid-item-card} {octicon}`code-square;2em` CLI reference
A reference guide to all cli functions, their parameters and options.

:::{button-ref} cli
:ref-type: myst
:expand:
:color: dark
:click-parent:

Check the CLI!
:::

:::

:::{grid-item-card} {octicon}`book;2em` Background
Read more about __multiomics__ data and MODO's features and structure!

:::{button-ref} intro/background
:ref-type: myst
:expand:
:color: dark
:click-parent:

Learn more about MODO!
:::

:::

:::{grid-item-card} {octicon}`beaker;2em` Modo-schema
Check our data model for advanced and custom usage of modo-api.

:::{button-link} https://sdsc-ordes.github.io/modo-schema/
:ref-type: myst
:expand:
:color: dark
:click-parent:

To the data model!
:::

:::

::::




:::{toctree}
:maxdepth: 1
:hidden:

intro/quickstart.md
intro/background.md
tutorials/tutorials.md
autoapi/index
cli.rst
changelog-link
:::