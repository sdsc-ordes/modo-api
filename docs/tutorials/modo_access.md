# Explore a MODO and access it's data

There are multiple ways to use a `MODO` and interact with it's elements.

(explore)=
## Explore a MODO and it's elements

`MODO` metadata can be easily accessed via CLI or using the python api:

::::{tab-set}

:::{tab-item} python
:sync: python
```{code-block} python
from modo.api import MODO

# Create modo object (see Create and modify MODO)
modo = MODO(path = "data/ex")

# Show metadata as dictionary
modo.show_contents()
# {'ex': {'@type': 'MODO', 'creation_date': '2024-02-19T00:00:00', 'description': ..}
```
:::

:::{tab-item} cli
:sync: cli
```{code-block} console
modo show "data/ex"
```
:::

::::

The objects structure can be visualized by displaying the internal hierachy:

::::{tab-set}

:::{tab-item} python
:sync: python
```{code-block} python
modo.list_arrays()
#/
# ├── assay
# │   └── assay1
# ├── data
# │   └── demo1
# ├── reference
# │   └── reference1
# └── sample
#     └── sample1
```
:::

:::{tab-item} cli
:sync: cli
```{code-block} console
modo show --zarr "data/ex"
#/
# ├── assay
# │   └── assay1
# ├── data
# │   └── demo1
# ├── reference
# │   └── reference1
# └── sample
#     └── sample1
```
:::

::::

:::{note}
`MODO` internally uses [zarr's hierachy groups](https://zarr.readthedocs.io/en/stable/api/hierarchy.html). Each subdirectory represents a new hierachy group. Any rray-like data can directly be stored within these hierachy groups, while other file formats are stored separately.
:::

All files part of a `MODO` can be listed:

::::{tab-set}

:::{tab-item} python
:sync: python
```{code-block} python
[fi for fi in modo.list_files()]
# [PosixPath('data/ex/reference1.fa'), PosixPath('data/ex/demo1.cram')]
```
:::

:::{tab-item} cli
:sync: cli
```{code-block} console
modo show --files "data/ex"
# data/ex/reference1.fa
# data/ex/demo1.cram
```
:::

::::

(publish)=
## Publish MODO as linked data

A semantic artifact can be created from the digital object and published as linked data.
In this process JSON metadata are converted to RDF and all relative paths are converted to URI's.

::::{tab-set}

:::{tab-item} python
:sync: python
```{code-block} python
modo.knowledge_graph(uri_prefix="http://demo-data")
# <Graph identifier=N1aa0d4f1f3d9428b9c01e703096c5c96 (<class 'rdflib.graph.Graph'>)>
```
:::

:::{tab-item} cli
:sync: cli
```{code-block} console
modo publish --base-uri "http://demo-data" "data/ex"
# @prefix EDAM: <http://edamontology.org/> .
# @prefix NCIT: <http://purl.obolibrary.org/obo/NCIT_> .
# @prefix modos: <https://w3id.org/sdsc-ordes/modo-schema/> .
# @prefix schema1: <http://schema.org/> .
# @prefix xsd: <http://www.w3.org/2001/XMLSchema#> .
#
# <http://demo-data/assay/assay1> a modos:Assay ;
#     schema1:description "Dummy assay for tests." ;
#     schema1:name "Assay 1" ;
#     modos:has_data <http://demo-data/demo1> ;
#     modos:has_sample <http://demo-data/sample1> ;
#     modos:omics_type NCIT:C84343 .
#
#     ...
```
:::

::::

## Enrich metadata
Genomic files such as cram files store relevant metadata, e.g. their reference sequences, in their header. These information can be automatically extracted and included into a `MODO`.

```{code-block} python
# Enrich modo
modo.enrich_metadata()

# Check the added elements
modo.list_arrays()
```