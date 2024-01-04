import json
from pathlib import Path
from typing import Any

from linkml_runtime.loaders import json_loader
import zarr

from .helpers import ElementType


# Initialize object's directory given the metadata graph
def init_zarr(root_directory: Path) -> zarr.Group:
    """Initialize object's directory given the metadata graph."""
    root_directory.mkdir(exist_ok=True)
    store = zarr.DirectoryStore(str(root_directory / "data.zarr"))
    data = zarr.group(store=store)

    return data


def add_metadata_group(parent_group: zarr.Group, metadata: dict) -> None:
    """Add input metadata dictionary to an existing zarr group."""
    # zarr groups cannot have slashes in their names
    group_name = metadata["id"].replace("/", "_")
    parent_group.create_group(group_name)
    # Fill attrs in the subject group for each predicate
    for key, value in metadata.items():
        if key == "id":
            continue
        parent_group[group_name].attrs[key] = value


def add_data(group: zarr.Group, data) -> None:
    """Add a numpy array to an existing zarr group."""
    group.create_dataset("data", data=data)


def attrs_to_instance(group: zarr.Group, id: str) -> Any:
    """Convert zarr group attributes to an instanced element.""" ""
    meta = dict(group.attrs)
    meta["id"] = id
    target_class = ElementType(meta["@type"].lower()).get_target_class()

    return json_loader.loads(json.dumps(meta), target_class=target_class)


def list_zarr_items(group: zarr.Group) -> list[zarr.Group | zarr.Array]:
    """Recursively list all zarr groups and arrays"""
    found = []

    def list_all(path: str, elem):
        found.append((path, elem))

    group.visititems(list_all)
    return found