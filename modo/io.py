from pathlib import Path
import re
from typing import Any, List

from linkml_runtime.loaders import (
    json_loader,
    yaml_loader,
    csv_loader,
    rdf_loader,
)
import modo_schema.datamodel as model
from .api import MODO
from .cram import slice_cram
from .helpers import dict_to_instance, update_haspart_id

ext2loader = {
    "json": json_loader,
    r"ya?ml": yaml_loader,
    "csv": csv_loader,
    r"rdf|ttl|nt(riples)?": rdf_loader,
}


def get_loader(path: Path):
    """Get a loader based on the file extension using regex."""
    ext = path.suffix[1:]
    for pattern, loader in ext2loader.items():
        if re.match(pattern, ext):
            return loader
    return None


def parse_instance(path: Path, target_class):
    """Load a model of target_class from a file."""
    loader = get_loader(path)
    if not loader:
        raise ValueError(f"Unsupported file format: {path}")
    return loader.load(str(path), target_class)


def parse_multiple_instances(path: Path) -> List:
    """Load one or more model from file. Model types must be specified as @type"""
    loader = get_loader(path)
    if not loader:
        raise ValueError(f"Unsupported file format: {path}")
    elems = loader.load_as_dict(str(path))
    if not isinstance(elems, list):
        elems = [elems]
    instances = []
    for elem in elems:
        instances.append(dict_to_instance(elem))
    return instances


def build_modo_from_file(path: Path, object_directory: Path) -> MODO:
    """build a modo from a yaml or json file"""
    instances = parse_multiple_instances(Path(path))
    # check for unique ids and fail early
    ids = [inst.id for inst in instances]
    if len(ids) > len(set(ids)):
        dup = {x for x in ids if ids.count(x) > 1}
        raise ValueError(
            f"Please specify a unique ID. Element(s) with ID(s) {dup} already exist."
        )
    # use full id for has_part attributes
    instances = [update_haspart_id(inst) for inst in instances]

    modo_inst = [
        instance for instance in instances if isinstance(instance, model.MODO)
    ]
    if len(modo_inst) != 1:
        raise ValueError(
            f"There must be exactly 1 MODO in the input file. Found {len(modo_inst)}"
        )
    modo_dict = modo_inst[0]._as_dict
    modo = MODO(path=object_directory, **modo_dict)
    for instance in instances:
        if not isinstance(instance, model.MODO):
            # copy data-path into modo
            if (
                isinstance(instance, model.DataEntity)
                and not modo.path in Path(instance.data_path).parents
            ):
                data_file = instance.data_path
                instance.data_path = Path(data_file).name
                modo.add_element(instance, data_file=data_file)
            else:
                modo.add_element(instance)
    return modo


def slice(data: model.DataEntity, region: str) -> Any:
    """Returns a slice of the input data for the requested region.

    Parameters
    ----------
    data
        A data object in any supported format.
    region
        The region string in UCSC format (i.e. chr:start-end).
    """

    match str(data.data_format):
        case "CRAM":
            return slice_cram(data.data_path, region)
        case _:
            raise ValueError(f"Unsupported data format: {data.data_format}")
