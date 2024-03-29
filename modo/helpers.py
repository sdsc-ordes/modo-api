from enum import Enum
from pathlib import Path
import re
import shutil
from typing import Any, Mapping, Optional
from urllib.parse import urlparse
import zarr

import modo_schema.datamodel as model

from .introspection import get_haspart_property, get_slot_range, load_schema


def class_from_name(name: str):
    class_names = list(load_schema().all_classes().keys())
    if name not in class_names:
        raise ValueError(f"Unknown class name: {name}")
    return getattr(model, name)


def dict_to_instance(element: Mapping[str, Any]) -> Any:
    elem_type = element.get("@type")
    target_class = class_from_name(elem_type)
    return target_class(
        **{k: v for k, v in element.items() if k not in "@type"}
    )


def copy_file_to_archive(
    data_file: Optional[str], base_path: Path, archive_path: Optional[Path]
):
    if data_file is not None:
        data_path = Path(data_file)
        shutil.copy(data_path, base_path / archive_path)


def set_haspart_relationship(
    child_class: str,
    child_path: str,
    parent_group: zarr.hierarchy.Group,
):
    """Add element to the hasPart attribute of a parent zarr group"""
    parent_type = getattr(
        model,
        parent_group.attrs.get("@type"),
    )

    has_prop = get_haspart_property(child_class)
    parent_slots = parent_type.__match_args__
    if has_prop not in parent_slots:
        raise ValueError(
            f"Cannot make {child_path} part of {parent_group.name}: {parent_type} does not have property {has_prop}"
        )
    # has_part is multivalued
    if has_prop not in parent_group.attrs:
        parent_group.attrs[has_prop] = []
    parent_group.attrs[has_prop] += [child_path]


def update_haspart_id(
    element: model.DataEntity
    | model.Sample
    | model.Assay
    | model.ReferenceGenome
    | model.MODO,
):
    """update the id of the has_part property of an element to use the full id including its type"""
    haspart_names = load_schema().slot_children("has_part")
    haspart_list = [
        haspart for haspart in haspart_names if haspart in vars(element).keys()
    ]
    if len(haspart_list) > 0:
        for has_part in haspart_list:
            haspart_type = get_slot_range(has_part)
            type_name = ElementType.from_model_name(haspart_type).value
            updated_ids = [
                f"{type_name}/{id}" for id in getattr(element, has_part)
            ]
            setattr(element, has_part, updated_ids)
    return element


class UserElementType(str, Enum):
    """Enumeration of element types exposed to the user."""

    SAMPLE = "sample"
    ASSAY = "assay"
    DATA_ENTITY = "data"
    REFERENCE_GENOME = "reference"

    def get_target_class(
        self,
    ) -> type:
        """Return the target class for the element type."""
        if self == UserElementType.SAMPLE:
            return model.Sample
        elif self == UserElementType.ASSAY:
            return model.Assay
        elif self == UserElementType.DATA_ENTITY:
            return model.DataEntity
        elif self == UserElementType.REFERENCE_GENOME:
            return model.ReferenceGenome
        else:
            raise ValueError(f"Unknown element type: {self}")

    @classmethod
    def from_object(cls, obj):
        """Return the element type from an object."""
        if isinstance(obj, model.Sample):
            return UserElementType.SAMPLE
        elif isinstance(obj, model.Assay):
            return UserElementType.ASSAY
        elif isinstance(obj, model.DataEntity):
            return UserElementType.DATA_ENTITY
        elif isinstance(obj, model.ReferenceGenome):
            return UserElementType.REFERENCE_GENOME
        else:
            raise ValueError(f"Unknown object type: {type(obj)}")


class ElementType(str, Enum):
    """Enumeration of all element types."""

    SAMPLE = "sample"
    ASSAY = "assay"
    DATA_ENTITY = "data"
    REFERENCE_GENOME = "reference"
    REFERENCE_SEQUENCE = "sequence"

    def get_target_class(
        self,
    ) -> type:
        """Return the target class for the element type."""
        if self == ElementType.SAMPLE:
            return model.Sample
        elif self == ElementType.ASSAY:
            return model.Assay
        elif self == ElementType.DATA_ENTITY:
            return model.DataEntity
        elif self == ElementType.REFERENCE_GENOME:
            return model.ReferenceGenome
        elif self == ElementType.REFERENCE_SEQUENCE:
            return model.ReferenceSequence
        else:
            raise ValueError(f"Unknown element type: {self}")

    @classmethod
    def from_object(cls, obj):
        """Return the element type from an object."""
        if isinstance(obj, model.Sample):
            return ElementType.SAMPLE
        elif isinstance(obj, model.Assay):
            return ElementType.ASSAY
        elif isinstance(obj, model.DataEntity):
            return ElementType.DATA_ENTITY
        elif isinstance(obj, model.ReferenceGenome):
            return ElementType.REFERENCE_GENOME
        elif isinstance(obj, model.ReferenceSequence):
            return ElementType.REFERENCE_SEQUENCE
        else:
            raise ValueError(f"Unknown object type: {type(obj)}")

    @classmethod
    def from_model_name(cls, obj):
        """Return the element type from an object."""
        if obj == "Sample":
            return ElementType.SAMPLE
        elif obj == "Assay":
            return ElementType.ASSAY
        elif obj == "DataEntity":
            return ElementType.DATA_ENTITY
        elif obj == "ReferenceGenome":
            return ElementType.REFERENCE_GENOME
        elif obj == "ReferenceSequence":
            return ElementType.REFERENCE_SEQUENCE
        else:
            raise ValueError(f"Unknown object type: {obj}")


def is_uri(text: str):
    """Checks if input is a valid URI."""
    try:
        result = urlparse(text)
        return all([result.scheme, result.netloc])
    except AttributeError:
        return False


def parse_region(region: str) -> tuple[str, int, int]:
    """Parses an input UCSC-format region string into
    (chrom, start, end).

    Examples
    --------
    >>> parse_region('chr1:10-320')
    ('chr1', 10, 320)
    >>> parse_region('chr-1ba:32-0100')
    ('chr-1ba', 32, 100)
    """

    if not re.match(r"[^:]+:[0-9]+-[0-9]+", region):
        raise ValueError(
            f"Invalid region format: {region}. Expected chr:start-end"
        )

    chrom, coords = region.split(":")
    start, end = coords.split("-")

    return (chrom, int(start), int(end))
