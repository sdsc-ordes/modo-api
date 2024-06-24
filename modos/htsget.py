import io
from pathlib import Path
import re
from typing import Iterator, Optional
from urllib.parse import urlparse, ParseResult

import htsget
from pysam import AlignedSegment, VariantRecord

from .helpers import (
    parse_region,
    GenomicFileSuffix,
    file_to_pysam_object,
    bytesio_to_iterator,
    iter_to_file,
)

HTSGET_FORMATS = ("CRAM", "VCF", "BCF")


def prepare_htsget_url(url: str) -> ParseResult:

    htsget_url = urlparse(url)
    htsget_url = htsget_url._replace(
        path=re.sub(r"\.vcf\.\w+$", ".vcf", htsget_url.path)
    )  # remove aditional extensions, e.g., .vcf.gz -> .vcf
    htsget_url = htsget_url._replace(
        path=str(Path(htsget_url.path).with_suffix(""))
    )
    return htsget_url


class HtsGetStream:
    buffer: io.BytesIO
    file_format: str
    reference_filename: Optional[str]

    def __init__(
        self,
        buffer: io.BytesIO,
        file_format: str,
        reference_filename: Optional[str] = None,
    ) -> None:
        self.buffer = buffer
        self.file_format = file_format
        self.reference_filename = reference_filename

    def __iter__(self) -> Iterator[AlignedSegment | VariantRecord]:
        return bytesio_to_iterator(
            self.buffer,
            file_format=self.file_format,
            reference_filename=self.reference_filename,
        )

    def to_file(self, output_filename: str) -> None:
        # To save remote slice to a local file without converting data type
        out_fileformat = GenomicFileSuffix.from_path(
            Path(output_filename)
        ).name
        if out_fileformat == self.file_format:
            with open(output_filename, "wb") as output:
                for chunk in self.buffer:
                    output.write(chunk)
        else:
            raise ValueError("Input and output formats do not match.")


def stream_htsget(
    url: str,
    region: Optional[str] = None,
    reference_filename: Optional[str] = None,
) -> HtsGetStream:
    """Stream or write to a local file a slice of a remote CRAM or VCF/BCF file"""

    htsget_url = prepare_htsget_url(url)
    format = GenomicFileSuffix.from_path(Path(htsget_url.path)).name
    if format not in HTSGET_FORMATS:
        raise ValueError(
            f"Unsupported format for htsget: {format}. "
            f"Supported formats: {', '.join(HTSGET_FORMATS)}"
        )
    reference_name, start, end = parse_region(region)
    htsget_response_buffer = io.BytesIO()

    htsget.get(
        url=htsget_url,
        output=htsget_response_buffer,
        reference_name=reference_name,
        start=start,
        end=end,
        data_format=format,
    )

    htsget_response_buffer.seek(0)

    return HtsGetStream(
        buffer=htsget_response_buffer,
        file_format=format,
        reference_filename=reference_filename,
    )
