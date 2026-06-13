"""
aiohttp server that serves genomic data files for IGV to consume.

Serves files under --data-dir over http://localhost:PORT/<relative-path> with
the byte-range support IGV needs. Robust by design (this is why it replaced the
stdlib http.server version, archived in .trash/serve_data_for_igv_stdlib.py):

  * A whole-file GET of a DATA file (BAM / bgzipped VCF / BED) returns 400.
    IGV must instead fetch the file's index (.bai/.tbi) and issue byte-Range
    requests for just the visible window. This makes it impossible to
    accidentally stream a genome-wide file in full -- the failure mode that
    saturated the SSH tunnel and corrupted concurrent track loads with the
    stdlib server.
  * Only index files and small text annotations (.bai/.tbi/.csi/.crai/.fai/
    .idx/.gtf/.gff/.gff3/.bed/.bedgraph) are returned whole.
  * Range responses stream in 1 MB chunks (bounded memory, clean EOF).

Pair with an IGV session whose Resource paths are
http://localhost:PORT/<relative-path> under DIR, and tunnel from your laptop
with `ssh -L PORT:localhost:PORT <host>`.

Usage:
    .venv/bin/python serve_data_for_igv.py [--port PORT] [--data-dir DIR]
"""

import argparse
import logging
import os

from aiohttp import web

CHUNK = 1024 * 1024  # 1 MB

# Files IGV legitimately fetches whole: indexes and small text annotations.
# Everything else (BAM, *.vcf.gz, *.bed.gz) must be range-requested via its index.
WHOLE_FILE_SUFFIXES = (
    ".bai", ".tbi", ".csi", ".crai", ".fai", ".idx",
    ".gtf", ".gff", ".gff3", ".bed", ".bedgraph",
)

DEFAULT_DATA_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "trio_dev_data")


def _resolve(data_dir, filename):
    """Resolve a request path under data_dir, rejecting path traversal."""
    path = os.path.normpath(os.path.join(data_dir, filename))
    if path != data_dir and not path.startswith(data_dir + os.sep):
        return None
    return path


async def handle(request):
    data_dir = request.app["data_dir"]
    filename = request.match_info.get("filename", "")
    path = _resolve(data_dir, filename)
    if path is None or not os.path.isfile(path):
        return web.Response(status=404, text=f"{filename} not found")

    file_size = os.path.getsize(path)
    range_header = request.headers.get("Range")

    if range_header:
        try:
            byte_range = range_header.strip().split("=")[1]
            start_s, end_s = byte_range.split("-")
            start = int(start_s) if start_s else 0
            end = int(end_s) if end_s else file_size - 1
        except (IndexError, ValueError):
            return web.Response(status=416, text="Invalid Range")
        end = min(end, file_size - 1)
        if start >= file_size or start > end:
            return web.Response(status=416, text="Requested Range Not Satisfiable")

        response = web.StreamResponse(status=206, headers={
            "Content-Type": "application/octet-stream",
            "Content-Range": f"bytes {start}-{end}/{file_size}",
            "Content-Length": str(end - start + 1),
            "Accept-Ranges": "bytes",
        })
        await response.prepare(request)
        with open(path, "rb") as f:
            f.seek(start)
            remaining = end - start + 1
            while remaining > 0:
                chunk = f.read(min(CHUNK, remaining))
                if not chunk:
                    break
                await response.write(chunk)
                remaining -= len(chunk)
        await response.write_eof()
        return response

    if filename.endswith(WHOLE_FILE_SUFFIXES):
        with open(path, "rb") as f:
            body = f.read()
        return web.Response(status=200, body=body, headers={
            "Content-Type": "application/octet-stream",
            "Content-Length": str(file_size),
            "Accept-Ranges": "bytes",
        })

    return web.Response(status=400, text=(
        "For range requests, please use the Range header.\n"
        "Only index/annotation files are served whole.\n"))


def main():
    parser = argparse.ArgumentParser(description="Serve genomic data files for IGV to consume")
    parser.add_argument("--port", type=int, default=8000)
    parser.add_argument("--data-dir", default=DEFAULT_DATA_DIR,
                        help="Directory to serve as the HTTP root (default: trio_dev_data)")
    args = parser.parse_args()

    # Configure logging so aiohttp's access logger (INFO) actually emits a line
    # per request; without a handler it is silent. Status codes show too, so a
    # refused whole-file data GET (400) or a range read (206) is visible.
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s",
                        datefmt="%Y-%m-%d %H:%M:%S")

    data_dir = os.path.abspath(args.data_dir)
    app = web.Application()
    app["data_dir"] = data_dir
    app.router.add_get("/{filename:.*}", handle)
    print(f"Serving {data_dir} on http://localhost:{args.port}")
    print("Press Ctrl+C to stop")
    web.run_app(app, host="localhost", port=args.port, print=None,
                access_log_format='%a "%r" %s %b')


if __name__ == "__main__":
    main()
