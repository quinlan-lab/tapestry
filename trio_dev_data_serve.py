"""
HTTP server with range-request support for serving trio_dev_data to IGV.

Usage:
    .venv/bin/python serve_trio_dev_data.py [--port PORT]

Then point IGV at e.g. http://localhost:8000/input/NA12878.GRCh38.haplotagged.bam
"""

import argparse
import mimetypes
import os
from http.server import HTTPServer, SimpleHTTPRequestHandler

DATA_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "trio_dev_data")

# Register MIME types for genomic file formats so IGV gets proper Content-Type
# headers, including for index files (.bai, .tbi, .csi, .fai).
GENOMIC_MIME_TYPES = {
    ".bam": "application/octet-stream",
    ".bai": "application/octet-stream",
    ".vcf": "text/plain",
    ".gz": "application/gzip",
    ".tbi": "application/octet-stream",
    ".csi": "application/octet-stream",
    ".fa": "text/plain",
    ".fasta": "text/plain",
    ".fai": "text/plain",
    ".bed": "text/plain",
    ".bw": "application/octet-stream",
    ".bigwig": "application/octet-stream",
    ".ped": "text/plain",
}
for ext, mime in GENOMIC_MIME_TYPES.items():
    mimetypes.add_type(mime, ext)


class RangeRequestHandler(SimpleHTTPRequestHandler):
    """HTTP handler that supports Range requests (required by IGV for BAM/VCF/etc)."""

    def __init__(self, *args, **kwargs):
        super().__init__(*args, directory=DATA_DIR, **kwargs)

    def end_headers(self):
        self.send_header("Access-Control-Allow-Origin", "*")
        self.send_header("Access-Control-Allow-Headers", "Range")
        self.send_header("Access-Control-Expose-Headers", "Content-Length, Content-Range, Accept-Ranges")
        super().end_headers()

    def do_OPTIONS(self):
        self.send_response(200)
        self.end_headers()

    def do_GET(self):
        range_header = self.headers.get("Range")
        if range_header is None:
            super().do_GET()
            return

        path = self.translate_path(self.path)
        if not os.path.isfile(path):
            self.send_error(404)
            return

        file_size = os.path.getsize(path)

        # Parse "bytes=start-end"
        try:
            range_spec = range_header.replace("bytes=", "")
            parts = range_spec.split("-")
            start = int(parts[0]) if parts[0] else 0
            end = int(parts[1]) if parts[1] else file_size - 1
        except (ValueError, IndexError):
            self.send_error(416, "Invalid range")
            return

        if start > end or start >= file_size:
            self.send_error(416, "Range not satisfiable")
            return

        end = min(end, file_size - 1)
        length = end - start + 1

        self.send_response(206)
        content_type = self.guess_type(path)
        self.send_header("Content-Type", content_type)
        self.send_header("Content-Range", f"bytes {start}-{end}/{file_size}")
        self.send_header("Content-Length", str(length))
        self.send_header("Accept-Ranges", "bytes")
        self.end_headers()

        with open(path, "rb") as f:
            f.seek(start)
            self.wfile.write(f.read(length))

    def do_HEAD(self):
        path = self.translate_path(self.path)
        if not os.path.isfile(path):
            self.send_error(404)
            return

        file_size = os.path.getsize(path)
        self.send_response(200)
        self.send_header("Content-Type", self.guess_type(path))
        self.send_header("Content-Length", str(file_size))
        self.send_header("Accept-Ranges", "bytes")
        self.end_headers()


def main():
    parser = argparse.ArgumentParser(description="Serve trio_dev_data with range request support for IGV")
    parser.add_argument("--port", type=int, default=8000)
    args = parser.parse_args()

    server = HTTPServer(("localhost", args.port), RangeRequestHandler)
    print(f"Serving {DATA_DIR} on http://localhost:{args.port}")
    print("Press Ctrl+C to stop")
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        print("\nStopped.")


if __name__ == "__main__":
    main()
