#!/usr/bin/env python3
"""HTTP subserver. POST /<name> -> dispatches to module <name>.run(data) -> dict.
Run with the venv Python that has the heavy deps (deepBE etc.) installed.

Usage: startSubServer.py PORT [ONLY]
  PORT  TCP port to listen on
  ONLY  optional module name; if set, only that module is served (other paths 404).
        Used to isolate one model per process."""

import sys, os, json, importlib, traceback
from http.server import BaseHTTPRequestHandler, HTTPServer
from subserversConf import SUBSERVERS

PORT = int(sys.argv[1]) if len(sys.argv) > 1 else 8001
ONLY = sys.argv[2] if len(sys.argv) > 2 else None
HOST = "127.0.0.1"
MAX_BODY = 1024 * 1024

class Handler(BaseHTTPRequestHandler):

    timeout = 30

    def log_message(self, fmt, *args):
        sys.stderr.write((fmt % args) + "\n")

    def do_POST(self):
        name = self.path.strip("/").removesuffix(".py")

        if name not in SUBSERVERS or (ONLY and name != ONLY):
            return self._reply(404, {"error": "unknown subserver"})
        try:
            length = int(self.headers.get("Content-Length", ""))
        except ValueError:
            return self._reply(411, {"error": "Content-Length required"})
        if length < 0 or length > MAX_BODY:
            return self._reply(413, {"error": "request body too large or invalid"})

        try:
            data = json.loads(self.rfile.read(length)) if length else {}
        except (ValueError, RecursionError):
            return self._reply(400, {"error": "invalid JSON"})
        try:
            result = importlib.import_module(name).run(data)
            self._reply(200, result)
        except Exception:
            # the traceback goes to the subserver log, never to the client
            traceback.print_exc()
            self._reply(500, {"error": "internal error, see the subserver log"})

    def _reply(self, code, payload):
        body = json.dumps(payload).encode("utf-8")
        self.send_response(code)
        self.send_header("Content-Type", "application/json")
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)


if __name__ == "__main__":
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    # Eager-load the model so the first request is fast.
    if ONLY:
        importlib.import_module(ONLY)

    print(f"Serving {ONLY or 'any module'} on {HOST}:{PORT}", file=sys.stderr)
    HTTPServer((HOST, PORT), Handler).serve_forever()
