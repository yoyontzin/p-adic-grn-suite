"""API mínima para demo institucional."""

from __future__ import annotations

from http.server import BaseHTTPRequestHandler, HTTPServer
import json


class Handler(BaseHTTPRequestHandler):
    def do_GET(self) -> None:  # noqa: N802
        if self.path == "/health":
            payload = {"status": "ok", "service": "web-demo-backend"}
            body = json.dumps(payload).encode("utf-8")
            self.send_response(200)
            self.send_header("Content-Type", "application/json; charset=utf-8")
            self.send_header("Content-Length", str(len(body)))
            self.end_headers()
            self.wfile.write(body)
            return
        self.send_response(404)
        self.end_headers()


def main() -> None:
    server = HTTPServer(("127.0.0.1", 8080), Handler)
    print("Backend demo escuchando en http://127.0.0.1:8080")
    server.serve_forever()


if __name__ == "__main__":
    main()
