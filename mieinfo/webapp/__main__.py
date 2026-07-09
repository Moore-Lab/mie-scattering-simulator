"""``python -m mieinfo.webapp`` launches the interactive pattern viewer."""
from __future__ import annotations

import argparse

from .app import main

if __name__ == "__main__":
    ap = argparse.ArgumentParser(prog="mieinfo.webapp")
    ap.add_argument("--host", default="127.0.0.1")
    ap.add_argument("--port", type=int, default=5000)
    a = ap.parse_args()
    main(host=a.host, port=a.port)
