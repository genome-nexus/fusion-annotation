"""Guards against a regression that took down the deployed MCP server: the
FUSION_ANNOTATION_ALLOWED_HOSTS env var was set to a full URL (scheme + path)
instead of a bare hostname, so the DNS-rebinding Host-header check never
matched the real Host header and rejected 100% of traffic with 421.
server.app.normalize_allowed_host() defensively extracts the bare host[:port]
from whatever is configured.
"""
import os
import sys

import pytest

# server/app.py pulls in the server extra (starlette, mcp, uvicorn), which isn't
# installed for the zero-dep core test matrix — skip cleanly there. The dedicated
# `test-server` CI job installs those deps so these assertions actually run.
pytest.importorskip("starlette")
pytest.importorskip("mcp")

sys.path.insert(0, os.path.join(os.path.dirname(__file__), os.pardir, "server"))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), os.pardir, "src"))

from app import normalize_allowed_host  # noqa: E402


@pytest.mark.parametrize("entry,expected", [
    ("my-service-abc123.a.run.app", "my-service-abc123.a.run.app"),
    ("https://my-service-abc123.a.run.app/mcp",
     "my-service-abc123.a.run.app"),
    ("http://foo.run.app", "foo.run.app"),
    ("foo.run.app/mcp/", "foo.run.app"),
    ("localhost:8080", "localhost:8080"),
    ("  foo.run.app  ", "foo.run.app"),
])
def test_normalize_allowed_host(entry, expected):
    assert normalize_allowed_host(entry) == expected
