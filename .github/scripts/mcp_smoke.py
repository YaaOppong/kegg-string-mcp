"""Drive a real MCP stdio handshake against the built container.

Two things this has to get right, both learned the hard way:

* Keep stdin open across requests. Piping them in with `printf` closes stdin
  before the server answers the second one, which reads as a broken server.
* Never block without a deadline. `readline()` on a live pipe blocks forever, so
  a container that starts but never answers would hang the job until GitHub's
  six-hour timeout instead of failing in a minute. stderr goes to a file rather
  than a pipe for the same reason -- an undrained pipe fills at 64KB and wedges
  the child.
"""

import json
import selectors
import subprocess
import sys
import tempfile
import time

EXPECTED = {"kegg_pathways", "string_partners", "pubmed_abstracts", "uniprot_protein"}
TIMEOUT = 60


def main(image: str) -> int:
    # Not a context manager: the handle must outlive the try/except so the failure
    # path can read back what the container wrote.
    errlog = tempfile.NamedTemporaryFile("w+", suffix=".err", delete=False)  # noqa: SIM115
    try:
        proc = subprocess.Popen(
            ["docker", "run", "-i", "--rm", image],
            stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=errlog,
            text=True, bufsize=1,
        )
    except (OSError, FileNotFoundError) as exc:
        print(f"could not start container: {exc}", file=sys.stderr)
        return 1
    selector = selectors.DefaultSelector()
    selector.register(proc.stdout, selectors.EVENT_READ)

    def fail(message: str) -> int:
        print(message, file=sys.stderr)
        proc.kill()
        errlog.flush()
        errlog.seek(0)
        print("--- container stderr ---", file=sys.stderr)
        print(errlog.read()[:4000], file=sys.stderr)
        return 1

    def send(payload: dict) -> None:
        proc.stdin.write(json.dumps(payload) + "\n")
        proc.stdin.flush()

    def read(want_id: int):
        deadline = time.monotonic() + TIMEOUT
        while time.monotonic() < deadline:
            if not selector.select(timeout=max(0.0, deadline - time.monotonic())):
                return None                      # timed out with nothing to read
            line = proc.stdout.readline()
            if not line:
                return None                      # pipe closed: the server exited
            try:
                message = json.loads(line)
            except ValueError:
                continue
            if message.get("id") == want_id:
                return message
        return None

    send({"jsonrpc": "2.0", "id": 1, "method": "initialize",
          "params": {"protocolVersion": "2025-06-18", "capabilities": {},
                     "clientInfo": {"name": "ci", "version": "1"}}})
    initialized = read(1)
    if not initialized:
        return fail(f"no response to initialize within {TIMEOUT}s")
    print("initialize:", initialized["result"]["serverInfo"])

    send({"jsonrpc": "2.0", "method": "notifications/initialized"})
    send({"jsonrpc": "2.0", "id": 2, "method": "tools/list"})
    listed = read(2)
    if not listed:
        return fail(f"no response to tools/list within {TIMEOUT}s")

    names = {t["name"] for t in listed["result"]["tools"]}
    print("tools/list:", sorted(names))

    proc.stdin.close()
    try:
        proc.wait(timeout=30)
    except subprocess.TimeoutExpired:
        # A clean handshake is what is being tested; a server that ignores EOF is
        # a separate (and lesser) problem, so say so rather than failing silently.
        print("warning: container did not exit within 30s of stdin close", file=sys.stderr)
        proc.kill()

    missing = EXPECTED - names
    if missing:
        return fail(f"missing tools: {sorted(missing)}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1]))
