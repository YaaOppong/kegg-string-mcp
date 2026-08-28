"""Drive a real MCP stdio handshake against the built container.

Keeps stdin open across requests: piping them in with `printf` closes stdin
before the server answers the second one, which reads as a broken server.
"""

import json
import subprocess
import sys
import time

EXPECTED = {"kegg_pathways", "string_partners", "pubmed_abstracts", "uniprot_protein"}


def main(image: str) -> int:
    proc = subprocess.Popen(
        ["docker", "run", "-i", "--rm", image],
        stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        text=True, bufsize=1,
    )

    def send(payload):
        proc.stdin.write(json.dumps(payload) + "\n")
        proc.stdin.flush()

    def read(want_id, timeout=60):
        deadline = time.time() + timeout
        while time.time() < deadline:
            line = proc.stdout.readline()
            if not line:
                return None
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
        print("no response to initialize", file=sys.stderr)
        print(proc.stderr.read()[:2000], file=sys.stderr)
        return 1
    print("initialize:", initialized["result"]["serverInfo"])

    send({"jsonrpc": "2.0", "method": "notifications/initialized"})
    send({"jsonrpc": "2.0", "id": 2, "method": "tools/list"})
    listed = read(2)
    if not listed:
        print("no response to tools/list", file=sys.stderr)
        print(proc.stderr.read()[:2000], file=sys.stderr)
        return 1

    names = {t["name"] for t in listed["result"]["tools"]}
    print("tools/list:", sorted(names))
    proc.stdin.close()
    proc.wait(timeout=30)

    missing = EXPECTED - names
    if missing:
        print(f"missing tools: {sorted(missing)}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1]))
