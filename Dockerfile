# MCP server over stdio. Run with:  docker run -i --rm kegg-string-mcp
#
# stdio transport means the container must not print anything to stdout that is
# not JSON-RPC, so nothing here logs to stdout.
FROM python:3.12-slim AS build

WORKDIR /src
COPY pyproject.toml README.md LICENSE ./
COPY src ./src
RUN pip install --no-cache-dir --prefix=/install .

FROM python:3.12-slim

# Non-root: the server only makes outbound HTTP requests and writes a cache.
RUN useradd --create-home --uid 10001 mcp
COPY --from=build /install /usr/local

ENV KEGG_STRING_MCP_CACHE=/home/mcp/.cache/kegg-string-mcp \
    PYTHONUNBUFFERED=1
USER mcp
WORKDIR /home/mcp
RUN mkdir -p "$KEGG_STRING_MCP_CACHE"

# NCBI and STRING ask callers to identify themselves; set these at run time:
#   docker run -i --rm -e NCBI_EMAIL=you@example.org -e STRING_CALLER_IDENTITY=you ...
ENTRYPOINT ["kegg-string-mcp"]
