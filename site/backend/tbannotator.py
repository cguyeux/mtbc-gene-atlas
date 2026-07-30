"""Resolve an SRA accession to its SPDI variant set via the TBannotator MCP server.

The deployed container has no local strain database, so to let a user type a public
accession (SRR/ERR/DRR...) and get a resistance profile, we query the read-only
TBannotator PostgreSQL through its MCP HTTP endpoint and read back the strain's
SPDI list. Results are cached per accession to avoid re-querying.

This is a soft dependency: if TBannotator is unreachable, the SRA mode fails
gracefully and the user can still paste variants directly.
"""
from __future__ import annotations

import ast
import csv
import io
import re
from functools import lru_cache

MCP_URL = "https://darthos.freeboxos.fr/mcp"
CHROM = "NC_000962.3"
_ACCESSION = re.compile(r"^[A-Za-z0-9_.-]{4,40}$")   # sanitise before string-formatting into SQL

_SQL = ("SELECT sp.spdi_variant_name FROM tb_report_strain st "
        "JOIN tb_report_strain_spdi ss ON ss.strain_id = st.strain_id "
        "JOIN tb_report_spdi sp ON sp.spdi_id = ss.spdi_id "
        "WHERE st.strain_name = '{sra}'")


async def _call(sql: str, timeout_s: int) -> str:
    from mcp import ClientSession
    from mcp.client.streamable_http import streamablehttp_client
    async with streamablehttp_client(MCP_URL) as (read, write, _):
        async with ClientSession(read, write) as s:
            await s.initialize()
            r = await s.call_tool("tool_query_postgres",
                                  {"query": sql, "max_rows": 100000, "timeout_seconds": timeout_s})
            obj = ast.literal_eval(r.content[0].text)   # server returns a python-repr dict
            if not obj.get("success"):
                raise RuntimeError(obj.get("error", "query failed"))
            return obj["data"]["csv"]


class TBannotatorError(RuntimeError):
    pass


@lru_cache(maxsize=512)
def fetch_spdi(accession: str) -> frozenset[str]:
    """SRA/ENA accession -> frozenset of SPDI. Raises TBannotatorError on failure."""
    acc = (accession or "").strip()
    if not _ACCESSION.match(acc):
        raise TBannotatorError("invalid accession (expected something like SRR1234567)")
    import anyio
    try:
        csv_txt = anyio.run(_call, _SQL.format(sra=acc), 30)
    except Exception as e:                              # network / MCP / server error
        raise TBannotatorError(f"TBannotator unreachable: {type(e).__name__}") from e
    rows = csv.DictReader(io.StringIO(csv_txt))
    spdis = frozenset(r["spdi_variant_name"] for r in rows
                      if (r.get("spdi_variant_name") or "").startswith(CHROM + ":"))
    return spdis
