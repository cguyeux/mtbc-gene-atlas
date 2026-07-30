"""In-process smoke test of the gene site (no port bind)."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from fastapi.testclient import TestClient  # noqa: E402

from backend.app import app  # noqa: E402


def run() -> None:
    with TestClient(app) as c:   # triggers startup -> init_db + auto_ingest
        r = c.get("/")
        assert r.status_code == 200, r.status_code
        assert "genes catalogued" in r.text
        print("home OK")

        r = c.get("/genes?hypo=1")
        assert r.status_code == 200 and "matched" in r.text
        print("genes list (hypothetical filter) OK")

        for rv, needle in [("Rv0004", "DciA"), ("Rv0025", "DUF4226"),
                           ("Rv0027", "EspC"), ("Rv0021c", "FMN_dh")]:
            r = c.get(f"/gene/{rv}")
            assert r.status_code == 200, (rv, r.status_code)
            assert needle in r.text, (rv, needle)
            assert "Sources" in r.text
            print(f"gene {rv} OK (contains {needle!r} + Sources)")

        r = c.get("/genes?verdict=requalified")
        assert r.status_code == 200
        print("verdict filter OK")

        r = c.get("/api/genes?q=Rv0004")
        assert r.status_code == 200 and r.json(), r.text
        print("api OK ->", r.json()[0])

        r = c.post("/feedback", data={"gene_rv": "Rv0025", "message": "test feedback"})
        assert r.status_code == 200  # followed redirect
        print("feedback OK")

        r = c.get("/gene/Rv9999")
        assert r.status_code == 200 and "not found" in r.text.lower()
        print("404-ish gene OK")

    print("\nALL SMOKE TESTS PASSED")


if __name__ == "__main__":
    run()
