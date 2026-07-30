from pathlib import Path

from pydantic_settings import BaseSettings, SettingsConfigDict

SITE_ROOT = Path(__file__).resolve().parent.parent          # annotation_mtbc/site/
PROJECT_ROOT = SITE_ROOT.parent                              # annotation_mtbc/
CONTENT = SITE_ROOT / "content"


class Settings(BaseSettings):
    model_config = SettingsConfigDict(env_file=".env", env_prefix="ANNOT_", extra="ignore")

    database_url: str = f"sqlite:///{SITE_ROOT / 'data' / 'db.sqlite'}"
    site_title: str = "MTBC Gene Atlas"
    auto_ingest: bool = True

    # Catalogue of all joined genes (tab-separated). Built by analyses/phase1.
    # In Docker a copy lives under content/; for local dev we fall back to the
    # pipeline output in the parent project.
    gene_xref: Path = CONTENT / "gene_xref.tsv"
    gene_xref_fallback: Path = PROJECT_ROOT / "data" / "gene_xref.tsv"
    # Enriched per-gene JSON records (built by analyses/phase4).
    genes_content: Path = CONTENT / "genes"
    # Separate microproteome track (built by analyses/phase83); may be absent.
    microproteins_content: Path = CONTENT / "microproteins.json"

    def resolved_xref(self) -> Path:
        return self.gene_xref if self.gene_xref.exists() else self.gene_xref_fallback


settings = Settings()
