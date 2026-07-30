from collections.abc import Iterator

from sqlalchemy import create_engine
from sqlalchemy.orm import DeclarativeBase, Session, sessionmaker

from .config import settings

engine = create_engine(
    settings.database_url,
    connect_args={"check_same_thread": False} if settings.database_url.startswith("sqlite") else {},
    future=True,
)
SessionLocal = sessionmaker(engine, expire_on_commit=False, autoflush=False)


class Base(DeclarativeBase):
    pass


def get_db() -> Iterator[Session]:
    db = SessionLocal()
    try:
        yield db
    finally:
        db.close()


def init_db() -> None:
    from . import models  # noqa: F401  registers the mappers

    # SQLite default path may not exist yet.
    if settings.database_url.startswith("sqlite"):
        (settings.database_url.split("///", 1)[-1])  # noqa
        from pathlib import Path
        Path(settings.database_url.split("///", 1)[-1]).parent.mkdir(parents=True, exist_ok=True)
    Base.metadata.create_all(engine)
