import sqlite3
import csv
from pathlib import Path

SCHEMA_SQL = """
PRAGMA foreign_keys = ON;

CREATE TABLE IF NOT EXISTS projects (
  project_id INTEGER PRIMARY KEY,
  project_code TEXT NOT NULL UNIQUE
);

CREATE TABLE IF NOT EXISTS subjects (
  subject_id INTEGER PRIMARY KEY,
  project_id INTEGER NOT NULL,
  subject_code TEXT NOT NULL,
  condition TEXT,
  age INTEGER,
  sex TEXT,
  treatment TEXT,
  response TEXT,
  UNIQUE(project_id, subject_code),
  FOREIGN KEY(project_id) REFERENCES projects(project_id)
);

CREATE TABLE IF NOT EXISTS samples (
  sample_id INTEGER PRIMARY KEY,
  subject_id INTEGER NOT NULL,
  sample_code TEXT NOT NULL UNIQUE,
  sample_type TEXT,
  time_from_treatment_start REAL,
  FOREIGN KEY(subject_id) REFERENCES subjects(subject_id)
);

CREATE TABLE IF NOT EXISTS cell_populations (
  population_id INTEGER PRIMARY KEY,
  population_name TEXT NOT NULL UNIQUE
);

CREATE TABLE IF NOT EXISTS cell_counts (
  sample_id INTEGER NOT NULL,
  population_id INTEGER NOT NULL,
  cell_count INTEGER NOT NULL,
  PRIMARY KEY(sample_id, population_id),
  FOREIGN KEY(sample_id) REFERENCES samples(sample_id),
  FOREIGN KEY(population_id) REFERENCES cell_populations(population_id)
);

CREATE INDEX IF NOT EXISTS idx_samples_subject_id ON samples(subject_id);
CREATE INDEX IF NOT EXISTS idx_cell_counts_population_id ON cell_counts(population_id);
"""

META_COLUMNS = [
    "project",
    "subject",
    "condition",
    "age",
    "sex",
    "treatment",
    "response",
    "sample",
    "sample_type",
    "time_from_treatment_start",
]


def init_db(conn: sqlite3.Connection) -> None:
    conn.executescript(SCHEMA_SQL)


def get_or_create_id(conn: sqlite3.Connection, insert_sql: str, select_sql: str, value_tuple):
    conn.execute(insert_sql, value_tuple)
    row = conn.execute(select_sql, value_tuple).fetchone()
    if row is None:
        raise RuntimeError("Failed to fetch id after insert")
    return row[0]


def load_csv(db_path: Path, csv_path: Path) -> None:
    with sqlite3.connect(db_path) as conn:
        init_db(conn)

        with csv_path.open(newline="", encoding="utf-8") as f:
            reader = csv.DictReader(f)
            if reader.fieldnames is None:
                raise ValueError("CSV missing header row")

            headers = [h.strip() for h in reader.fieldnames]
            population_cols = [h for h in headers if h not in META_COLUMNS]

            for pop in population_cols:
                conn.execute(
                    "INSERT OR IGNORE INTO cell_populations(population_name) VALUES (?);",
                    (pop,),
                )

            for row in reader:
                project = row["project"].strip()
                subject = row["subject"].strip()
                sample = row["sample"].strip()

                project_id = get_or_create_id(
                    conn,
                    "INSERT OR IGNORE INTO projects(project_code) VALUES (?);",
                    "SELECT project_id FROM projects WHERE project_code = ?;",
                    (project,),
                )

                conn.execute(
                    """
                    INSERT OR IGNORE INTO subjects(
                      project_id, subject_code, condition, age, sex, treatment, response
                    )
                    VALUES (?, ?, ?, ?, ?, ?, ?);
                    """,
                    (
                        project_id,
                        subject,
                        row.get("condition"),
                        int(row["age"]) if row.get("age") else None,
                        row.get("sex"),
                        row.get("treatment"),
                        row.get("response"),
                    ),
                )
                subject_id = conn.execute(
                    "SELECT subject_id FROM subjects WHERE project_id = ? AND subject_code = ?;",
                    (project_id, subject),
                ).fetchone()[0]

                conn.execute(
                    """
                    INSERT OR IGNORE INTO samples(
                      subject_id, sample_code, sample_type, time_from_treatment_start
                    )
                    VALUES (?, ?, ?, ?);
                    """,
                    (
                        subject_id,
                        sample,
                        row.get("sample_type"),
                        float(row["time_from_treatment_start"]) if row.get("time_from_treatment_start") else None,
                    ),
                )
                sample_id = conn.execute(
                    "SELECT sample_id FROM samples WHERE sample_code = ?;",
                    (sample,),
                ).fetchone()[0]

                for pop in population_cols:
                    raw = row.get(pop)
                    if raw is None or raw == "":
                        continue
                    cell_count = int(float(raw))

                    population_id = conn.execute(
                        "SELECT population_id FROM cell_populations WHERE population_name = ?;",
                        (pop,),
                    ).fetchone()[0]

                    conn.execute(
                        """
                        INSERT INTO cell_counts(sample_id, population_id, cell_count)
                        VALUES (?, ?, ?)
                        ON CONFLICT(sample_id, population_id) DO UPDATE SET
                          cell_count = excluded.cell_count;
                        """,
                        (sample_id, population_id, cell_count),
                    )


def main():
    db_path = Path("cell_counts.sqlite3")
    csv_path = Path("data/cell-count.csv")
    load_csv(db_path, csv_path)
    print("Done. Created", db_path)


if __name__ == "__main__":
    main()

