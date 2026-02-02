import sqlite3
import pandas as pd

db_path = "cell_counts.sqlite3"  # change if needed
params = ("melanoma", "PBMC", "miraclib", 0)

def print_df(df, n=20):
    if df.empty:
        print("No rows")
        return
    print(df.head(n).to_string(index=False))

with sqlite3.connect(db_path) as conn:
    q_samples = """
    SELECT
      p.project_code,
      s.subject_code,
      sa.sample_code,
      sa.sample_type,
      sa.time_from_treatment_start,
      s.treatment,
      s.response,
      s.sex
    FROM samples sa
    JOIN subjects s ON s.subject_id = sa.subject_id
    JOIN projects p ON p.project_id = s.project_id
    WHERE s.condition = ?
      AND sa.sample_type = ?
      AND s.treatment = ?
      AND sa.time_from_treatment_start = ?
    ORDER BY p.project_code, s.subject_code, sa.sample_code;
    """
    baseline_samples = pd.read_sql_query(q_samples, conn, params=params)
    print("Baseline samples matching filters:", len(baseline_samples))
    print_df(baseline_samples, n=20)
    print()

    q_samples_per_project = """
    SELECT
      p.project_code AS project,
      COUNT(DISTINCT sa.sample_code) AS n_samples
    FROM samples sa
    JOIN subjects s ON s.subject_id = sa.subject_id
    JOIN projects p ON p.project_id = s.project_id
    WHERE s.condition = ?
      AND sa.sample_type = ?
      AND s.treatment = ?
      AND sa.time_from_treatment_start = ?
    GROUP BY p.project_code
    ORDER BY n_samples DESC;
    """
    samples_per_project = pd.read_sql_query(q_samples_per_project, conn, params=params)
    print("Samples per project")
    print_df(samples_per_project, n=100)
    print()

    q_subjects_by_response = """
    SELECT
      COALESCE(s.response, 'unknown') AS response,
      COUNT(DISTINCT s.subject_id) AS n_subjects
    FROM samples sa
    JOIN subjects s ON s.subject_id = sa.subject_id
    WHERE s.condition = ?
      AND sa.sample_type = ?
      AND s.treatment = ?
      AND sa.time_from_treatment_start = ?
    GROUP BY COALESCE(s.response, 'unknown')
    ORDER BY n_subjects DESC;
    """
    subjects_by_response = pd.read_sql_query(q_subjects_by_response, conn, params=params)
    print("Subjects by response")
    print_df(subjects_by_response, n=100)
    print()

    q_subjects_by_sex = """
    SELECT
      COALESCE(s.sex, 'unknown') AS sex,
      COUNT(DISTINCT s.subject_id) AS n_subjects
    FROM samples sa
    JOIN subjects s ON s.subject_id = sa.subject_id
    WHERE s.condition = ?
      AND sa.sample_type = ?
      AND s.treatment = ?
      AND sa.time_from_treatment_start = ?
    GROUP BY COALESCE(s.sex, 'unknown')
    ORDER BY n_subjects DESC;
    """
    subjects_by_sex = pd.read_sql_query(q_subjects_by_sex, conn, params=params)
    print("Subjects by sex")
    print_df(subjects_by_sex, n=100)

