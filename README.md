# Teiko Teiknical

## Instructions for Running Code

- Parts 2 and 3 of Analysis are completed in a Jupyter Notebook called "Initial Data Analysis.ipynb"
    - Just run the notebook top down
- After ensuring that "cell_counts.csv" is located in a subfolder called "data", run "load.py" and 
subsequently, "baseline_query.py"
- "app.py" contains the code for the dashboard

## Data Schema

The schema consists of 5 tables:
1. Projects table
    Consisting of columns:
    a. project
    b. index project_id
2. Subjects table
    Consisting of columns:
    a. project
    b. subject
    c. condition
    d. age
    e. sex
    f. treatment
    g. response
    h. index subject_id
    i. project_id
3. Samples table
    Consisting of columns:
    a. sample
    b. sample_type
    c. time_from_treatment_start
    d. index sample_id
    e. subject_id
4. Populations table
    Consisting of columns:
    a. population
    b. index population_id
5. Cell Counts table
    Consisting of columns:
    a. Cell Count
    b. population_id
    c. sample_id

This type of schema enables different components of data to be kept track of without needing
to keep track of everything together. The individual tables are designed to reference each
other so that if information on an entity such as the subject needs to be accessed, it can
be referenced.

This scales well because with hundreds of different projects to keep track of you only need
to make updates to the projects table indicating what the new project is and the subjects
table indicating which project the subject is tied to. For example another sample from the
same subject of the same project only needs to be added to the samples table and you don't
need to worry about adding a row with the same identifier such as the subject's age or sex
since those don't change.

## Code Structure
I used Jupyter Notebook for the analytics sections in part 2 and 3 because it was convenient
to generate figures and extract/manipulate data at the same time.

For parts 1 and 4 I used SQL to make the data schema and query from the relational database
to generate insights since the directions specified making a database and querying from it.

## Dashboard
http://127.0.0.1:8050/

## LLM Usage
I used OpenAI's ChatGPT to generate code for the data schema and Dashboard.