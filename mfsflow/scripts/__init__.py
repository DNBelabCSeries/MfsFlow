"""
Standalone scripts executed as subprocesses by the pipeline.

Each script runs in its own Python process via ``subprocess``. Heavy logic is
kept in importable functions where practical, while ``main()`` handles argv and
configuration loading.
"""
