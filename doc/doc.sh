export PYTHONPATH=..
find .. -maxdepth 1 -name "*.py" -exec pdoc --html --overwrite "{}" ";"
