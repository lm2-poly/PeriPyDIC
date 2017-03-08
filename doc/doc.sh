export PYTHONPATH=..
find .. -maxdepth 2 -name "*.py" -exec pdoc --html --overwrite "{}" ";"
