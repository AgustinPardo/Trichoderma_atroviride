find . -type f -name *.py | cut -c 3- | sed 's/^/!/' > .gitignore
