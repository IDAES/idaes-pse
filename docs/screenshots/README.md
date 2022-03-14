# Screenshot tools for documentation

## Shot-scraper

This directory has configuration for the "shot-scraper" tool.

See https://pypi.org/project/shot-scraper/

The main script to automate the screenshots is: `sshot.py`.
Information on usage is in the docstring at the top of this script.

The script depends on the existence of Python scripts that print to standard output
one line in the following format (they may print other things, too):

    [scraper] <URL>

The `<URL>` represents the URL that is given to shot-scraper to convert into an image.
