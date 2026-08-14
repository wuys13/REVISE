"""Render the repository README inside Sphinx without duplicating its text."""

from __future__ import annotations

import re

from myst_parser.parsers.sphinx_ import MystParser


class Parser(MystParser):
    """Resolve repository-relative README assets for the documentation site."""

    def parse(self, inputstring, document):
        source = inputstring.replace(
            'src="./logo/',
            'src="./_static/',
        )
        source = source.replace("](png/", "](../png/")

        def repository_url(match):
            target = match.group(1)
            if target.startswith(
                ("http://", "https://", "mailto:", "#", "../png/")
            ):
                return match.group(0)
            route = "tree" if target.endswith("/") else "blob"
            return f"](https://github.com/wuys13/REVISE/{route}/main/{target})"

        source = re.sub(r"\]\(([^)]+)\)", repository_url, source)
        super().parse(source, document)
