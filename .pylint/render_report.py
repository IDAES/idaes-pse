import argparse
from dataclasses import dataclass
import itertools
import json
from pathlib import Path
import sys
from typing import Iterable, Optional


@dataclass
class PylintRecord:
    type: str
    module: str
    obj: str
    line: int
    column: int
    path: Path
    symbol: str
    message: str
    message_id: str

    path_url: str
    line_url: str = None

    @classmethod
    def from_pylint_output(cls, path: Path, **kwargs) -> Iterable["PylintRecord"]:
        with Path(path).open() as f:
            dicts = json.load(f)
        for d in dicts:
            yield cls.from_dict(d, **kwargs)

    @classmethod
    def from_dict(cls, d: dict, base_url: str = "file://") -> "PylintRecord":
        path = Path(d.pop("path"))
        line = int(d.pop("line"))
        path_url = f"{base_url}/{path}"
        line_url = f"{path_url}#L{line}"
        return cls(
            type=d.pop("type"),
            module=d.pop("module"),
            obj=d.pop("obj"),
            message_id=d.pop("message-id"),
            path=path,
            line=line,
            path_url=path_url,
            line_url=line_url,
            column=int(d.pop("column")),
            symbol=d.pop("symbol"),
            message=d.pop("message"),
        )


@dataclass
class RecordsByModule:
    module: str
    url: str
    records: Iterable[PylintRecord]

    def __iter__(self):
        return iter(self.records)

    @classmethod
    def from_records(
        cls, all_records: Iterable[PylintRecord]
    ) -> Iterable["RecordsByModule"]:
        grouped = itertools.groupby(all_records, key=lambda r: r.module)
        for groupby_key, records_for_module in grouped:
            records = list(records_for_module)
            first = records[0]
            yield cls(module=first.module, url=first.path_url, records=records)


class RecordRenderer:
    def __call__(self, records: Iterable[PylintRecord]) -> Iterable[str]:
        yield from []


class MarkdownChecklistRenderer(RecordRenderer):
    def render_record(self, record: PylintRecord) -> str:
        location = f"At line {record.line}"
        if record.obj:
            location += f", in `{record.obj}`"
        descr = f"(`{record.message_id}`) `{record.message}` (`{record.symbol}`)"
        return f"{location}: ([link]({record.line_url})) {descr}"

    def render_record_group(self, record_group: RecordsByModule) -> Iterable[str]:
        group_header = f"[`{record_group.module}`]({record_group.url})"
        yield f"\n#### {group_header}"
        for rec in record_group:
            yield f"- [ ] {self.render_record(rec)}"

    def __call__(self, records: Iterable[PylintRecord]) -> Iterable[str]:
        for record_group in RecordsByModule.from_records(records):
            yield from self.render_record_group(record_group)


@dataclass
class Options:
    pylint_data_path: Path
    base_url: str
    markdown_output_path: Optional[str] = None

    @classmethod
    def from_cli_args(cls, args: Optional[Iterable[str]] = None) -> "Options":
        parser = argparse.ArgumentParser()
        parser.add_argument("--pylint-data-path", default="./pylint.json")
        parser.add_argument(
            "--base-url",
            default="https://github.com/IDAES/idaes-pse/tree/main",
        )
        parser.add_argument(
            "--markdown-output-path",
            default=None,
        )
        parsed = parser.parse_args(args)

        return cls(
            pylint_data_path=Path(parsed.pylint_data_path),
            base_url=str(parsed.base_url),
            markdown_output_path=Path(parsed.markdown_output_path)
            if parsed.markdown_output_path
            else None,
        )


def main(args=None):
    options = Options.from_cli_args(args=args)

    try:
        records = PylintRecord.from_pylint_output(
            options.pylint_data_path,
            base_url=options.base_url,
        )
    except Exception as e:
        raise e

    # TODO add other renderers as needed
    render = MarkdownChecklistRenderer()

    file_to_print_to = sys.stdout
    if options.markdown_output_path:
        file_to_print_to = open(options.markdown_output_path, "a")
    for line in render(records):
        print(line, file=file_to_print_to)


if __name__ == "__main__":
    main()
