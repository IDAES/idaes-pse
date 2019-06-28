# coding: utf-8

import pathlib, importlib, time
import idaes

def importr(root, max_sec=10):
    root = pathlib.Path(root)
    base = root.parent
    failures, total = {}, 0
    for path in root.rglob("*.py"):
        if path.parts[-1] == '__init__.py':
            continue
        module_path = path.relative_to(base).with_suffix('')
        module_name = '.'.join(module_path.parts)
        try:
            start = time.time()
            importlib.import_module(module_name)
            sec = time.time() - start
            if sec > max_sec:
                raise ImportError(f"Import took too long ({sec:.1f}s)")
            print(f"{sec:3.1f}s :: {module_name}")
        except ImportError as e:
            failures[module_name] = str(e)
        total += 1
    return failures, total

def test_import():
    root_dir = pathlib.Path(idaes.__file__).parent
    failures, total = importr(root_dir)
    n = len(failures)
    assert n ==0, f"{n:d} failures in {total:d} tests: {failures}"

