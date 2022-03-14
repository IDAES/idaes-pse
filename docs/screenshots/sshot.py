"""
Automate screenshots.

Depends on 'shot-scraper' (pip package).

Uses a YAML configuration file, e.g., `sshot.yaml`, that is in a given
directory.

The YAML file is a mapping using the following schema:

    desc: "Describe the app/tool for this subdirectory"
    screens:
        <screen_name_1>:
            script: <script_to_run>.py
            output: <image_filename>.png
        <screen_name_2>:
                ...

The program will run the `script` and look for a line of output to standard
output that looks like:

    [scraper] <URL>

The program will use the value of `<URL>` as the URl to feed to shot-scraper, which
will fire up a web page at that URL and then convert the screen to an image, which
will be saved in the `output` image filename.

Additional options for the shot-scraper program (e.g., width, height, but also
any advanced option) can be given in an `options` section for the screen; they will
be copied into the shot-scraper configuration. See shot-scraper documentation for
details.

Basic usage:

    python sshot.py

This will look for subdirectories containing the file "sshot.yaml".
For each subdirectory, it will 'chdir' to it, execute the commands in the configuration, and stop.
There is no guarantee as to the order in which the directories are visited.

You can provide a starting (and possibly terminal) directory:

    python sshot.py my_app_dir

You can also select screens (or patterns of screens) from the command line.
You can provide either a simple string name or a comma-searated list of regular expressions as a filter for which section(s) or screen(s) to run.

Selecting the "hda" screen:

    python sshot.py --screen hda

"""
import argparse
import logging
import os
from pathlib import Path
import re
from subprocess import Popen, PIPE
import sys
import time
from typing import List
import yaml

# Logging setup
_log = logging.getLogger("screenshot")
_h = logging.StreamHandler()
_h.setFormatter(
    logging.Formatter(fmt="({name}) {asctime} [{levelname}] {message}", style="{")
)
_log.addHandler(_h)


class ConfigError(Exception):
    pass


class ExecError(Exception):
    pass


class Screenshot:
    encoding = "utf-8"

    def __init__(self, conf):
        if isinstance(conf, dict):
            self._conf = conf
        elif hasattr(conf, "read"):
            try:
                self._conf = yaml.safe_load(conf)
            except Exception as err:
                raise ConfigError(f"Cannot load configuration file: {err}")
        # Regex filters, see filter_match() functions
        self._screens = set()

    @property
    def screens(self):
        return self._screens.copy()

    @screens.setter
    def screens(self, value):
        for v in value:
            self._screens.add(re.compile(v, flags=re.I))

    def run(self):
        desc = self._conf.get("desc", "No description")
        try:
            screens = self._conf["screens"]
        except KeyError:
            raise ConfigError(f"Missing 'screens' configuration")
        # Create and scrape each 'screen'
        for label, info in screens.items():
            if not filter_match(label, self._screens):
                _log.debug(f"Skip screen. name={label}")
                continue
            _log.info(f"Run screen. name={label}")
            try:
                output_file = info["output"]
            except KeyError:
                raise ConfigError(
                    f"Screen '{label}' is missing output file field 'output'"
                )
            try:
                script = info["script"]
            except KeyError:
                raise ConfigError(
                    f"Screen '{label}' is missing script file field 'script'"
                )
            if not Path(script).exists():
                raise ConfigError(
                    f"Screen '{label}': script file '{script}' does not exist"
                )
            _log.info(
                f"Run script. name={script}"
            )
            url, proc = self._run_script(script)
            _log.info(
                f"Scrape screen. name={label}, output={output_file}"
            )
            options = info.get("options", {})
            try:
                self._scrape_screen(url, output_file, options)
            except ExecError as err:
                _log.error(f"Screen scraping error. Error={err}")
            proc.kill()
            time.sleep(1)

    def _run_script(self, script):
        if script.endswith(".py"):
            script_args = ["python", script]
        else:
            script_args = [script]
        proc = Popen(script_args, stdout=PIPE)
        while 1:
            proc.poll()
            if proc.returncode is not None:
                raise ExecError(f"Script had an error. returncode={proc.returncode}")
            line = proc.stdout.readline().strip().decode(self.encoding)
            if line.startswith("[scraper]"):
                url = line.split(" ")[1]
                break
        return url, proc

    def _scrape_screen(self, url, ofile, options):
        # Create a YAML config for the screen scraper tool
        opath = Path(ofile)
        conf_file = change_suffix(opath, ".yaml")
        with open(conf_file, "w", encoding=self.encoding) as f:
            f.write(f"- url: {url}\n")
            f.write(f"  output: {ofile}\n")
            f.write(f"  wait: 1000\n")
            for okey, ovalue in options.items():
                f.write(f"  {okey}: {ovalue}\n")
        # Run the screen scraper tool
        try:
            opath.unlink(missing_ok=True)
        except Exception as err:
            _log.warning(f"Failed to remove old output file. Error={err}")
        proc = Popen(["shot-scraper", "multi", str(conf_file)])
        _log.info(f"Scraping url. input={url}, output={ofile}")
        proc.wait(120)
        try:
            conf_file.unlink()
        except Exception:
            _log.error(f"Unable to remove temporary config file. path={conf_file}")
        if proc.returncode != 0:
            raise ExecError(f"Scraping screen failed. url={url} code={proc.returncode}")


def filter_match(s: str, expressions: List[re.Pattern]):
    if not expressions:
        is_match = True
    else:
        is_match = False
        for e in expressions:
            if e.match(s):
                is_match = True
                break
    return is_match


def change_suffix(p, suffix):
    new_name = Path(f"{p.stem}{suffix}")
    file_parts = list(p.parts)
    file_parts[-1] = new_name
    return Path(*file_parts)

g_conf_filename = "sshot.yaml"

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("dirname", help="Directory root for files", nargs="?", default=".")
    ap.add_argument("--screen", "-s", help="Regex to select screen(s)", default="")
    ap.add_argument(
        "--verbose", "-v", action="count", help="Increase verbosity", default=0
    )
    args = ap.parse_args()
    if args.verbose > 1:
        _log.setLevel(logging.DEBUG)
    elif args.verbose > 0:
        _log.setLevel(logging.INFO)
    else:
        _log.setLevel(logging.WARNING)
    for root, dirs, files in os.walk(args.dirname):
        _log.debug(f"Examine directory. path={root}")
        if g_conf_filename in files:
            _log.info(f"Generating screenshots. path={root}")
            cwd = os.getcwd()
            os.chdir(root)
            status = run_subdir(args)
            _log.info(f"Generated screenshots. path={root}, code={status}")
            os.chdir(cwd)
        else:
            _log.debug(f"Skipping directory. path={root}")

def run_subdir(args):
    with Path(g_conf_filename).open() as f:
        shot = Screenshot(f)
    if args.screen:
        shot.screens = args.screen.split(",")
    try:
        shot.run()
    except ConfigError as err:
        _log.error(f"Invalid configuration: {err}")
        return -1
    except ExecError as err:
        _log.error(f"Execution failed: {err}")
        return -1
    return 0


if __name__ == "__main__":
    sys.exit(main())


