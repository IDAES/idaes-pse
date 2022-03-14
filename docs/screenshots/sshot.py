"""
Automate screenshots.
"""
import argparse
from pathlib import Path
from subprocess import Popen, PIPE
import sys
import time
import yaml


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

    def run(self):
        for section, contents in self._conf.items():
            desc = contents.get("desc", section)
            print(f"Section: {desc}")
            self._run_section(section, contents)

    def _run_section(self, section, contents):
        try:
            screens = contents["screens"]
        except KeyError:
            raise ConfigError(f"Section '{section}' is missing 'screens' configuration")
        for label, info in screens.items():
            try:
                output_file = info["output"]
            except KeyError:
                raise ConfigError(f"Screen '{label}' in section '{section}' is "
                                  f"missing output file field 'output'")
            try:
                script = info["script"]
            except KeyError:
                raise ConfigError(f"Screen '{label}' in section '{section}' is "
                                  f"missing script file field 'script'")
            if not Path(script).exists():
                raise ConfigError(f"Screen '{label}' in section '{section}': script "
                                  f"file '{script}' does not exist")
            print(f"Screen {label}: {script} -> {output_file}")
            url, proc = self._run_script(script)
            options = info.get("options", {})
            self._scrape_screen(url, output_file, options)
            proc.kill()
            time.sleep(1)

    def _run_script(self, script):
        if script.endswith(".py"):
            script_args = ["python", script]
        else:
            script_args = [script]
        proc = Popen(script_args, stdout=PIPE)
        while 1:
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
            print(f"Failed to remove old output file: {err}")
        with Popen(["shot-scraper", "multi", str(conf_file)]):
            print(f"Scraping url {url} -> {ofile}")


def change_suffix(p, suffix):
    new_name = Path(f"{p.stem}{suffix}")
    file_parts = list(p.parts)
    file_parts[-1] = new_name
    return Path(*file_parts)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("conf", help="Configuration file", type=argparse.FileType("r"))
    args = ap.parse_args()
    shot = Screenshot(args.conf)
    try:
        shot.run()
    except ConfigError as err:
        print(f"Configuration error! {err}")
        return -1
    return 0


if __name__ == "__main__":
    sys.exit(main())
