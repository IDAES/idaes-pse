import logging
import os


def level_num(name, default_level):
    """Get level num from name"""
    level_mapping = dict(ERROR=logging.ERROR, INFO=logging.INFO,
                         WARNING=logging.WARNING, DEBUG=logging.DEBUG)
    return level_mapping.get(name.upper(), default_level)


# set up test logger
_log = logging.getLogger(__name__)

level = logging.WARNING
if 'IDAES_TEST_LOG_LEVEL' in os.environ:
    env_level = os.environ['IDAES_TEST_LOG_LEVEL']
    level = level_num(env_level, level)

_h = logging.StreamHandler()
_h.setFormatter(logging.Formatter('%(asctime)s [%(levelname)s] '
                                  '%(filename)s:%(lineno)d :: %(message)s'))
_log.addHandler(_h)
_log.setLevel(level)
_log.propagate = False