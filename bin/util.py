
import argparse
import logging
import os
from pathlib import Path

_log_name = None


def get_main_logger(name):
    """Create the top-level logger."""
    global _log_name
    _log_name = name
    logging.basicConfig(
        format='[%(asctime)s - %(name)s] %(message)s',
        datefmt='%H:%M:%S', level=logging.INFO)
    return logging.getLogger(name)


def get_named_logger(name):
    """Create a logger with a name.

    :param name: name of logger.
    """
    name = name.ljust(10)[:10]  # so logging is aligned
    logger = logging.getLogger('{}.{}'.format(_log_name, name))
    return logger


def my_parser(name):
    """Make an argument parser for a workflow command."""
    return argparse.ArgumentParser(
        name,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        add_help=False)


def _log_level():
    """Parser to set logging level and acquire software version/commit."""
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter, add_help=False)

    modify_log_level = parser.add_mutually_exclusive_group()
    modify_log_level.add_argument(
        '--debug', action='store_const',
        dest='log_level', const=logging.DEBUG, default=logging.INFO,
        help='Verbose logging of debug information.')
    modify_log_level.add_argument(
        '--quiet', action='store_const',
        dest='log_level', const=logging.WARNING, default=logging.INFO,
        help='Minimal logging; warnings only.')

    return parser



def check_and_create_folder(folder):
    """
    Check if the specified folder exists; if not, create it.
    """
    if not os.path.exists(folder):
        os.makedirs(folder)



def check_folder_for_files(folder_path):
    """
    Check if the specified folder contains any files using pathlib.
    """
    folder = Path(folder_path)

    if not folder.exists():
        print(f"The folder '{folder_path}' does not exist.")
        return False

    # Check if there are any files in the folder
    if any(folder.iterdir()):  # Check if there's any item in the folder (file or directory)
        if any(item.is_file() for item in folder.iterdir()):
            print(f"The folder '{folder_path}' contains files.")
            return True
        else:
            print(f"The folder '{folder_path}' does not contain any files, only directories.")
            return False
    else:
        print(f"The folder '{folder_path}' is empty.")
        return False



