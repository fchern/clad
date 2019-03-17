import os
import glob
import systemd
import traceback
import os
import datetime
import sys
import subprocess
import traceback
import argparse
import getpass
import configparser


class ZergEnvironment(object):
    def __init__(self, args):
        pass

    def get_benchmark_dir(self):
        return dir

    def get_source_dir(self):
        return dir

    def get_fitting_dir(self):
        return dir

    def get_hive_dir(self):
        return dir


def is_identical(env):
    fs = glob(env.get_fitting_dir())

    for filename in fs:
        new_frame = pd.read_csv(filename)

        benchmark_frame = pd.read_csv(env.get_benchmark_dir() + filename.split('/')[-1])

        if not new_frame.equals(benchmark_frame):
            return False

    return True


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("-g", "--graph", help="the directory that stores the node AND edge file ")

    args = parser.parse_args()

    # config = configparser.ConfigParser()
    source_dir_list = ['python1', 'python2', 'cpp1']

    env = ZergEnvironment(args)

    os.system("compile")
    os.system("run from date to date")
    os.system("run the math")
    os.system("create model file")

    if is_identical(env):
        commit_message = input("please input git commit message: ")

        for target_dir in source_dir_list:

            os.system('gitcommit.sh ' + target_dir + ' ' + commit_message)



if __name__ == "__main__":
    if getpass.getuser() != 'ubuntu':
        raise NameError('operational scripts can only be run by ubuntu')

    try:
        main()
    except Exception as e:
        print('\n')
        traceback.print_exc()
        print("\n")
        sys.exit(1)
