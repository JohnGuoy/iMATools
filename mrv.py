# -*- coding: utf-8 -*-

"""
mrv
~~~~~~~~~~~

An integrated methylation pattern region identification and annotation platform iMATools based on long-read sequencing
data is used to process the coding, identification and feature annotation of methylation pattern regions in
ultra-large-scale long-read methylation information, with a view to providing methylation Methylation researchers
provide specialized methylation analysis and visualization tools for long-read sequencing to precisely reveal DNA
methylation patterns at the cellular and read levels.
"""

import sys
import os
import os.path
import shutil
import argparse
import textwrap
import re
import configparser
import hashlib
import sqlite3
import pickle
import portion
import matplotlib
import matplotlib.pyplot as plt
import tqdm

DEBUG = False


CHROMOSOMES_CpG_RANGES = {}

DATA_FILE = ""

DATA_FILE_ROW_COUNT = 0

DATA_FILE_SHA256SUM = ""

OUT_PUTDIR = "."

PREPROCESS_DIR = ""


def parse_args(args):
    prompt = """usage: python mrv.py <--data-file DATAFILE> {<--chromosome CHROMOSOME> <--cpg-range RANGE \
[RANGE ...]> | <--to-visualize-file TOVISUALIZEFILE>} [--output-dir OUTPUTDIR] 
        | --version 
        | --help

mrv is a visualization tool used to visualize whether the CpG site of a certain read data obtained by long-read \
sequencing is methylated.

optional arguments:
  -h, --help            show this help message and exit
  --data-file DATAFILE  file path, specify the path of a text file containing long-read sequencing reads data
  --chromosome CHROMOSOME
                        text, specify chromosome name
  --output-dir OUTPUTDIR
                        directory path, specify the path of the output directory default: current directory
  --version, -V         show mrv version
  --cpg-range RANGE [RANGE ...]
                        text, specify the range of CpG sites. The syntax format of RANGE is: [start,end].
                        You can give several RANGEs at once, such as
                        --cpg-range [1085,1200] [1220,1280] [2300,2395]
  --to-visualize-file TOVISUALIZEFILE
                        file path, specify the path of a text file. You can put the names of multiple
                        chromosomes and their CpG site intervals to be visualized in this text file, and
                        hand them to mrv to calculate and output multiple visual files. An example of
                        TOVISUALIZEFILE:
                        [Y]
                        5431,9587
                        15680,17506
                        12003,12210
                        80,3327
    
                        [KI270580.1]
                        1154,1669
                        756,1321
                        800,1154

examples:
$ python mrv.py --data-file /home/someone/data.txt --chromosome Y --cpg-range [10802025,10861195]

$ python mrv.py --data-file /home/someone/data.txt --chromosome Y --cpg-range [1085,1200] [1220,1280] [2300,2395]

$ python mrv.py --data-file /home/someone/data.txt --to-visualize-file /home/someone/chromosomes_CpGs.txt"""

    if "mrv.py" in sys.argv:
        sys.argv.remove("mrv.py")

    if len(sys.argv) == 0:
        # print("error: you did not specify any parameters.\n")
        print(prompt)
        sys.exit(1)

    prog = "mrv.py"

    usage = """python mrv.py <--data-file DATAFILE> {<--chromosome CHROMOSOME> <--cpg-range RANGE [RANGE ...]> | \
<--to-visualize-file TOVISUALIZEFILE>} [--output-dir OUTPUTDIR] 
        | --version 
        | --help"""

    description = """mrv is a visualization tool used to visualize whether the CpG site of a certain read data \
obtained by long-read sequencing is methylated."""

    epilog = """examples:
$ python mrv.py --data-file /home/someone/data.txt --chromosome Y --cpg-range [10802025,10861195]

$ python mrv.py --data-file /home/someone/data.txt --chromosome Y --cpg-range [1085,1200] [1220,1280] [2300,2395]

$ python mrv.py --data-file /home/someone/data.txt --to-visualize-file /home/someone/chromosomes_CpGs.txt"""

    parser = argparse.ArgumentParser(prog=prog,
                                     usage=usage,
                                     description=description,
                                     epilog=epilog,
                                     formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument("--data-file",
                        dest="datafile",
                        required=True,
                        type=str,
                        help="file path, specify the path of a text file containing long-read sequencing reads data")

    parser.add_argument("--chromosome",
                        dest="chromosome",
                        type=str,
                        help="text, specify chromosome name")

    parser.add_argument("--output-dir",
                        dest="outputdir",
                        default=".",
                        type=str,
                        help="directory path, specify the path of the output directory default: current directory")

    parser.add_argument("--version",
                        "-V",
                        dest="version",
                        action="version",
                        version="%(prog)s 1.0.0",
                        help="show mrv version")

    exclusive_group = parser.add_mutually_exclusive_group()

    exclusive_group.add_argument("--cpg-range",
                                 dest="range",
                                 action="extend",
                                 nargs="+",
                                 type=str,
                                 help=textwrap.dedent("""\
                                text, specify the range of CpG sites. The syntax format of RANGE is: [start,end].
                                You can give several RANGEs at once, such as
                                --cpg-range [1085,1200] [1220,1280] [2300,2395]"""))

    exclusive_group.add_argument("--to-visualize-file",
                                 dest="tovisualizefile",
                                 type=str,
                                 help=textwrap.dedent("""\
                                file path, specify the path of a text file. You can put the names of multiple
                                chromosomes and their CpG site intervals to be visualized in this text file, and
                                hand them to mrv to calculate and output multiple visual files. An example of
                                TOVISUALIZEFILE:
                                [Y]
                                5431,9587
                                15680,17506
                                12003,12210
                                80,3327
            
                                [KI270580.1]
                                1154,1669
                                756,1321
                                800,1154"""))

    parsed_args = parser.parse_args(args)

    if "--to-visualize-file" in args:
        parsed_args.chromosome = None
    else:
        if "--chromosome" not in args and "--cpg-range" not in args:
            print("""usage: python mrv.py <--data-file DATAFILE> {<--chromosome CHROMOSOME> <--cpg-range RANGE \
[RANGE ...]> | <--to-visualize-file TOVISUALIZEFILE>} [--output-dir OUTPUTDIR] 
        | --version 
        | --help
mrv.py: error: the following arguments are required: --chromosome and --cpg-range""")
            return False
        elif "--chromosome" not in args:
            print("""usage: python mrv.py <--data-file DATAFILE> {<--chromosome CHROMOSOME> <--cpg-range RANGE \
[RANGE ...]> | <--to-visualize-file TOVISUALIZEFILE>} [--output-dir OUTPUTDIR] 
        | --version 
        | --help
mrv.py: error: the following arguments are required: --chromosome""")
            return False
        elif "--cpg-range" not in args:
            print("""usage: python mrv.py <--data-file DATAFILE> {<--chromosome CHROMOSOME> <--cpg-range RANGE \
[RANGE ...]> | <--to-visualize-file TOVISUALIZEFILE>} [--output-dir OUTPUTDIR] 
        | --version 
        | --help
mrv.py: error: the following arguments are required: --cpg-range""")
            return False

    if DEBUG:
        print("------DEBUG: parsed_args------")
        print(parsed_args)
        print("------parsed_args :DEBUG------")

    if not os.path.exists(parsed_args.datafile):
        print("""usage: python mrv.py <--data-file DATAFILE> {<--chromosome CHROMOSOME> <--cpg-range RANGE \
[RANGE ...]> | <--to-visualize-file TOVISUALIZEFILE>} [--output-dir OUTPUTDIR] 
        | --version 
        | --help
mrv.py: error: %s does not exist""" % parsed_args.datafile)
        return False
    elif not os.path.isfile(parsed_args.datafile):
        print("""usage: python mrv.py <--data-file DATAFILE> {<--chromosome CHROMOSOME> <--cpg-range RANGE \
[RANGE ...]> | <--to-visualize-file TOVISUALIZEFILE>} [--output-dir OUTPUTDIR] 
        | --version 
        | --help
mrv.py: error: %s is not a file""" % parsed_args.datafile)
        return False

    global DATA_FILE
    DATA_FILE = parsed_args.datafile

    def add_to_chromosomes_cpg_ranges(ch, ra):
        global CHROMOSOMES_CpG_RANGES
        CHROMOSOMES_CpG_RANGES[ch] = ra

    if parsed_args.range is not None:
        ranges = []
        rc = re.compile(r"^\[\d+,\d+\]$")
        for range_str in parsed_args.range:
            matched_obj = rc.match(range_str)
            if matched_obj is not None:
                li_str = matched_obj.group()
                li = eval(li_str)
                if li[0] == li[1]:
                    print("""usage: python mrv.py <--data-file DATAFILE> {<--chromosome CHROMOSOME> <--cpg-range \
RANGE [RANGE ...]> | <--to-visualize-file TOVISUALIZEFILE>} [--output-dir OUTPUTDIR] 
        | --version 
        | --help
mrv.py: error: %s is not a valid RANGE, the left and right endpoints of an interval should not be equal.""" %
                          range_str)
                    return False
                elif li[0] > li[1]:
                    print("""usage: python mrv.py <--data-file DATAFILE> {<--chromosome CHROMOSOME> <--cpg-range \
RANGE [RANGE ...]> | <--to-visualize-file TOVISUALIZEFILE>} [--output-dir OUTPUTDIR] 
        | --version 
        | --help
mrv.py: error: %s is not a valid RANGE, the left endpoints of the interval should be less than the right endpoints.""" %
                          range_str)
                    return False
                ranges.append(li)
            else:
                print("""usage: python mrv.py <--data-file DATAFILE> {<--chromosome CHROMOSOME> <--cpg-range RANGE \
[RANGE ...]> | <--to-visualize-file TOVISUALIZEFILE>} [--output-dir OUTPUTDIR] 
        | --version 
        | --help
mrv.py: error: %s is not a valid RANGE, the correct format of RANGE is such as [5431,9587]. Notice: there is no \
space between [ and ], and the numbers should be integer.""" % range_str)
                return False
        add_to_chromosomes_cpg_ranges(parsed_args.chromosome, ranges)

    if parsed_args.tovisualizefile is not None:
        if not os.path.exists(parsed_args.tovisualizefile):
            print("""usage: python mrv.py <--data-file DATAFILE> {<--chromosome CHROMOSOME> <--cpg-range RANGE \
[RANGE ...]> | <--to-visualize-file TOVISUALIZEFILE>} [--output-dir OUTPUTDIR] 
        | --version 
        | --help
mrv.py: error: %s does not exist""" % parsed_args.tovisualizefile)
            return False
        elif not os.path.isfile(parsed_args.tovisualizefile):
            print("""usage: python mrv.py <--data-file DATAFILE> {<--chromosome CHROMOSOME> <--cpg-range RANGE \
[RANGE ...]> | <--to-visualize-file TOVISUALIZEFILE>} [--output-dir OUTPUTDIR] 
        | --version 
        | --help
mrv.py: error: %s is not a file""" % parsed_args.tovisualizefile)
            return False

        tovisualizefile_example = textwrap.dedent("""\
                                        [Y]
                                        5431,9587
                                        15680,17506
                                        12003,12210
                                        80,3327
                                        
                                        [KI270580.1]
                                        1154,1669
                                        756,1321
                                        800,1154""")

        config = configparser.ConfigParser(delimiters=",")

        try:
            config.read(parsed_args.tovisualizefile, encoding="UTF-8")
        except configparser.ParsingError as e:
            print("""usage: python mrv.py <--data-file DATAFILE> {<--chromosome CHROMOSOME> <--cpg-range RANGE \
[RANGE ...]> | <--to-visualize-file TOVISUALIZEFILE>} [--output-dir OUTPUTDIR] 
        | --version 
        | --help
mrv.py: error: %s syntax error. Here is a correct syntax example of TOVISUALIZEFILE：\n%s""" %
                  (parsed_args.tovisualizefile, tovisualizefile_example))
            return False

        chromosomes = config.sections()
        if len(chromosomes) == 0:
            print("""usage: python mrv.py <--data-file DATAFILE> {<--chromosome CHROMOSOME> <--cpg-range RANGE \
[RANGE ...]> | <--to-visualize-file TOVISUALIZEFILE>} [--output-dir OUTPUTDIR] 
        | --version 
        | --help
mrv.py: error: %s syntax error: %s is an empty file. Here is a correct syntax example of TOVISUALIZEFILE：\n%s""" %
                  (parsed_args.tovisualizefile, parsed_args.tovisualizefile, tovisualizefile_example))
            return False

        if DEBUG:
            print("------DEBUG: chromosomes------")
            print(chromosomes)
            print("------chromosomes :DEBUG------")

        for chromosome in chromosomes:
            ranges = []

            ranges_list = config.items(chromosome)

            for element in ranges_list:
                try:
                    e1 = int(element[0])
                    e2 = int(element[1])
                except ValueError as e:
                    print("""usage: python mrv.py <--data-file DATAFILE> {<--chromosome CHROMOSOME> <--cpg-range \
RANGE [RANGE ...]> | <--to-visualize-file TOVISUALIZEFILE>} [--output-dir OUTPUTDIR] 
        | --version 
        | --help
mrv.py: error: %s syntax error. The left and right endpoints of the interval should be integers. Here is \
a correct syntax example of TOVISUALIZEFILE：\n%s""" % (parsed_args.tovisualizefile, tovisualizefile_example))
                    return False
                if e1 < e2:
                    li = [e1, e2]
                elif e1 > e2:
                    li = [e2, e1]
                else:
                    print("""usage: python mrv.py <--data-file DATAFILE> {<--chromosome CHROMOSOME> <--cpg-range \
RANGE [RANGE ...]> | <--to-visualize-file TOVISUALIZEFILE>} [--output-dir OUTPUTDIR] 
        | --version 
        | --help
mrv.py: error: %s. The left and right endpoints of the interval should not be equal. Here is \
a correct syntax example of TOVISUALIZEFILE：\n%s""" % (parsed_args.tovisualizefile, tovisualizefile_example))
                    return False

                ranges.append(li)

            if len(ranges_list) == 0:
                continue

            add_to_chromosomes_cpg_ranges(chromosome, ranges)

    global OUT_PUTDIR
    OUT_PUTDIR = parsed_args.outputdir

    return True


def get_file_row_count(file_path):
    with open(file_path, 'r', encoding="UTF-8") as f:
        count = 0
        for i in f:
            count += 1
        return count


def get_file_sha256sum(file_path):
    with open(file_path, "rb") as f:
        sha256sum = hashlib.new("sha256", b"")
        while True:
            data = f.read(64 * 1024)
            if not data:
                break
            sha256sum.update(data)
        return sha256sum.hexdigest()


def is_preprocessed(file_path, dir_path="."):
    mrv_output_dir_path = "%s/%s/" % (dir_path.rstrip("/"), "mrv_output")

    global DATA_FILE_ROW_COUNT
    DATA_FILE_ROW_COUNT = get_file_row_count(file_path)

    global DATA_FILE_SHA256SUM
    DATA_FILE_SHA256SUM = get_file_sha256sum(file_path)

    if DEBUG:
        print("------DEBUG: DATA_FILE_SHA256SUM------")
        print(DATA_FILE_SHA256SUM)
        print("------DATA_FILE_SHA256SUM :DEBUG------")

    global PREPROCESS_DIR
    PREPROCESS_DIR = "%s/%s/" % (mrv_output_dir_path.rstrip("/"), DATA_FILE_SHA256SUM)
    if not os.path.exists(PREPROCESS_DIR):
        return False

    if DATA_FILE_SHA256SUM not in os.listdir(mrv_output_dir_path):
        return False
    else:
        meta_data_file_path = "%s/%s/%s" % (PREPROCESS_DIR.rstrip("/"), "data", "meta_data")

        if not os.path.exists(meta_data_file_path):
            return False

        meta_data_file = open(meta_data_file_path, 'rb')
        count = pickle.load(meta_data_file)
        meta_data_file.close()

        data_file_row_count = DATA_FILE_ROW_COUNT

        if count == data_file_row_count:
            return True
        else:
            return False


def create_output_directory(file_path, dir_path="."):
    if dir_path != ".":
        if not os.path.exists(dir_path):
            try:
                os.makedirs(dir_path)
            except OSError as e:
                print("""usage: python mrv.py <--data-file DATAFILE> {<--chromosome CHROMOSOME> <--cpg-range RANGE \
[RANGE ...]> | <--to-visualize-file TOVISUALIZEFILE>} [--output-dir OUTPUTDIR] 
        | --version 
        | --help
mrv.py: error: can not create directory %s""" % dir_path)
                return False

    mrv_output_dir_path = "%s/%s/" % (dir_path.rstrip("/"), "mrv_output")
    if "mrv_output" not in os.listdir(dir_path):
        try:
            os.makedirs(mrv_output_dir_path)
        except OSError as e:
            print("""usage: python mrv.py <--data-file DATAFILE> {<--chromosome CHROMOSOME> <--cpg-range RANGE \
[RANGE ...]> | <--to-visualize-file TOVISUALIZEFILE>} [--output-dir OUTPUTDIR] 
        | --version 
        | --help
mrv.py: error: can not create directory %s""" % mrv_output_dir_path)
            return False

    sha256sum = DATA_FILE_SHA256SUM

    global PREPROCESS_DIR
    PREPROCESS_DIR = "%s/%s/" % (mrv_output_dir_path.rstrip("/"), sha256sum)
    if sha256sum not in os.listdir(mrv_output_dir_path):
        os.makedirs(PREPROCESS_DIR)
        os.makedirs("%s/%s/" % (PREPROCESS_DIR.rstrip("/"), "data/"))
        os.makedirs("%s/%s/" % (PREPROCESS_DIR.rstrip("/"), "visualization/"))
    else:
        shutil.rmtree(PREPROCESS_DIR)
        os.makedirs(PREPROCESS_DIR)
        os.makedirs("%s/%s/" % (PREPROCESS_DIR.rstrip("/"), "data/"))
        os.makedirs("%s/%s/" % (PREPROCESS_DIR.rstrip("/"), "visualization/"))

    return True


def preprocess_file(file_path):
    if DEBUG:
        print("------DEBUG: PREPROCESS_DIR------")
        print(PREPROCESS_DIR)
        print("------PREPROCESS_DIR :DEBUG------")

    with open(file_path, 'r', encoding="UTF-8") as f:
        chromosomes = []
        sqlite3_conns = {}

        count = 0
        # progress_bar = 0
        pbar = tqdm.tqdm(total=DATA_FILE_ROW_COUNT, desc='Processing', unit="rows", colour="GREEN")
        for each_row in f:
            row = each_row.strip()
            row_split = row.split('\t')

            if len(row_split) < 6:
                print("Data format error! At least 6 columns of data, with each column separated by the tab key.")
                pbar.close()
                return False

            count += 1

            if count == 1:
                continue

            chromosome = row_split[0]
            if chromosome not in chromosomes:
                chromosomes.append(chromosome)

                db = "%s/%s/%s.db" % (PREPROCESS_DIR.rstrip("/"), "data", chromosome)
                conn = sqlite3.connect(db)
                if DEBUG:
                    print("------DEBUG: SQLite db------")
                    print("creating %s..." % db)
                    print("------SQLite db :DEBUG------")

                sqlite3_conns[chromosome] = conn

                cur = conn.cursor()
                cur.execute("""
                    CREATE TABLE '%s'(
                        start int,
                        read_name nvarchar(40),
                        is_methylated boolean DEFAULT 0
                        );""" % chromosome)

                if DEBUG:
                    print("------DEBUG: SQLite table------")
                    print("creating table %s..." % chromosome)
                    print("------SQLite table :DEBUG------")

                is_methylated = 0
                if not row_split[5].startswith('-'):
                    is_methylated = 1

                start = row_split[2]
                read_name = row_split[4]
                cur.execute("insert into '%s'(start,read_name,is_methylated) values(%s,'%s',%s);" %
                            (chromosome, start, read_name, is_methylated))
                conn.commit()
            else:
                is_methylated = 0
                if not row_split[5].startswith('-'):
                    is_methylated = 1

                start = row_split[2]
                read_name = row_split[4]
                sqlite3_conns[chromosome].execute(
                    "insert into '%s'(start,read_name,is_methylated) values(%s,'%s',%s);" % (
                        chromosome, start, read_name, is_methylated))
                sqlite3_conns[chromosome].commit()

            pbar.update(1)
        pbar.close()

        if DEBUG:
            print("------DEBUG: count------")
            print("processed %s rows totally..." % count)
            print("------count :DEBUG------")

        meta_data_file_path = "%s/%s/%s" % (PREPROCESS_DIR.rstrip("/"), "data", "meta_data")
        meta_data_file = open(meta_data_file_path, 'wb')
        pickle.dump(count, meta_data_file)
        meta_data_file.close()

        for key in sqlite3_conns.keys():
            chromosome = key
            conn = sqlite3_conns[key]
            conn.commit()

            if DEBUG:
                print("------DEBUG: SQLite index------")
                print("creating index of %s table..." % chromosome)
                print("------SQLite index :DEBUG------")

            conn.execute("CREATE INDEX start_index ON '%s'(start);" % chromosome)
            conn.commit()
            conn.execute("CREATE INDEX read_name_index ON '%s'(read_name);" % chromosome)
            conn.commit()
            conn.close()

    return True


def preprocess_chromosomes_cpg_ranges(chromosomes_cpg_ranges):
    def preprocess_chromosome_cpg_ranges(chromosome, cpg_ranges):
        if DEBUG:
            print("------DEBUG: chromosome and cpg_ranges------")
            print("%s:%s" % (chromosome, cpg_ranges))
            print("------chromosome and cpg_ranges :DEBUG------")

        old_intervals = portion.empty()
        for e in cpg_ranges:
            interval = portion.closed(e[0], e[1])
            old_intervals = old_intervals.union(interval)

        global CHROMOSOMES_CpG_RANGES
        new_ranges = []
        for e in portion.to_data(old_intervals):
            new_ranges.append([e[1], e[2]])

        if len(new_ranges) == 0:
            CHROMOSOMES_CpG_RANGES[chromosome] = None
            print("There is no CpGs' information of chromosome %s in given ranges." % chromosome)
        else:
            CHROMOSOMES_CpG_RANGES[chromosome] = new_ranges

    global CHROMOSOMES_CpG_RANGES
    for key in CHROMOSOMES_CpG_RANGES.keys():
        dbs_dir = "%s/%s/" % (PREPROCESS_DIR.rstrip("/"), "data")
        if (key + ".db") not in os.listdir(dbs_dir):
            CHROMOSOMES_CpG_RANGES[key] = None
            print("The data of chromosome %s you specified does not exist in the data file." % key)
            continue

        preprocess_chromosome_cpg_ranges(key, CHROMOSOMES_CpG_RANGES[key])

    if DEBUG:
        print("------DEBUG: new CHROMOSOMES_CpG_RANGES------")
        print(CHROMOSOMES_CpG_RANGES)
        print("------new CHROMOSOMES_CpG_RANGES :DEBUG------")

    return True


def visualize(chromosomes_cpg_ranges):
    def visualize_one(chromosome, cpg_ranges):
        print("For chromosome %s and the given ranges:" % chromosome)

        read_names = []
        cpg_positions = []

        db = "%s/%s/%s.db" % (PREPROCESS_DIR.rstrip("/"), "data", chromosome)
        if not os.path.exists(db):
            print("the data of chromosome %s you specified does not exist in the data file." % chromosome)
            return False

        conn = sqlite3.connect(db)
        cur = conn.cursor()

        cpg_ranges_len = len(cpg_ranges)
        if cpg_ranges_len == 0:
            print("there is no CpGs' information of chromosome %s in given ranges." % chromosome)
            return False

        if cpg_ranges_len == 1:
            interval = cpg_ranges[0]
            dql = "SELECT DISTINCT read_name FROM (SELECT DISTINCT read_name,start FROM '%s' WHERE start BETWEEN %s\
                  AND %s ORDER BY start);" % (chromosome, interval[0], interval[1])

            cur.execute(dql)
            for row in cur:
                read_names.append(row[0])
        elif cpg_ranges_len <= 300:
            dql = "SELECT DISTINCT read_name FROM "
            count = 1
            for interval in cpg_ranges:
                if count == 1:
                    dql += "(SELECT DISTINCT read_name,start FROM '%s' WHERE start BETWEEN %s AND %s" % (chromosome,
                                                                                                         interval[0],
                                                                                                         interval[1])
                elif count == cpg_ranges_len:
                    dql += " UNION SELECT DISTINCT read_name,start FROM '%s' WHERE start BETWEEN %s AND %s\
                     ORDER BY start);" % (chromosome, interval[0], interval[1])
                else:
                    dql += " UNION SELECT DISTINCT read_name,start FROM '%s' WHERE start BETWEEN %s AND %s" % (
                        chromosome, interval[0], interval[1])
                count += 1

            cur.execute(dql)
            for row in cur:
                read_names.append(row[0])
        else:
            try:
                cur.execute("DROP TABLE read_names;")
            except sqlite3.OperationalError:
                pass
            cur.execute("CREATE TABLE read_names(start int, read_name nvarchar(40));")
            conn.commit()

            if DEBUG:
                print("------DEBUG: SQLite table------")
                print("creating temporary table read_names...")
                print("------SQLite table :DEBUG------")

            dql = "SELECT DISTINCT read_name FROM "
            count = 1
            times300 = 1
            for interval in cpg_ranges:
                if count == 1:
                    dql += "(SELECT DISTINCT read_name,start FROM '%s' WHERE start BETWEEN %s AND %s" % (chromosome,
                                                                                                         interval[0],
                                                                                                         interval[1])
                elif count == cpg_ranges_len:
                    dql += " UNION SELECT DISTINCT read_name,start FROM '%s' WHERE start BETWEEN %s AND %s\
                         ORDER BY start);" % (chromosome, interval[0], interval[1])
                    cur.execute(dql)
                    for row in cur:
                        cur_temp = conn.cursor()
                        cur_temp.execute("insert into read_names(read_name) values('%s');" % row[0])
                        conn.commit()
                    break
                else:
                    dql += " UNION SELECT DISTINCT read_name,start FROM '%s' WHERE start BETWEEN %s AND %s" % (
                        chromosome, interval[0], interval[1])

                if count >= 300 * times300 - 1:
                    dql += " UNION SELECT DISTINCT read_name,start FROM '%s' WHERE start BETWEEN %s AND %s\
                     ORDER BY start);" % (chromosome, interval[0], interval[1])
                    cur.execute(dql)
                    for row in cur:
                        dql_temp = "SELECT start FROM '%s' WHERE read_name='%s' ORDER BY start LIMIT 1;" % (
                            chromosome, row[0])
                        cur_temp = conn.cursor()
                        cur_temp.execute(dql_temp)
                        for row_temp in cur_temp:
                            the_first_CpG_position = row_temp[0]

                        cur_temp.execute("insert into read_names(start, read_name) values(%d,'%s');" % (
                            the_first_CpG_position, row[0]))
                        conn.commit()

                    times300 += 1

                    dql = "SELECT DISTINCT read_name FROM (SELECT DISTINCT read_name,start FROM '%s' WHERE 1=0 " % chromosome

                count += 1

            cur.execute("CREATE INDEX read_names_table_read_name_index ON read_names (read_name);")
            cur.execute("CREATE INDEX read_names_table_start_index ON read_names (start);")
            conn.commit()

            dql = "SELECT DISTINCT read_name FROM read_names ORDER BY start;"
            cur.execute(dql)
            for row in cur:
                read_names.append(row[0])

        if DEBUG:
            print("------DEBUG: read_names------")
            # print(read_names)
            print("length of read_names is %d" % len(read_names))
            print("------read_names :DEBUG------")

        read_names_len = len(read_names)
        if read_names_len == 0:
            print("there is no CpGs' information of chromosome %s in given ranges." % chromosome)
            return False

        where_clause_part = ""
        count = 1
        for r in cpg_ranges:
            if count == 1:
                where_clause_part += " start between %s and %s " % (r[0], r[1])
            else:
                where_clause_part += "or start between %s and %s " % (r[0], r[1])
            count += 1

        if read_names_len == 1:
            read_name = read_names[0]
            dql = "SELECT DISTINCT start FROM '%s' WHERE read_name='%s' and %s ORDER BY start;" % (
                chromosome, read_name, where_clause_part)
            cur.execute(dql)
            for row in cur:
                cpg_positions.append(row[0])
        elif read_names_len <= 300:
            dql = ""
            count = 1
            for read_name in read_names:
                if count == 1:
                    dql += "SELECT DISTINCT start FROM '%s' WHERE read_name='%s' and %s " % (
                        chromosome, read_name, where_clause_part)
                elif count == read_names_len:
                    dql += " UNION SELECT DISTINCT start FROM '%s' WHERE read_name='%s' and %s ORDER BY start;" % (
                        chromosome, read_name, where_clause_part)
                else:
                    dql += " UNION SELECT DISTINCT start FROM '%s' WHERE read_name='%s' and %s " % (
                        chromosome, read_name, where_clause_part)

                count += 1

            cur.execute(dql)
            for row in cur:
                cpg_positions.append(row[0])
        else:
            try:
                cur.execute("DROP TABLE cpg_positions;")
            except sqlite3.OperationalError:
                pass
            cur.execute("CREATE TABLE cpg_positions(start int);")
            conn.commit()

            dql = ""
            count = 1
            times300 = 1
            for read_name in read_names:
                if count == 1:
                    dql += "SELECT DISTINCT start FROM '%s' WHERE read_name='%s' and %s " % (
                        chromosome, read_name, where_clause_part)
                elif count == read_names_len:
                    dql += " UNION SELECT DISTINCT start FROM '%s' WHERE read_name='%s' and %s  ORDER BY start;" % (
                        chromosome, read_name, where_clause_part)
                    for row in cur:
                        cur_temp = conn.cursor()
                        cur_temp.execute("insert into cpg_positions(start) values(%d);" % row[0])
                        conn.commit()
                    break
                else:
                    dql += " UNION SELECT DISTINCT start FROM '%s' WHERE read_name='%s' and %s " % (
                        chromosome, read_name, where_clause_part)

                if count >= 300 * times300 - 1:
                    dql += " UNION SELECT DISTINCT start FROM '%s' WHERE read_name='%s' and %s  ORDER BY start;" % (
                        chromosome, read_name, where_clause_part)
                    cur.execute(dql)
                    for row in cur:
                        cur_temp = conn.cursor()
                        cur_temp.execute("insert into cpg_positions(start) values(%d);" % row[0])
                        conn.commit()

                    times300 += 1

                    dql = "SELECT DISTINCT start FROM '%s' where 1=0 " % chromosome

                count += 1

            cur.execute("CREATE INDEX cpg_positions_table_start_index ON cpg_positions (start);")
            conn.commit()

            dql = "SELECT DISTINCT start FROM cpg_positions ORDER BY start;"
            cur.execute(dql)
            for row in cur:
                cpg_positions.append(row[0])

        if DEBUG:
            print("------DEBUG: cpg_positions------")
            print("length of cpg_positions is %d" % len(cpg_positions))
            print("------cpg_positions :DEBUG------")

        matrix = [[-1] * len(cpg_positions) for _ in range(len(read_names))]

        for i in range(len(read_names)):
            the_first_CpG_position = 0
            dql = "SELECT start FROM '%s' WHERE read_name='%s' and %s ORDER BY start LIMIT 1;" % (
                chromosome, read_names[i], where_clause_part)
            cur.execute(dql)
            for row in cur:
                the_first_CpG_position = row[0]

            the_first_CpG_position_index = cpg_positions.index(the_first_CpG_position)

            dql = "SELECT DISTINCT start, is_methylated FROM '%s' WHERE read_name='%s' and %s  ORDER BY start;" % (
                chromosome, read_names[i], where_clause_part)

            cur.execute(dql)
            j = the_first_CpG_position_index
            for row in cur:
                if j >= len(cpg_positions):
                    if DEBUG:
                        print("------Error:j >= len(cpg_positions)------")
                        print("read_name=%s, row[0]=%d, j=%d" % (read_names[i], row[0], j))
                        print("------j >= len(cpg_positions) :Error------")
                        return False
                    continue

                if row[0] != cpg_positions[j]:
                    pass
                else:
                    if row[1] == 1:
                        matrix[i][j] = 1
                    else:
                        matrix[i][j] = 0
                j += 1

        conn.close()

        indexes = []
        for i in range(len(matrix)):
            if matrix[i].count(1) == 1 and matrix[i].count(1) + matrix[i].count(-1) == len(matrix[i]) or \
                    matrix[i].count(0) == 1 and matrix[i].count(0) + matrix[i].count(-1) == len(matrix[i]):
                indexes.append(i)

        read_names_to_remove = [read_names[i] for i in indexes]
        rows_to_remove = [matrix[i] for i in indexes]

        for obj in read_names_to_remove:
            read_names.remove(obj)
        for obj in rows_to_remove:
            matrix.remove(obj)

        finished_flag = 0
        for j in reversed(list(range(len(cpg_positions)))):
            del_flag = 1
            for i in range(len(read_names)):
                if matrix[i][j] != -1:
                    del_flag = 0
                    finished_flag = 1
                    break

            if finished_flag:
                break

            for i in range(len(read_names)):
                if del_flag == 1:
                    matrix[i].pop()

            cpg_positions.pop()

        if finished_flag == 0:
            print("there is no CpGs' information of chromosome %s in given ranges." % chromosome)
            return False

        visualization_txt_file_path = "%s/%s/%s_%s_%s_visualization.txt" % (
            PREPROCESS_DIR.rstrip("/"), "visualization", chromosome, cpg_positions[0], cpg_positions[-1])
        with open(visualization_txt_file_path, 'w', encoding="UTF-8") as f:
            blank = ""
            for i in range(len(read_names[0])):
                blank += " "
            f.write(blank)
            f.write('\t')
            for cpg in cpg_positions:
                f.write(str(cpg))
                f.write('\t')
            f.write(os.linesep)

            for r in range(len(read_names)):
                f.write(read_names[r])
                f.write('\t')

                i = 0
                for c in matrix[r]:
                    for j in range(len(str(cpg_positions[i])) - 1):
                        f.write(' ')

                    if c == -1:
                        f.write(' ')
                    elif c == 0:
                        f.write('0')
                    else:
                        f.write('1')

                    f.write('\t')

                    i += 1

                f.write(os.linesep)
        print("the visualization txt file is at %s" % visualization_txt_file_path)

        import matplotlib.ticker as ticker
        matplotlib.use('svg')
        if len(cpg_positions) * 2 <= len(read_names):
            fig, ax = plt.subplots(figsize=(len(cpg_positions) * 4 * 0.6, len(read_names) * 0.6),
                                   constrained_layout=True)
        else:
            fig, ax = plt.subplots(figsize=(len(cpg_positions) * 2.5 * 0.6, len(read_names) * 0.6),
                                   constrained_layout=True)

        plt.ylim(0, len(read_names) + 1)
        plt.xlim(0, len(cpg_positions) + 1)

        xticks = range(1, len(cpg_positions) + 1)
        xlabels = cpg_positions
        ax.set_xticks(xticks)
        ax.set_xticklabels(labels=xlabels, rotation=30, ha="center")

        yticks = range(1, len(read_names) + 1)
        ylabels = read_names
        ax.set_yticks(yticks)
        ax.set_yticklabels(labels=ylabels)

        ax.xaxis.set_major_locator(ticker.MultipleLocator(1))
        ax.yaxis.set_major_locator(ticker.MultipleLocator(1))

        def draw_a_read(i):
            start_index = 0
            end_index = -1
            flag = 0
            for j in range(len(cpg_positions)):
                if matrix[i][j] != -1 and flag == 0:
                    start_index = j
                    flag = 1
                    continue

                if matrix[i][j] == -1 and flag == 1:
                    k = j
                    while k < len(cpg_positions) - 1 and matrix[i][k] == -1:
                        k = k + 1
                    if not k < len(cpg_positions) - 1:
                        end_index = j - 1
                        break
                    else:
                        continue

                if j == len(cpg_positions) - 1 and flag == 1:
                    end_index = j

            x = list(range(start_index, end_index + 1))

            markers0 = []
            markers1 = []
            for j in x:
                if matrix[i][j] == -1:
                    pass
                elif matrix[i][j] == 0:
                    markers0.append(j)
                else:
                    markers1.append(j)

            x = [_ + 1 for _ in x]
            markers0 = [_ + 1 for _ in markers0]
            markers1 = [_ + 1 for _ in markers1]

            y = [i + 1 for _ in range(len(x))]

            ax.plot(x, y, linestyle='-', color="black", linewidth=2)

            for maker in markers0:
                ax.plot(maker, y[0], marker='o', markersize=8, markerfacecolor='white', markeredgecolor='black')
            for maker in markers1:
                ax.plot(maker, y[0], marker='o', markersize=8, markerfacecolor='black', markeredgecolor='black')

        for i in range(len(read_names)):
            draw_a_read(i)

        ax.set_title("CpG ranges [%s, %s] of %s Chromosome" % (cpg_positions[0], cpg_positions[-1], chromosome),
                     fontsize=12)
        ax.set_xlabel("CpG sites", fontsize=12)
        ax.set_ylabel("reads", fontsize=12)

        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

        visualization_svg_file_path = "%s/%s/%s_%s_%s_visualization.svg" % (
            PREPROCESS_DIR.rstrip("/"), "visualization", chromosome, cpg_positions[0], cpg_positions[-1])
        fig.savefig(visualization_svg_file_path, dpi=600, format='svg')
        print("the visualization svg file is at %s" % visualization_svg_file_path)

    global CHROMOSOMES_CpG_RANGES
    for key in CHROMOSOMES_CpG_RANGES.keys():
        if CHROMOSOMES_CpG_RANGES[key] is None:
            continue
        visualize_one(key, CHROMOSOMES_CpG_RANGES[key])

    return True


def main():
    if sys.version_info.major < 3:
        print("You current Python version is %d.%d, the mrv require Python 3, Python 3.6+ is better." %
              (sys.version_info.major, sys.version_info.minor))
        sys.exit(1)

    # print("Parsing options...")
    if not parse_args(sys.argv):
        print("Error parsing command line options.")
        sys.exit(1)
    print("Parsing command line options...")
    print("Done.")

    if DEBUG:
        print("------DEBUG: global variables------")
        print("DATA_FILE:")
        print(DATA_FILE)

        print("CHROMOSOMES_CpG_RANGES:")
        print(CHROMOSOMES_CpG_RANGES)

        print("OUT_PUTDIR:")
        print(OUT_PUTDIR)

        print("------global variables :DEBUG------")

    if not is_preprocessed(DATA_FILE, OUT_PUTDIR):
        print("\nCreating output directory...")
        if not create_output_directory(DATA_FILE, OUT_PUTDIR):
            sys.exit(1)
        print("Done.")

        print("\nPreprocessing %s..." % DATA_FILE)
        if not preprocess_file(DATA_FILE):
            sys.exit(1)
        print("Done.")

    print("\nPreprocessing Chromosomes and their CpG ranges...")
    if not preprocess_chromosomes_cpg_ranges(CHROMOSOMES_CpG_RANGES):
        sys.exit(1)
    print("Done.")

    print("\nVisualizing...")
    if not visualize(CHROMOSOMES_CpG_RANGES):
        sys.exit(1)
    print("\nAll have done.")
    sys.exit(0)


if __name__ == "__main__":
    main()
