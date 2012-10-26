#!/usr/bin/env python

import sys


def process_line (line):
    """
    Parses a line from the input file, extracting times for
    the load balancing analysis. Returns a dictionary.-
    """
    ret_value = dict ( )
    columns = line.split (':')
    ret_value['end_time'] = float (columns[1])
    ret_value['action']   = columns[2].split ('\t')[1]
    ret_value['elapsed_time'] = columns[2].split ('\t')[2]
    ret_value['elapsed_time'] = float (ret_value['elapsed_time'].split (' ')[0])
    ret_value['start_time']   = ret_value['end_time'] - ret_value['elapsed_time']
    return ret_value


def load_balancing (worker_id):
    worker_id         = '(%d)' % int (worker_id)
    input_file        = sys.argv[2]
    line_data         = None
    previous_end_time = None
    for line in open (input_file):
        if 'TIME' in line:
            if worker_id in line:
                if line_data is not None:
                    previous_end_time = line_data['end_time']
                line_data = process_line (line)
                if previous_end_time is not None:
                    idle_time = line_data['start_time'] - previous_end_time
                    print ('%s %s\t%.9f\t%.9f\t%.9f' % ('Idle',
                                                        worker_id,
                                                        previous_end_time,
                                                        idle_time,
                                                        line_data['start_time']))
                print ('%(action)s\t%(start_time).10f\t%(elapsed_time).10f\t%(end_time).10f' % line_data)


if __name__ == "__main__":
    if len (sys.argv) > 1:
        load_balancing (sys.argv[1])
    else:
        print "Usage:\t$0 [worker id] [log files ...]"
        print "Extracts the load-balancing time from a bunch of experiment log files."

