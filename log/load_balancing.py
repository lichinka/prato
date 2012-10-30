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
    """
    Prints out the detailed running time for worker 'worker_id',
    returning its total running time.-
    """
    ret_value         = 0.0
    worker_id         = '(%d)' % int (worker_id)
    input_file        = sys.argv[2]
    line_data         = None
    previous_end_time = None
    with open (input_file, 'r') as f:
        for line in f:
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
                    ret_value += line_data['elapsed_time']
                    print ('%(action)s\t%(start_time).10f\t%(elapsed_time).10f\t%(end_time).10f' % line_data)
    return ret_value



if __name__ == "__main__":
    if len (sys.argv) > 1:
        #
        # output the headers of the extracted data
        #
        print ("# action\tstart time\telapsed time\tend time")
        #
        # start extracting the data ...
        #
        NP = int (sys.argv[1])
        fastest_worker = {'id': 0,
                          'time': float ('+inf')}
        for np in range (NP):
            running_time = load_balancing (np)
            if fastest_worker['time'] > running_time:
                fastest_worker['id']   = np
                fastest_worker['time'] = running_time
            print ("---> total running time for worker (%d) is\t%.9f" % (np,
                                                                         running_time))
        print ('---> fastest worker (%(id)d) time is\t%(time).9f' % fastest_worker)
    else:
        print ("Usage:\t%s [number of processes] [log file]" % sys.argv[0])
        print ("Extracts per-worker and fastest running times from an experimental log file.")

