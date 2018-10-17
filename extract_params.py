#!python3
# extract_params.py

from pprint import pprint
import datetime as dt
import os
import numpy as np
import scipy.io


def get_final_params(logfile):
    def is_param_line(l):
        return ':' in l and 'a href' not in l

    def is_param_block_start(l):
        return 'param =' in l

    def get_open_file():
        return open(logfile, 'r')

    with get_open_file() as fle:
        param_start_line = 0
        for i, line in enumerate(fle):
            if is_param_block_start(line):
                param_start_line = i

        did_find_param_block = False
        did_begin_getting_lines = False
        did_finish_getting_lines = False
        param_line_cache = {}

    with get_open_file() as fle:
        for i, line in enumerate(fle):
            if i == param_start_line:
                did_find_param_block = True

            if did_find_param_block and is_param_line(line) \
                    and not did_begin_getting_lines:
                did_begin_getting_lines = True

            if did_begin_getting_lines:
                if is_param_line(line):
                    lst = line.split(':')
                    key = lst[0].strip()
                    val = ''.join(lst[1:]).strip()
                    param_line_cache[key] = val
                else:
                    did_finish_getting_lines = True

            if did_finish_getting_lines:
                break

    return param_line_cache


def save_transform_params(params, save_dir):
    selected_params = dict(offset='offset',
                           trans='trans',
                           rot='rot',
                           centroid='centroid')
    numpy_vars_out = {}
    for save_name, param_key in selected_params.items():
        raw_val = list(params[param_key])
        striped_val = raw_val[1:] if raw_val[0] in ['[', ] else raw_val
        striped_val = striped_val[:-1] if striped_val[-1] in [']', ] \
            else striped_val
        striped_val = ''.join(striped_val)
        split_val = striped_val.split(' ')
        numpy_vars_out[save_name] = np.array([float(s) for s in split_val])

    numpy_vars_out['metaCreationDate'] = str(dt.datetime.now())
    numpy_vars_out['metaFileSource'] = 'Python Extraction'

    save_path = os.path.join(save_dir, 'register_params.mat')
    scipy.io.savemat(save_path, numpy_vars_out, oned_as='row')
    print('Saved extracted params to \n\t"%s"' % save_path)


if __name__ == '__main__':
    cwd = os.getcwd()
    base_dir = os.path.join(cwd,
                            '..',
                            '..',
                            'data',
                            'pipeline_sim',
                            'run_1_(initial_tests)_on_181005_at_1830')

    for subdir in next(os.walk(base_dir))[1]:
        register_dir = os.path.join(base_dir, subdir, 'reg_output')
        logfile_name = [f for f in os.listdir(register_dir)
                        if f.endswith('.log')]
        if len(logfile_name) > 1:
            print('WARNING, too many logfiles in the directory:\n\t "%s"' %
                  register_dir)
        logfile_name = logfile_name[0]
        logfile_path = os.path.join(register_dir, logfile_name)
        param_dict = get_final_params(logfile_path)
        save_transform_params(param_dict, register_dir)


