#!python3
# split_to_function_files.m

# Corban Swain, 2018

import re
import os
import pprint as pp


if __name__ == "__main__":
    # Inputs
    in_filename = os.path.join(os.getcwd(), 'register.m')
    out_dir = os.path.join(os.getcwd())

    # Script
    function_keyword = 'function'
    block_keywords = {'for', 'if', 'try', 'switch', 'parfor', 'while',
                      function_keyword}
    end_keyword = 'end'

    discard_out_file_key = '__discard__'
    out_files = {discard_out_file_key: ['%% DISCARDED LINES', ]}

    func_name_pattern = re.compile(
        r'function\s+(.*=\s*)?(?P<function_name>\w+)\s*?(\(.*\))?$')

    print('Beginning automatic MATLAB function splitter: ')
    print('\t Loading from MATLAB file: %s' % in_filename)
    with open(in_filename, 'r') as fle:
        lines = fle.readlines()

        # FIXME - implement peek_next_full_line which will incorporate the `...`
        #         syntax into its logic.

        def peek_next():
            return lines[0]

        def pop_next():
            return lines.pop(0)

        current_function = None
        first_function = None
        level = 0
        while lines:
            if current_function:
                line = pop_next()
                out_files[current_function].append(line)

                if line.strip():
                    first_word = line.split()[0].strip()
                    if first_word == end_keyword:
                        if level is 0:
                            current_function = None
                        else:
                            level -= 1
                    elif first_word in block_keywords:
                        # FIXME - should check if the line ends with an `end`
                        level += 1
            else:
                if peek_next().strip().startswith(function_keyword):
                    function_line = peek_next().strip()
                    match = func_name_pattern.match(function_line)
                    if match:
                        current_function = match.group('function_name')
                        out_files[current_function] = [pop_next(), ]
                        if not first_function:
                            first_function = current_function
                    else:
                        out_files[discard_out_file_key].append(pop_next())
                else:
                    out_files[discard_out_file_key].append(pop_next())

    print('Identified functions:')
    pp.pprint(list(out_files.keys()))
    for function_name, file_lines in out_files.items():
        new_file_name = os.path.join(out_dir, function_name + '.m')
        with open(new_file_name, 'w') as fle:
            lines = file_lines
            if function_name == first_function:
                lines += ['\n', ] * 2 + ['% The following functions were split '
                                         'into separate files\n', ]
                all_funcs = [f for f in out_files.keys()
                             if f not in [first_function, discard_out_file_key]]
                lines += ['%   ' + f + '\n' for f in all_funcs]
                lines += ['\n', ] * 2 + ['% Discarded lines were placed into a '
                                         'file named \"__discard__.m\"\n', ]
            fle.writelines(lines)
