#!/usr/bin/env python3
import string
import os
import sys
def check_input(args):
    fh = open(args[0], 'r')
    return (fh)

def reletter_chains(fh):
    chainlist = list(
        string.digits[::-1] + string.ascii_lowercase[::-1] + string.ascii_uppercase[::-1]
    )  # 987...zyx...cbaZYX...BCA.
    max_chains = len(chainlist)
    
    prev_chain='oogabooga' #there will never be a first chain with this name 
    curchain = chainlist.pop()
    terline='TER'+' '*77
    for line in fh:
        if line.startswith('ATOM'):
            if prev_chain != 'oogabooga':
                if prev_chain != line[21]:
                    curchain = chainlist.pop()
            prev_chain = line[21]
            line = line[:21] + curchain + line[22:]
            yield line

        elif line.startswith('TER'):
            yield terline+'\n'

def main():
    pdbfh = check_input(sys.argv[1:])
    new_pdb = reletter_chains(pdbfh)

    try:
        _buffer = []
        _buffer_size = 5000  # write N lines at a time
        for lineno, line in enumerate(new_pdb):
            if not (lineno % _buffer_size):
                sys.stdout.write(''.join(_buffer))
                _buffer = []
            _buffer.append(line)

        sys.stdout.write(''.join(_buffer))
        sys.stdout.flush()
    except IOError:
        pass
    pdbfh.close()
    sys.exit(0)


if __name__ == '__main__':
    main()
