def write_header(filename, header):
    with open(filename, 'w') as f:
        f.write('\n'.join(header) + '\n')
