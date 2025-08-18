def get_header(filename):
    with open(filename) as fh: 
        lines = fh.readlines()
        lines = [line.strip() for line in lines]
    return lines
